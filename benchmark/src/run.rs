//! Presolve every instance of a corpus under a named settings profile.

use crate::corpus::Instance;
use crate::mps;
use crate::snapshot::{Presolved, Record, Run};
use presolve::{Presolver, Settings};
use rayon::prelude::*;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::time::{Duration, Instant};

/// Named settings. Every profile runs presolve on one thread with no time
/// limit, so its records depend only on the input and the code.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Profile {
    Default,
    Aggressive,
}
impl Profile {
    pub fn settings(self) -> Settings {
        let settings = match self {
            Self::Default => Settings {
                time_limit: Duration::MAX,
                ..Settings::default()
            },
            Self::Aggressive => Settings::aggressive(Duration::MAX),
        };
        Settings {
            threads: 1,
            ..settings
        }
    }
}

/// One instance's record and the time its presolve call took, excluding
/// loading. Failed instances report zero time.
pub struct Outcome {
    pub record: Record,
    pub elapsed: Duration,
}

/// Run `instances` on `jobs` worker threads, zero meaning one per core.
/// Each instance's presolve call is itself serial.
pub fn run(instances: &[Instance], profile: Profile, jobs: usize) -> Result<Vec<Outcome>, String> {
    let presolver = Presolver::new(profile.settings()).map_err(|e| e.to_string())?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .map_err(|e| e.to_string())?;
    Ok(pool.install(|| {
        instances
            .par_iter()
            .map(|instance| run_one(&presolver, instance))
            .collect()
    }))
}

fn run_one(presolver: &Presolver, instance: &Instance) -> Outcome {
    let (run, elapsed) = match mps::read(&instance.path) {
        Err(e) => (
            Run::LoadError {
                message: e.to_string(),
            },
            Duration::ZERO,
        ),
        Ok(problem) => {
            let cones = problem.cones.len();
            let start = Instant::now();
            let result = catch_unwind(AssertUnwindSafe(|| presolver.presolve(problem)));
            let elapsed = start.elapsed();
            match result {
                Ok(result) => (
                    Run::Presolved(Box::new(Presolved::new(&result, cones))),
                    elapsed,
                ),
                Err(payload) => (
                    Run::Panic {
                        message: panic_message(payload.as_ref()),
                    },
                    Duration::ZERO,
                ),
            }
        }
    };
    Outcome {
        record: Record {
            instance: instance.id(),
            run,
        },
        elapsed,
    }
}

fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        s.to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "non-string panic payload".into()
    }
}
