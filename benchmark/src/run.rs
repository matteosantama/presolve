use crate::results::{Case, Measurement, Metadata, RunWriter};
use crate::{Kind, PoolMode, Preset, Result, Selection, SparsificationMode, Tuning, data};
use presolve::{
    Presolver,
    problem::ProblemData,
    result::Outcome,
    settings::{Rules, Settings},
};
use std::hint::black_box;
use std::path::Path;
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

// Exhaustive initialization makes a new library rule require an entry here.
macro_rules! rules {
    ($($field:ident),+ $(,)?) => {
        [
            ("all", Rules { $($field: true),+ }),
            $((stringify!($field), Rules { $field: true, ..Rules::none() })),+
        ]
    };
}

const RULES: [(&str, Rules); 15] = rules!(
    fixed_variables,
    empty_columns,
    empty_rows,
    dual_fixing,
    singleton_rows,
    singleton_columns,
    doubleton_equalities,
    short_equalities,
    bound_propagation,
    redundant_bounds,
    parallel_rows,
    parallel_columns,
    sparsification,
    cones,
);

pub fn settings(rule: &str, threads: usize, tuning: &Tuning) -> Result<Settings> {
    let rules = RULES
        .iter()
        .find(|(name, _)| *name == rule)
        .ok_or_else(|| {
            format!(
                "unknown rule '{rule}'; choose {}",
                RULES.map(|(name, _)| name).join(", ")
            )
        })?
        .1;
    let mut settings = match tuning.preset {
        Preset::Default | Preset::Fill => Settings::default(),
        Preset::Aggressive | Preset::Unrestricted => {
            Settings::aggressive(std::time::Duration::from_secs(2))
        }
    };
    settings.rules = rules;
    settings.threads = threads;
    if matches!(tuning.preset, Preset::Fill) {
        settings.substitution_fill = usize::MAX;
        settings.allow_hessian_growth = true;
    }
    if matches!(tuning.preset, Preset::Unrestricted) {
        settings.equalities.max_row_length = usize::MAX;
        settings.equalities.max_column_length = usize::MAX;
    }
    if let Some(ms) = tuning.time_limit_ms {
        settings.time_limit = std::time::Duration::from_millis(ms);
    }
    if let Some(n) = tuning.equality_row_limit {
        settings.equalities.max_row_length = n;
    }
    if let Some(n) = tuning.equality_column_limit {
        settings.equalities.max_column_length = n;
    }
    if let Some(mode) = tuning.sparsification {
        settings.rules.sparsification &= !matches!(mode, SparsificationMode::Off);
        settings.sparsification.allow_auxiliary_variables = matches!(mode, SparsificationMode::All);
    }
    Ok(settings)
}

fn bound_sides(bounds: impl IntoIterator<Item = presolve::problem::Bounds>) -> usize {
    bounds
        .into_iter()
        .map(|b| usize::from(b.lower.is_finite()) + usize::from(b.upper.is_finite()))
        .sum()
}

pub fn measure(
    input: ProblemData,
    settings: &Settings,
    timed: bool,
    mode: PoolMode,
) -> Result<Measurement> {
    let before_bound_sides = Some(bound_sides(input.variable_bounds.iter().copied()));
    let ready = if mode == PoolMode::Reused {
        let presolver = Presolver::new(settings.clone())?;
        if timed {
            drop(presolver.presolve(input.clone()));
        }
        Some(presolver)
    } else {
        None
    };
    // Loading, warm-up, cloning, and destruction of both the returned model and
    // executor stay outside timing. Cold mode includes Presolver construction.
    let input = black_box(input);
    let start = timed.then(Instant::now);
    let presolver = match ready {
        Some(presolver) => presolver,
        None => Presolver::new(black_box(settings).clone())?,
    };
    let result = presolver.presolve(input);
    let elapsed_ns = start.map(|start| start.elapsed().as_nanos().min(u64::MAX as u128) as u64);
    let after_bound_sides = match &result.outcome {
        Outcome::Unchanged(p) => Some(bound_sides(
            (0..p.variable_count()).map(|j| p.variable_bounds(j)),
        )),
        Outcome::Reduced(r) => Some(bound_sides(
            (0..r.problem.variable_count()).map(|j| r.problem.variable_bounds(j)),
        )),
        Outcome::Solved(_) => Some(0),
        _ => None,
    };
    let outcome = match &result.outcome {
        Outcome::Unchanged(_) => "unchanged",
        Outcome::Reduced(_) => "reduced",
        Outcome::Solved(_) => "solved",
        Outcome::Infeasible(_) => "infeasible",
        Outcome::Unbounded(_) => "unbounded",
    };
    Ok(Measurement {
        elapsed_ns,
        before_bound_sides,
        after_bound_sides,
        outcome: outcome.into(),
        before: result.stats.before.into(),
        after: result.stats.after.map(Into::into),
        time_limit_reached: result.stats.time_limit_reached,
    })
}

fn trial(
    executable: &Path,
    path: &Path,
    rule: &str,
    threads: usize,
    mode: PoolMode,
    tuning: &Tuning,
) -> Result<Measurement> {
    let mut command = Command::new(executable);
    command.arg("worker").arg(path).arg(rule);
    command.arg("--threads").arg(threads.to_string());
    command.arg("--pool-mode").arg(match mode {
        PoolMode::Cold => "cold",
        PoolMode::Reused => "reused",
    });
    command.arg("--preset").arg(tuning.preset.as_str());
    for (flag, value) in [
        (
            "--time-limit-ms",
            tuning.time_limit_ms.map(|v| v.to_string()),
        ),
        (
            "--equality-row-limit",
            tuning.equality_row_limit.map(|v| v.to_string()),
        ),
        (
            "--equality-column-limit",
            tuning.equality_column_limit.map(|v| v.to_string()),
        ),
        (
            "--sparsification",
            tuning.sparsification.map(|v| v.as_str().to_owned()),
        ),
    ] {
        if let Some(value) = value {
            command.arg(flag).arg(value);
        }
    }
    let output = command.output()?;
    if !output.status.success() {
        return Err(format!(
            "worker {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        )
        .into());
    }
    serde_json::from_slice(&output.stdout).map_err(|e| format!("invalid worker output: {e}").into())
}

// Stable fingerprint to detect changed input files, not a security checksum.
fn fingerprint(path: &Path) -> Result<String> {
    let hash = std::fs::read(path)?
        .into_iter()
        .fold(0xcbf29ce484222325u64, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
        });
    Ok(format!("fnv1a64:{hash:016x}"))
}

pub fn run(root: &Path, kind: Kind, selection: Selection, trials: usize) -> Result<()> {
    if let Some(rule) = &selection.rule {
        settings(rule, selection.threads, &selection.tuning)?;
    }
    let rules: Vec<_> = RULES
        .iter()
        .filter(|(name, _)| selection.rule.as_deref().is_none_or(|rule| rule == *name))
        .collect();
    let mut problems = Vec::new();
    for suite in data::SUITES {
        for path in data::problems(suite)? {
            let name = path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .trim_end_matches(".mps.gz");
            let id = format!("{suite}/{name}");
            if selection
                .problem
                .as_deref()
                .is_none_or(|filter| filter == name || filter == id)
            {
                problems.push((id, path));
            }
        }
    }
    if problems.is_empty() {
        return Err("no matching problem; use an exact name or suite/name".into());
    }
    if selection.problem.is_some() && problems.len() != 1 {
        return Err("problem name is ambiguous; specify suite/name".into());
    }
    let host = Command::new("hostname")
        .output()
        .ok()
        .filter(|out| out.status.success())
        .map(|out| String::from_utf8_lossy(&out.stdout).trim().to_owned())
        .unwrap_or_default();
    let config = settings("all", selection.threads, &selection.tuning)?;
    let metadata = Metadata {
        version: 1,
        name: selection.name,
        kind,
        trials,
        created_unix_seconds: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs(),
        machine: format!("{}/{}/{host}", std::env::consts::OS, std::env::consts::ARCH),
        // Execution width may differ in a valid scaling comparison. Keep it
        // separate from the algorithm settings, which must still match.
        settings: format!(
            "{:?}",
            Settings {
                threads: 1,
                ..config
            }
        ),
        threads: selection.threads,
        pool_mode: selection.pool_mode,
    };
    let mut writer = RunWriter::create(root, metadata)?;
    let executable = std::env::current_exe()?;
    let total = problems.len() * rules.len();
    let mut completed = 0;
    let mut failed = 0;
    for (problem, path) in problems {
        let hash = fingerprint(&path)?;
        // Size measurements can share parsed input; only the timing path needs
        // a fresh process and freshly prepared input for each call.
        let input = (kind == Kind::Size).then(|| data::read(&path));
        for &(rule, _) in &rules {
            eprintln!(
                "[{}/{}] {problem}/{rule} ({trials} trial{})",
                completed + 1,
                total,
                if trials == 1 { "" } else { "s" }
            );
            let mut case = Case {
                input_hash: hash.clone(),
                measurements: Vec::with_capacity(trials),
                errors: Vec::new(),
            };
            for at in 0..trials {
                let result = match &input {
                    Some(Ok(input)) => measure(
                        input.clone(),
                        &settings(rule, selection.threads, &selection.tuning)?,
                        false,
                        selection.pool_mode,
                    ),
                    Some(Err(e)) => Err(e.to_string().into()),
                    None => trial(
                        &executable,
                        &path,
                        rule,
                        selection.threads,
                        selection.pool_mode,
                        &selection.tuning,
                    ),
                };
                match result {
                    Ok(measurement) => {
                        if measurement.time_limit_reached {
                            case.errors
                                .push(format!("trial {}: presolve time limit reached", at + 1));
                        }
                        case.measurements.push(measurement);
                    }
                    Err(e) => case.errors.push(format!("trial {}: {e}", at + 1)),
                }
            }
            if !case.errors.is_empty() {
                failed += 1;
                for error in &case.errors {
                    eprintln!("  {error}");
                }
            }
            writer.case(format!("{problem}/{rule}"), case)?;
            completed += 1;
        }
    }
    writer.finish()?;
    println!(
        "Saved {completed} cases to {} ({failed} failed cases).",
        writer.path.display()
    );
    if failed > 0 {
        Err("some cases failed; their errors were saved with the results".into())
    } else {
        Ok(())
    }
}
