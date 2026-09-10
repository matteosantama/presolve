mod data;
mod report;
mod results;
mod run;

use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Kind {
    Time,
    Size,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum PoolMode {
    /// Include Presolver construction in each measurement.
    #[default]
    Cold,
    /// Reuse a Presolver after an untimed call on the same problem.
    Reused,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum Preset {
    #[default]
    Default,
    /// Lift the substitution fill cap and allow Hessian growth.
    Fill,
    /// Use Settings::aggressive with a two-second budget unless overridden.
    Aggressive,
    /// Also remove the aggressive preset's equality row and column length caps.
    Unrestricted,
}
impl Preset {
    fn as_str(self) -> &'static str {
        match self {
            Self::Default => "default",
            Self::Fill => "fill",
            Self::Aggressive => "aggressive",
            Self::Unrestricted => "unrestricted",
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SparsificationMode {
    All,
    Equalities,
    Off,
}
impl SparsificationMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Equalities => "equalities",
            Self::Off => "off",
        }
    }
}

#[derive(Clone, Debug, Default, Args)]
struct Tuning {
    #[arg(long, value_enum, default_value_t = Preset::Default)]
    preset: Preset,
    /// Override the per-call soft time budget in milliseconds.
    #[arg(long)]
    time_limit_ms: Option<u64>,
    #[arg(long)]
    equality_row_limit: Option<usize>,
    #[arg(long)]
    equality_column_limit: Option<usize>,
    #[arg(long)]
    equality_pivot_relative: Option<f64>,
    #[arg(long)]
    equality_pivot_attempts: Option<usize>,
    #[arg(long, action = clap::ArgAction::Set)]
    equality_cost_aware: Option<bool>,
    #[arg(long)]
    propagation_relative_gain: Option<f64>,
    #[arg(long)]
    propagation_gain_factor: Option<f64>,
    #[arg(long, value_enum)]
    sparsification: Option<SparsificationMode>,
}

impl Kind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Time => "time",
            Self::Size => "size",
        }
    }
}

#[derive(Parser)]
#[command(
    name = "benchmark",
    about = "Measure presolve time and model size",
    version
)]
struct Cli {
    /// Directory containing named time/size runs.
    #[arg(long, global = true, default_value = concat!(env!("CARGO_MANIFEST_DIR"), "/results"))]
    results_dir: PathBuf,
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Run every selected problem/rule and save results. Requires a release build.
    Run {
        #[command(subcommand)]
        mode: RunCommand,
    },
    /// Print differences between two saved runs. Negative deltas mean less time/size.
    Compare {
        kind: Kind,
        name_1: String,
        name_2: String,
    },
    /// Execute a single timing trial; used by the parent process.
    #[command(hide = true)]
    Worker {
        path: PathBuf,
        rule: String,
        #[command(flatten)]
        tuning: Tuning,
        #[arg(long, default_value_t = 1)]
        threads: usize,
        #[arg(long, value_enum, default_value_t = PoolMode::Cold)]
        pool_mode: PoolMode,
    },
}

#[derive(Subcommand)]
enum RunCommand {
    /// Measure one presolve call per fresh process; excludes loading and process startup.
    Time {
        #[command(flatten)]
        selection: Selection,
        #[arg(long, default_value_t = 100, value_parser = clap::value_parser!(u32).range(1..))]
        trials: u32,
    },
    /// Record the outcome and dimensions/nonzeros before and after presolve, once per case.
    Size {
        #[command(flatten)]
        selection: Selection,
    },
}

#[derive(Args)]
struct Selection {
    #[command(flatten)]
    tuning: Tuning,
    /// Unique run name; existing results are never overwritten.
    #[arg(long)]
    name: String,
    /// Exact problem name (AFIRO) or suite/name (netlib/AFIRO). Default: all problems.
    #[arg(long)]
    problem: Option<String>,
    /// Exact rule name (sparsification, singleton_rows, ...); 'all' is the full pipeline.
    #[arg(long)]
    rule: Option<String>,
    /// Presolve threads: 1 is serial; 0 lets Rayon select automatically.
    #[arg(long, default_value_t = 1)]
    threads: usize,
    /// Include pool setup, or measure a reused pool after an untimed warm-up.
    #[arg(long, value_enum, default_value_t = PoolMode::Cold)]
    pool_mode: PoolMode,
}

fn execute(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Run { mode } => {
            if cfg!(debug_assertions) {
                return Err("measurement requires a release build: cargo run --release -p benchmark -- run ...".into());
            }
            match mode {
                RunCommand::Time { selection, trials } => {
                    run::run(&cli.results_dir, Kind::Time, selection, trials as usize)
                }
                RunCommand::Size { selection } => {
                    run::run(&cli.results_dir, Kind::Size, selection, 1)
                }
            }
        }
        Command::Compare {
            kind,
            name_1,
            name_2,
        } => {
            let old = results::load(&cli.results_dir, kind, &name_1)?;
            let new = results::load(&cli.results_dir, kind, &name_2)?;
            if report::compare(&old, &new) {
                Err("comparison has incompatible, incomplete, or failed cases; see report".into())
            } else {
                Ok(())
            }
        }
        Command::Worker {
            path,
            rule,
            tuning,
            threads,
            pool_mode,
        } => {
            if cfg!(debug_assertions) {
                return Err("worker requires a release build".into());
            }
            let input = data::read(&path)?;
            let settings = run::settings(&rule, threads, &tuning)?;
            let result = run::measure(input, &settings, true, pool_mode)?;
            serde_json::to_writer(std::io::stdout().lock(), &result)?;
            Ok(())
        }
    }
}

fn main() -> std::process::ExitCode {
    match execute(Cli::parse()) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}
