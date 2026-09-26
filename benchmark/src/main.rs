//! Deterministic presolve snapshots of the benchmark corpus, and their
//! comparison.
//!
//! `snapshot` presolves every instance under a profile and writes one JSON
//! record per line. `compare` diffs two snapshots and exits with status 1
//! when statuses, outcomes, sizes, or reductions differ, 0 when only work
//! counters or nothing differ, and 2 on error.

use benchmark::compare::{self, Comparison};
use benchmark::corpus;
use benchmark::run::{self, Profile};
use benchmark::snapshot::{self, Run};
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

#[derive(Parser)]
#[command(about = "Presolve reduction snapshots of the benchmark corpus")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Presolve the corpus and write a snapshot.
    Snapshot {
        #[arg(long, value_enum, default_value = "default")]
        profile: Profile,
        /// Directory whose subdirectories are instance families.
        #[arg(long, default_value = concat!(env!("CARGO_MANIFEST_DIR"), "/data"))]
        data: PathBuf,
        /// Restrict to a family; repeat for several. Default: all.
        #[arg(long = "family")]
        families: Vec<String>,
        /// Skip a family's files above this many MiB on disk, as FAMILY=MIB.
        /// Repeat for several families. Giving any limit replaces the default.
        #[arg(long = "max-size-mib", value_name = "FAMILY=MIB", default_values = ["miplib=4"])]
        limits: Vec<corpus::SizeLimit>,
        /// Ignore every size limit and run the whole corpus.
        #[arg(long, conflicts_with = "limits")]
        all_sizes: bool,
        /// Instances presolved concurrently; 0 means one per core.
        #[arg(long, default_value_t = 0)]
        jobs: usize,
        /// Instances to list in the timing summary.
        #[arg(long, default_value_t = 10)]
        slowest: usize,
        #[arg(long)]
        out: PathBuf,
    },
    /// Compare two snapshots.
    Compare {
        base: PathBuf,
        head: PathBuf,
        #[arg(long, value_enum, default_value = "text")]
        format: compare::Format,
        /// Rows shown in the per-instance table.
        #[arg(long, default_value_t = 100)]
        max_rows: usize,
    },
}

fn main() -> ExitCode {
    match Cli::parse().command {
        Command::Snapshot {
            profile,
            data,
            families,
            limits,
            all_sizes,
            jobs,
            slowest,
            out,
        } => match take_snapshot(
            profile,
            &data,
            &families,
            if all_sizes { &[] } else { &limits },
            jobs,
            slowest,
            &out,
        ) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::from(2)
            }
        },
        Command::Compare {
            base,
            head,
            format,
            max_rows,
        } => {
            let read = |path: &PathBuf| {
                snapshot::read(path).map_err(|e| format!("{}: {e}", path.display()))
            };
            match read(&base).and_then(|b| Ok((b, read(&head)?))) {
                Ok((base, head)) => {
                    let comparison = Comparison::new(&base, &head);
                    print!("{}", comparison.render(format, max_rows));
                    if comparison.has_changes() {
                        ExitCode::FAILURE
                    } else {
                        ExitCode::SUCCESS
                    }
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    ExitCode::from(2)
                }
            }
        }
    }
}

fn take_snapshot(
    profile: Profile,
    data: &std::path::Path,
    families: &[String],
    limits: &[corpus::SizeLimit],
    jobs: usize,
    slowest: usize,
    out: &std::path::Path,
) -> Result<(), String> {
    let instances =
        corpus::discover(data, families).map_err(|e| format!("{}: {e}", data.display()))?;
    // With --family, a limit for a family not being run does not apply.
    // Otherwise every limit must name a family, so a typo is an error.
    let limits: Vec<_> = limits
        .iter()
        .filter(|l| families.is_empty() || families.contains(&l.family))
        .cloned()
        .collect();
    let (instances, excluded) = corpus::within_limits(instances, &limits)?;
    if !excluded.is_empty() {
        let limits: Vec<String> = limits
            .iter()
            .map(|l| format!("{}={} MiB", l.family, l.mib))
            .collect();
        eprintln!(
            "Skipping {} instances above the size limits {}; pass --all-sizes to include them.",
            excluded.len(),
            limits.join(", ")
        );
    }
    if instances.is_empty() {
        return Err(format!("no instances under {}", data.display()));
    }
    let start = Instant::now();
    let outcomes = run::run(&instances, profile, jobs)?;
    let wall = start.elapsed();
    let mut timed: Vec<_> = outcomes
        .iter()
        .map(|o| (o.elapsed, o.record.instance.as_str()))
        .collect();
    let presolve_total: f64 = timed.iter().map(|(t, _)| t.as_secs_f64()).sum();
    timed.sort_by(|a, b| b.cmp(a));
    let failures: Vec<String> = outcomes
        .iter()
        .filter_map(|o| match &o.record.run {
            Run::Presolved(_) => None,
            Run::LoadError { message } => {
                Some(format!("  {}: load error: {message}", o.record.instance))
            }
            Run::Panic { message } => Some(format!("  {}: panic: {message}", o.record.instance)),
        })
        .collect();
    eprintln!(
        "{} instances under the {profile:?} profile in {:.1} s wall; presolve {:.1} s summed over instances.",
        outcomes.len(),
        wall.as_secs_f64(),
        presolve_total
    );
    if slowest > 0 {
        eprintln!("Slowest presolve calls:");
        for (t, instance) in timed.iter().take(slowest) {
            eprintln!("  {:>8.3} s  {instance}", t.as_secs_f64());
        }
    }
    if !failures.is_empty() {
        eprintln!(
            "{} instances failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
    let mut records: Vec<_> = outcomes.into_iter().map(|o| o.record).collect();
    snapshot::write(out, &mut records).map_err(|e| format!("{}: {e}", out.display()))
}
