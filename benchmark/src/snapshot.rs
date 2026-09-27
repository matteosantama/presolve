//! Deterministic presolve records, one JSON object per line and instance.
//!
//! Records hold only what is a function of the input and settings when the
//! time limit is not reached: outcome, sizes, reductions by rule and kind,
//! and work counters. Timings are excluded. Snapshots are written from the
//! typed records below but read back as generic JSON, so a snapshot written
//! by an older or newer binary still compares field by field.

use presolve::Outcome;
use presolve::result::{PresolveResult, Size as LibrarySize};
use serde::Serialize;
use serde_json::Value;
use std::collections::BTreeMap;
use std::io::{self, BufRead, Write};
use std::path::Path;

#[derive(Clone, Debug, Serialize)]
pub struct Record {
    /// `family/name`.
    pub instance: String,
    #[serde(flatten)]
    pub run: Run,
}

#[derive(Clone, Debug, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum Run {
    Presolved(Box<Presolved>),
    LoadError { message: String },
    Panic { message: String },
}

#[derive(Clone, Debug, Serialize)]
pub struct Presolved {
    /// `unchanged`, `reduced`, `solved`, `infeasible`, or `unbounded`.
    pub outcome: &'static str,
    /// When true, the rest of the record depends on machine speed.
    pub time_limit_reached: bool,
    pub quadratic_changed: bool,
    pub before: Size,
    /// Absent when presolve stops on a certificate.
    pub after: Option<Size>,
    /// Nonzero counts as `{rule: {kind: count}}`.
    pub reductions: BTreeMap<&'static str, BTreeMap<&'static str, usize>>,
    pub work: Work,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct Size {
    pub variables: usize,
    pub linear_rows: usize,
    pub conic_rows: usize,
    pub cone_blocks: usize,
    pub a_nonzeros: usize,
    pub g_nonzeros: usize,
    pub p_nonzeros: usize,
}
impl Size {
    fn new(size: LibrarySize, cone_blocks: usize) -> Self {
        Self {
            variables: size.variables,
            linear_rows: size.linear_rows,
            conic_rows: size.conic_rows,
            cone_blocks,
            a_nonzeros: size.a_nonzeros,
            g_nonzeros: size.g_nonzeros,
            p_nonzeros: size.p_nonzeros,
        }
    }
}

/// Effort counters. They change when an algorithm does more or less work
/// for the same reductions.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct Work {
    pub fast_phases: usize,
    pub medium_phases: usize,
    pub parallel_comparisons: usize,
    pub equalities: Equalities,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct Equalities {
    pub rows_examined: usize,
    pub structural_rejections: usize,
    pub pivot_rejections: usize,
    pub work_rejections: usize,
    pub attempts: usize,
    pub rejected_updates: usize,
    pub numerical_rejections: usize,
    pub constraint_fill_rejections: usize,
    pub quadratic_fill_rejections: usize,
    pub hessian_growth_rejections: usize,
    pub deadline_rejections: usize,
    pub accepted: usize,
    pub estimated_work: usize,
}

impl Presolved {
    /// Summarize a presolve call whose input had `cone_blocks` cone blocks.
    pub fn new(result: &PresolveResult, cone_blocks: usize) -> Self {
        let stats = &result.stats;
        let (outcome, cones_after) = match &result.outcome {
            Outcome::Unchanged(problem) => ("unchanged", problem.cones.len()),
            Outcome::Reduced(reduced) => ("reduced", reduced.problem.cones.len()),
            Outcome::Solved(_) => ("solved", 0),
            Outcome::Infeasible(_) => ("infeasible", 0),
            Outcome::Unbounded(_) => ("unbounded", 0),
        };
        let mut reductions: BTreeMap<_, BTreeMap<_, _>> = BTreeMap::new();
        for (rule, kind, count) in stats.reductions.iter() {
            reductions
                .entry(rule.name())
                .or_default()
                .insert(kind.name(), count);
        }
        let e = stats.equalities;
        Self {
            outcome,
            time_limit_reached: stats.time_limit_reached,
            quadratic_changed: stats.quadratic_changed,
            before: Size::new(stats.before, cone_blocks),
            after: stats.after.map(|size| Size::new(size, cones_after)),
            reductions,
            work: Work {
                fast_phases: stats.phases.fast,
                medium_phases: stats.phases.medium,
                parallel_comparisons: stats.parallel_comparisons,
                equalities: Equalities {
                    rows_examined: e.rows_examined,
                    structural_rejections: e.structural_rejections,
                    pivot_rejections: e.pivot_rejections,
                    work_rejections: e.work_rejections,
                    attempts: e.attempts,
                    rejected_updates: e.rejected_updates,
                    numerical_rejections: e.numerical_rejections,
                    constraint_fill_rejections: e.constraint_fill_rejections,
                    quadratic_fill_rejections: e.quadratic_fill_rejections,
                    hessian_growth_rejections: e.hessian_growth_rejections,
                    deadline_rejections: e.deadline_rejections,
                    accepted: e.accepted,
                    estimated_work: e.estimated_work,
                },
            },
        }
    }
}

/// Write records sorted by instance, one JSON object per line, creating the
/// parent directory. Duplicate instances are an error.
pub fn write(path: &Path, records: &mut [Record]) -> io::Result<()> {
    records.sort_by(|a, b| a.instance.cmp(&b.instance));
    if let Some(pair) = records.windows(2).find(|p| p[0].instance == p[1].instance) {
        return Err(io::Error::other(format!(
            "instance {} recorded twice",
            pair[0].instance
        )));
    }
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let mut out = io::BufWriter::new(std::fs::File::create(path)?);
    for record in records.iter() {
        serde_json::to_writer(&mut out, record)?;
        out.write_all(b"\n")?;
    }
    out.flush()
}

/// Read a snapshot as generic JSON keyed by instance, with the `instance`
/// field removed from each value.
pub fn read(path: &Path) -> io::Result<BTreeMap<String, Value>> {
    let file = std::fs::File::open(path)?;
    let mut records = BTreeMap::new();
    for (index, line) in io::BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fail = |message: String| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{}:{}: {message}", path.display(), index + 1),
            )
        };
        let mut value: Value = serde_json::from_str(&line).map_err(|e| fail(e.to_string()))?;
        let instance = match value.as_object_mut().and_then(|o| o.remove("instance")) {
            Some(Value::String(instance)) => instance,
            _ => return Err(fail("record without a string instance".into())),
        };
        if records.contains_key(&instance) {
            return Err(fail(format!("instance {instance} recorded twice")));
        }
        records.insert(instance, value);
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use presolve::{Presolver, Problem};

    fn presolved() -> Presolved {
        let text = "\
NAME
ROWS
 N obj
 E r
COLUMNS
 x obj 1 r 1
 y obj 1 r 1
RHS
 rhs r 2
BOUNDS
 FX b x 1
ENDATA
";
        let problem: Problem = crate::mps::parse(text).unwrap();
        let cones = problem.cones.len();
        Presolved::new(&Presolver::default().presolve(problem), cones)
    }

    #[test]
    fn records_round_trip_as_sorted_json_lines() {
        let path =
            std::env::temp_dir().join(format!("benchmark-snapshot-{}.jsonl", std::process::id()));
        let mut records = vec![
            Record {
                instance: "b/second".into(),
                run: Run::LoadError {
                    message: "line 3: unknown row \"q\"".into(),
                },
            },
            Record {
                instance: "a/first".into(),
                run: Run::Presolved(Box::new(presolved())),
            },
        ];
        write(&path, &mut records).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(
            lines[0]
                .starts_with(r#"{"instance":"a/first","status":"presolved","outcome":"solved""#),
            "{}",
            lines[0]
        );
        assert_eq!(
            lines[1],
            r#"{"instance":"b/second","status":"load_error","message":"line 3: unknown row \"q\""}"#
        );
        let read_back = read(&path).unwrap();
        assert_eq!(
            read_back.keys().collect::<Vec<_>>(),
            ["a/first", "b/second"]
        );
        let first = &read_back["a/first"];
        assert_eq!(first["before"]["variables"], 2);
        assert_eq!(first["before"]["cone_blocks"], 0);
        assert_eq!(first["after"]["variables"], 0);
        // Fixing x leaves the singleton equality y = 1, which bounds y on
        // both sides; the fixed-variables rule then removes y.
        assert_eq!(
            first["reductions"],
            serde_json::json!({
                "fixed_variables": {"fixed": 2},
                "singleton_rows": {"deleted_row": 1, "tightened_bound": 2},
            })
        );
        assert_eq!(first["time_limit_reached"], false);
        assert!(first["work"]["fast_phases"].as_u64().unwrap() >= 1);
        assert!(first.get("instance").is_none());
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn duplicate_instances_are_rejected() {
        let path = std::env::temp_dir().join(format!(
            "benchmark-snapshot-dup-{}.jsonl",
            std::process::id()
        ));
        let record = Record {
            instance: "a/x".into(),
            run: Run::Panic {
                message: "boom".into(),
            },
        };
        let err = write(&path, &mut [record.clone(), record]).unwrap_err();
        assert_eq!(err.to_string(), "instance a/x recorded twice");
        std::fs::write(&path, "{\"instance\":\"a/x\"}\n{\"instance\":\"a/x\"}\n").unwrap();
        let err = read(&path).unwrap_err();
        assert!(
            err.to_string().ends_with(":2: instance a/x recorded twice"),
            "{err}"
        );
        std::fs::write(&path, "{\"status\":\"panic\"}\n").unwrap();
        let err = read(&path).unwrap_err();
        assert!(
            err.to_string()
                .ends_with(":1: record without a string instance"),
            "{err}"
        );
        std::fs::remove_file(&path).unwrap();
    }
}
