//! Presolve checked against Clarabel on the original problem.
//!
//! Each instance is presolved; the original is solved with Clarabel; a
//! reduced problem is solved too and its answer mapped back through
//! postsolve. The answer in original coordinates, whether recovered or
//! produced by presolve itself, must satisfy the measures of
//! [`crate::checks`] and agree with the reference status and objective.
//!
//! Clarabel's tolerances, 1e-8 on gap and feasibility by default, hold in
//! its scaled space and leave elementwise relative violations up to about
//! 1e-4 in original coordinates. So the same measures of Clarabel's own
//! answers, to the original (`reference_checks`) and to the reduced problem
//! (`reduced_checks`), calibrate the tolerance: the largest of them is the
//! solver's accuracy on the instance, and every measure of a recovered
//! answer, and its objective, may reach [`FACTOR`] times that accuracy,
//! never less than [`BASE`]. Measures of an answer presolve computes itself
//! use `BASE`, except its objective, which is calibrated on the reference.
//! `certificate_value` and `ray_objective` must stay below `-BASE`.
//!
//! An instance fails when a measure exceeds its tolerance, presolve or
//! postsolve panics, a reduced problem is malformed, or presolve and
//! Clarabel contradict each other (say optimal against infeasible). It is
//! inconclusive when nothing failed but Clarabel stopped without a
//! full-accuracy answer on either problem, returned an infeasibility
//! certificate that fails its own checks, returned an optimum that the
//! checked point beats while feasible, or was so inaccurate that the
//! tolerance exceeds [`LIMIT`].
//!
//! Postsolve can amplify the solver's residuals, for example through a
//! bound relaxed because rows imply it, so an instance that does not pass
//! is verified again with Clarabel at [`TIGHT`], unless it hit the time
//! limit or failed in a way no solve can change; the second verdict stands
//! unless it is inconclusive.

use crate::checks::{Measures, Model};
use crate::corpus::Instance;
use crate::mps;
use crate::reference::{self, Answer, Reference};
use crate::run::{Profile, panic_message};
use presolve::matrix::{CscMatrix, CscMatrixRef};
use presolve::postsolve::Solution;
use presolve::problem::{Constraint, Problem};
use presolve::{Outcome, Presolver};
use rayon::prelude::*;
use serde::Serialize;
use std::fmt::Write as _;
use std::io::{self, Write};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::time::{Duration, Instant};

/// The smallest relative tolerance.
pub const BASE: f64 = 1e-6;
/// How much less accurate than Clarabel's answers a recovered answer may be.
pub const FACTOR: f64 = 10.0;
/// The largest meaningful tolerance.
pub const LIMIT: f64 = 1e-3;
/// Clarabel's gap, feasibility and infeasibility tolerance, first as by
/// default, then for instances that do not pass.
pub const DEFAULT: f64 = 1e-8;
pub const TIGHT: f64 = 1e-10;

/// Measures that must be negative rather than small.
fn signed(check: &str) -> bool {
    matches!(check, "certificate_value" | "ray_objective")
}

/// Measures that calibrate the tolerance but never fail an answer.
fn informational(check: &str) -> bool {
    check == "duality_gap"
}

/// The tolerance calibrated on Clarabel's measures of its own answers.
fn calibrated(measures: &[&Measures]) -> f64 {
    let accuracy = measures
        .iter()
        .flat_map(|m| m.iter())
        .filter(|(check, v)| !signed(check) && v.is_finite())
        .fold(0.0_f64, |a, (_, &v)| a.max(v));
    (FACTOR * accuracy).max(BASE)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Verdict {
    Pass,
    Fail,
    Inconclusive,
}

#[derive(Clone, Debug, Serialize)]
pub struct Failure {
    pub check: &'static str,
    /// Absent for failures without a magnitude.
    pub value: Option<f64>,
    pub tolerance: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

/// One Clarabel solve.
#[derive(Clone, Debug, Serialize)]
pub struct Solve {
    pub status: String,
    pub iterations: u32,
    pub seconds: f64,
    /// Objective of the solver's point, including the constant.
    pub objective: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Record {
    /// `family/name`.
    pub instance: String,
    pub verdict: Verdict,
    /// Presolve's outcome; absent when it did not return one.
    pub outcome: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub failures: Vec<Failure>,
    /// Failures at the default solver tolerance that the tight one cleared
    /// or confirmed.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub default_failures: Vec<Failure>,
    /// Clarabel's gap and feasibility tolerance for this record.
    pub solver_tolerance: f64,
    /// The largest tolerance applied to a measure.
    pub tolerance: f64,
    /// Measures of presolve's answer in original coordinates.
    pub checks: Measures,
    /// The same measures of Clarabel's answer to the original.
    pub reference_checks: Measures,
    /// The same measures of Clarabel's answer to the reduced problem, on it.
    pub reduced_checks: Measures,
    pub original: Option<Solve>,
    pub reduced: Option<Solve>,
    pub presolve_seconds: f64,
}

impl Record {
    fn new(instance: String, solver_tolerance: f64) -> Self {
        Self {
            instance,
            verdict: Verdict::Pass,
            outcome: None,
            reason: None,
            failures: Vec::new(),
            default_failures: Vec::new(),
            solver_tolerance,
            tolerance: 0.0,
            checks: Measures::new(),
            reference_checks: Measures::new(),
            reduced_checks: Measures::new(),
            original: None,
            reduced: None,
            presolve_seconds: 0.0,
        }
    }

    fn measure(&mut self, measures: Measures, tolerance: f64) {
        self.tolerance = self.tolerance.max(tolerance);
        for (check, value) in measures {
            let tolerance = if signed(check) { -BASE } else { tolerance };
            if !informational(check) && (value.is_nan() || value > tolerance) {
                self.failures.push(Failure {
                    check,
                    value: Some(value),
                    tolerance: Some(tolerance),
                    detail: None,
                });
            }
            self.checks.insert(check, value);
        }
    }

    fn fail(&mut self, check: &'static str, detail: String) {
        self.failures.push(Failure {
            check,
            value: None,
            tolerance: None,
            detail: Some(detail),
        });
    }

    fn inconclusive(&mut self, reason: String) {
        match &mut self.reason {
            Some(r) => {
                r.push_str("; ");
                r.push_str(&reason);
            }
            None => self.reason = Some(reason),
        }
    }

    fn finish(mut self) -> Self {
        if self.failures.is_empty() && self.tolerance > LIMIT {
            let reason = format!("Clarabel accuracy sets tolerance {:.0e}", self.tolerance);
            self.inconclusive(reason);
        }
        self.verdict = if !self.failures.is_empty() {
            Verdict::Fail
        } else if self.reason.is_some() {
            Verdict::Inconclusive
        } else {
            Verdict::Pass
        };
        self
    }

    /// Whether a more accurate solve could change the verdict.
    fn retryable(&self) -> bool {
        self.verdict != Verdict::Pass
            && !self
                .reason
                .as_deref()
                .is_some_and(|r| r.contains("MaxTime"))
            && self.failures.iter().all(|f| {
                !matches!(
                    f.check,
                    "presolve_panic"
                        | "postsolve_panic"
                        | "reduced_structure"
                        | "unchanged_identity"
                )
            })
    }
}

/// What an answer says about a problem.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Optimal,
    Infeasible,
    Unbounded,
}
impl Kind {
    fn of(answer: &Answer) -> Option<Self> {
        match answer {
            Answer::Optimal(_) => Some(Self::Optimal),
            Answer::Infeasible(_) => Some(Self::Infeasible),
            Answer::Unbounded(_) => Some(Self::Unbounded),
            Answer::Unknown => None,
        }
    }
}

fn summary(reference: &Reference, model: &Model) -> Solve {
    Solve {
        status: reference.status.clone(),
        iterations: reference.iterations,
        seconds: reference.seconds,
        objective: match &reference.result {
            Answer::Optimal(s) => Some(model.objective(&s.x)),
            _ => None,
        },
    }
}

/// Clarabel's own measures of an answer to `model`'s problem, and whether
/// the answer can be used: an infeasibility certificate must pass its own
/// checks within `FACTOR * BASE`, since Clarabel can stop on a direction
/// that only looks like one in its scaled space.
fn own_checks(model: &Model, answer: &Answer) -> (Measures, bool) {
    let certificate = |m: Measures| {
        let trusted = m.iter().all(|(check, &v)| {
            if signed(check) {
                v <= -BASE
            } else {
                v <= FACTOR * BASE
            }
        });
        (m, trusted)
    };
    match answer {
        Answer::Optimal(s) => (model.optimality(s), true),
        Answer::Infeasible(c) => certificate(model.farkas(c)),
        Answer::Unbounded(d) => certificate(model.ray(d)),
        Answer::Unknown => (Measures::new(), false),
    }
}

fn unusable(problem: &str, solve: &Reference) -> String {
    match solve.result {
        Answer::Unknown => format!("{problem}: {}", solve.status),
        _ => format!("{problem}: {} certificate fails its checks", solve.status),
    }
}

fn relative(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1.0)
}

/// Presolve `problem` and check the outcome against Clarabel, whose solves
/// stop after `time_limit` seconds. The record's instance name is empty.
pub fn verify(problem: Problem, presolver: &Presolver, time_limit: f64) -> Record {
    verify_with(problem, presolver, time_limit, &|_| {})
}

/// As `verify`, with `tamper` applied to the solution in original
/// coordinates before it is checked.
fn verify_with(
    problem: Problem,
    presolver: &Presolver,
    time_limit: f64,
    tamper: &dyn Fn(&mut Solution),
) -> Record {
    let first = verify_at(problem.clone(), presolver, DEFAULT, time_limit, tamper);
    if !first.retryable() {
        return first;
    }
    let mut second = verify_at(problem, presolver, TIGHT, time_limit, tamper);
    if second.verdict == Verdict::Inconclusive {
        return first;
    }
    second.default_failures = first.failures;
    second
}

fn verify_at(
    problem: Problem,
    presolver: &Presolver,
    accuracy: f64,
    time_limit: f64,
    tamper: &dyn Fn(&mut Solution),
) -> Record {
    let mut r = Record::new(String::new(), accuracy);
    let original = problem.clone();
    let model = Model::new(&original);
    let start = Instant::now();
    let result = catch_unwind(AssertUnwindSafe(|| presolver.presolve(problem)));
    r.presolve_seconds = start.elapsed().as_secs_f64();
    let outcome = match result {
        Ok(result) => result.outcome,
        Err(payload) => {
            r.fail("presolve_panic", panic_message(payload.as_ref()));
            return r.finish();
        }
    };
    r.outcome = Some(match &outcome {
        Outcome::Unchanged(_) => "unchanged",
        Outcome::Reduced(_) => "reduced",
        Outcome::Solved(_) => "solved",
        Outcome::Infeasible(_) => "infeasible",
        Outcome::Unbounded(_) => "unbounded",
    });
    if let Outcome::Unchanged(returned) = &outcome {
        if !identical(returned, &original) {
            r.fail(
                "unchanged_identity",
                "returned problem differs from the input".into(),
            );
        }
        return r.finish();
    }

    let (reference_kind, reference_objective) =
        match reference::solve(&original, accuracy, time_limit) {
            Err(e) => {
                r.inconclusive(format!("original: {e}"));
                (None, None)
            }
            Ok(reference) => {
                r.original = Some(summary(&reference, &model));
                let (checks, usable) = own_checks(&model, &reference.result);
                r.reference_checks = checks;
                if !usable {
                    r.inconclusive(unusable("original", &reference));
                }
                let objective = r.original.as_ref().and_then(|o| o.objective);
                (Kind::of(&reference.result).filter(|_| usable), objective)
            }
        };
    // Compare what presolve's side concluded with the reference.
    let compare = |r: &mut Record, kind: Kind, source: &str| match reference_kind {
        None => {}
        Some(reference) if reference == kind => {}
        Some(reference @ (Kind::Infeasible | Kind::Unbounded))
            if kind != Kind::Optimal && reference != Kind::Optimal =>
        {
            r.inconclusive(format!("{source} {kind:?}, original {reference:?}"));
        }
        Some(reference) => r.fail(
            "status",
            format!("{source} {kind:?}, original {reference:?}"),
        ),
    };
    // `tolerance` applies to the answer, `objective` to its objective.
    // A feasible point better than Clarabel's shows Clarabel's is not optimal.
    let check_solution = |r: &mut Record, mut solution: Solution, tolerance, objective| {
        tamper(&mut solution);
        r.measure(model.optimality(&solution), tolerance);
        if let Some(f) = reference_objective {
            let value = model.objective(&solution.x);
            let error = relative(value, f);
            let feasible = ["primal_rows", "primal_bounds", "primal_cones"]
                .iter()
                .all(|k| r.checks.get(k).is_none_or(|&v| v <= tolerance));
            if error > objective && value < f && feasible {
                r.checks.insert("objective", error);
                r.inconclusive(format!(
                    "checked point is feasible and better than Clarabel's by {error:.1e}"
                ));
            } else {
                r.measure(Measures::from([("objective", error)]), objective);
            }
        }
        solution
    };

    match outcome {
        Outcome::Unchanged(_) => unreachable!(),
        Outcome::Solved(solution) => {
            let objective = calibrated(&[&r.reference_checks]);
            check_solution(&mut r, solution, BASE, objective);
            compare(&mut r, Kind::Optimal, "presolve");
        }
        Outcome::Infeasible(certificate) => {
            r.measure(model.farkas(&certificate), BASE);
            compare(&mut r, Kind::Infeasible, "presolve");
        }
        Outcome::Unbounded(certificate) => {
            r.measure(model.primal(&certificate.point), BASE);
            r.measure(model.ray(&certificate.ray), BASE);
            compare(&mut r, Kind::Unbounded, "presolve");
        }
        Outcome::Reduced(reduced) => {
            if let Err(e) = well_formed(&reduced.problem) {
                r.fail("reduced_structure", e);
                return r.finish();
            }
            let reduced_model = Model::new(&reduced.problem);
            let solved = match reference::solve(&reduced.problem, accuracy, time_limit) {
                Ok(solved) => solved,
                Err(e) => {
                    r.inconclusive(format!("reduced: {e}"));
                    return r.finish();
                }
            };
            r.reduced = Some(summary(&solved, &reduced_model));
            let (checks, usable) = own_checks(&reduced_model, &solved.result);
            r.reduced_checks = checks;
            if !usable {
                r.inconclusive(unusable("reduced", &solved));
                return r.finish();
            }
            let tolerance = calibrated(&[&r.reference_checks, &r.reduced_checks]);
            let postsolve = &reduced.postsolve;
            let recovered = catch_unwind(AssertUnwindSafe(|| match &solved.result {
                Answer::Optimal(s) => Some((
                    Kind::Optimal,
                    Recovered::Solution(postsolve.recover_solution(s.as_ref())),
                )),
                Answer::Infeasible(c) => Some((
                    Kind::Infeasible,
                    Recovered::Certificate(postsolve.recover_primal_certificate(c.as_ref())),
                )),
                Answer::Unbounded(d) => Some((
                    Kind::Unbounded,
                    Recovered::Ray(postsolve.recover_primal_ray(d)),
                )),
                Answer::Unknown => None,
            }));
            match recovered {
                Err(payload) => r.fail("postsolve_panic", panic_message(payload.as_ref())),
                Ok(None) => unreachable!(),
                Ok(Some((kind, recovered))) => {
                    match recovered {
                        Recovered::Solution(solution) => {
                            let solution = check_solution(&mut r, solution, tolerance, tolerance);
                            if let Some(f) = r.reduced.as_ref().and_then(|s| s.objective) {
                                let value = relative(model.objective(&solution.x), f);
                                let measure = Measures::from([("objective_recovery", value)]);
                                r.measure(measure, tolerance);
                            }
                        }
                        Recovered::Certificate(c) => r.measure(model.farkas(&c), tolerance),
                        Recovered::Ray(d) => r.measure(model.ray(&d), tolerance),
                    }
                    compare(&mut r, kind, "reduced");
                }
            }
        }
    }
    r.finish()
}

enum Recovered {
    Solution(Solution),
    Certificate(presolve::postsolve::PrimalCertificate),
    Ray(Vec<f64>),
}

fn identical(a: &Problem, b: &Problem) -> bool {
    a.p == b.p
        && a.c == b.c
        && a.c0 == b.c0
        && a.a == b.a
        && a.rows == b.rows
        && a.variable_bounds == b.variable_bounds
        && a.cones == b.cones
}

fn valid_csc(m: &CscMatrix) -> Result<(), String> {
    CscMatrixRef::new(
        m.rows(),
        m.columns(),
        m.column_pointers(),
        m.row_indices(),
        m.values(),
    )
    .map(|_| ())
    .map_err(|e| e.to_string())
}

/// Dimensions, sparse structure, finite data, an upper-triangular `P`, and
/// cone blocks that match their rows.
fn well_formed(p: &Problem) -> Result<(), String> {
    let n = p.c.len();
    if p.a.columns() != n || p.a.rows() != p.rows.len() {
        return Err(format!(
            "A is {}x{} for {} rows and {n} variables",
            p.a.rows(),
            p.a.columns(),
            p.rows.len()
        ));
    }
    if !p.variable_bounds.is_empty() && p.variable_bounds.len() != n {
        return Err(format!(
            "{} variable bounds for {n} variables",
            p.variable_bounds.len()
        ));
    }
    valid_csc(&p.a).map_err(|e| format!("A: {e}"))?;
    if let Some(q) = &p.p {
        if q.rows() != n || q.columns() != n {
            return Err(format!(
                "P is {}x{} for {n} variables",
                q.rows(),
                q.columns()
            ));
        }
        valid_csc(q).map_err(|e| format!("P: {e}"))?;
        if (0..n).any(|j| q.as_ref().column(j).any(|(i, _)| i > j)) {
            return Err("P has entries below the diagonal".into());
        }
    }
    if !p.c.iter().all(|v| v.is_finite()) || !p.c0.is_finite() {
        return Err("objective has non-finite values".into());
    }
    let bad = |b: presolve::problem::Bounds| {
        b.lower.is_nan()
            || b.upper.is_nan()
            || b.lower == f64::INFINITY
            || b.upper == f64::NEG_INFINITY
    };
    if p.variable_bounds.iter().any(|&b| bad(b)) {
        return Err("variable bounds are NaN or empty at infinity".into());
    }
    let mut counts = vec![0; p.cones.len()];
    for row in &p.rows {
        match *row {
            Constraint::Linear(b) if bad(b) => {
                return Err("row bounds are NaN or empty at infinity".into());
            }
            Constraint::Linear(_) => {}
            Constraint::Cone { rhs, block } => {
                if !rhs.is_finite() || block >= counts.len() {
                    return Err(format!("cone row with rhs {rhs} in block {block}"));
                }
                counts[block] += 1;
            }
        }
    }
    if let Some(k) = (0..counts.len()).find(|&k| counts[k] != p.cones[k].dimension()) {
        return Err(format!(
            "cone block {k} has {} rows for {:?}",
            counts[k], p.cones[k]
        ));
    }
    Ok(())
}

/// Verify `instances` on `jobs` worker threads, zero meaning one per core.
pub fn run(
    instances: &[Instance],
    profile: Profile,
    jobs: usize,
    time_limit: f64,
) -> Result<Vec<Record>, String> {
    let presolver = Presolver::new(profile.settings()).map_err(|e| e.to_string())?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .map_err(|e| e.to_string())?;
    Ok(pool.install(|| {
        instances
            .par_iter()
            .map(|instance| {
                let mut record = match mps::read(&instance.path) {
                    Ok(problem) => verify(problem, &presolver, time_limit),
                    Err(e) => {
                        let mut r = Record::new(String::new(), DEFAULT);
                        r.inconclusive(format!("load error: {e}"));
                        r.finish()
                    }
                };
                record.instance = instance.id();
                record
            })
            .collect()
    }))
}

/// One JSON record per line, sorted by instance.
pub fn write(path: &Path, records: &mut [Record]) -> io::Result<()> {
    records.sort_by(|a, b| a.instance.cmp(&b.instance));
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let mut out = io::BufWriter::new(std::fs::File::create(path)?);
    for record in records.iter() {
        serde_json::to_writer(&mut out, record)?;
        out.write_all(b"\n")?;
    }
    out.flush()
}

/// Counts by verdict, then every failure and inconclusive instance.
pub fn summarize(records: &[Record], profile: Profile, wall: Duration) -> String {
    let count = |v: Verdict| records.iter().filter(|r| r.verdict == v).count();
    let mut s = format!(
        "{} instances under the {profile:?} profile in {:.1} s wall: {} pass, {} fail, {} inconclusive.\n",
        records.len(),
        wall.as_secs_f64(),
        count(Verdict::Pass),
        count(Verdict::Fail),
        count(Verdict::Inconclusive),
    );
    let cleared: Vec<&str> = records
        .iter()
        .filter(|r| r.verdict == Verdict::Pass && r.solver_tolerance == TIGHT)
        .map(|r| r.instance.as_str())
        .collect();
    if !cleared.is_empty() {
        let _ = writeln!(
            s,
            "  {} pass only with Clarabel at {TIGHT:.0e}: {}",
            cleared.len(),
            cleared.join(", ")
        );
    }
    for r in records.iter().filter(|r| r.verdict == Verdict::Fail) {
        let failures: Vec<String> = r
            .failures
            .iter()
            .map(|f| match (f.value, &f.detail) {
                (Some(v), _) => format!("{} {v:.2e} > {:.0e}", f.check, f.tolerance.unwrap_or(0.0)),
                (None, Some(d)) => format!("{}: {d}", f.check),
                (None, None) => f.check.to_string(),
            })
            .collect();
        let outcome = r.outcome.unwrap_or("-");
        let _ = writeln!(
            s,
            "  FAIL {} ({outcome}, Clarabel at {:.0e}): {}",
            r.instance,
            r.solver_tolerance,
            failures.join(", ")
        );
    }
    for r in records
        .iter()
        .filter(|r| r.verdict == Verdict::Inconclusive)
    {
        let outcome = r.outcome.unwrap_or("-");
        let reason = r.reason.as_deref().unwrap_or("");
        let _ = writeln!(s, "  inconclusive {} ({outcome}): {reason}", r.instance);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use presolve::problem::Bounds;

    fn bounds(lower: f64, upper: f64) -> Bounds {
        Bounds { lower, upper }
    }

    /// A problem from dense rows and an optional dense upper triangle of P.
    fn problem(
        c: &[f64],
        p: Option<&[&[f64]]>,
        rows: &[(&[f64], Bounds)],
        variables: &[Bounds],
    ) -> Problem {
        let n = c.len();
        let triplets = |dense: &[&[f64]]| {
            let (mut i, mut j, mut v) = (Vec::new(), Vec::new(), Vec::new());
            for (r, row) in dense.iter().enumerate() {
                for (k, &value) in row.iter().enumerate() {
                    if value != 0.0 {
                        i.push(r);
                        j.push(k);
                        v.push(value);
                    }
                }
            }
            CscMatrix::from_triplets(dense.len(), n, i, j, v).unwrap()
        };
        let dense: Vec<&[f64]> = rows.iter().map(|(row, _)| *row).collect();
        Problem {
            p: p.map(triplets),
            c: c.to_vec(),
            c0: 1.5,
            a: triplets(&dense),
            rows: rows.iter().map(|&(_, b)| Constraint::Linear(b)).collect(),
            variable_bounds: variables.to_vec(),
            cones: Vec::new(),
        }
    }

    fn presolver() -> Presolver {
        Presolver::new(Profile::Default.settings()).unwrap()
    }

    /// Fixed x4 and x5, with rows long enough to leave a reduced LP.
    fn lp() -> Problem {
        let inf = f64::INFINITY;
        problem(
            &[1.0, -2.0, 1.0, 3.0, 1.0, -1.0],
            None,
            &[
                (&[1.0, 1.0, 1.0, 1.0, 1.0, 0.0], bounds(2.0, 8.0)),
                (&[1.0, -1.0, 2.0, 0.0, 0.0, 1.0], bounds(-1.0, 3.0)),
                (&[0.0, 1.0, 1.0, -1.0, 2.0, 1.0], bounds(-inf, 4.0)),
                (&[2.0, 1.0, 0.0, 1.0, 0.0, 1.0], bounds(1.0, inf)),
            ],
            &[
                bounds(0.0, 5.0),
                bounds(0.0, 5.0),
                bounds(-1.0, 5.0),
                bounds(0.0, inf),
                bounds(0.5, 0.5),
                bounds(1.0, 1.0),
            ],
        )
    }

    #[test]
    fn lp_with_fixed_variables_passes() {
        let record = verify(lp(), &presolver(), 10.0);
        assert_eq!(record.outcome, Some("reduced"), "{record:?}");
        assert_eq!(record.verdict, Verdict::Pass, "{record:?}");
        assert!(record.checks.contains_key("stationarity"));
        assert!(record.checks["objective"] <= BASE);
    }

    #[test]
    fn qp_passes() {
        let inf = f64::INFINITY;
        let qp = problem(
            &[-1.0, -2.0, 0.5, 1.0],
            Some(&[
                &[2.0, 0.5, 0.0, 0.0],
                &[0.0, 1.0, 0.0, 0.0],
                &[0.0, 0.0, 3.0, 1.0],
                &[0.0, 0.0, 0.0, 0.0],
            ]),
            &[
                (&[1.0, 1.0, 1.0, 1.0], bounds(1.0, 1.0)),
                (&[1.0, -1.0, 0.0, 2.0], bounds(-inf, 0.5)),
            ],
            &[
                bounds(0.0, inf),
                bounds(0.0, inf),
                bounds(-1.0, 1.0),
                bounds(0.25, 0.25),
            ],
        );
        let record = verify(qp, &presolver(), 10.0);
        assert_ne!(record.outcome, Some("unchanged"), "{record:?}");
        assert_eq!(record.verdict, Verdict::Pass, "{record:?}");
        assert!(record.checks.contains_key("duality_gap"));
    }

    #[test]
    fn corrupted_recovered_solutions_fail() {
        let shifted = verify_with(lp(), &presolver(), 10.0, &|s| s.x[0] += 0.1);
        assert_eq!(shifted.verdict, Verdict::Fail, "{shifted:?}");
        assert!(!shifted.default_failures.is_empty());
        let flipped = verify_with(lp(), &presolver(), 10.0, &|s| s.y[0] = -s.y[0] - 1.0);
        assert_eq!(flipped.verdict, Verdict::Fail, "{flipped:?}");
        let failed = |check| flipped.failures.iter().any(|f| f.check == check);
        assert!(failed("stationarity"), "{flipped:?}");
    }

    #[test]
    fn infeasible_and_unbounded_problems_pass() {
        let inf = f64::INFINITY;
        let infeasible = problem(
            &[1.0, 1.0],
            None,
            &[(&[1.0, 1.0], bounds(3.0, inf))],
            &[bounds(0.0, 1.0), bounds(0.0, 1.0)],
        );
        let record = verify(infeasible, &presolver(), 10.0);
        assert_eq!(record.verdict, Verdict::Pass, "{record:?}");
        assert!(record.checks["certificate_value"] < 0.0);

        let unbounded = problem(
            &[-1.0, -1.0],
            None,
            &[(&[1.0, -1.0], bounds(-inf, 1.0))],
            &[bounds(0.0, inf), bounds(0.0, inf)],
        );
        let record = verify(unbounded, &presolver(), 10.0);
        assert_eq!(record.verdict, Verdict::Pass, "{record:?}");
        assert!(
            record.checks.get("ray_objective").is_some_and(|&v| v < 0.0),
            "{record:?}"
        );
    }

    #[test]
    fn malformed_reduced_problems_are_detected() {
        let mut p = lp();
        assert!(well_formed(&p).is_ok());
        p.p = Some(CscMatrix::from_triplets(6, 6, vec![1], vec![0], vec![1.0]).unwrap());
        assert!(well_formed(&p).unwrap_err().contains("below the diagonal"));
    }
}
