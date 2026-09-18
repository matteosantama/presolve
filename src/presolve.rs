//! Presolve entry point, working-model preparation, and result packing.
use crate::{
    core::{
        execution::Executor,
        model::{Model, RowDomain},
    },
    matrix::quadratic::Quadratic,
    postsolve::{
        Coordinates, OriginalMap, Postsolve, PrimalCertificate, SolutionRef, tape::Recovery,
    },
    problem::{Bounds, Cone, Constraint, Matrix, Problem, row_indices},
    result::{Outcome, PresolveResult, ReducedProblem, Size, Stats, UnboundednessCertificate},
    settings::Settings,
};
use std::sync::Arc;
use std::time::{Duration, Instant};

/// Failure to create the requested execution resources.
#[derive(Debug, thiserror::Error)]
#[error("failed to initialize presolver threads: {0}")]
pub struct InitError(#[from] rayon::ThreadPoolBuildError);

/// Reusable presolve configuration and its dedicated execution resources.
/// Each call owns a separate working model and recovery tape. Results do not
/// borrow this object, and concurrent calls may share it through `&self`.
#[derive(Debug)]
pub struct Presolver {
    settings: Settings,
    executor: Executor,
}

/// Default settings run serially, so construction cannot fail.
impl Default for Presolver {
    fn default() -> Self {
        Self {
            settings: Settings::default(),
            executor: Executor::Serial,
        }
    }
}

impl Presolver {
    /// Configure execution once. A thread count of 1 creates no worker threads;
    /// 0 lets Rayon choose automatically. Pool creation errors are returned.
    pub fn new(settings: Settings) -> Result<Self, InitError> {
        let executor = Executor::new(settings.threads)?;
        Ok(Self { settings, executor })
    }

    /// Immutable settings used for every call.
    pub fn settings(&self) -> &Settings {
        &self.settings
    }

    /// Effective execution thread count, including automatic selection.
    pub fn threads(&self) -> usize {
        self.executor.threads()
    }

    /// Consume a problem, reusing this presolver's threads. The time budget and
    /// elapsed statistic start fresh for each call and exclude `Self::new`.
    /// No constraint output matrices are built until explicitly exported.
    /// Input validity is the caller's responsibility.
    pub fn presolve(&self, problem: Problem) -> PresolveResult {
        presolve_owned(problem, &self.settings, &self.executor, Instant::now())
    }
}

fn presolve_owned(
    mut problem: Problem,
    settings: &Settings,
    executor: &Executor,
    start: Instant,
) -> PresolveResult {
    let n = problem.c.len();
    let free_bounds = problem.variable_bounds.is_empty();
    let mut model = working_model(&mut problem, settings, start);
    let before = size(&model);
    let mut stats = Stats {
        equalities: Default::default(),
        elapsed: Duration::ZERO,
        time_limit_reached: false,
        before,
        after: Some(before),
        parallel_comparisons: 0,
        quadratic_changed: false,
    };
    let phases = model.run(
        settings.time_limit.saturating_sub(start.elapsed()),
        executor,
    );
    stats.equalities = model.equality_stats;
    stats.quadratic_changed = model.objective.p.revision != 0;
    let outcome = match phases {
        Err(certificate) => {
            // Domains can have been deleted; use the caller's original row tags.
            let original = OriginalMap {
                columns: n,
                linear: &problem.linear_rows,
                conic: &problem.conic_rows,
            };
            match certificate.mode {
                Recovery::PrimalInfeasibility => {
                    stats.after = None;
                    let (y, z, conic_dual) = original.gather_dual(&certificate.point);
                    Outcome::Infeasible(PrimalCertificate { y, z, conic_dual })
                }
                Recovery::DualInfeasibility => match feasible_point(&model) {
                    Some(mut point) => {
                        stats.after = None;
                        point.truncate(n);
                        let ray = original.gather(certificate.point, Vec::new()).x;
                        Outcome::Unbounded(UnboundednessCertificate { point, ray })
                    }
                    None => finish(model, problem, free_bounds, before, &mut stats),
                },
                Recovery::Solution => unreachable!(),
            }
        }
        Ok(phases) => {
            stats.parallel_comparisons = phases.parallel_comparisons;
            stats.time_limit_reached |= phases.time_limit;
            finish(model, problem, free_bounds, before, &mut stats)
        }
    };
    stats.elapsed = start.elapsed();
    stats.time_limit_reached |= stats.elapsed >= settings.time_limit;
    PresolveResult { outcome, stats }
}

/// Build the working model, taking the input's editable storage and objective.
fn working_model(problem: &mut Problem, settings: &Settings, start: Instant) -> Model {
    let n = problem.c.len();
    let rows = problem
        .rows
        .iter()
        .map(|&row| match row {
            Constraint::Linear(b) => RowDomain::Linear(b),
            Constraint::Cone { rhs, block } => RowDomain::Cone { rhs, block },
        })
        .collect();
    let bounds = if problem.variable_bounds.is_empty() {
        vec![Bounds::FREE; n]
    } else {
        std::mem::take(&mut problem.variable_bounds)
    };
    let mut model = Model::from_parts(
        problem.working_matrix(),
        problem.take_objective(),
        rows,
        bounds,
    );
    model.objective.constant = problem.objective_constant;
    model.configure(settings, start.checked_add(settings.time_limit));
    model.set_cones(problem.cones.clone());
    model
}

/// Hand the input back untouched, or pack the reduced model.
fn finish(
    model: Model,
    problem: Problem,
    free_bounds: bool,
    before: Size,
    stats: &mut Stats,
) -> Outcome {
    // With no edits this equals `before`, so both paths report the same size.
    stats.after = Some(size(&model));
    if model.revision == 0 {
        unchanged(model, problem, free_bounds)
    } else {
        pack(model, problem, before)
    }
}

/// Return the input allocations: storage lent to the model comes back, and
/// empty input bounds stay empty.
fn unchanged(model: Model, mut problem: Problem, free_bounds: bool) -> Outcome {
    problem.restore_working_matrix(model.a);
    problem.c = model.objective.c;
    if !free_bounds {
        problem.variable_bounds = model.bounds;
    }
    Outcome::Unchanged(problem)
}

/// Verify a cheap candidate only when a recession ray has been discovered.
/// Failure to find a witness simply returns a reduced or unchanged problem.
fn feasible_point(model: &Model) -> Option<Vec<f64>> {
    let mut p = model.point();
    for (j, b) in model.bounds.iter().enumerate() {
        if model.alive[j] {
            p.x[j] = 0.0_f64.max(b.lower).min(b.upper);
        }
    }
    for (i, row) in model.rows.iter().enumerate() {
        if let RowDomain::Linear(b) = row {
            let value = model.a.row(i).iter().map(|(j, a)| a * p.x[j]).sum::<f64>();
            if !b.contains(value) {
                return None;
            }
        }
    }
    for (block, rows) in model.cone_rows.iter().enumerate() {
        let rhs: Vec<_> = rows
            .iter()
            .filter_map(|&i| match model.rows[i] {
                RowDomain::Cone { rhs, block: b } if b == block => {
                    Some(rhs - model.a.row(i).iter().map(|(j, a)| a * p.x[j]).sum::<f64>())
                }
                _ => None,
            })
            .collect();
        if !rhs.is_empty()
            && !matches!(
                model.cones[block].classify(&rhs, 0.),
                crate::problem::Membership::Inside
            )
        {
            return None;
        }
    }
    model.postsolve.recover(&mut p, Recovery::Solution);
    // Auxiliary columns are appended; the caller truncates to input dimensions.
    Some(p.x)
}

fn size(model: &Model) -> Size {
    let mut size = Size {
        variables: model.alive.iter().filter(|&&v| v).count(),
        p_nonzeros: model.objective.p.nnz(),
        ..Size::default()
    };
    for (i, row) in model.rows.iter().enumerate() {
        match row {
            RowDomain::Linear(_) => {
                size.linear_rows += 1;
                size.a_nonzeros += model.a.row(i).len();
            }
            RowDomain::Cone { .. } => {
                size.conic_rows += 1;
                size.g_nonzeros += model.a.row(i).len();
            }
            RowDomain::Deleted => (),
        }
    }
    size
}
/// Surviving columns and rows of the working model, in stable indices, with
/// the cone blocks that keep at least one coordinate.
struct Survivors {
    columns: Vec<usize>,
    linear_rows: Vec<usize>,
    conic_rows: Vec<usize>,
    cones: Vec<Cone>,
}

fn survivors(model: &Model) -> Survivors {
    let columns = (0..model.alive.len()).filter(|&j| model.alive[j]).collect();
    let linear_rows = (0..model.rows.len())
        .filter(|&i| matches!(model.rows[i], RowDomain::Linear(_)))
        .collect();
    let mut conic_rows = Vec::new();
    let mut cones = Vec::new();
    for (block, block_rows) in model.cone_rows.iter().enumerate() {
        let surviving: Vec<_> = block_rows
            .iter()
            .copied()
            .filter(|&i| matches!(model.rows[i], RowDomain::Cone { block: b, .. } if b == block))
            .collect();
        if !surviving.is_empty() {
            cones.push(model.cones[block]);
            conic_rows.extend(surviving);
        }
    }
    Survivors {
        columns,
        linear_rows,
        conic_rows,
        cones,
    }
}

/// Inverse of a compact-to-stable map; unmapped stable indices hold `usize::MAX`.
fn inverse(map: &[usize], len: usize) -> Vec<usize> {
    let mut out = vec![usize::MAX; len];
    for (at, &j) in map.iter().enumerate() {
        out[j] = at;
    }
    out
}

/// The recovery map takes the survivor index vectors and the model's tape.
fn build_postsolve(
    model: &mut Model,
    original: &mut Problem,
    survivors: Survivors,
    before: Size,
) -> (Postsolve, Vec<Cone>) {
    let input_linear = std::mem::take(&mut original.linear_rows);
    let input_conic = std::mem::take(&mut original.conic_rows);
    let coordinates = Coordinates {
        compact_to_stable_columns: Arc::new(survivors.columns),
        compact_to_stable_linear_rows: survivors.linear_rows,
        compact_to_stable_conic_rows: survivors.conic_rows,
    };
    // Surviving cone coordinates keep their input position for slack copies.
    let original_positions = if input_conic.is_empty() {
        Vec::new()
    } else {
        inverse(&input_conic, model.rows.len())
    };
    let direct_slacks = coordinates
        .compact_to_stable_conic_rows
        .iter()
        .enumerate()
        .map(|(at, &i)| (original_positions[i], at))
        .collect();
    let postsolve = Postsolve {
        total_columns: model.bounds.len(),
        total_rows: model.rows.len(),
        original_columns: before.variables,
        tape: std::mem::take(&mut model.postsolve),
        direct_slacks,
        coordinates,
        input_linear,
        input_conic,
    };
    (postsolve, survivors.cones)
}

/// Compact the working model's objective, bounds, and domains into a problem
/// whose constraint storage is the model's, addressed through the survivor maps.
fn compact_problem(
    model: Model,
    original: Problem,
    coordinates: &Coordinates,
    cones: Vec<Cone>,
    before: Size,
) -> Problem {
    let compact_to_stable_columns = &coordinates.compact_to_stable_columns;
    let stable_to_compact_columns = Arc::new(inverse(compact_to_stable_columns, model.alive.len()));
    let n = compact_to_stable_columns.len();
    let p = if model.objective.p.revision == 0
        && compact_to_stable_columns
            .iter()
            .copied()
            .eq(0..before.variables)
    {
        original.p
    } else if original.p.is_some() || model.objective.p.nnz() != 0 {
        Some(Quadratic::from_sparse(
            model.objective.p,
            Arc::clone(compact_to_stable_columns),
            Arc::clone(&stable_to_compact_columns),
        ))
    } else {
        None
    };
    let mut c = model.objective.c;
    let mut bounds = model.bounds;
    for (j, &old) in compact_to_stable_columns.iter().enumerate() {
        c[j] = c[old];
        bounds[j] = bounds[old];
    }
    c.truncate(n);
    bounds.truncate(n);
    let mut domains = Vec::new();
    for &i in &coordinates.compact_to_stable_linear_rows {
        let RowDomain::Linear(b) = model.rows[i] else {
            unreachable!()
        };
        domains.push(Constraint::Linear(b));
    }
    let mut block_index = 0;
    let mut last_block = None;
    for &i in &coordinates.compact_to_stable_conic_rows {
        let RowDomain::Cone { rhs, block } = model.rows[i] else {
            unreachable!()
        };
        if last_block.is_some_and(|b| b != block) {
            block_index += 1;
        }
        last_block = Some(block);
        domains.push(Constraint::Cone {
            rhs,
            block: block_index,
        });
    }
    let (linear_rows, conic_rows) = row_indices(&domains);
    Problem {
        p,
        c,
        objective_constant: model.objective.constant,
        variable_bounds: bounds,
        rows: domains,
        cones,
        linear_rows,
        conic_rows,
        a: Matrix::Linked {
            matrix: Some(model.a),
            compact_to_stable_rows: coordinates
                .compact_to_stable_linear_rows
                .iter()
                .chain(&coordinates.compact_to_stable_conic_rows)
                .copied()
                .collect(),
            stable_to_compact_columns,
        },
    }
}

fn pack(mut model: Model, mut original: Problem, before: Size) -> Outcome {
    let survivors = survivors(&model);
    let (postsolve, cones) = build_postsolve(&mut model, &mut original, survivors, before);
    let coordinates = &postsolve.coordinates;
    if coordinates.compact_to_stable_columns.is_empty()
        && coordinates.compact_to_stable_linear_rows.is_empty()
        && coordinates.compact_to_stable_conic_rows.is_empty()
    {
        return Outcome::Solved(postsolve.recover_solution(SolutionRef {
            x: &[],
            y: &[],
            z: &[],
            conic_dual: &[],
            conic_slack: &[],
        }));
    }
    let problem = compact_problem(model, original, &postsolve.coordinates, cones, before);
    Outcome::Reduced(Box::new(ReducedProblem { problem, postsolve }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{matrix::CscMatrix, problem::ProblemData};
    #[test]
    fn auxiliary_replacing_an_original_column_rebuilds_the_hessian() {
        let p = CscMatrix::from_triplets(2, 2, vec![0, 1], vec![0, 1], vec![2., 3.]).unwrap();
        let mut original = Problem::from(ProblemData {
            p: Some(p),
            c: vec![0.; 2],
            objective_constant: 0.,
            a: CscMatrix::from_parts(1, 2, vec![0, 1, 2], vec![0, 0], vec![1., 1.]),
            rows: vec![Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 2.,
            })],
            variable_bounds: vec![],
            cones: vec![],
        });
        let mut model = Model::from_parts(
            original.working_matrix(),
            original.take_objective(),
            vec![RowDomain::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 2.,
            })],
            vec![Bounds::FREE; 2],
        );
        let before = size(&model);
        assert!(model.fix(0, 1.));
        model.add_variable(Bounds::FREE);
        let Outcome::Reduced(r) = pack(model, original, before) else {
            panic!()
        };
        let p = r.problem.p().unwrap();
        assert!(p.as_csc().is_none());
        assert_eq!(p.column(0).collect::<Vec<_>>(), [(0, 3.)]);
        assert_eq!(p.column(1).count(), 0);
        assert_eq!(r.problem.objective_constant(), 1.);
    }
}
