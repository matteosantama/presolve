//! Presolve entry point, working-model preparation, and result packing.
use crate::{
    core::{
        execution::Executor,
        model::{Model, RowDomain},
        schedule::Limits,
    },
    matrix::quadratic::Quadratic,
    postsolve::{
        Coordinates, Postsolve, PrimalCertificate, SolutionRef, original_point,
        tape::{Point, Recovery},
    },
    problem::{Bounds, Constraint, Matrix, Problem, row_indices},
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
    pub fn presolve(&self, problem: impl Into<Problem>) -> PresolveResult {
        let start = Instant::now();
        presolve_owned(problem.into(), &self.settings, &self.executor, start)
    }
}

/// Presolve one problem with temporary execution resources.
/// Use `Presolver` to reuse the pool across calls. Initialization errors are
/// distinct from optimization outcomes such as infeasibility or unboundedness.
pub fn presolve(
    problem: impl Into<Problem>,
    settings: &Settings,
) -> Result<PresolveResult, InitError> {
    Ok(Presolver::new(settings.clone())?.presolve(problem))
}
fn presolve_owned(
    mut problem: Problem,
    settings: &Settings,
    executor: &Executor,
    start: Instant,
) -> PresolveResult {
    let n = problem.c.len();
    let rows: Vec<_> = problem
        .rows
        .iter()
        .map(|&row| match row {
            Constraint::Linear(b) => RowDomain::Linear(b),
            Constraint::Cone { rhs, block } => RowDomain::Cone { rhs, block },
        })
        .collect();
    let free_bounds = problem.variable_bounds.is_empty();
    let mut model = Model::from_parts(
        problem.working_matrix(),
        problem.take_objective(),
        rows,
        if problem.variable_bounds.is_empty() {
            vec![Bounds::FREE; n]
        } else {
            std::mem::take(&mut problem.variable_bounds)
        },
    );
    model.objective.constant = problem.objective_constant;
    model.rules = settings.rules;
    model.numerics = settings.numerics;
    model.set_cones(problem.cones.clone());
    let before = size(&model);
    let mut stats = Stats {
        elapsed: Duration::ZERO,
        time_limit_reached: false,
        before,
        after: Some(before),
        parallel_comparisons: 0,
        quadratic_changed: false,
    };
    let phases = model.run(
        Limits {
            time: settings.time_limit.saturating_sub(start.elapsed()),
            fill: settings.substitution_fill,
            sparsify: settings.rules.sparsification,
        },
        executor,
    );
    stats.quadratic_changed = model.objective.p.revision != 0;
    let outcome = match phases {
        Err(certificate) => {
            stats.after = None;
            // Domains can have been deleted; iterate the caller's original tags.
            let (linear, conic) = original_rows(&problem.rows);
            let p = original_point(certificate.point, n, &linear, &conic, Vec::new());
            match certificate.mode {
                Recovery::PrimalInfeasibility => Outcome::Infeasible(PrimalCertificate {
                    y: p.y,
                    z: p.z,
                    conic_dual: p.conic_dual,
                }),
                Recovery::DualInfeasibility => {
                    if let Some(mut point) = feasible_point(&model) {
                        point.truncate(n);
                        Outcome::Unbounded(UnboundednessCertificate { point, ray: p.x })
                    } else {
                        stats.after = Some(size(&model));
                        if model.revision == 0 {
                            problem.restore_working_matrix(model.a);
                            problem.c = model.objective.c;
                            if !free_bounds {
                                problem.variable_bounds = model.bounds;
                            }
                            Outcome::Unchanged(problem)
                        } else {
                            pack(model, problem, before)
                        }
                    }
                }
                Recovery::Solution => unreachable!(),
            }
        }
        Ok(phases) => {
            stats.parallel_comparisons = phases.parallel_comparisons;
            stats.time_limit_reached |= phases.time_limit;
            if model.revision == 0 {
                problem.restore_working_matrix(model.a);
                problem.c = model.objective.c;
                if !free_bounds {
                    problem.variable_bounds = model.bounds;
                }
                Outcome::Unchanged(problem)
            } else {
                stats.after = Some(size(&model));
                pack(model, problem, before)
            }
        }
    };
    stats.elapsed = start.elapsed();
    stats.time_limit_reached |= stats.elapsed >= settings.time_limit;
    PresolveResult { outcome, stats }
}

/// Verify a cheap candidate only when a recession ray has been discovered.
/// Failure to find a witness simply returns a reduced or unchanged problem.
fn feasible_point(model: &Model) -> Option<Vec<f64>> {
    let mut p = Point::zeros(model.bounds.len(), model.rows.len());
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
fn original_rows(rows: &[Constraint]) -> (Vec<usize>, Vec<usize>) {
    let mut linear = Vec::new();
    let mut conic = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        match row {
            Constraint::Linear(_) => linear.push(i),
            Constraint::Cone { .. } => conic.push(i),
        }
    }
    (linear, conic)
}

fn pack(model: Model, mut original: Problem, before: Size) -> Outcome {
    let columns: Vec<_> = (0..model.alive.len()).filter(|&j| model.alive[j]).collect();
    let mut stable_to_compact_columns = vec![usize::MAX; model.alive.len()];
    for (at, &j) in columns.iter().enumerate() {
        stable_to_compact_columns[j] = at;
    }
    let stable_to_compact_columns = Arc::new(stable_to_compact_columns);
    let mut rows = Vec::new();
    let mut conic_rows = Vec::new();
    let mut cones = Vec::new();
    let mut domains = Vec::new();
    for (i, domain) in model.rows.iter().enumerate() {
        if matches!(domain, RowDomain::Linear(_)) {
            rows.push(i);
        }
    }
    for (block, block_rows) in model.cone_rows.iter().enumerate() {
        let surviving: Vec<_> = block_rows
            .iter()
            .copied()
            .filter(|&i| matches!(model.rows[i],RowDomain::Cone {block:b,..} if b==block))
            .collect();
        if !surviving.is_empty() {
            cones.push(model.cones[block]);
            conic_rows.extend(surviving);
        }
    }
    let input_linear = std::mem::take(&mut original.linear_rows);
    let input_conic = std::mem::take(&mut original.conic_rows);
    let solved = columns.is_empty() && rows.is_empty() && conic_rows.is_empty();
    let coordinates = Coordinates {
        compact_to_stable_columns: Arc::new(columns),
        compact_to_stable_linear_rows: rows,
        compact_to_stable_conic_rows: conic_rows,
    };
    let mut direct_slacks = Vec::new();
    let mut original_positions = vec![
        usize::MAX;
        if input_conic.is_empty() {
            0
        } else {
            model.rows.len()
        }
    ];
    for (pos, &i) in input_conic.iter().enumerate() {
        original_positions[i] = pos;
    }
    for (at, &i) in coordinates.compact_to_stable_conic_rows.iter().enumerate() {
        let pos = original_positions[i];
        direct_slacks.push((pos, at));
    }
    let postsolve = Postsolve {
        total_columns: model.bounds.len(),
        total_rows: model.rows.len(),
        original_columns: before.variables,
        tape: model.postsolve,
        direct_slacks,
        coordinates,
        input_linear,
        input_conic,
    };
    let coordinates = &postsolve.coordinates;
    if solved {
        return Outcome::Solved(postsolve.recover_solution(SolutionRef {
            x: &[],
            y: &[],
            z: &[],
            conic_dual: &[],
            conic_slack: &[],
        }));
    }
    let objective_constant = model.objective.constant;
    let n = coordinates.compact_to_stable_columns.len();
    let p = if model.objective.p.revision == 0
        && coordinates
            .compact_to_stable_columns
            .iter()
            .copied()
            .eq(0..before.variables)
    {
        original.p
    } else if original.p.is_some() || model.objective.p.nnz() != 0 {
        Some(Quadratic::from_sparse(
            model.objective.p,
            Arc::clone(&coordinates.compact_to_stable_columns),
            Arc::clone(&stable_to_compact_columns),
        ))
    } else {
        None
    };
    let mut c = model.objective.c;
    let mut bounds = model.bounds;
    for (j, &old) in coordinates.compact_to_stable_columns.iter().enumerate() {
        c[j] = c[old];
        bounds[j] = bounds[old];
    }
    c.truncate(n);
    bounds.truncate(n);
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
    let problem = Problem {
        p,
        c,
        objective_constant,
        variable_bounds: bounds,
        rows: domains,
        cones,
        linear_rows,
        conic_rows,
        a: Matrix::Linked {
            matrix: model.a,
            compact_to_stable_rows: coordinates
                .compact_to_stable_linear_rows
                .iter()
                .chain(&coordinates.compact_to_stable_conic_rows)
                .copied()
                .collect(),
            stable_to_compact_columns,
        },
    };
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
