use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, ProblemData},
    settings::Rules,
};

fn only_dual_propagation() -> Settings {
    Settings {
        rules: Rules {
            dual_propagation: true,
            ..Rules::none()
        },
        ..Settings::default()
    }
}

/// min -x subject to x - y <= 0, x free, 0 <= y <= 5.
/// The free column forces y_0 = -1, so the row is tight at every optimum.
fn tight_row() -> ProblemData {
    ProblemData {
        p: None,
        c: vec![-1., 0.],
        objective_constant: 0.,
        a: CscMatrix::from_triplets(1, 2, vec![0, 0], vec![0, 1], vec![1., -1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: 0.,
        })],
        variable_bounds: vec![
            Bounds::FREE,
            Bounds {
                lower: 0.,
                upper: 5.,
            },
        ],
        cones: vec![],
    }
}

/// min x + w/2 subject to x + w >= 1, x >= 0, w free.
/// The free column forces y_0 = 1/2, so z_x = 1/2 > 0 fixes x at its bound
/// although the row locks x downward.
fn fixed_column() -> ProblemData {
    ProblemData {
        p: None,
        c: vec![1., 0.5],
        objective_constant: 0.,
        a: CscMatrix::from_triplets(1, 2, vec![0, 0], vec![0, 1], vec![1., 1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: 1.,
            upper: f64::INFINITY,
        })],
        variable_bounds: vec![
            Bounds {
                lower: 0.,
                upper: f64::INFINITY,
            },
            Bounds::FREE,
        ],
        cones: vec![],
    }
}

#[test]
fn strictly_negative_multiplier_makes_the_row_tight_and_fixes_its_bound() {
    let result = Presolver::new(only_dual_propagation())
        .unwrap()
        .presolve(tight_row());
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    // y_0 = -1 makes the row tight; z_y = y_0 < 0 also fixes y at 5.
    assert_eq!(r.problem.variable_count(), 1);
    assert_eq!(r.problem.row_bounds(0), Bounds::fixed(5.));
    // The solver's equality multiplier is the original upper-side multiplier.
    let reduced = Solution {
        x: vec![5.],
        y: vec![-1.],
        z: vec![0.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, [5., 5.]);
    assert_eq!(recovered.y, [-1.]);
    assert_eq!(recovered.z, [0., -1.]);
}

#[test]
fn default_pipeline_solves_the_tight_row_problem() {
    let result = Presolver::new(Settings::default())
        .unwrap()
        .presolve(tight_row());
    let Outcome::Solved(solution) = result.outcome else {
        panic!("expected a solved problem")
    };
    assert_eq!(solution.x, [5., 5.]);
    // Stationarity: c - A^T y - z = 0 in the original problem.
    assert_eq!(-1. - solution.y[0], solution.z[0]);
    assert_eq!(solution.y[0], solution.z[1]);
    assert!(solution.y[0] <= 0.);
}

#[test]
fn strictly_positive_reduced_cost_fixes_a_locked_column() {
    let result = Presolver::new(only_dual_propagation())
        .unwrap()
        .presolve(fixed_column());
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    // y_0 = 1/2 > 0 makes the row tight and z_x = 1/2 > 0 fixes x at zero.
    assert_eq!(r.problem.variable_count(), 1);
    assert_eq!(r.problem.row_bounds(0), Bounds::fixed(1.));
    let reduced = Solution {
        x: vec![1.],
        y: vec![0.5],
        z: vec![0.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, [0., 1.]);
    assert_eq!(recovered.y, [0.5]);
    assert_eq!(recovered.z, [0.5, 0.]);
}

#[test]
fn dual_infeasible_systems_produce_no_reductions() {
    // min -x - y subject to x - y <= 0, x, y >= 0 is unbounded. The
    // multiplier bounds contradict, so the pass must not tighten anything.
    let problem = ProblemData {
        p: None,
        c: vec![-1., -1.],
        objective_constant: 0.,
        a: CscMatrix::from_triplets(1, 2, vec![0, 0], vec![0, 1], vec![1., -1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: 0.,
        })],
        variable_bounds: vec![
            Bounds {
                lower: 0.,
                upper: f64::INFINITY
            };
            2
        ],
        cones: vec![],
    };
    let result = Presolver::new(only_dual_propagation())
        .unwrap()
        .presolve(problem);
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
}

#[test]
fn quadratic_columns_do_not_contribute_dual_rows() {
    // The same data as the tight row, but x carries curvature. Its reduced
    // cost depends on x, so nothing can be proved and the problem is unchanged.
    let mut problem = tight_row();
    problem.p = Some(CscMatrix::from_triplets(2, 2, vec![0], vec![0], vec![1.]).unwrap());
    let result = Presolver::new(only_dual_propagation())
        .unwrap()
        .presolve(problem);
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
}
