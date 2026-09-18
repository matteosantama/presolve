use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, ProblemData},
    settings::Rules,
};

/// min ½(x² + 2xy + 3y²) + x with y >= 1 and x free. Only y has a row.
/// Minimizing over x gives x = -(1 + y); the reduced objective is y² - y - ½.
fn problem() -> ProblemData {
    ProblemData {
        p: Some(
            CscMatrix::from_triplets(2, 2, vec![0, 0, 1], vec![0, 1, 1], vec![1., 1., 3.]).unwrap(),
        ),
        c: vec![1., 0.],
        objective_constant: 0.,
        a: CscMatrix::from_triplets(1, 2, vec![0], vec![1], vec![1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: 1.,
            upper: f64::INFINITY,
        })],
        variable_bounds: vec![Bounds::FREE, Bounds::FREE],
        cones: vec![],
    }
}

fn only(rules: Rules) -> Settings {
    Settings {
        rules,
        ..Settings::default()
    }
}

#[test]
fn coupled_free_column_is_minimized_out_of_the_objective() {
    let settings = only(Rules {
        empty_columns: true,
        quadratic_elimination: true,
        ..Rules::none()
    });
    let result = Presolver::new(settings).unwrap().presolve(problem());
    assert!(result.stats.quadratic_changed);
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    assert_eq!(r.problem.variable_count(), 1);
    assert_eq!(r.problem.c(), [-1.]);
    assert_eq!(r.problem.objective_constant(), -0.5);
    let p = r.problem.p().unwrap();
    assert_eq!(p.column(0).collect::<Vec<_>>(), [(0, 2.)]);
    // Reduced optimum y = 1 with multiplier 1 on the row.
    let reduced = Solution {
        x: vec![1.],
        y: vec![1.],
        z: vec![0.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, [-2., 1.]);
    assert_eq!(recovered.z, [0., 0.]);
    assert_eq!(recovered.y, [1.]);
    // Stationarity in the original: P x + c - A^T y - z = 0.
    assert_eq!(recovered.x[0] + recovered.x[1] + 1., 0.);
    assert_eq!(recovered.x[0] + 3. * recovered.x[1] - recovered.y[0], 0.);
}

#[test]
fn bounded_or_flat_coupled_columns_are_retained() {
    let settings = only(Rules {
        empty_columns: true,
        quadratic_elimination: true,
        ..Rules::none()
    });
    let mut bounded = problem();
    bounded.variable_bounds[0] = Bounds {
        lower: -1.,
        upper: 1.,
    };
    let result = Presolver::new(settings.clone()).unwrap().presolve(bounded);
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
    let mut flat = problem();
    flat.p = Some(CscMatrix::from_triplets(2, 2, vec![0, 1], vec![1, 1], vec![1., 3.]).unwrap());
    let result = Presolver::new(settings).unwrap().presolve(flat);
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
    let off = only(Rules {
        empty_columns: true,
        ..Rules::none()
    });
    let result = Presolver::new(off).unwrap().presolve(problem());
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
}

#[test]
fn default_pipeline_solves_the_example() {
    let result = Presolver::new(Settings::default())
        .unwrap()
        .presolve(problem());
    let Outcome::Solved(solution) = result.outcome else {
        panic!("expected a solved problem")
    };
    assert_eq!(solution.x, [-2., 1.]);
}
