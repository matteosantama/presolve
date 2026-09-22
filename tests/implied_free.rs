use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, Problem},
    settings::{Rules, WorkLimit},
};

fn settings() -> Settings {
    Settings {
        equalities: presolve::settings::EqualitySettings {
            min_column_length: 1,
            ..Default::default()
        },
        rules: Rules {
            implied_free_equalities: true,
            ..Rules::none()
        },
        ..Settings::default()
    }
}

fn nonnegative() -> Bounds {
    Bounds {
        lower: 0.,
        upper: f64::INFINITY,
    }
}

fn long_row(sign: f64) -> Problem {
    // x = sum(y), x + t >= 2. Only x has its lower bound implied by
    // the equality. Its degree is two although the equation is long.
    let n = 12;
    let mut rows = vec![0; n - 1];
    rows.extend([1, 1]);
    let mut columns: Vec<_> = (0..n - 1).collect();
    columns.extend([0, n - 1]);
    let mut values = vec![-sign; n - 1];
    values[0] = sign;
    values.extend([1., 1.]);
    let mut c = vec![0.; n];
    c[0] = 1.;
    Problem {
        p: None,
        c,
        c0: 0.,
        a: CscMatrix::from_triplets(2, n, rows, columns, values).unwrap(),
        rows: vec![
            Constraint::Linear(Bounds::fixed(0.)),
            Constraint::Linear(Bounds {
                lower: 2.,
                upper: f64::INFINITY,
            }),
        ],
        variable_bounds: vec![nonnegative(); n],
        cones: vec![],
    }
}

#[test]
fn long_degree_two_pivot_removes_row_and_recovers_primal_and_dual() {
    for sign in [-1., 1.] {
        let result = Presolver::new(settings()).unwrap().presolve(long_row(sign));
        let Outcome::Reduced(r) = result.outcome else {
            panic!("expected reduction")
        };
        assert_eq!(r.problem.variable_count(), 11);
        assert_eq!(r.problem.rows.len(), 1);
        assert_eq!(r.problem.a.values().len(), 11);
        let mut x = vec![0.; 11];
        x[10] = 2.;
        let mut z = vec![1.; 11];
        z[10] = 0.;
        let reduced = Solution {
            x,
            y: vec![0.],
            z,
            conic_dual: vec![],
            conic_slack: vec![],
        };
        let recovered = r.postsolve.recover_solution(reduced.as_ref());
        assert_eq!(recovered.x[0], 0.);
        assert_eq!(recovered.x[11], 2.);
        assert_eq!(recovered.y, [sign, 0.]);
        assert_eq!(recovered.z[0], 0.);
        assert_eq!(&recovered.z[1..11], &[1.; 10]);
    }
}

#[test]
fn optional_rule_is_off_by_default() {
    assert!(!Rules::default().implied_free_equalities);
    let mut s = settings();
    s.rules.implied_free_equalities = false;
    assert!(matches!(
        Presolver::new(s).unwrap().presolve(long_row(1.)).outcome,
        Outcome::Unchanged(_)
    ));
}

#[test]
fn nonredundant_pivot_bound_is_preserved() {
    let p = Problem {
        p: None,
        c: vec![1.; 3],
        c0: 0.,
        a: CscMatrix::from_triplets(1, 3, vec![0; 3], vec![0, 1, 2], vec![1.; 3]).unwrap(),
        rows: vec![Constraint::Linear(Bounds::fixed(1.))],
        variable_bounds: vec![nonnegative(); 3],
        cones: vec![],
    };
    assert!(matches!(
        Presolver::new(settings()).unwrap().presolve(p).outcome,
        Outcome::Unchanged(_)
    ));
}

#[test]
fn both_finite_bounds_can_be_implied() {
    let p = Problem {
        p: None,
        c: vec![0., 0., 1.],
        c0: 0.,
        a: CscMatrix::from_triplets(1, 3, vec![0; 3], vec![0, 1, 2], vec![-1., -1., 1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds::fixed(0.))],
        variable_bounds: vec![
            Bounds {
                lower: 0.,
                upper: 1.,
            },
            Bounds {
                lower: 0.,
                upper: 1.,
            },
            Bounds {
                lower: 0.,
                upper: 2.,
            },
        ],
        cones: vec![],
    };
    let Outcome::Reduced(r) = Presolver::new(settings()).unwrap().presolve(p).outcome else {
        panic!("expected reduction")
    };
    assert_eq!(r.problem.variable_count(), 2);
    assert!(r.problem.rows.is_empty());
    assert_eq!(r.problem.c, [1., 1.]);
    let s = Solution {
        x: vec![0., 0.],
        y: vec![],
        z: vec![1., 1.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(s.as_ref());
    assert_eq!(recovered.x, [0., 0., 0.]);
    assert_eq!(recovered.y, [1.]);
    assert_eq!(recovered.z, [1., 1., 0.]);
}

#[test]
fn fill_and_work_limits_are_respected() {
    let mut fill = settings();
    fill.substitution_fill = 0;
    let mut work = settings();
    work.equalities.work_limit = WorkLimit::Entries(0);
    for s in [fill, work] {
        assert!(matches!(
            Presolver::new(s).unwrap().presolve(long_row(1.)).outcome,
            Outcome::Unchanged(_)
        ));
    }
}

#[test]
fn quadratic_pivot_obeys_hessian_growth_and_recovers_solution() {
    let p = Problem {
        p: Some(CscMatrix::from_triplets(3, 3, vec![0], vec![0], vec![1.]).unwrap()),
        c: vec![-1., 0., 0.],
        c0: 0.,
        a: CscMatrix::from_triplets(1, 3, vec![0; 3], vec![0, 1, 2], vec![1., -1., -1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds::fixed(0.))],
        variable_bounds: vec![nonnegative(); 3],
        cones: vec![],
    };
    let mut s = settings();
    s.equalities.require_linear_variable = false;
    assert!(matches!(
        Presolver::new(s.clone())
            .unwrap()
            .presolve(p.clone())
            .outcome,
        Outcome::Unchanged(_)
    ));
    s.allow_hessian_growth = true;
    s.equalities.preserve_nonzeros = false;
    let Outcome::Reduced(r) = Presolver::new(s).unwrap().presolve(p).outcome else {
        panic!("expected quadratic reduction")
    };
    assert_eq!(r.problem.p.as_ref().unwrap().values(), &[1., 1., 1.]);
    assert_eq!(r.problem.c, [-1., -1.]);
    let reduced = Solution {
        x: vec![0.5, 0.5],
        y: vec![],
        z: vec![0., 0.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, [1., 0.5, 0.5]);
    assert_eq!(recovered.y, [0.]);
    assert_eq!(recovered.z, [0., 0., 0.]);
}

#[test]
fn zero_pivot_attempts_disable_aggregation() {
    let mut s = settings();
    s.equalities.max_pivot_attempts = 0;
    assert!(matches!(
        Presolver::new(s).unwrap().presolve(long_row(1.)).outcome,
        Outcome::Unchanged(_)
    ));
}
