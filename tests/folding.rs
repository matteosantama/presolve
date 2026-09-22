use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::{CertificateRef, Solution},
    problem::{Bounds, Constraint, Problem},
    result::ReducedProblem,
    settings::{Rules, WorkLimit},
};

fn problem() -> Problem {
    Problem {
        p: None,
        c: vec![1.0; 4],
        c0: 3.0,
        a: CscMatrix::from_triplets(2, 4, vec![0, 0, 1, 1], vec![0, 1, 2, 3], vec![1.0; 4])
            .unwrap(),
        rows: vec![
            Constraint::Linear(Bounds {
                lower: 2.0,
                upper: 2.0
            });
            2
        ],
        variable_bounds: vec![
            Bounds {
                lower: 0.0,
                upper: f64::INFINITY
            };
            4
        ],
        cones: vec![],
    }
}
fn settings() -> Settings {
    Settings {
        rules: Rules {
            lp_folding: true,
            ..Rules::none()
        },
        ..Settings::default()
    }
}
fn fold(p: Problem) -> Box<ReducedProblem> {
    match Presolver::new(settings()).unwrap().presolve(p).outcome {
        Outcome::Reduced(r) => r,
        other => panic!("expected folded problem, got {other:?}"),
    }
}
fn point(x: Vec<f64>, y: Vec<f64>, z: Vec<f64>) -> Solution {
    Solution {
        x,
        y,
        z,
        conic_dual: vec![],
        conic_slack: vec![],
    }
}

#[test]
fn folding_projects_and_lifts_primal_dual_coordinates_with_different_class_sizes() {
    let r = fold(problem());
    assert_eq!(r.problem.c, vec![4.0]);
    assert_eq!(r.problem.c0, 3.0);
    assert_eq!(
        r.problem.a,
        CscMatrix::from_triplets(1, 1, vec![0], vec![0], vec![2.0]).unwrap()
    );
    let reduced = point(vec![1.0], vec![2.0], vec![0.0]);
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, vec![1.0; 4]);
    assert_eq!(recovered.y, vec![1.0; 2]);
    assert_eq!(recovered.z, vec![0.0; 4]);
    // Feasible but asymmetric points project onto a class-constant optimum.
    let warm = point(
        vec![0.5, 1.5, 0.25, 1.75],
        vec![0.5, 1.5],
        vec![1.0, 2.0, 3.0, 4.0],
    );
    let mapped = r.postsolve.reduce_warm_start(warm.as_ref());
    assert_eq!(mapped.x, vec![1.0]);
    assert_eq!(mapped.y, vec![2.0]);
    assert_eq!(mapped.z, vec![10.0]);
    let lifted = r.postsolve.recover_solution(mapped.as_ref());
    assert_eq!(lifted.z, vec![2.5; 4]);
}

#[test]
fn folding_invalidates_the_original_empty_hessian_dimensions() {
    let mut p = problem();
    p.p = Some(CscMatrix::zeros(4, 4).unwrap());
    let result = Presolver::new(settings()).unwrap().presolve(p);
    assert!(result.stats.quadratic_changed);
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected folding");
    };
    assert_eq!(r.problem.variable_count(), 1);
    if let Some(p) = &r.problem.p {
        assert_eq!(p.columns(), 1);
        assert_eq!(p.rows(), 1);
    }
}

#[test]
fn folding_lifts_farkas_and_recession_certificates() {
    let mut p = problem();
    p.rows = vec![
        Constraint::Linear(Bounds {
            lower: -2.0,
            upper: -2.0
        });
        2
    ];
    let r = fold(p);
    let c = r.postsolve.recover_primal_certificate(CertificateRef {
        y: &[-1.0],
        z: &[2.0],
        conic_dual: &[],
    });
    assert_eq!(c.y, vec![-0.5; 2]);
    assert_eq!(c.z, vec![0.5; 4]);
    for z in &c.z {
        assert_eq!(c.y[0] + z, 0.0);
    }
    assert!(-2.0 * c.y.iter().sum::<f64>() > 0.0);
    let mut p = problem();
    p.c = vec![-1.0; 4];
    p.rows = vec![
        Constraint::Linear(Bounds {
            lower: 0.0,
            upper: f64::INFINITY
        });
        2
    ];
    let r = fold(p);
    let ray = r.postsolve.recover_primal_ray(&[1.0]);
    assert_eq!(ray, vec![1.0; 4]);
    assert!(-ray.iter().sum::<f64>() < 0.0);
}

#[test]
fn folding_is_opt_in_and_rejects_quadratic_and_unfinished_partitions() {
    let mut s = settings();
    s.rules.lp_folding = false;
    assert!(matches!(
        Presolver::new(s).unwrap().presolve(problem()).outcome,
        Outcome::Unchanged(_)
    ));
    let mut p = problem();
    p.p = Some(CscMatrix::from_triplets(4, 4, vec![0], vec![0], vec![1.0]).unwrap());
    assert!(matches!(
        Presolver::new(settings()).unwrap().presolve(p).outcome,
        Outcome::Unchanged(_)
    ));
    let mut s = settings();
    s.folding.work_limit = WorkLimit::Entries(0);
    assert!(matches!(
        Presolver::new(s).unwrap().presolve(problem()).outcome,
        Outcome::Unchanged(_)
    ));
    // Initially identical row bounds need multiple rounds to distinguish these rows.
    let mut p = problem();
    p.c[0] = 2.0;
    let mut s = settings();
    s.folding.max_rounds = 1;
    assert!(matches!(
        Presolver::new(s).unwrap().presolve(p).outcome,
        Outcome::Unchanged(_)
    ));
}

#[test]
fn folding_preserves_distinct_bounds_and_rounds_scaled_objectives_once() {
    let mut p = problem();
    p.variable_bounds[0].upper = 1.0;
    p.variable_bounds[1].upper = 2.0;
    let r = fold(p);
    // Singleton classes and a folded pair must retain distinct, consecutive
    // indices for both the reduced matrix and the recovery mapping.
    assert_eq!(r.problem.c, vec![1.0, 1.0, 2.0]);
    assert_eq!(r.problem.rows.len(), 2);
    assert_eq!(
        r.problem.a,
        CscMatrix::from_triplets(2, 3, vec![0, 0, 1], vec![0, 1, 2], vec![1.0, 1.0, 2.0]).unwrap()
    );
    let reduced = point(vec![0.5, 1.5, 1.0], vec![2.0, 3.0], vec![4.0, 5.0, 6.0]);
    let lifted = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(lifted.x, vec![0.5, 1.5, 1.0, 1.0]);
    assert_eq!(lifted.y, vec![2.0, 3.0]);
    assert_eq!(lifted.z, vec![4.0, 5.0, 3.0, 3.0]);
    let mut p = problem();
    p.c = vec![0.1; 4];
    // Four times 0.1 is representable exactly although sequential addition is not.
    assert_eq!(fold(p).problem.c, vec![0.4]);
    let p = Problem {
        p: None,
        c: vec![0.1; 3],
        c0: 0.0,
        a: CscMatrix::from_triplets(1, 3, vec![0; 3], vec![0, 1, 2], vec![1.0; 3]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: 1.0,
            upper: 1.0,
        })],
        variable_bounds: vec![
            Bounds {
                lower: 0.0,
                upper: f64::INFINITY
            };
            3
        ],
        cones: vec![],
    };
    assert_eq!(fold(p).problem.c, vec![3.0 * 0.1]);
}

#[test]
fn folding_warm_start_mean_stays_finite_when_the_sum_overflows() {
    let p = Problem {
        p: None,
        c: vec![0.; 4],
        c0: 0.,
        a: CscMatrix::from_triplets(1, 4, vec![0; 4], vec![0, 1, 2, 3], vec![0.25; 4]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: 1.5e308,
        })],
        variable_bounds: vec![Bounds::FREE; 4],
        cones: vec![],
    };
    let r = fold(p);
    for value in [1e308, -1e308, f64::from_bits(1)] {
        let warm = point(vec![value; 4], vec![0.], vec![0.; 4]);
        let mapped = r.postsolve.reduce_warm_start(warm.as_ref());
        assert_eq!(mapped.x, vec![value]);
        assert_eq!(
            r.postsolve.recover_solution(mapped.as_ref()).x,
            vec![value; 4]
        );
    }
    let warm = point(vec![1e308, 1e308, -1e308, -1e308], vec![0.], vec![0.; 4]);
    let mapped = r.postsolve.reduce_warm_start(warm.as_ref());
    assert_eq!(mapped.x, vec![0.]);
}
