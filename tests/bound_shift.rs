use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::{CertificateRef, Solution},
    problem::{Bounds, Constraint, Problem},
    settings::{Rules, WorkLimit},
};

fn settings() -> Settings {
    Settings {
        rules: Rules {
            bound_shift: true,
            ..Rules::none()
        },
        ..Settings::default()
    }
}
fn fixture(lower: bool, pivot_sign: f64) -> Problem {
    let sign = if lower { 1. } else { -1. };
    let pivot = 2. * pivot_sign;
    let slope = 0.5 * sign;
    let offset = 2. * sign;
    let rhs = pivot * offset;
    let domain = if (pivot > 0.) == lower {
        Bounds {
            lower: rhs,
            upper: f64::INFINITY,
        }
    } else {
        Bounds {
            lower: f64::NEG_INFINITY,
            upper: rhs,
        }
    };
    Problem {
        p: None,
        c: vec![2. * sign, 1.],
        c0: 3.,
        a: CscMatrix::from_triplets(
            2,
            2,
            vec![0, 1, 0, 1],
            vec![0, 0, 1, 1],
            vec![pivot, sign, -pivot * slope, 1.],
        )
        .unwrap(),
        rows: vec![
            Constraint::Linear(domain),
            Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 8.,
            }),
        ],
        variable_bounds: vec![
            if lower {
                Bounds {
                    lower: offset,
                    upper: f64::INFINITY,
                }
            } else {
                Bounds {
                    lower: f64::NEG_INFINITY,
                    upper: offset,
                }
            },
            Bounds {
                lower: 0.,
                upper: 4.,
            },
        ],
        cones: vec![],
    }
}
fn stationary(p: &Problem, s: &Solution) {
    for j in 0..p.c.len() {
        let ay: f64 = (p.a.column_pointers()[j]..p.a.column_pointers()[j + 1])
            .map(|k| p.a.values()[k] * s.y[p.a.row_indices()[k]])
            .sum();
        assert!((p.c[j] - ay - s.z[j]).abs() < 1e-12, "column {j}");
    }
}
fn objective(p: &Problem, x: &[f64]) -> f64 {
    p.c0 + p.c.iter().zip(x).map(|(c, x)| c * x).sum::<f64>()
}

#[test]
fn shift_preserves_objective_stationarity_and_warm_starts_for_all_signs() {
    for lower in [false, true] {
        for pivot_sign in [-1., 1.] {
            let input = fixture(lower, pivot_sign);
            let sign = if lower { 1. } else { -1. };
            let result = Presolver::new(settings()).unwrap().presolve(input.clone());
            let Outcome::Reduced(reduced) = result.outcome else {
                panic!("expected shift")
            };
            assert_eq!(
                (reduced.problem.c.len(), reduced.problem.rows.len()),
                (2, 1)
            );
            let point = Solution {
                x: vec![2. * sign, 0.],
                y: vec![sign / pivot_sign, 0.],
                z: vec![0., 2.],
                conic_dual: vec![],
                conic_slack: vec![],
            };
            stationary(&input, &point);
            let warm = reduced.postsolve.reduce_warm_start(point.as_ref());
            assert_eq!(warm.x, vec![0., 0.]);
            stationary(&reduced.problem, &warm);
            let recovered = reduced.postsolve.recover_solution(warm.as_ref());
            assert_eq!(point.x, recovered.x);
            assert_eq!(point.y, recovered.y);
            assert_eq!(point.z, recovered.z);
            assert_eq!(
                objective(&input, &point.x),
                objective(&reduced.problem, &warm.x)
            );
            // Also project an active multiplier on the redundant original bound.
            let alternate = Solution {
                x: point.x.clone(),
                y: vec![0., 0.],
                z: vec![2. * sign, 1.],
                conic_dual: vec![],
                conic_slack: vec![],
            };
            let warm = reduced.postsolve.reduce_warm_start(alternate.as_ref());
            stationary(&reduced.problem, &warm);
            stationary(&input, &reduced.postsolve.recover_solution(warm.as_ref()));
        }
    }
}

#[test]
fn shift_recovers_farkas_multipliers_and_recession_directions() {
    let mut input = fixture(true, 1.);
    input.a = CscMatrix::from_triplets(
        2,
        2,
        vec![0, 1, 0, 1],
        vec![0, 0, 1, 1],
        vec![2., 1., -1., -0.5],
    )
    .unwrap();
    input.rows[1] = Constraint::Linear(Bounds {
        lower: f64::NEG_INFINITY,
        upper: 1.,
    });
    let result = Presolver::new(settings()).unwrap().presolve(input.clone());
    let Outcome::Reduced(reduced) = result.outcome else {
        panic!("expected reduced")
    };
    let certificate = reduced
        .postsolve
        .recover_primal_certificate(CertificateRef {
            y: &[-1.],
            z: &[1., 0.],
            conic_dual: &[],
        });
    assert_eq!(certificate.y, vec![0.5, -1.]);
    assert_eq!(certificate.z, vec![0., 0.]);
    assert_eq!(2. * certificate.y[0] + certificate.y[1], 0.);
    assert_eq!(-certificate.y[0] - 0.5 * certificate.y[1], 0.);
    assert!(4. * certificate.y[0] + certificate.y[1] > 0.);
    // The offset 2 is omitted when lifting a direction.
    assert_eq!(
        reduced.postsolve.recover_primal_ray(&[1., 0.]),
        vec![1., 0.]
    );
    assert_eq!(
        reduced.postsolve.recover_primal_ray(&[0., 2.]),
        vec![1., 2.]
    );
}

#[test]
fn shift_limits_and_quadratic_pivots_reject_without_partial_changes() {
    let original = fixture(true, 1.);
    let mut variants = vec![];
    let mut s = settings();
    s.bound_shift.max_column_length = 1;
    variants.push(s);
    let mut s = settings();
    s.bound_shift.max_ratio = 0.25;
    variants.push(s);
    let mut s = settings();
    s.bound_shift.work_limit = WorkLimit::Entries(1);
    variants.push(s);
    for s in variants {
        assert!(matches!(
            Presolver::new(s)
                .unwrap()
                .presolve(original.clone())
                .outcome,
            Outcome::Unchanged(_)
        ));
    }
    let mut input = original.clone();
    input.p = Some(CscMatrix::from_triplets(2, 2, vec![0], vec![0], vec![1.]).unwrap());
    assert!(matches!(
        Presolver::new(settings()).unwrap().presolve(input).outcome,
        Outcome::Unchanged(_)
    ));
    let mut input = original;
    input.variable_bounds[0].upper = 100.;
    assert!(matches!(
        Presolver::new(settings()).unwrap().presolve(input).outcome,
        Outcome::Unchanged(_)
    ));
}

#[test]
fn fill_limit_is_transactional_and_chained_shifts_reverse_in_order() {
    let mut input = fixture(true, 1.);
    input.a = CscMatrix::from_triplets(
        2,
        3,
        vec![0, 1, 0, 1],
        vec![0, 0, 1, 2],
        vec![2., 1., -1., 1.],
    )
    .unwrap();
    input.c.push(1.);
    input.variable_bounds.push(Bounds {
        lower: 0.,
        upper: 1.,
    });
    let mut limited = settings();
    limited.bound_shift.max_fill = 0;
    assert!(matches!(
        Presolver::new(limited)
            .unwrap()
            .presolve(input.clone())
            .outcome,
        Outcome::Unchanged(_)
    ));
    let result = Presolver::new(settings()).unwrap().presolve(input);
    assert_eq!(result.stats.after.unwrap().linear_rows, 1);

    let input = Problem {
        p: None,
        c: vec![1.; 3],
        c0: 0.,
        a: CscMatrix::from_triplets(
            3,
            3,
            vec![0, 2, 0, 1, 2, 1, 2],
            vec![0, 0, 1, 1, 1, 2, 2],
            vec![1., 1., -2., 1., 1., -2., 1.],
        )
        .unwrap(),
        rows: vec![
            Constraint::Linear(Bounds {
                lower: 1.,
                upper: f64::INFINITY,
            }),
            Constraint::Linear(Bounds {
                lower: 1.,
                upper: f64::INFINITY,
            }),
            Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 100.,
            }),
        ],
        variable_bounds: vec![
            Bounds {
                lower: 0.,
                upper: f64::INFINITY
            };
            3
        ],
        cones: vec![],
    };
    let result = Presolver::new(settings()).unwrap().presolve(input.clone());
    let Outcome::Reduced(reduced) = result.outcome else {
        panic!("expected shifts")
    };
    assert_eq!(reduced.problem.rows.len(), 1);
    let point = Solution {
        x: vec![3., 1., 0.],
        y: vec![1., 3., 0.],
        z: vec![0., 0., 7.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    stationary(&input, &point);
    let warm = reduced.postsolve.reduce_warm_start(point.as_ref());
    assert_eq!(warm.x, vec![0.; 3]);
    stationary(&reduced.problem, &warm);
    assert_eq!(
        objective(&input, &point.x),
        objective(&reduced.problem, &warm.x)
    );
    let recovered = reduced.postsolve.recover_solution(warm.as_ref());
    assert_eq!(point.x, recovered.x);
    assert_eq!(point.y, recovered.y);
    assert_eq!(point.z, recovered.z);
}

#[test]
fn invalid_ratio_limits_and_underflowing_updates_skip_the_transaction() {
    for ratio in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1., 0.] {
        let mut s = settings();
        s.bound_shift.max_ratio = ratio;
        assert!(matches!(
            Presolver::new(s)
                .unwrap()
                .presolve(fixture(true, 1.))
                .outcome,
            Outcome::Unchanged(_)
        ));
    }
    let mut input = fixture(true, 1.);
    // x >= .5*y + 2 would add .5*min_subnormal to the last column.
    // That product must not silently disappear during the coordinate change.
    input.a = CscMatrix::from_triplets(
        2,
        3,
        vec![0, 1, 0, 1],
        vec![0, 0, 1, 2],
        vec![2., f64::from_bits(1), -1., 1.],
    )
    .unwrap();
    input.c.push(1.);
    input.variable_bounds.push(Bounds {
        lower: 0.,
        upper: 1.,
    });
    assert!(matches!(
        Presolver::new(settings()).unwrap().presolve(input).outcome,
        Outcome::Unchanged(_)
    ));
}
