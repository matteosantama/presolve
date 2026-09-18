use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::{CertificateRef, Solution},
    problem::{Bounds, Cone, Constraint, Problem},
    settings::Rules,
};

fn settings() -> Settings {
    Settings {
        rules: Rules {
            cones: true,
            ..Rules::none()
        },
        ..Settings::default()
    }
}

// min t, 3x+4y+12z >= 1, (t,x,-y,z,v) in Q_5.
fn problem() -> Problem {
    Problem {
        p: None,
        c: vec![1., 0., 0., 0., 0.],
        c0: 0.,
        a: CscMatrix::from_triplets(
            6,
            5,
            vec![0, 0, 0, 1, 2, 3, 4, 5],
            vec![1, 2, 3, 0, 1, 2, 3, 4],
            vec![3., 4., 12., -1., -1., 1., -1., -1.],
        )
        .unwrap(),
        rows: std::iter::once(Constraint::Linear(Bounds {
            lower: 1.,
            upper: f64::INFINITY,
        }))
        .chain((0..5).map(|_| Constraint::Cone { rhs: 0., block: 0 }))
        .collect(),
        variable_bounds: vec![Bounds::FREE; 5],
        cones: vec![Cone::SecondOrder(5)],
    }
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-12 * (1. + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}
fn equal(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        close(a, b);
    }
}
fn stationary(p: &Problem, y: &[f64], z: &[f64], w: &[f64], objective: bool) {
    for (j, &zj) in z.iter().enumerate() {
        let mut residual = if objective { p.c[j] } else { 0. } - zj;
        for (i, a) in p.a.as_ref().column(j) {
            residual += if i == 0 { -a * y[0] } else { a * w[i - 1] };
        }
        close(residual, 0.);
    }
}
fn in_soc(s: &[f64]) {
    assert!(s[0] + 1e-12 >= s[1..].iter().fold(0.0_f64, |n, &v| n.hypot(v)));
}

#[test]
fn aggregates_soc_direction_and_recovers_primal_dual_solution() {
    let p = problem();
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(p.clone())
        .outcome
    else {
        panic!()
    };
    assert_eq!(r.problem.variable_count(), 3);
    assert_eq!(r.problem.cones, [Cone::SecondOrder(3)]);
    assert_eq!(
        r.problem.a.as_ref().column(1).collect::<Vec<_>>(),
        [(0, 13.), (2, -1.)]
    );
    assert_eq!(r.postsolve.surviving_conic_coordinates(), [(0, 0), (4, 2)]);
    let reduced = Solution {
        x: vec![1. / 13., 1. / 13., 0.],
        y: vec![1. / 13.],
        z: vec![0.; 3],
        conic_dual: vec![1., -1., 0.],
        conic_slack: vec![1. / 13., 1. / 13., 0.],
    };
    let original = r.postsolve.recover_solution(reduced.as_ref());
    equal(
        &original.x,
        &[1. / 13., 3. / 169., 4. / 169., 12. / 169., 0.],
    );
    equal(
        &original.conic_dual,
        &[1., -3. / 13., 4. / 13., -12. / 13., 0.],
    );
    equal(
        &original.conic_slack,
        &[1. / 13., 3. / 169., -4. / 169., 12. / 169., 0.],
    );
    stationary(&p, &original.y, &original.z, &original.conic_dual, true);
    in_soc(&original.conic_slack);
    in_soc(&original.conic_dual);
    let warm = r.postsolve.reduce_warm_start(original.as_ref());
    equal(&warm.x, &reduced.x);
    equal(&warm.y, &reduced.y);
    equal(&warm.z, &reduced.z);
    equal(&warm.conic_dual, &reduced.conic_dual);
    equal(&warm.conic_slack, &reduced.conic_slack);
}

#[test]
fn warm_start_projects_off_optimal_primal_dual_and_independent_slacks() {
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(problem())
        .outcome
    else {
        panic!()
    };
    // Slacks deliberately differ from Gx: preserve their own projection.
    let original = Solution {
        x: vec![5., 1., 2., -1., 0.5],
        y: vec![0.7],
        z: vec![0.; 5],
        conic_dual: vec![5., 2., -1., 3., 0.25],
        conic_slack: vec![6., 2., -3., 1., 0.75],
    };
    let warm = r.postsolve.reduce_warm_start(original.as_ref());
    equal(&warm.x, &[5., -1. / 13., 0.5]);
    equal(&warm.conic_dual, &[5., 46. / 13., 0.25]);
    equal(&warm.conic_slack, &[6., 30. / 13., 0.75]);
    in_soc(&warm.conic_dual);
    in_soc(&warm.conic_slack);
}

#[test]
fn recovers_soc_farkas_certificate() {
    let mut p = problem();
    p.variable_bounds[0] = Bounds {
        lower: f64::NEG_INFINITY,
        upper: 0.,
    };
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(p.clone())
        .outcome
    else {
        panic!()
    };
    let cert = r.postsolve.recover_primal_certificate(CertificateRef {
        y: &[1.],
        z: &[-13., 0., 0.],
        conic_dual: &[13., -13., 0.],
    });
    equal(&cert.conic_dual, &[13., -3., 4., -12., 0.]);
    stationary(&p, &cert.y, &cert.z, &cert.conic_dual, false);
    in_soc(&cert.conic_dual);
    close(cert.y[0] + cert.z[0] * p.variable_bounds[0].upper, 1.);
}

#[test]
fn recovers_soc_recession_ray() {
    let mut p = problem();
    p.c[0] = -1.;
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(p.clone())
        .outcome
    else {
        panic!()
    };
    let ray = r.postsolve.recover_primal_ray(&[1., 0.5, 0.]);
    equal(&ray, &[1., 1.5 / 13., 2. / 13., 6. / 13., 0.]);
    in_soc(&[ray[0], ray[1], -ray[2], ray[3], ray[4]]);
    assert!(3. * ray[1] + 4. * ray[2] + 12. * ray[3] >= 0.);
    assert!(p.c.iter().zip(&ray).map(|(c, x)| c * x).sum::<f64>() < 0.);
}

#[test]
fn skips_bounded_curved_shifted_and_nonzero_objective_coordinates() {
    for mode in 0..4 {
        let mut p = problem();
        for j in 1..4 {
            match mode {
                0 => {
                    p.variable_bounds[j] = Bounds {
                        lower: 0.,
                        upper: f64::INFINITY,
                    }
                }
                1 => p.c[j] = 1.,
                2 => p.rows[j + 1] = Constraint::Cone { rhs: 1., block: 0 },
                _ => (),
            }
        }
        if mode == 3 {
            p.p = Some(
                CscMatrix::from_triplets(5, 5, vec![1, 2, 3], vec![1, 2, 3], vec![1.; 3]).unwrap(),
            );
        }
        let result = Presolver::new(settings()).unwrap().presolve(p);
        assert!(matches!(result.outcome, Outcome::Unchanged(_)));
    }
}

#[test]
fn aggregation_composes_with_soc_linearization() {
    let mut p = problem();
    p.c.pop();
    p.variable_bounds.pop();
    p.a = CscMatrix::from_triplets(
        6,
        4,
        vec![0, 0, 0, 1, 2, 3, 4],
        vec![1, 2, 3, 0, 1, 2, 3],
        vec![3., 4., 12., -1., -1., 1., -1.],
    )
    .unwrap();
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(p.clone())
        .outcome
    else {
        panic!()
    };
    assert_eq!(r.problem.variable_count(), 2);
    assert!(r.problem.cones.is_empty());
    let reduced = Solution {
        x: vec![1. / 13., 1. / 13.],
        y: vec![1. / 13., 0., -1.],
        z: vec![0.; 2],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let original = r.postsolve.recover_solution(reduced.as_ref());
    equal(&original.x, &[1. / 13., 3. / 169., 4. / 169., 12. / 169.]);
    equal(
        &original.conic_dual,
        &[1., -3. / 13., 4. / 13., -12. / 13., 0.],
    );
    equal(
        &original.conic_slack,
        &[1. / 13., 3. / 169., -4. / 169., 12. / 169., 0.],
    );
    stationary(&p, &original.y, &original.z, &original.conic_dual, true);
    let warm = r.postsolve.reduce_warm_start(original.as_ref());
    equal(&warm.y, &reduced.y);
    equal(&warm.x, &reduced.x);
    assert!(warm.conic_slack.is_empty());
}

#[test]
fn projected_slack_round_trip_is_independent_of_primal_iterate() {
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(problem())
        .outcome
    else {
        panic!()
    };
    let reduced = Solution {
        x: vec![5., 1., 0.],
        y: vec![0.],
        z: vec![0.; 3],
        conic_dual: vec![3., 1., 0.],
        conic_slack: vec![4., 2., 0.5],
    };
    let original = r.postsolve.recover_solution(reduced.as_ref());
    equal(
        &original.conic_slack,
        &[4., 6. / 13., -8. / 13., 24. / 13., 0.5],
    );
    let warm = r.postsolve.reduce_warm_start(original.as_ref());
    equal(&warm.conic_slack, &reduced.conic_slack);
    equal(&warm.conic_dual, &reduced.conic_dual);
}

#[test]
fn skips_nonunit_coordinates_and_unrepresentable_normalization() {
    for (coefficients, cone_coefficients) in [
        ([3., 4., 12.], [-2., 2., -2.]),
        ([1.1e308; 3], [-1., 1., -1.]),
        ([f64::from_bits(1); 3], [-1., 1., -1.]),
        ([f64::MIN_POSITIVE, 1e308, 1e308], [-1., 1., -1.]),
    ] {
        let mut p = problem();
        p.a = CscMatrix::from_triplets(
            6,
            5,
            vec![0, 0, 0, 1, 2, 3, 4, 5],
            vec![1, 2, 3, 0, 1, 2, 3, 4],
            vec![
                coefficients[0],
                coefficients[1],
                coefficients[2],
                -1.,
                cone_coefficients[0],
                cone_coefficients[1],
                cone_coefficients[2],
                -1.,
            ],
        )
        .unwrap();
        let result = Presolver::new(settings()).unwrap().presolve(p);
        assert!(matches!(result.outcome, Outcome::Unchanged(_)));
    }
}

#[test]
fn negative_external_coefficients_and_irrational_norm_preserve_stationarity() {
    let mut p = problem();
    p.a = CscMatrix::from_triplets(
        6,
        5,
        vec![0, 0, 0, 1, 2, 3, 4, 5],
        vec![1, 2, 3, 0, 1, 2, 3, 4],
        vec![1., -2., 3., -1., -1., 1., -1., -1.],
    )
    .unwrap();
    let Outcome::Reduced(r) = Presolver::new(settings())
        .unwrap()
        .presolve(p.clone())
        .outcome
    else {
        panic!()
    };
    let norm = 14_f64.sqrt();
    let reduced = Solution {
        x: vec![1. / norm, 1. / norm, 0.],
        y: vec![1. / norm],
        z: vec![0.; 3],
        conic_dual: vec![1., -1., 0.],
        conic_slack: vec![1. / norm, 1. / norm, 0.],
    };
    let original = r.postsolve.recover_solution(reduced.as_ref());
    equal(&original.x, &[1. / norm, 1. / 14., -2. / 14., 3. / 14., 0.]);
    stationary(&p, &original.y, &original.z, &original.conic_dual, true);
    in_soc(&original.conic_slack);
    in_soc(&original.conic_dual);
    let warm = r.postsolve.reduce_warm_start(original.as_ref());
    equal(&warm.x, &reduced.x);
    equal(&warm.conic_dual, &reduced.conic_dual);
}
