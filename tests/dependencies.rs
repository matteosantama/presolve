use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, ProblemData},
    settings::{Rules, WorkLimit},
};

fn fixture(quadratic: bool, perturbation: f64, inconsistent: bool) -> ProblemData {
    // Third equality is the sum of two nonparallel, dense equalities.
    let rows = [
        [1., 2., 1., 2., 1., 2.],
        [2., 1., 2., 1., 2., 1.],
        [3., 3., 3., 3., 3., 3. + perturbation],
    ];
    let rhs = if inconsistent {
        vec![9., 9., 19.]
    } else {
        vec![9., 9., 18. + perturbation]
    };
    let c = (0..6)
        .map(|j| rows[0][j] + 2. * rows[1][j] + 3. * rows[2][j] - if quadratic { 1. } else { 0. })
        .collect();
    ProblemData {
        p: quadratic.then(|| {
            CscMatrix::from_triplets(6, 6, (0..6).collect(), (0..6).collect(), vec![1.; 6]).unwrap()
        }),
        c,
        objective_constant: 2.,
        a: CscMatrix::from_triplets(
            3,
            6,
            (0..3).flat_map(|i| [i; 6]).collect(),
            (0..3).flat_map(|_| 0..6).collect(),
            rows.into_iter().flatten().collect(),
        )
        .unwrap(),
        rows: rhs
            .into_iter()
            .map(|v| Constraint::Linear(Bounds::fixed(v)))
            .collect(),
        variable_bounds: vec![
            Bounds {
                lower: -10.,
                upper: 10.
            };
            6
        ],
        cones: vec![],
    }
}
fn settings() -> Settings {
    let mut settings = Settings {
        rules: Rules {
            equality_dependencies: true,
            ..Rules::none()
        },
        ..Settings::default()
    };
    settings.dependencies.work_limit = WorkLimit::Unlimited;
    settings
}
fn stationarity(data: &ProblemData, point: &Solution) {
    for j in 0..data.c.len() {
        let a = (data.a.column_pointers()[j]..data.a.column_pointers()[j + 1])
            .map(|k| data.a.values()[k] * point.y[data.a.row_indices()[k]])
            .sum::<f64>();
        let gradient = data.c[j] + if data.p.is_some() { point.x[j] } else { 0. };
        assert!((gradient - a - point.z[j]).abs() < 1e-12);
    }
}
#[test]
fn dependent_equalities_preserve_quadratic_structure_and_dual_warm_starts() {
    for quadratic in [false, true] {
        let input = fixture(quadratic, 0., false);
        let result = Presolver::new(settings()).unwrap().presolve(input.clone());
        let after = result.stats.after.unwrap();
        assert_eq!(
            (after.variables, after.linear_rows, after.a_nonzeros),
            (6, 2, 12)
        );
        assert!(!result.stats.quadratic_changed);
        let Outcome::Reduced(reduced) = result.outcome else {
            panic!("expected reduction")
        };
        let point = Solution {
            x: vec![1.; 6],
            y: vec![1., 2., 3.],
            z: vec![0.; 6],
            conic_dual: vec![],
            conic_slack: vec![],
        };
        stationarity(&input, &point);
        let warm = reduced.postsolve.reduce_warm_start(point.as_ref());
        let recovered = reduced.postsolve.recover_solution(warm.as_ref());
        assert_eq!(recovered.x, point.x);
        stationarity(&input, &recovered);
        let output = reduced.problem.into_csc();
        assert_eq!(input.p, output.p);
        assert_eq!(input.c, output.c);
        assert_eq!(input.variable_bounds, output.variable_bounds);
        stationarity(&output, &warm);
    }
}
#[test]
fn dependencies_detect_a_contradiction_but_retain_near_dependencies() {
    let input = fixture(true, 0., true);
    let result = Presolver::new(settings()).unwrap().presolve(input.clone());
    let Outcome::Infeasible(certificate) = result.outcome else {
        panic!("expected certificate")
    };
    for j in 0..6 {
        let value: f64 = (input.a.column_pointers()[j]..input.a.column_pointers()[j + 1])
            .map(|k| input.a.values()[k] * certificate.y[input.a.row_indices()[k]])
            .sum();
        assert_eq!(value, 0.);
    }
    let contradiction: f64 = input
        .rows
        .iter()
        .zip(&certificate.y)
        .map(|(r, y)| {
            let Constraint::Linear(b) = r else {
                unreachable!()
            };
            b.lower * y
        })
        .sum();
    assert!(contradiction > 0.);
    let result = Presolver::new(settings())
        .unwrap()
        .presolve(fixture(true, 1e-10, false));
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
}
#[test]
fn scratch_and_work_limits_do_not_partially_transform_the_problem() {
    let mut variants = Vec::new();
    let mut s = settings();
    s.dependencies.max_row_length = 5;
    variants.push(s);
    let mut s = settings();
    s.dependencies.max_basis_rows = 1;
    variants.push(s);
    let mut s = settings();
    s.dependencies.work_limit = WorkLimit::Entries(1);
    variants.push(s);
    let mut s = settings();
    s.time_limit = std::time::Duration::ZERO;
    variants.push(s);
    for s in variants {
        let result = Presolver::new(s)
            .unwrap()
            .presolve(fixture(true, 0., false));
        assert!(matches!(result.outcome, Outcome::Unchanged(_)));
    }
}

#[test]
fn multirow_proofs_survive_prior_dependency_deletions() {
    let (n, rank, m) = (8, 4, 12);
    let mut rows: Vec<Vec<f64>> = (0..rank)
        .map(|i| {
            (0usize..n)
                .map(|j| {
                    if (i & j).count_ones() % 2 == 0 {
                        1.
                    } else {
                        -1.
                    }
                })
                .collect()
        })
        .collect();
    for i in rank..m {
        rows.push(
            (0..n)
                .map(|j| {
                    (0..rank)
                        .map(|k| (1 + (i + k) % 3) as f64 * rows[k][j])
                        .sum()
                })
                .collect(),
        );
    }
    let y: Vec<_> = (0..m).map(|i| (i as f64 - 4.) / 2.).collect();
    let input = ProblemData {
        p: Some(
            CscMatrix::from_triplets(n, n, (0..n).collect(), (0..n).collect(), vec![1.; n])
                .unwrap(),
        ),
        c: (0..n)
            .map(|j| -1. + (0..m).map(|i| rows[i][j] * y[i]).sum::<f64>())
            .collect(),
        objective_constant: 0.,
        a: CscMatrix::from_triplets(
            m,
            n,
            (0..m).flat_map(|i| [i; 8]).collect(),
            (0..m).flat_map(|_| 0..n).collect(),
            rows.iter().flatten().copied().collect(),
        )
        .unwrap(),
        rows: rows
            .iter()
            .map(|r| Constraint::Linear(Bounds::fixed(r.iter().sum())))
            .collect(),
        variable_bounds: vec![
            Bounds {
                lower: -10.,
                upper: 10.
            };
            n
        ],
        cones: vec![],
    };
    let mut options = settings();
    options.dependencies.max_basis_rows = rank;
    let result = Presolver::new(options.clone())
        .unwrap()
        .presolve(input.clone());
    assert_eq!(result.stats.after.unwrap().linear_rows, rank);
    let Outcome::Reduced(reduced) = result.outcome else {
        panic!("expected reduction")
    };
    let point = Solution {
        x: vec![1.; n],
        y,
        z: vec![0.; n],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let warm = reduced.postsolve.reduce_warm_start(point.as_ref());
    stationarity(&input, &reduced.postsolve.recover_solution(warm.as_ref()));
    stationarity(&reduced.problem.into_csc(), &warm);
    let mut inconsistent = input;
    let Constraint::Linear(ref mut b) = inconsistent.rows[m - 1] else {
        unreachable!()
    };
    b.lower += 1.;
    b.upper += 1.;
    let result = Presolver::new(options)
        .unwrap()
        .presolve(inconsistent.clone());
    let Outcome::Infeasible(certificate) = result.outcome else {
        panic!("expected contradiction")
    };
    for j in 0..n {
        let residual: f64 = (inconsistent.a.column_pointers()[j]
            ..inconsistent.a.column_pointers()[j + 1])
            .map(|k| inconsistent.a.values()[k] * certificate.y[inconsistent.a.row_indices()[k]])
            .sum();
        assert_eq!(residual, 0.);
    }
}
