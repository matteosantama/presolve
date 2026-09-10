use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, ProblemData},
    settings::{EqualitySettings, Progress, Rules, WorkLimit},
};
use std::time::Duration;

fn settings() -> Settings {
    Settings {
        substitution_fill: usize::MAX,
        allow_hessian_growth: true,
        equalities: EqualitySettings {
            max_row_length: usize::MAX,
            min_column_length: 1,
            max_column_length: usize::MAX,
            require_free_variable: false,
            require_linear_variable: false,
            preserve_nonzeros: false,
            work_limit: WorkLimit::Unlimited,
            ..EqualitySettings::default()
        },
        progress: Progress::AnyChange,
        rules: Rules {
            short_equalities: true,
            ..Rules::none()
        },
        ..Settings::default()
    }
}

fn equation(n: usize, bounded: bool) -> ProblemData {
    ProblemData {
        p: Some(
            CscMatrix::from_triplets(n, n, (0..n).collect(), (0..n).collect(), vec![1.; n])
                .unwrap(),
        ),
        // x = 1, y = 1/4, z = 0 satisfies stationarity and feasibility.
        c: vec![-0.75; n],
        objective_constant: 3.,
        a: CscMatrix::from_triplets(1, n, vec![0; n], (0..n).collect(), vec![1.; n]).unwrap(),
        rows: vec![Constraint::Linear(Bounds::fixed(n as f64))],
        variable_bounds: vec![
            if bounded {
                Bounds {
                    lower: 0.,
                    upper: 2.,
                }
            } else {
                Bounds::FREE
            };
            n
        ],
        cones: vec![],
    }
}

fn column(a: &CscMatrix, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
    (a.column_pointers()[j]..a.column_pointers()[j + 1])
        .map(|k| (a.row_indices()[k], a.values()[k]))
}

fn check_kkt(p: &ProblemData, s: &Solution) -> f64 {
    let n = p.c.len();
    let mut gradient = p.c.clone();
    let mut value = p.objective_constant + p.c.iter().zip(&s.x).map(|(c, x)| c * x).sum::<f64>();
    if let Some(h) = &p.p {
        for j in 0..n {
            for (i, a) in column(h, j) {
                gradient[i] += a * s.x[j];
                if i != j {
                    gradient[j] += a * s.x[i];
                }
                value += if i == j { 0.5 } else { 1. } * a * s.x[i] * s.x[j];
            }
        }
    }
    let mut activity = vec![0.; p.rows.len()];
    for (j, gradient) in gradient.iter_mut().enumerate() {
        for (i, a) in column(&p.a, j) {
            activity[i] += a * s.x[j];
            *gradient -= a * s.y[i];
        }
        assert!((*gradient - s.z[j]).abs() < 1e-8);
        let b = p.variable_bounds.get(j).copied().unwrap_or(Bounds::FREE);
        assert!(s.x[j] >= b.lower - 1e-8 && s.x[j] <= b.upper + 1e-8);
        if s.z[j] != 0. {
            let side = if s.z[j] > 0. { b.lower } else { b.upper };
            assert!((s.z[j] * (s.x[j] - side)).abs() < 1e-8);
        }
    }
    for (i, row) in p.rows.iter().enumerate() {
        let Constraint::Linear(b) = row else {
            unreachable!()
        };
        assert!(activity[i] >= b.lower - 1e-8 && activity[i] <= b.upper + 1e-8);
        if s.y[i] != 0. {
            let side = if s.y[i] > 0. { b.lower } else { b.upper };
            assert!((s.y[i] * (activity[i] - side)).abs() < 1e-8);
        }
    }
    value
}

#[test]
fn wide_quadratic_substitution_preserves_objective_kkt_and_retained_bounds() {
    for (bounded, sign) in [(false, 1.), (true, 1.), (false, -1.), (true, -1.)] {
        let mut p = equation(12, bounded);
        for value in p.a.values_mut() {
            *value *= sign;
        }
        p.rows[0] = Constraint::Linear(Bounds::fixed(sign * 12.));
        let result = Presolver::new(settings()).unwrap().presolve(p.clone());
        let size = result.stats.after.unwrap();
        assert_eq!(size.variables, 11);
        assert_eq!(size.linear_rows, usize::from(bounded));
        assert_eq!(size.p_nonzeros, 121);
        let Outcome::Reduced(r) = result.outcome else {
            panic!("expected reduction")
        };
        let original = Solution {
            x: vec![1.; 12],
            y: vec![sign * 0.25],
            z: vec![0.; 12],
            conic_dual: vec![],
            conic_slack: vec![],
        };
        let warm = r.postsolve.reduce_warm_start(original.as_ref());
        // The known reduced optimum has no active retained bound row.
        let reduced = Solution {
            y: vec![0.; warm.y.len()],
            z: vec![0.; warm.z.len()],
            ..warm
        };
        let recovered = r.postsolve.recover_solution(reduced.as_ref());
        let output = r.problem.into_csc();
        let objective = check_kkt(&p, &original);
        assert!((check_kkt(&output, &reduced) - objective).abs() < 1e-8);
        assert!((check_kkt(&p, &recovered) - objective).abs() < 1e-8);
    }
}

#[test]
fn structural_and_fill_policies_independently_block_the_same_substitution() {
    let p = equation(12, true);
    let mut cases = Vec::new();
    let mut s = settings();
    s.allow_hessian_growth = false;
    cases.push(s);
    let mut s = settings();
    s.substitution_fill = 64;
    cases.push(s);
    let mut s = settings();
    s.equalities.max_row_length = 8;
    cases.push(s);
    let mut s = settings();
    s.equalities.min_column_length = 2;
    cases.push(s);
    let mut s = settings();
    s.equalities.max_column_length = 0;
    cases.push(s);
    let mut s = settings();
    s.equalities.require_free_variable = true;
    cases.push(s);
    let mut s = settings();
    s.equalities.require_linear_variable = true;
    cases.push(s);
    let mut s = settings();
    s.equalities.preserve_nonzeros = true;
    cases.push(s);
    let mut s = settings();
    s.equalities.work_limit = WorkLimit::Entries(0);
    cases.push(s);
    for s in cases {
        let r = Presolver::new(s).unwrap().presolve(p.clone());
        assert!(matches!(r.outcome, Outcome::Unchanged(_)));
        assert_eq!(r.stats.before, r.stats.after.unwrap());
    }
}

#[test]
fn zero_time_budget_leaves_a_valid_unchanged_problem() {
    let mut s = settings();
    s.time_limit = Duration::ZERO;
    let r = Presolver::new(s).unwrap().presolve(equation(12, true));
    assert!(r.stats.time_limit_reached);
    assert!(matches!(r.outcome, Outcome::Unchanged(_)));
    assert_eq!(r.stats.before, r.stats.after.unwrap());
}

#[test]
fn propagation_rounds_work_and_progress_control_bound_only_chains() {
    let n = 12;
    let mut p = equation(n, false);
    p.p = None;
    p.c.fill(0.);
    p.a = CscMatrix::from_triplets(
        n - 1,
        n,
        (0..n - 1).flat_map(|i| [i, i]).collect(),
        (0..n - 1).flat_map(|i| [i, i + 1]).collect(),
        (0..n - 1).flat_map(|_| [1., -1.]).collect(),
    )
    .unwrap();
    p.rows = vec![
        Constraint::Linear(Bounds {
            lower: 0.,
            upper: f64::INFINITY
        });
        n - 1
    ];
    p.variable_bounds = vec![
        Bounds {
            lower: 0.,
            upper: f64::INFINITY
        };
        n
    ];
    p.variable_bounds[n - 1].lower = 1.;
    let mut base = Settings {
        rules: Rules {
            bound_propagation: true,
            ..Rules::none()
        },
        ..Settings::default()
    };
    base.propagation.additional_rounds = 0;
    let lower = |s| {
        let r = Presolver::new(s).unwrap().presolve(p.clone());
        let Outcome::Reduced(r) = r.outcome else {
            panic!("expected bound changes")
        };
        r.problem.variable_bounds(0).lower
    };
    assert_eq!(lower(base.clone()), 0.);
    let mut exhaustive = base.clone();
    exhaustive.progress = Progress::AnyChange;
    assert_eq!(lower(exhaustive), 1.);
    base.propagation.additional_rounds = usize::MAX;
    base.propagation.work_limit = WorkLimit::Unlimited;
    assert_eq!(lower(base.clone()), 1.);
    base.propagation.work_limit = WorkLimit::Entries(0);
    assert_eq!(lower(base), 0.);
}

#[test]
fn sparsification_can_exclude_auxiliary_variables_and_limit_work() {
    let n = 21;
    let mut p = equation(n, false);
    p.p = None;
    p.c.fill(0.);
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut values = Vec::new();
    for i in 0..3 {
        for j in 0..20 {
            rows.push(i);
            cols.push(j);
            values.push(1.);
        }
        if i != 0 {
            rows.push(i);
            cols.push(20);
            values.push(i as f64);
        }
    }
    p.a = CscMatrix::from_triplets(3, n, rows, cols, values).unwrap();
    p.rows = vec![
        Constraint::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: 20.
        });
        3
    ];
    let base = Settings {
        rules: Rules {
            sparsification: true,
            ..Rules::none()
        },
        ..Settings::default()
    };
    let run = |s| Presolver::new(s).unwrap().presolve(p.clone());
    assert_eq!(run(base.clone()).stats.after.unwrap().variables, n + 1);
    let mut equality_only = base.clone();
    equality_only.sparsification.allow_auxiliary_variables = false;
    assert!(matches!(run(equality_only).outcome, Outcome::Unchanged(_)));
    let mut no_work = base;
    no_work.sparsification.work_limit = WorkLimit::Entries(0);
    assert!(matches!(run(no_work).outcome, Outcome::Unchanged(_)));
}

#[test]
fn near_dependent_equalities_do_not_create_roundoff_pivots() {
    let mut p = equation(3, false);
    p.p = None;
    p.c.fill(0.);
    p.a = CscMatrix::from_triplets(
        2,
        3,
        vec![0, 0, 0, 1, 1, 1],
        vec![0, 1, 2, 0, 1, 2],
        vec![1., 1., 1., 1., 1. + 1e-12, 1. + 2e-12],
    )
    .unwrap();
    p.rows = vec![
        Constraint::Linear(Bounds::fixed(3.)),
        Constraint::Linear(Bounds::fixed(3. + 3e-12)),
    ];
    let result = Presolver::new(settings()).unwrap().presolve(p.clone());
    let Outcome::Unchanged(output) = result.outcome else {
        panic!("ill-conditioned substitutions should be rejected")
    };
    let output = output.into_csc();
    assert_eq!(output.a, p.a);
    assert_eq!(output.rows, p.rows);
}

#[test]
fn exhaustive_cycles_revisit_equalities_after_pivot_degrees_change() {
    let mut p = equation(7, false);
    p.p = None;
    p.c.fill(0.);
    // Row 0's only largest pivot starts at degree 3, outside the configured cap.
    // Eliminating row 1's free singleton makes that pivot eligible at degree 2.
    p.a = CscMatrix::from_triplets(
        3,
        7,
        vec![0, 0, 0, 1, 1, 1, 2, 2, 2],
        vec![0, 1, 2, 0, 3, 4, 0, 5, 6],
        vec![2., 1., 1., 1., 1., 2., 1., 1., 1.],
    )
    .unwrap();
    p.rows = vec![Constraint::Linear(Bounds::fixed(0.)); 3];
    let mut s = settings();
    s.equalities.max_column_length = 2;
    let r = Presolver::new(s).unwrap().presolve(p);
    assert_eq!(r.stats.after.unwrap().linear_rows, 0);
    assert_eq!(r.stats.after.unwrap().variables, 4);
}

#[test]
fn alternative_pivots_recover_after_a_fill_rejection() {
    let input = ProblemData {
        p: None,
        c: vec![0.; 4],
        objective_constant: 0.,
        a: CscMatrix::from_triplets(
            3,
            4,
            vec![0, 0, 0, 1, 1, 2, 2, 2],
            vec![0, 1, 2, 0, 3, 1, 2, 3],
            vec![1.; 8],
        )
        .unwrap(),
        rows: vec![],
        variable_bounds: vec![Bounds::FREE; 4],
        cones: vec![],
    };
    let mut input = input;
    input.rows = vec![
        Constraint::Linear(Bounds::fixed(0.)),
        Constraint::Linear(Bounds {
            lower: -1.,
            upper: 1.,
        }),
        Constraint::Linear(Bounds {
            lower: -1.,
            upper: 1.,
        }),
    ];
    let mut options = settings();
    options.substitution_fill = 1;
    options.equalities.max_pivot_attempts = 1;
    let first = Presolver::new(options.clone())
        .unwrap()
        .presolve(input.clone());
    assert_eq!(first.stats.after.unwrap().variables, 4);
    assert_eq!(first.stats.equalities.rejected_updates, 1);
    assert_eq!(first.stats.equalities.constraint_fill_rejections, 1);
    options.equalities.max_pivot_attempts = 2;
    let next = Presolver::new(options).unwrap().presolve(input.clone());
    assert_eq!(next.stats.after.unwrap().variables, 3);
    assert_eq!(next.stats.equalities.attempts, 2);
    assert_eq!(next.stats.equalities.accepted, 1);
    let Outcome::Reduced(reduced) = next.outcome else {
        panic!("expected reduction")
    };
    let point = Solution {
        x: vec![0.; 4],
        y: vec![0.; 3],
        z: vec![0.; 4],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let warm = reduced.postsolve.reduce_warm_start(point.as_ref());
    check_kkt(&input, &reduced.postsolve.recover_solution(warm.as_ref()));
}

#[test]
fn candidate_scoring_counts_quadratic_work_and_pivot_screening_is_configurable() {
    let mut input = equation(3, false);
    input.p = Some(CscMatrix::from_triplets(3, 3, vec![0], vec![0], vec![1.]).unwrap());
    let mut options = settings();
    options.allow_hessian_growth = false;
    options.equalities.max_pivot_attempts = 1;
    options.equalities.cost_aware = false;
    assert_eq!(
        Presolver::new(options.clone())
            .unwrap()
            .presolve(input.clone())
            .stats
            .after
            .unwrap()
            .variables,
        3
    );
    options.equalities.cost_aware = true;
    let result = Presolver::new(options.clone())
        .unwrap()
        .presolve(input.clone());
    assert_eq!(result.stats.after.unwrap().variables, 2);
    assert_eq!(result.stats.after.unwrap().p_nonzeros, 1);
    input.a.values_mut()[0] = 2.;
    options.equalities.require_linear_variable = true;
    for relative in [1., 0., -1., f64::NAN, 2.] {
        options.equalities.relative_pivot = relative;
        assert_eq!(
            Presolver::new(options.clone())
                .unwrap()
                .presolve(input.clone())
                .stats
                .after
                .unwrap()
                .variables,
            3
        );
    }
    options.equalities.relative_pivot = 0.5;
    assert_eq!(
        Presolver::new(options)
            .unwrap()
            .presolve(input)
            .stats
            .after
            .unwrap()
            .variables,
        2
    );
}

#[test]
fn pivot_policies_preserve_known_quadratic_optima() {
    for seed in 0..24 {
        let (m, n) = (4, 8);
        let mut state = seed + 1u64;
        let mut values = Vec::new();
        let mut c = vec![-1.; n];
        let mut rhs = vec![0.; m];
        let y = [0.25, -0.5, 0.75, -1.];
        for coefficient in &mut c {
            for i in 0..m {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let value = ((state >> 32) % 9) as f64 - 4.;
                values.push(value);
                *coefficient += y[i] * value;
                rhs[i] += value;
            }
        }
        let input = ProblemData {
            p: Some(
                CscMatrix::from_triplets(n, n, (0..n).collect(), (0..n).collect(), vec![1.; n])
                    .unwrap(),
            ),
            c,
            objective_constant: 3.,
            a: CscMatrix::from_triplets(
                m,
                n,
                (0..n).flat_map(|_| 0..m).collect(),
                (0..n).flat_map(|j| std::iter::repeat_n(j, m)).collect(),
                values,
            )
            .unwrap(),
            rows: rhs
                .into_iter()
                .map(|v| Constraint::Linear(Bounds::fixed(v)))
                .collect(),
            variable_bounds: vec![Bounds::FREE; n],
            cones: vec![],
        };
        let optimum = Solution {
            x: vec![1.; n],
            y: y.to_vec(),
            z: vec![0.; n],
            conic_dual: vec![],
            conic_slack: vec![],
        };
        let objective = check_kkt(&input, &optimum);
        for relative in [1., 0.5, 0.1] {
            for aware in [false, true] {
                let mut options = settings();
                options.equalities.relative_pivot = relative;
                options.equalities.cost_aware = aware;
                options.equalities.max_pivot_attempts = 4;
                let result = Presolver::new(options).unwrap().presolve(input.clone());
                match result.outcome {
                    Outcome::Reduced(reduced) => {
                        let warm = reduced.postsolve.reduce_warm_start(optimum.as_ref());
                        let recovered = reduced.postsolve.recover_solution(warm.as_ref());
                        assert!((check_kkt(&input, &recovered) - objective).abs() < 1e-8);
                        assert!(
                            (check_kkt(&reduced.problem.into_csc(), &warm) - objective).abs()
                                < 1e-8
                        );
                    }
                    Outcome::Unchanged(_) => {}
                    outcome => panic!("unexpected {outcome:?}"),
                }
            }
        }
    }
}

#[test]
fn propagation_gain_controls_are_independent_and_preserve_dual_recovery() {
    for (rhs, old_upper) in [(1000., 1005.), (1., 1. + 1e-6)] {
        let input = ProblemData {
            p: None,
            c: vec![-1., 0.],
            objective_constant: 0.,
            a: CscMatrix::from_triplets(1, 2, vec![0, 0], vec![0, 1], vec![1., 1.]).unwrap(),
            rows: vec![Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: rhs,
            })],
            variable_bounds: vec![
                Bounds {
                    lower: 0.,
                    upper: old_upper
                };
                2
            ],
            cones: vec![],
        };
        let mut options = Settings {
            rules: Rules {
                bound_propagation: true,
                ..Rules::none()
            },
            ..Settings::default()
        };
        assert!(matches!(
            Presolver::new(options.clone())
                .unwrap()
                .presolve(input.clone())
                .outcome,
            Outcome::Unchanged(_)
        ));
        options.propagation.minimum_relative_gain = 0.;
        if rhs == 1. {
            assert!(matches!(
                Presolver::new(options.clone())
                    .unwrap()
                    .presolve(input.clone())
                    .outcome,
                Outcome::Unchanged(_)
            ));
        }
        options.propagation.minimum_gain_factor = 1.;
        let result = Presolver::new(options).unwrap().presolve(input.clone());
        let Outcome::Reduced(reduced) = result.outcome else {
            panic!("expected tighter bounds")
        };
        let point = Solution {
            x: vec![rhs, 0.],
            y: vec![0.],
            z: vec![-1., 0.],
            conic_dual: vec![],
            conic_slack: vec![],
        };
        let recovered = reduced.postsolve.recover_solution(point.as_ref());
        assert_eq!(check_kkt(&input, &recovered), -rhs);
        let output = reduced.problem.into_csc();
        assert!(output.variable_bounds.iter().all(|b| b.upper == rhs));
        assert_eq!(check_kkt(&output, &point), -rhs);
    }
}
