use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, ProblemData},
    result::PresolveResult,
    settings::Rules,
};

const WIDTH: usize = 8;
// 32,768 constraint entries exercise the parallel path, including row scans.
const BLOCKS: usize = 2048;

fn problem(blocks: usize, quadratic: bool) -> ProblemData {
    let n = WIDTH * blocks;
    let mut ai = Vec::new();
    let mut aj = Vec::new();
    let mut av = Vec::new();
    let mut pi = Vec::new();
    let mut pj = Vec::new();
    for block in 0..blocks {
        for j in 0..WIDTH {
            // Negative proportional rows test sign handling and tie ordering.
            ai.extend([2 * block, 2 * block + 1]);
            aj.extend([WIDTH * block + j; 2]);
            av.extend([1., -2.]);
            for i in 0..=j {
                pi.push(WIDTH * block + i);
                pj.push(WIDTH * block + j);
            }
        }
    }
    let pv = vec![1.; pi.len()];
    ProblemData {
        p: quadratic.then(|| CscMatrix::from_triplets(n, n, pi, pj, pv).unwrap()),
        c: vec![-2.; n],
        objective_constant: 2. * blocks as f64,
        a: CscMatrix::from_triplets(2 * blocks, n, ai, aj, av).unwrap(),
        rows: (0..blocks)
            .flat_map(|_| {
                [
                    Constraint::Linear(Bounds {
                        lower: 1.,
                        upper: 4.,
                    }),
                    Constraint::Linear(Bounds {
                        lower: -6.,
                        upper: -2.,
                    }),
                ]
            })
            .collect(),
        variable_bounds: vec![
            Bounds {
                lower: 0.,
                upper: 1.
            };
            n
        ],
        cones: vec![],
    }
}

fn run(input: &ProblemData, rules: Rules, threads: usize) -> PresolveResult {
    presolve::presolve(
        input.clone(),
        &Settings {
            rules,
            threads,
            ..Settings::default()
        },
    )
    .unwrap()
}

fn same_solution(a: &Solution, b: &Solution) {
    assert_eq!(a.x, b.x);
    assert_eq!(a.y, b.y);
    assert_eq!(a.z, b.z);
    assert_eq!(a.conic_dual, b.conic_dual);
    assert_eq!(a.conic_slack, b.conic_slack);
}

fn same_result(a: PresolveResult, b: PresolveResult, input: &ProblemData) {
    assert!(!a.stats.time_limit_reached && !b.stats.time_limit_reached);
    assert_eq!(a.stats.before, b.stats.before);
    assert_eq!(a.stats.after, b.stats.after);
    assert_eq!(a.stats.parallel_comparisons, b.stats.parallel_comparisons);
    assert_eq!(a.stats.quadratic_changed, b.stats.quadratic_changed);
    match (a.outcome, b.outcome) {
        (Outcome::Reduced(a), Outcome::Reduced(b)) => {
            let point = Solution {
                x: vec![0.25; input.c.len()],
                y: vec![0.; input.rows.len()],
                z: vec![0.; input.c.len()],
                conic_dual: vec![],
                conic_slack: vec![],
            };
            let ar = a.postsolve.reduce_warm_start(point.as_ref());
            let br = b.postsolve.reduce_warm_start(point.as_ref());
            same_solution(&ar, &br);
            let recovered = a.postsolve.recover_solution(ar.as_ref());
            same_solution(&recovered, &b.postsolve.recover_solution(br.as_ref()));
            check_optimum(&recovered);
            let a = a.problem.into_csc();
            let b = b.problem.into_csc();
            assert_eq!(a.a, b.a);
            assert_eq!(a.p, b.p);
            assert_eq!(a.c, b.c);
            assert_eq!(a.objective_constant, b.objective_constant);
            assert_eq!(a.rows, b.rows);
            assert_eq!(a.variable_bounds, b.variable_bounds);
            assert_eq!(a.cones, b.cones);
        }
        (Outcome::Solved(a), Outcome::Solved(b)) => {
            same_solution(&a, &b);
            check_optimum(&a);
        }
        (Outcome::Infeasible(a), Outcome::Infeasible(b)) => {
            assert_eq!(a.y, b.y);
            assert_eq!(a.z, b.z);
            assert_eq!(a.conic_dual, b.conic_dual);
            // Original-coordinate Farkas stationarity and positive contradiction.
            for j in 0..input.c.len() {
                let value: f64 = column(&input.a, j).map(|(i, v)| v * a.y[i]).sum();
                assert_eq!(value + a.z[j], 0.);
            }
            let contradiction: f64 = input
                .rows
                .iter()
                .zip(&a.y)
                .map(|(row, y)| {
                    let Constraint::Linear(bounds) = row else {
                        unreachable!()
                    };
                    y * if *y >= 0. { bounds.lower } else { bounds.upper }
                })
                .sum();
            assert!(contradiction > 0.);
        }
        (Outcome::Unbounded(a), Outcome::Unbounded(b)) => {
            assert_eq!(a.point, b.point);
            assert_eq!(a.ray, b.ray);
            assert!(input.c.iter().zip(&a.ray).map(|(c, x)| c * x).sum::<f64>() < 0.);
            let mut activity = vec![0.; input.rows.len()];
            for (j, x) in a.ray.iter().enumerate() {
                for (i, v) in column(&input.a, j) {
                    activity[i] += v * x;
                }
            }
            assert!(activity.iter().all(|&x| x == 0.));
        }
        (a, b) => panic!("unexpected outcomes: {a:?}, {b:?}"),
    }
}

fn check_optimum(point: &Solution) {
    // Each independent block minimizes 0.5 * (sum(x) - 2)^2.
    assert!(point.x.iter().all(|x| (0.0..=1.0).contains(x)));
    for block in point.x.chunks_exact(WIDTH) {
        assert_eq!(block.iter().sum::<f64>(), 2.);
    }
    assert!(point.y.iter().chain(&point.z).all(|&v| v == 0.));
}

#[test]
fn thread_counts_preserve_reductions_and_postsolve() {
    assert_eq!(Settings::default().threads, 1);
    for blocks in [1, BLOCKS] {
        let input = problem(blocks, true);
        for rules in [
            Rules {
                parallel_rows: true,
                ..Rules::none()
            },
            Rules {
                parallel_columns: true,
                ..Rules::none()
            },
            Rules {
                parallel_rows: true,
                parallel_columns: true,
                ..Rules::none()
            },
            Rules::default(),
        ] {
            for threads in [0, 1, 2, 4] {
                let parallel = run(&input, rules, threads);
                same_result(run(&input, rules, 1), parallel, &input);
            }
        }
    }
}

#[test]
fn parallel_discovery_preserves_certificates() {
    let mut infeasible = problem(BLOCKS, false);
    // sum(x) >= 1 and -2 sum(x) >= 2 are incompatible.
    infeasible.rows[1] = Constraint::Linear(Bounds {
        lower: 2.,
        upper: 4.,
    });
    let rows = Rules {
        parallel_rows: true,
        ..Rules::none()
    };
    let serial = run(&infeasible, rows, 1);
    assert!(matches!(serial.outcome, Outcome::Infeasible(_)));
    same_result(serial, run(&infeasible, rows, 4), &infeasible);

    let mut unbounded = problem(BLOCKS, false);
    unbounded.variable_bounds.fill(Bounds::FREE);
    unbounded.rows.fill(Constraint::Linear(Bounds {
        lower: -1.,
        upper: 1.,
    }));
    unbounded.c = (0..unbounded.c.len())
        .map(|j| (j % WIDTH != 0) as u8 as f64)
        .collect();
    let columns = Rules {
        parallel_columns: true,
        ..Rules::none()
    };
    let serial = run(&unbounded, columns, 1);
    assert!(matches!(serial.outcome, Outcome::Unbounded(_)));
    same_result(serial, run(&unbounded, columns, 4), &unbounded);
}

#[test]
fn reusable_presolver_supports_concurrent_calls_and_independent_results() {
    let rules = Rules {
        parallel_rows: true,
        parallel_columns: true,
        ..Rules::none()
    };
    let presolver = Presolver::new(Settings {
        rules,
        threads: 4,
        ..Settings::default()
    })
    .unwrap();
    assert_eq!(presolver.threads(), 4);
    assert_eq!(presolver.settings().threads, 4);
    let large = problem(BLOCKS, true);
    let small = problem(2, true);
    let (first, second) = std::thread::scope(|scope| {
        let a = scope.spawn(|| presolver.presolve(large.clone()));
        let b = scope.spawn(|| presolver.presolve(small.clone()));
        (a.join().unwrap(), b.join().unwrap())
    });
    let third = presolver.presolve(large.clone());
    drop(presolver);
    // Postsolve and reduced storage outlive the executor and all later calls.
    same_result(run(&large, rules, 1), first, &large);
    same_result(run(&small, rules, 1), second, &small);
    same_result(run(&large, rules, 1), third, &large);
}

fn column(matrix: &CscMatrix, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
    let start = matrix.column_pointers()[j];
    let end = matrix.column_pointers()[j + 1];
    matrix.row_indices()[start..end]
        .iter()
        .copied()
        .zip(matrix.values()[start..end].iter().copied())
}

fn collided_rows(classes: usize, width: usize, factor: f64) -> ProblemData {
    let m = 1 + 2 * classes;
    let mut ai = Vec::with_capacity(m * width);
    let mut aj = Vec::with_capacity(m * width);
    let mut av = Vec::with_capacity(m * width);
    for i in 0..m {
        let class = i.div_ceil(2);
        let scale = if i > 0 && i % 2 == 0 { factor } else { 1.0 };
        for j in 0..width {
            ai.push(i);
            aj.push(j);
            av.push(
                scale
                    * if j + 1 == width {
                        1.0 + class as f64 * 2e-10
                    } else {
                        1.0
                    },
            );
        }
    }
    ProblemData {
        p: None,
        c: vec![0.; width],
        objective_constant: 0.,
        a: CscMatrix::from_triplets(m, width, ai, aj, av).unwrap(),
        rows: vec![
            Constraint::Linear(Bounds {
                lower: -1.,
                upper: 1.
            });
            m
        ],
        variable_bounds: vec![Bounds::FREE; width],
        cones: vec![],
    }
}

#[test]
fn hash_collisions_preserve_hidden_classes_with_bounded_comparisons() {
    let rules = Rules {
        parallel_rows: true,
        ..Rules::none()
    };
    for (classes, width) in [(1, 2), (1024, 16)] {
        for factor in [2., -2., 1e-250, -1e250] {
            let input = collided_rows(classes, width, factor);
            let mut reference = None;
            for threads in [1, 4] {
                let result = run(&input, rules, threads);
                assert_eq!(result.stats.after.unwrap().linear_rows, classes + 1);
                // The scheduler repeats discovery after the first reduction.
                assert!(result.stats.parallel_comparisons <= 4 * input.rows.len());
                let Outcome::Reduced(reduced) = result.outcome else {
                    panic!("expected reduction")
                };
                let original = Solution {
                    x: vec![0.; width],
                    y: vec![0.; input.rows.len()],
                    z: vec![0.; width],
                    conic_dual: vec![],
                    conic_slack: vec![],
                };
                let warm = reduced.postsolve.reduce_warm_start(original.as_ref());
                same_solution(
                    &original,
                    &reduced.postsolve.recover_solution(warm.as_ref()),
                );
                let data = reduced.problem.into_csc();
                if let Some((a, rows)) = reference.as_ref() {
                    assert_eq!(&data.a, a);
                    assert_eq!(&data.rows, rows);
                } else {
                    reference = Some((data.a, data.rows));
                }
            }
        }
    }
}

#[test]
fn hash_collision_does_not_hide_an_infeasibility_certificate() {
    let mut input = collided_rows(1, 2, -2.);
    input.rows[1] = Constraint::Linear(Bounds {
        lower: 1.,
        upper: 2.,
    });
    input.rows[2] = Constraint::Linear(Bounds {
        lower: 1.,
        upper: 2.,
    });
    let rules = Rules {
        parallel_rows: true,
        ..Rules::none()
    };
    let result = run(&input, rules, 1);
    assert!(matches!(result.outcome, Outcome::Infeasible(_)));
    same_result(result, run(&input, rules, 4), &input);
}

#[test]
fn redundant_parallel_row_keeps_warm_start_stationarity_without_tightening() {
    for scale in [1., -2.] {
        let input = ProblemData {
            p: None,
            c: vec![2.; 2],
            objective_constant: 0.,
            a: CscMatrix::from_triplets(
                2,
                2,
                vec![0, 0, 1, 1],
                vec![0, 1, 0, 1],
                vec![1., 1., scale, scale],
            )
            .unwrap(),
            rows: vec![
                Constraint::Linear(Bounds {
                    lower: 1.,
                    upper: f64::INFINITY,
                }),
                Constraint::Linear(if scale > 0. {
                    Bounds {
                        lower: scale,
                        upper: f64::INFINITY,
                    }
                } else {
                    Bounds {
                        lower: f64::NEG_INFINITY,
                        upper: scale,
                    }
                }),
            ],
            variable_bounds: vec![Bounds::FREE; 2],
            cones: vec![],
        };
        let result = run(
            &input,
            Rules {
                parallel_rows: true,
                ..Rules::none()
            },
            1,
        );
        let Outcome::Reduced(r) = result.outcome else {
            panic!("expected merged row");
        };
        let original = Solution {
            x: vec![0.5; 2],
            y: vec![0., 2. / scale],
            z: vec![0.; 2],
            conic_dual: vec![],
            conic_slack: vec![],
        };
        let warm = r.postsolve.reduce_warm_start(original.as_ref());
        assert_eq!(warm.y, vec![2.]);
        let recovered = r.postsolve.recover_solution(warm.as_ref());
        assert_eq!(recovered.y, vec![2., 0.]);
        assert_eq!(recovered.x, original.x);
        assert_eq!(recovered.z, original.z);
    }
}
