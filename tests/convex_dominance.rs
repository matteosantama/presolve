use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, Problem},
    settings::Rules,
};

fn settings(enabled: bool) -> Settings {
    Settings {
        rules: Rules {
            convex_dominance: enabled,
            ..Rules::none()
        },
        ..Settings::default()
    }
}
fn problem() -> Problem {
    Problem {
        p: None,
        c: vec![1., 3., 3.],
        c0: 0.,
        a: CscMatrix::from_triplets(
            2,
            3,
            vec![0, 1, 0, 1, 0, 1],
            vec![0, 0, 1, 1, 2, 2],
            vec![1., 1., 1., 2., 1., 3.],
        )
        .unwrap(),
        rows: vec![
            Constraint::Linear(Bounds::fixed(1.)),
            Constraint::Linear(Bounds::fixed(2.)),
        ],
        variable_bounds: vec![
            Bounds {
                lower: 0.,
                upper: f64::INFINITY
            };
            3
        ],
        cones: vec![],
    }
}
#[test]
fn mixture_dominance_preserves_primal_dual_optimum() {
    let result = Presolver::new(settings(true)).unwrap().presolve(problem());
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected reduction")
    };
    assert_eq!(r.problem.c, [1., 3.]);
    let solution = Solution {
        x: vec![0.5, 0.5],
        y: vec![0., 1.],
        z: vec![0., 0.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let restored = r.postsolve.recover_solution(solution.as_ref());
    assert_eq!(restored.x, [0.5, 0., 0.5]);
    assert_eq!(restored.y, [0., 1.]);
    assert_eq!(restored.z, [0., 1., 0.]);
    let warm = r.postsolve.reduce_warm_start(restored.as_ref());
    assert_eq!(warm.x, solution.x);
    assert_eq!(warm.y, solution.y);
    assert_eq!(warm.z, solution.z);
    assert!(matches!(
        Presolver::new(settings(false))
            .unwrap()
            .presolve(problem())
            .outcome,
        Outcome::Unchanged(_)
    ));
}
#[test]
fn uncertain_or_invalid_mixtures_are_retained() {
    for kind in 0..6 {
        let mut p = problem();
        match kind {
            0 => p.c[1] = 2., // collinear: numerical proof conservatively declines
            1 => p.variable_bounds[0].upper = 0.1, // witness lacks capacity
            2 => p.variable_bounds[1].lower = -1.,
            3 => p.p = Some(CscMatrix::from_triplets(3, 3, vec![1], vec![1], vec![1.]).unwrap()),
            4 => {
                p.a = CscMatrix::from_triplets(
                    2,
                    3,
                    vec![0, 1, 0, 1, 0, 1],
                    vec![0, 0, 1, 1, 2, 2],
                    vec![1., 1., 1.0_f64.next_up(), 2., 1., 3.],
                )
                .unwrap()
            }
            _ => p.c[1] = 1., // below the chord
        }
        assert!(
            matches!(
                Presolver::new(settings(true)).unwrap().presolve(p).outcome,
                Outcome::Unchanged(_)
            ),
            "case {kind}"
        );
    }
}
#[test]
fn equal_abscissae_keep_the_cheapest_column() {
    let mut p = problem();
    p.a = CscMatrix::from_triplets(
        2,
        3,
        vec![0, 1, 0, 1, 0, 1],
        vec![0, 0, 1, 1, 2, 2],
        vec![1., 1., 1., 1., 1., 3.],
    )
    .unwrap();
    let Outcome::Reduced(r) = Presolver::new(settings(true)).unwrap().presolve(p).outcome else {
        panic!("expected duplicate removal")
    };
    assert_eq!(r.problem.c, [1., 3.]);
}

#[test]
fn configurable_gate_accounts_for_quadratic_work_and_is_independent_of_progress() {
    use presolve::settings::Progress;
    let mut input = problem();
    // The same three-column hull sits beside fifty quadratic coordinates.
    // Only two A entries can be removed, while all fifty P entries remain.
    input.c.resize(53, 0.0);
    input.variable_bounds.resize(
        53,
        Bounds {
            lower: 0.0,
            upper: f64::INFINITY,
        },
    );
    input.a = CscMatrix::from_triplets(
        2,
        53,
        vec![0, 1, 0, 1, 0, 1],
        vec![0, 0, 1, 1, 2, 2],
        vec![1., 1., 1., 2., 1., 3.],
    )
    .unwrap();
    input.p = Some(
        CscMatrix::from_triplets(53, 53, (3..53).collect(), (3..53).collect(), vec![1.; 50])
            .unwrap(),
    );
    assert_eq!(Settings::default().convex_dominance.minimum_reduction, 0.05);
    let aggressive = Settings::aggressive(std::time::Duration::from_secs(1));
    assert_eq!(aggressive.convex_dominance.minimum_reduction, 0.0);
    // Repeating exploration after any edit must not bypass the separate gate.
    let mut options = settings(true);
    options.progress = Progress::AnyChange;
    assert!(matches!(
        Presolver::new(options)
            .unwrap()
            .presolve(input.clone())
            .outcome,
        Outcome::Unchanged(_)
    ));
    for (threshold, accepted) in [
        (aggressive.convex_dominance.minimum_reduction, true),
        (0.01, true),
        (0.05, false),
        (f64::NAN, false),
        (-1.0, true),
        (1.0, false),
        (2.0, false),
    ] {
        let mut options = settings(true);
        options.convex_dominance.minimum_reduction = threshold;
        // Disabling exploration repeats must not block an accepted hull batch.
        options.progress = Progress::Nonzeros {
            minimum_reduction: 1.0,
        };
        let result = Presolver::new(options).unwrap().presolve(input.clone());
        match result.outcome {
            Outcome::Reduced(reduced) if accepted => {
                assert_eq!(reduced.problem.c.len(), 52);
                assert_eq!(reduced.problem.p.unwrap().values().len(), 50);
            }
            Outcome::Unchanged(_) if !accepted => {}
            other => panic!("unexpected outcome for threshold {threshold}: {other:?}"),
        }
    }
}
