use presolve::{
    Outcome, Presolver, Settings,
    matrix::CscMatrix,
    postsolve::Solution,
    problem::{Bounds, Constraint, Problem},
    settings::{DominatedColumnSettings, Rules},
};

fn only_dominated_columns() -> Settings {
    Settings {
        rules: Rules {
            dominated_columns: true,
            ..Rules::none()
        },
        dominated_columns: DominatedColumnSettings {
            general_search: true,
            ..DominatedColumnSettings::default()
        },
        ..Settings::default()
    }
}

/// min x + 2y subject to x + y >= 1, y + w <= 5, with the given bounds on
/// x and y and w >= 0. Column x dominates column y on both rows.
fn problem(x: Bounds, y: Bounds) -> Problem {
    Problem {
        p: None,
        c: vec![1., 2., 0.],
        c0: 0.,
        a: CscMatrix::from_triplets(2, 3, vec![0, 0, 1, 1], vec![0, 1, 1, 2], vec![1.; 4]).unwrap(),
        rows: vec![
            Constraint::Linear(Bounds {
                lower: 1.,
                upper: f64::INFINITY,
            }),
            Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 5.,
            }),
        ],
        variable_bounds: vec![
            x,
            y,
            Bounds {
                lower: 0.,
                upper: f64::INFINITY,
            },
        ],
        cones: vec![],
    }
}

const NONNEGATIVE: Bounds = Bounds {
    lower: 0.,
    upper: f64::INFINITY,
};

#[test]
fn dominated_column_is_fixed_at_its_lower_bound() {
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(problem(NONNEGATIVE, NONNEGATIVE));
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    assert_eq!(r.problem.variable_count(), 2);
    assert_eq!(r.problem.row_count(), 2);
    let reduced = Solution {
        x: vec![1., 0.],
        y: vec![1., 0.],
        z: vec![0., 0.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, [1., 0., 0.]);
    // z_y = c_y - y_0 - y_1 = 1 >= 0 without any dual transformation.
    assert_eq!(recovered.z, [0., 1., 0.]);
    assert_eq!(recovered.y, [1., 0.]);
}

#[test]
fn dominating_column_is_fixed_at_its_upper_bound_when_the_other_is_free_below() {
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(problem(
            Bounds {
                lower: 0.,
                upper: 3.,
            },
            Bounds {
                lower: f64::NEG_INFINITY,
                upper: f64::INFINITY,
            },
        ));
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    assert_eq!(r.problem.variable_count(), 2);
    assert_eq!(
        r.problem.rows[0],
        Constraint::Linear(Bounds {
            lower: -2.,
            upper: f64::INFINITY
        })
    );
    assert_eq!(r.problem.c0, 3.);
}

#[test]
fn unbounded_shift_direction_is_reported() {
    // The origin must be feasible for the entry point to confirm the ray.
    let mut unbounded = problem(
        NONNEGATIVE,
        Bounds {
            lower: f64::NEG_INFINITY,
            upper: f64::INFINITY,
        },
    );
    unbounded.rows[0] = Constraint::Linear(Bounds {
        lower: -1.,
        upper: f64::INFINITY,
    });
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(unbounded);
    let Outcome::Unbounded(certificate) = result.outcome else {
        panic!("expected an unbounded certificate")
    };
    assert_eq!(certificate.ray, [1., -1., 0.]);
}

#[test]
fn identical_support_pairs_are_found_inside_the_parallel_column_scan() {
    // min x + 2y subject to x + y >= 1, x + y + w <= 5: x and y share their
    // support without being proportional to w, so the default search finds them.
    let problem = Problem {
        p: None,
        c: vec![1., 2., 0.],
        c0: 0.,
        a: CscMatrix::from_triplets(
            2,
            3,
            vec![0, 1, 0, 1, 1],
            vec![0, 0, 1, 1, 2],
            vec![1., 1., 1., 1., 1.],
        )
        .unwrap(),
        rows: vec![
            Constraint::Linear(Bounds {
                lower: 1.,
                upper: f64::INFINITY,
            }),
            Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 5.,
            }),
        ],
        variable_bounds: vec![NONNEGATIVE; 3],
        cones: vec![],
    };
    let settings = Settings {
        rules: Rules {
            parallel_columns: true,
            dominated_columns: true,
            ..Rules::none()
        },
        ..Settings::default()
    };
    let result = Presolver::new(settings.clone())
        .unwrap()
        .presolve(problem.clone());
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    assert_eq!(r.problem.variable_count(), 2);
    let recovered = r.postsolve.recover_solution(
        Solution {
            x: vec![1., 0.],
            y: vec![1., 0.],
            z: vec![0., 0.],
            conic_dual: vec![],
            conic_slack: vec![],
        }
        .as_ref(),
    );
    assert_eq!(recovered.x, [1., 0., 0.]);
    assert_eq!(recovered.z, [0., 1., 0.]);
    // Without the general search, the nested-support pairs stay untouched.
    let nested = Presolver::new(settings).unwrap().presolve(super_problem());
    assert!(matches!(nested.outcome, Outcome::Unchanged(_)));
}

fn super_problem() -> Problem {
    problem(NONNEGATIVE, NONNEGATIVE)
}

#[test]
fn quadratic_and_boxed_columns_are_left_alone() {
    let mut curved = problem(NONNEGATIVE, NONNEGATIVE);
    curved.p = Some(CscMatrix::from_triplets(3, 3, vec![1], vec![1], vec![1.]).unwrap());
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(curved);
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
    let boxed = problem(
        Bounds {
            lower: 0.,
            upper: 3.,
        },
        Bounds {
            lower: 0.,
            upper: 3.,
        },
    );
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(boxed);
    assert!(matches!(result.outcome, Outcome::Unchanged(_)));
}

#[test]
fn multiplier_on_an_implied_side_moves_to_the_implying_row() {
    // min x + 2y subject to x + y >= 1 with x, y in [0, 1]. The row implies
    // y >= 0, so x is fixed at its upper bound. A reduced solution may put
    // y's whole reduced cost on that implied bound; recovery must move it to
    // the row, or x's reduced cost gets the lower-side sign at its upper bound.
    let unit = Bounds {
        lower: 0.,
        upper: 1.,
    };
    let problem = Problem {
        p: None,
        c: vec![1., 2.],
        c0: 0.,
        a: CscMatrix::from_triplets(1, 2, vec![0, 0], vec![0, 1], vec![1., 1.]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: 1.,
            upper: f64::INFINITY,
        })],
        variable_bounds: vec![unit; 2],
        cones: vec![],
    };
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(problem);
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    assert_eq!(r.problem.variable_count(), 1);
    assert_eq!(r.problem.c0, 1.);
    let reduced = Solution {
        x: vec![0.],
        y: vec![0.],
        z: vec![2.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x, [1., 0.]);
    // z = c - Aᵀy with x at its upper bound needs z_x <= 0.
    assert_eq!(recovered.y, [2.]);
    assert_eq!(recovered.z, [-1., 0.]);
}

#[test]
fn fixes_citing_one_implying_row_recover_a_valid_dual() {
    // min x1 + 2y1 + x2 + 2y2 subject to x1 + y1 + x2 + y2 >= 3 in [0, 1]^4.
    // The row implies every lower bound, so three columns are fixed at their
    // upper bounds on the strength of the same row, whose saved copy still
    // holds the columns fixed before it is cited again.
    let unit = Bounds {
        lower: 0.,
        upper: 1.,
    };
    let c = [1., 2., 1., 2.];
    let problem = Problem {
        p: None,
        c: c.to_vec(),
        c0: 0.,
        a: CscMatrix::from_triplets(1, 4, vec![0; 4], vec![0, 1, 2, 3], vec![1.; 4]).unwrap(),
        rows: vec![Constraint::Linear(Bounds {
            lower: 3.,
            upper: f64::INFINITY,
        })],
        variable_bounds: vec![unit; 4],
        cones: vec![],
    };
    let result = Presolver::new(only_dominated_columns())
        .unwrap()
        .presolve(problem);
    let Outcome::Reduced(r) = result.outcome else {
        panic!("expected a reduced problem")
    };
    assert_eq!(r.problem.variable_count(), 1);
    assert_eq!(r.problem.c, [2.]);
    assert_eq!(r.problem.c0, 4.);
    // The remaining y sits at its bound, which carries its whole reduced cost.
    let reduced = Solution {
        x: vec![0.],
        y: vec![0.],
        z: vec![2.],
        conic_dual: vec![],
        conic_slack: vec![],
    };
    let recovered = r.postsolve.recover_solution(reduced.as_ref());
    assert_eq!(recovered.x.iter().sum::<f64>(), 3.);
    assert_eq!(recovered.y, [2.]);
    for (j, (&x, &z)) in recovered.x.iter().zip(&recovered.z).enumerate() {
        assert_eq!(z, c[j] - recovered.y[0]);
        if x == 1. {
            assert!(z <= 0., "column {j}: {:?}", recovered.z);
        } else {
            assert_eq!(x, 0.);
            assert!(z >= 0., "column {j}: {:?}", recovered.z);
        }
    }
}
