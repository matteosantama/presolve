//! Per-rule reduction counts and phase statistics.
use presolve::matrix::CscMatrix;
use presolve::problem::{Bounds, Constraint};
use presolve::result::{ReductionKind, RuleId};
use presolve::settings::Rules;
use presolve::{Outcome, Presolver, Problem, Settings};
use std::time::Duration;

/// Minimize x0 + x1 with x0 fixed at 3 and x1 in [2, 4], both outside every row.
fn two_columns() -> Problem {
    Problem {
        p: None,
        c: vec![1.0, 1.0],
        c0: 0.0,
        a: CscMatrix::zeros(0, 2).unwrap(),
        rows: vec![],
        variable_bounds: vec![
            Bounds::fixed(3.0),
            Bounds {
                lower: 2.0,
                upper: 4.0,
            },
        ],
        cones: vec![],
    }
}

#[test]
fn reductions_are_credited_to_the_rule_that_applied_them() {
    let result = Presolver::default().presolve(two_columns());
    assert!(matches!(result.outcome, Outcome::Solved(_)));
    let reductions = &result.stats.reductions;
    assert_eq!(
        reductions.count(RuleId::FixedVariables, ReductionKind::Fixed),
        1
    );
    assert_eq!(
        reductions.count(RuleId::EmptyColumns, ReductionKind::Fixed),
        1
    );
    assert_eq!(reductions.total(), 2);
    assert_eq!(reductions.kind_total(ReductionKind::Fixed), 2);
    assert_eq!(reductions.rule_total(RuleId::EmptyColumns), 1);
    let listed: Vec<_> = reductions.iter().collect();
    assert_eq!(
        listed,
        [
            (RuleId::FixedVariables, ReductionKind::Fixed, 1),
            (RuleId::EmptyColumns, ReductionKind::Fixed, 1),
        ]
    );
    assert_eq!(
        format!("{reductions:?}"),
        "{fixed_variables.fixed: 1, empty_columns.fixed: 1}"
    );
    assert!(result.stats.phases.fast >= 1);
}

#[test]
fn relaxed_bounds_are_counted_without_a_tape_record() {
    // Minimize x0 - x1 subject to x0 + x1 <= 2 and x0 - x1 >= -5 with x0 in
    // [0, 2] and x1 in [0, 3]. Only redundant-bound removal is enabled. The
    // first row and the lower bounds imply both upper bounds, so the rule
    // drops both, and each removal is counted although nothing is recorded
    // for recovery.
    let problem = Problem {
        p: None,
        c: vec![1.0, -1.0],
        c0: 0.0,
        a: CscMatrix::from_triplets(
            2,
            2,
            vec![0, 1, 0, 1],
            vec![0, 0, 1, 1],
            vec![1.0, 1.0, 1.0, -1.0],
        )
        .unwrap(),
        rows: vec![
            Constraint::Linear(Bounds {
                lower: f64::NEG_INFINITY,
                upper: 2.0,
            }),
            Constraint::Linear(Bounds {
                lower: -5.0,
                upper: f64::INFINITY,
            }),
        ],
        variable_bounds: vec![
            Bounds {
                lower: 0.0,
                upper: 2.0,
            },
            Bounds {
                lower: 0.0,
                upper: 3.0,
            },
        ],
        cones: vec![],
    };
    let settings = Settings {
        rules: Rules {
            redundant_bounds: true,
            ..Rules::none()
        },
        ..Settings::default()
    };
    let result = Presolver::new(settings).unwrap().presolve(problem);
    let reductions = &result.stats.reductions;
    assert_eq!(
        reductions.iter().collect::<Vec<_>>(),
        [(RuleId::RedundantBounds, ReductionKind::RelaxedBound, 2)]
    );
}

#[test]
fn an_unlimited_time_budget_does_not_overflow() {
    let settings = Settings {
        time_limit: Duration::MAX,
        ..Settings::default()
    };
    let result = Presolver::new(settings).unwrap().presolve(two_columns());
    assert!(matches!(result.outcome, Outcome::Solved(_)));
    assert!(!result.stats.time_limit_reached);
    assert_eq!(result.stats.reductions.total(), 2);
}

#[test]
fn every_rule_and_kind_has_a_distinct_name() {
    let rules: std::collections::BTreeSet<_> = RuleId::ALL.iter().map(|r| r.name()).collect();
    assert_eq!(rules.len(), RuleId::ALL.len());
    let kinds: std::collections::BTreeSet<_> =
        ReductionKind::ALL.iter().map(|k| k.name()).collect();
    assert_eq!(kinds.len(), ReductionKind::ALL.len());
    assert_eq!(
        RuleId::ALL.iter().position(|&r| r == RuleId::Cones),
        Some(RuleId::Cones as usize)
    );
    assert_eq!(
        ReductionKind::ALL
            .iter()
            .position(|&k| k == ReductionKind::Eliminated),
        Some(ReductionKind::Eliminated as usize)
    );
}
