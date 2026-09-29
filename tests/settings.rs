use presolve::{Settings, settings::WorkLimit};
use std::time::Duration;
#[test]
fn defaults_enable_structural_reductions_without_unbounded_or_quadratic_fill() {
    let s = Settings::default();
    assert_eq!(s.time_limit, Duration::from_secs(60));
    assert!(s.rules.lp_folding && s.rules.implied_free_equalities);
    assert!(s.allow_unit_doubleton_fill && s.final_parallel_scan);
    assert!(!s.allow_hessian_growth);
    assert_eq!(s.substitution_fill, 64);
    assert_eq!(s.equalities.relative_pivot, 1.0);
    assert!(matches!(s.equalities.work_limit, WorkLimit::Default));
    assert!(matches!(s.folding.work_limit, WorkLimit::Default));
    let aggressive = Settings::aggressive(Duration::from_secs(7));
    assert_eq!(aggressive.time_limit, Duration::from_secs(7));
    assert!(aggressive.rules.lp_folding && aggressive.rules.implied_free_equalities);
    assert!(aggressive.allow_unit_doubleton_fill && aggressive.final_parallel_scan);
    assert!(aggressive.allow_hessian_growth);
    assert_eq!(aggressive.substitution_fill, usize::MAX);
}
