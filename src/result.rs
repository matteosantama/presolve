//! Presolve outcomes, model sizes, and execution statistics.
use crate::postsolve::{Postsolve, PrimalCertificate, Solution};
use crate::problem::Problem;
use std::time::Duration;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Size {
    pub variables: usize,
    pub linear_rows: usize,
    pub conic_rows: usize,
    /// Nonzero values, excluding explicitly stored zeros.
    pub a_nonzeros: usize,
    pub g_nonzeros: usize,
    /// Nonzeros of the full symmetric matrix: off-diagonal pairs count twice.
    pub p_nonzeros: usize,
}
/// Short-equality candidate decisions. Counts include revisits across passes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EqualityStats {
    pub rows_examined: usize,
    /// Rows outside row limits, plus individual candidates outside column/domain limits.
    pub structural_rejections: usize,
    pub pivot_rejections: usize,
    pub work_rejections: usize,
    pub attempts: usize,
    /// Sum of the five transactional rejection categories below.
    pub rejected_updates: usize,
    pub numerical_rejections: usize,
    pub constraint_fill_rejections: usize,
    pub quadratic_fill_rejections: usize,
    pub hessian_growth_rejections: usize,
    pub deadline_rejections: usize,
    pub accepted: usize,
    pub estimated_work: usize,
}
#[derive(Clone, Debug)]
pub struct Stats {
    pub equalities: EqualityStats,
    pub elapsed: Duration,
    pub time_limit_reached: bool,
    pub before: Size,
    /// None when presolve stops on a certificate before producing a problem.
    pub after: Option<Size>,
    pub parallel_comparisons: usize,
    /// Whether Hessian coefficients or variable coordinates changed. A caller
    /// retaining its original matrix can reuse it when this is false.
    pub quadratic_changed: bool,
    pub reductions: Reductions,
    pub phases: Phases,
}

#[derive(Debug)]
pub struct PresolveResult {
    pub outcome: Outcome,
    pub stats: Stats,
}
// Returning unchanged data inline avoids an allocation on that path.
#[allow(clippy::large_enum_variant)]
#[derive(Debug)]
pub enum Outcome {
    /// Return the input allocations when no rule was applied.
    Unchanged(Problem),
    Reduced(Box<ReducedProblem>),
    /// Original-coordinate solution. Never emitted with remaining opaque cone blocks:
    /// the caller must check their feasibility even if all variables are fixed.
    Solved(Solution),
    Infeasible(PrimalCertificate),
    /// A feasible point and an objective-decreasing recession direction.
    Unbounded(UnboundednessCertificate),
}
#[derive(Debug)]
pub struct ReducedProblem {
    pub problem: Problem,
    pub postsolve: Postsolve,
}

#[derive(Clone, Debug)]
pub struct UnboundednessCertificate {
    pub point: Vec<f64>,
    pub ray: Vec<f64>,
}

/// Rule families that attribute reductions, one per switch in `settings::Rules`.
/// Names match the switch fields.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RuleId {
    FixedVariables,
    EmptyColumns,
    QuadraticElimination,
    EmptyRows,
    DualFixing,
    DualPropagation,
    SingletonRows,
    SingletonColumns,
    DoubletonEqualities,
    ShortEqualities,
    ImpliedFreeEqualities,
    BoundShift,
    LpFolding,
    EqualityDependencies,
    BoundPropagation,
    RedundantBounds,
    ParallelRows,
    ParallelColumns,
    DominatedColumns,
    Sparsification,
    Cones,
}
impl RuleId {
    pub const ALL: [Self; 21] = [
        Self::FixedVariables,
        Self::EmptyColumns,
        Self::QuadraticElimination,
        Self::EmptyRows,
        Self::DualFixing,
        Self::DualPropagation,
        Self::SingletonRows,
        Self::SingletonColumns,
        Self::DoubletonEqualities,
        Self::ShortEqualities,
        Self::ImpliedFreeEqualities,
        Self::BoundShift,
        Self::LpFolding,
        Self::EqualityDependencies,
        Self::BoundPropagation,
        Self::RedundantBounds,
        Self::ParallelRows,
        Self::ParallelColumns,
        Self::DominatedColumns,
        Self::Sparsification,
        Self::Cones,
    ];
    /// The `settings::Rules` field name.
    pub fn name(self) -> &'static str {
        match self {
            Self::FixedVariables => "fixed_variables",
            Self::EmptyColumns => "empty_columns",
            Self::QuadraticElimination => "quadratic_elimination",
            Self::EmptyRows => "empty_rows",
            Self::DualFixing => "dual_fixing",
            Self::DualPropagation => "dual_propagation",
            Self::SingletonRows => "singleton_rows",
            Self::SingletonColumns => "singleton_columns",
            Self::DoubletonEqualities => "doubleton_equalities",
            Self::ShortEqualities => "short_equalities",
            Self::ImpliedFreeEqualities => "implied_free_equalities",
            Self::BoundShift => "bound_shift",
            Self::LpFolding => "lp_folding",
            Self::EqualityDependencies => "equality_dependencies",
            Self::BoundPropagation => "bound_propagation",
            Self::RedundantBounds => "redundant_bounds",
            Self::ParallelRows => "parallel_rows",
            Self::ParallelColumns => "parallel_columns",
            Self::DominatedColumns => "dominated_columns",
            Self::Sparsification => "sparsification",
            Self::Cones => "cones",
        }
    }
}

/// Kinds of model transformation. All but `RelaxedBound` correspond to one
/// recovery-tape record; relaxing an implied bound needs no recovery.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReductionKind {
    /// One equitable LP quotient, however many rows and columns it merges.
    LpFold,
    BoundShift,
    /// One chain of doubleton substitutions applied as a batch.
    DoubletonChain,
    SocAggregated,
    ConeSlack,
    SocToLinear,
    SocFace,
    PsdZeroFace,
    /// One reference row combined into its targets.
    RowCombination,
    Fixed,
    Substituted,
    DependentRow,
    MergedRow,
    DeletedRow,
    TightenedBound,
    RelaxedBound,
    TightenedRow,
    ParallelColumns,
    Unlocked,
    Eliminated,
}
impl ReductionKind {
    pub const ALL: [Self; 20] = [
        Self::LpFold,
        Self::BoundShift,
        Self::DoubletonChain,
        Self::SocAggregated,
        Self::ConeSlack,
        Self::SocToLinear,
        Self::SocFace,
        Self::PsdZeroFace,
        Self::RowCombination,
        Self::Fixed,
        Self::Substituted,
        Self::DependentRow,
        Self::MergedRow,
        Self::DeletedRow,
        Self::TightenedBound,
        Self::RelaxedBound,
        Self::TightenedRow,
        Self::ParallelColumns,
        Self::Unlocked,
        Self::Eliminated,
    ];
    pub fn name(self) -> &'static str {
        match self {
            Self::LpFold => "lp_fold",
            Self::BoundShift => "bound_shift",
            Self::DoubletonChain => "doubleton_chain",
            Self::SocAggregated => "soc_aggregated",
            Self::ConeSlack => "cone_slack",
            Self::SocToLinear => "soc_to_linear",
            Self::SocFace => "soc_face",
            Self::PsdZeroFace => "psd_zero_face",
            Self::RowCombination => "row_combination",
            Self::Fixed => "fixed",
            Self::Substituted => "substituted",
            Self::DependentRow => "dependent_row",
            Self::MergedRow => "merged_row",
            Self::DeletedRow => "deleted_row",
            Self::TightenedBound => "tightened_bound",
            Self::RelaxedBound => "relaxed_bound",
            Self::TightenedRow => "tightened_row",
            Self::ParallelColumns => "parallel_columns",
            Self::Unlocked => "unlocked",
            Self::Eliminated => "eliminated",
        }
    }
}

/// Reductions applied during one presolve call, counted by the rule that
/// applied them and by kind. The count is a pure function of the input and
/// the settings when no time limit is reached.
#[derive(Clone, Default, PartialEq, Eq)]
pub struct Reductions {
    counts: [[usize; ReductionKind::ALL.len()]; RuleId::ALL.len()],
}
impl Reductions {
    pub fn count(&self, rule: RuleId, kind: ReductionKind) -> usize {
        self.counts[rule as usize][kind as usize]
    }
    pub fn rule_total(&self, rule: RuleId) -> usize {
        self.counts[rule as usize].iter().sum()
    }
    pub fn kind_total(&self, kind: ReductionKind) -> usize {
        self.counts.iter().map(|row| row[kind as usize]).sum()
    }
    pub fn total(&self) -> usize {
        self.counts.iter().flatten().sum()
    }
    /// Nonzero entries in `RuleId::ALL` then `ReductionKind::ALL` order.
    pub fn iter(&self) -> impl Iterator<Item = (RuleId, ReductionKind, usize)> + '_ {
        RuleId::ALL.into_iter().flat_map(move |rule| {
            ReductionKind::ALL
                .into_iter()
                .map(move |kind| (rule, kind, self.count(rule, kind)))
                .filter(|&(_, _, n)| n != 0)
        })
    }
    pub(crate) fn add(&mut self, rule: RuleId, kind: ReductionKind) {
        self.counts[rule as usize][kind as usize] += 1;
    }
}
impl std::fmt::Debug for Reductions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut map = f.debug_map();
        for (rule, kind, n) in self.iter() {
            map.entry(&format_args!("{}.{}", rule.name(), kind.name()), &n);
        }
        map.finish()
    }
}

/// Scheduler phases completed, including a phase cut short by a certificate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Phases {
    pub fast: usize,
    pub medium: usize,
}
