//! Rule selection, numerical tolerances, and work limits.
use std::time::Duration;

#[derive(Clone, Debug)]
pub struct Settings {
    /// Soft budget per presolve call, including model construction and export
    /// preparation, excluding Presolver initialization. A transaction already
    /// in progress finishes before the time limit is observed.
    pub time_limit: Duration,
    pub rules: Rules,
    pub numerics: Numerics,
    /// Maximum newly introduced A and P coefficients per substitution.
    /// Use usize::MAX for unrestricted fill; Hessian growth is a separate policy.
    pub substitution_fill: usize,
    /// Allow substitutions to increase the total number of Hessian nonzeros.
    pub allow_hessian_growth: bool,
    pub equalities: EqualitySettings,
    pub propagation: PropagationSettings,
    pub progress: Progress,
    pub sparsification: SparsificationSettings,
    /// Maximum thread count for parallel fingerprinting and sorting: 1 is serial
    /// (default), 0 lets Rayon select automatically (honoring RAYON_NUM_THREADS).
    /// Presolver owns the pool, initialized once on construction; small scans
    /// still run serially. Reductions and postsolve records are applied in
    /// deterministic serial order. Pool creation failures return InitError.
    pub threads: usize,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            time_limit: Duration::from_secs(60),
            rules: Rules::default(),
            numerics: Numerics::default(),
            substitution_fill: 64,
            allow_hessian_growth: false,
            equalities: EqualitySettings::default(),
            propagation: PropagationSettings::default(),
            progress: Progress::default(),
            sparsification: SparsificationSettings::default(),
            threads: 1,
        }
    }
}
impl Settings {
    /// Favor dimensional reduction within the caller's soft time budget.
    /// Allow unlimited fill and Hessian growth, bounded and quadratic equality
    /// pivots, and exploration after any model edit. Candidate row and column
    /// lengths are capped at 16 to avoid spending the budget on very expensive
    /// substitutions. These limits and every other field remain configurable.
    /// Sparsification uses only equality references, adding no activity variables.
    /// The default execution width remains one thread.
    pub fn aggressive(time_limit: Duration) -> Self {
        Self {
            time_limit,
            substitution_fill: usize::MAX,
            allow_hessian_growth: true,
            equalities: EqualitySettings {
                max_row_length: 16,
                min_column_length: 1,
                max_column_length: 16,
                require_free_variable: false,
                require_linear_variable: false,
                preserve_nonzeros: false,
                work_limit: WorkLimit::Unlimited,
            },
            propagation: PropagationSettings {
                additional_rounds: usize::MAX,
                work_limit: WorkLimit::Unlimited,
            },
            progress: Progress::AnyChange,
            sparsification: SparsificationSettings {
                allow_auxiliary_variables: false,
                ..SparsificationSettings::default()
            },
            ..Self::default()
        }
    }
}

/// Independent switches for rule families.
#[derive(Clone, Copy, Debug)]
pub struct Rules {
    pub fixed_variables: bool,
    pub empty_columns: bool,
    pub empty_rows: bool,
    pub dual_fixing: bool,
    pub singleton_rows: bool,
    pub singleton_columns: bool,
    pub doubleton_equalities: bool,
    pub short_equalities: bool,
    pub bound_propagation: bool,
    pub redundant_bounds: bool,
    pub parallel_rows: bool,
    pub parallel_columns: bool,
    pub sparsification: bool,
    /// Structural and constant-block cone rules with finite dual recovery.
    /// Partial PSD faces that can destroy dual attainment are retained.
    pub cones: bool,
}
impl Rules {
    pub const fn none() -> Self {
        Self {
            fixed_variables: false,
            empty_columns: false,
            empty_rows: false,
            dual_fixing: false,
            singleton_rows: false,
            singleton_columns: false,
            doubleton_equalities: false,
            short_equalities: false,
            bound_propagation: false,
            redundant_bounds: false,
            parallel_rows: false,
            parallel_columns: false,
            sparsification: false,
            cones: false,
        }
    }
}
impl Default for Rules {
    fn default() -> Self {
        Self {
            fixed_variables: true,
            empty_columns: true,
            empty_rows: true,
            dual_fixing: true,
            singleton_rows: true,
            singleton_columns: true,
            doubleton_equalities: true,
            short_equalities: true,
            bound_propagation: true,
            redundant_bounds: true,
            parallel_rows: true,
            parallel_columns: true,
            sparsification: true,
            cones: true,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Numerics {
    /// Relative feasibility margin. Does not round intervals or cone coordinates
    /// to zero; geometric rules require exact structural zeros.
    pub feasibility: f64,
    pub parallel: f64,
    pub huge_bound: f64,
}
impl Default for Numerics {
    fn default() -> Self {
        Self {
            feasibility: 1e-9,
            parallel: 1e-12,
            huge_bound: 1e7,
        }
    }
}

/// Candidate restrictions for the `short_equalities` rule. The name is retained
/// for compatibility; raising these limits also enables long equalities.
/// The pivot must still have maximum coefficient magnitude in its row.
#[derive(Clone, Copy, Debug)]
pub struct EqualitySettings {
    pub max_row_length: usize,
    pub min_column_length: usize,
    pub max_column_length: usize,
    pub require_free_variable: bool,
    pub require_linear_variable: bool,
    /// Limit fill to the entries removed with the equation and pivot column.
    pub preserve_nonzeros: bool,
    /// Default allowance: twice the current constraint nonzeros per pass.
    pub work_limit: WorkLimit,
}
impl Default for EqualitySettings {
    fn default() -> Self {
        Self {
            max_row_length: 8,
            min_column_length: 2,
            max_column_length: 8,
            require_free_variable: true,
            require_linear_variable: true,
            preserve_nonzeros: true,
            work_limit: WorkLimit::Default,
        }
    }
}

/// Work allowances count estimated sparse-entry visits, not elapsed time.
/// The independent time limit applies even with Unlimited.
#[derive(Clone, Copy, Debug, Default)]
pub enum WorkLimit {
    #[default]
    Default,
    Entries(usize),
    Unlimited,
}
impl WorkLimit {
    pub(crate) fn resolve(self, default: usize) -> usize {
        match self {
            Self::Default => default,
            Self::Entries(n) => n,
            Self::Unlimited => usize::MAX,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct PropagationSettings {
    /// Extra rounds after the initial propagation pass. usize::MAX removes the cap.
    pub additional_rounds: usize,
    /// Default allowance: max(A nonzeros / 4, 256) across the extra rounds.
    pub work_limit: WorkLimit,
}
impl Default for PropagationSettings {
    fn default() -> Self {
        Self {
            additional_rounds: 3,
            work_limit: WorkLimit::Default,
        }
    }
}

/// Decide whether to repeat the fast phases and full exploration cycles.
#[derive(Clone, Copy, Debug)]
pub enum Progress {
    /// Require a strictly larger fractional reduction in A + P nonzeros.
    /// Values are clamped to [0, 1]; NaN is treated as the default 0.05.
    Nonzeros { minimum_reduction: f64 },
    /// Continue after any model edit, including bound changes and substitutions
    /// that increase nonzeros. This does not assert a mathematical fixed point:
    /// rule-specific candidate, numerical, work, and time limits still apply.
    AnyChange,
}
impl Default for Progress {
    fn default() -> Self {
        Self::Nonzeros {
            minimum_reduction: 0.05,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SparsificationSettings {
    /// Allow inequality references, which introduce an auxiliary activity variable.
    /// False restricts the pass to equality references.
    pub allow_auxiliary_variables: bool,
    /// Default allowance: max(8 * (A nonzeros + bound entries), 1024).
    pub work_limit: WorkLimit,
}
impl Default for SparsificationSettings {
    fn default() -> Self {
        Self {
            allow_auxiliary_variables: true,
            work_limit: WorkLimit::Default,
        }
    }
}
