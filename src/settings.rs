//! Rule selection, numerical tolerances, and work limits.
use std::time::Duration;

#[derive(Clone, Debug)]
pub struct Settings {
    /// Soft budget including construction and export preparation. A transaction
    /// already in progress finishes before the time limit is observed.
    pub time_limit: Duration,
    pub rules: Rules,
    pub numerics: Numerics,
    pub substitution_fill: usize,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            time_limit: Duration::from_secs(60),
            rules: Rules::default(),
            numerics: Numerics::default(),
            substitution_fill: 64,
        }
    }
}
/// Independent switches for rule families. Fixed-column substitution also
/// observes the substitution fill policy; Hessian nnz can never increase.
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
