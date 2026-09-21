// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Row activity bounds and counts of constraints that lock variable directions.

use crate::{model::RowDomain, model::tape::Side, problem::Bounds};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Locks {
    pub up: usize,
    pub down: usize,
}

impl Locks {
    #[inline]
    pub fn contribution(coefficient: f64, domain: RowDomain) -> Self {
        match domain {
            RowDomain::Deleted => Self::default(),
            RowDomain::Cone { .. } => Self { up: 1, down: 1 },
            RowDomain::Linear(b) => {
                let lower = usize::from(b.lower.is_finite());
                let upper = usize::from(b.upper.is_finite());
                if coefficient > 0.0 {
                    Self {
                        up: upper,
                        down: lower,
                    }
                } else {
                    Self {
                        up: lower,
                        down: upper,
                    }
                }
            }
        }
    }

    #[inline]
    pub fn add(&mut self, other: Self) {
        self.up += other.up;
        self.down += other.down;
    }
    #[inline]
    pub fn remove(&mut self, other: Self) {
        self.up -= other.up;
        self.down -= other.down;
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Extreme {
    /// Sum of finite terms, retained even when another term is infinite.
    pub sum: f64,
    pub infinite: usize,
}

impl Extreme {
    /// Incremental bound updates, with a direct-recomputation fallback when
    /// subtracting or adding a term would lose the remaining finite sum.
    fn replace(&mut self, old: f64, new: f64) -> bool {
        if old == new {
            return true;
        }
        if old.is_finite() {
            self.sum -= old;
            if old != 0.0 && self.sum.abs() < 32.0 * f64::EPSILON * old.abs() {
                return false;
            }
        } else {
            self.infinite -= 1;
        }
        if new.is_finite() {
            self.sum += new;
            if new != 0.0 && self.sum.abs() < 32.0 * f64::EPSILON * new.abs() {
                return false;
            }
        } else {
            self.infinite += 1;
        }
        self.sum.is_finite()
    }
    /// Branch free: `sum` starts at +0.0 and only finite terms are added, so
    /// it is never -0.0 and adding 0.0 leaves it bit-identical.
    #[inline]
    fn add(&mut self, term: f64) {
        let finite = term.is_finite();
        self.sum += if finite { term } else { 0.0 };
        self.infinite += usize::from(!finite);
    }

    #[inline]
    pub fn value(self) -> Option<f64> {
        (self.infinite == 0 && self.sum.is_finite()).then_some(self.sum)
    }

    /// Bound implied for the excluded term's variable by `rhs`:
    /// `(rhs - residual) / a`. When the cached sum cannot exclude the term
    /// because of cancellation, `recompute` supplies the residual directly,
    /// but only when this term is the sole infinite one: any other infinite
    /// contribution keeps the residual infinite.
    #[inline]
    pub fn implied(
        self,
        term: f64,
        rhs: f64,
        a: f64,
        recompute: impl FnOnce() -> Option<f64>,
    ) -> Option<f64> {
        self.excluding(term)
            .or_else(|| {
                (self.infinite == usize::from(!term.is_finite()))
                    .then(recompute)
                    .flatten()
            })
            .map(|v| (rhs - v) / a)
    }

    #[inline]
    pub fn excluding(self, term: f64) -> Option<f64> {
        if term.is_finite() {
            let residual = self.sum - term;
            // The caller recomputes a cancellation-sensitive residual directly.
            let cancellation = term != 0.0 && residual.abs() < 32.0 * f64::EPSILON * term.abs();
            (self.infinite == 0 && residual.is_finite() && !cancellation).then_some(residual)
        } else {
            (self.infinite == 1 && self.sum.is_finite()).then_some(self.sum)
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Activity {
    pub min: Extreme,
    pub max: Extreme,
}

impl Activity {
    /// An activity that yields no finite extreme, for callers that decline to
    /// recompute a residual.
    pub const UNKNOWN: Self = Self {
        min: Extreme {
            sum: 0.0,
            infinite: 1,
        },
        max: Extreme {
            sum: 0.0,
            infinite: 1,
        },
    };

    /// Both bounds one row implies for the variable with coefficient `a` and
    /// bounds `b`, as (from the row's lower side, from its upper side), each
    /// only when that side is finite. `residual` recomputes the activity
    /// without the variable when the cached extremes cannot exclude its term.
    #[inline]
    pub fn implied(
        self,
        a: f64,
        b: Bounds,
        rhs: Bounds,
        residual: impl Fn() -> Activity,
    ) -> (Option<f64>, Option<f64>) {
        let (min_term, max_term) = Self::terms(a, b);
        let lower = rhs
            .lower
            .is_finite()
            .then(|| {
                self.max
                    .implied(max_term, rhs.lower, a, || residual().max.value())
            })
            .flatten();
        let upper = rhs
            .upper
            .is_finite()
            .then(|| {
                self.min
                    .implied(min_term, rhs.upper, a, || residual().min.value())
            })
            .flatten();
        (lower, upper)
    }

    /// The bound implied for the variable's `side` by one row, when that
    /// row's relevant extreme excludes the variable's own term.
    #[inline]
    pub fn implied_side(
        self,
        a: f64,
        b: Bounds,
        rhs: Bounds,
        side: Side,
        residual: impl Fn() -> Activity,
    ) -> Option<f64> {
        let (min_term, max_term) = Self::terms(a, b);
        // A positive coefficient takes the variable's lower bound from the
        // row's lower side; a negative one flips the sides.
        if (side == Side::Lower) == (a > 0.0) {
            rhs.lower
                .is_finite()
                .then(|| {
                    self.max
                        .implied(max_term, rhs.lower, a, || residual().max.value())
                })
                .flatten()
        } else {
            rhs.upper
                .is_finite()
                .then(|| {
                    self.min
                        .implied(min_term, rhs.upper, a, || residual().min.value())
                })
                .flatten()
        }
    }

    #[inline]
    pub fn replace_bound(&mut self, a: f64, old: Bounds, new: Bounds) -> bool {
        if old.lower == new.lower {
            let extreme = if a > 0.0 {
                &mut self.max
            } else {
                &mut self.min
            };
            return extreme.replace(a * old.upper, a * new.upper);
        }
        if old.upper == new.upper {
            let extreme = if a > 0.0 {
                &mut self.min
            } else {
                &mut self.max
            };
            return extreme.replace(a * old.lower, a * new.lower);
        }
        let (old_min, old_max) = Self::terms(a, old);
        let (new_min, new_max) = Self::terms(a, new);
        self.min.replace(old_min, new_min) && self.max.replace(old_max, new_max)
    }
    #[inline]
    pub fn terms(a: f64, b: Bounds) -> (f64, f64) {
        if a > 0.0 {
            (a * b.lower, a * b.upper)
        } else {
            (a * b.upper, a * b.lower)
        }
    }

    pub fn compute(row: impl IntoIterator<Item = (usize, f64)>, bounds: &[Bounds]) -> Self {
        use wide::{f64x2, u64x2};
        // One lane per extreme preserves the scalar accumulation order. Explicit
        // masks keep nonfinite classification free of per-term constant setup.
        const EXPONENT: u64x2 = u64x2::new([0x7ff0000000000000; 2]);
        const ONE: u64x2 = u64x2::new([1; 2]);
        let mut sums = f64x2::ZERO;
        let mut infinite = u64x2::ZERO;
        for (j, a) in row {
            let b = bounds[j];
            let coefficient = f64x2::splat(a);
            let sides = coefficient.simd_gt(f64x2::ZERO).bitselect(
                f64x2::new([b.lower, b.upper]),
                f64x2::new([b.upper, b.lower]),
            );
            let terms = sides * coefficient;
            let bits = u64x2::new(terms.to_array().map(f64::to_bits));
            let nonfinite = (bits & EXPONENT).simd_eq(EXPONENT);
            let finite = f64x2::new((!nonfinite).to_array().map(f64::from_bits));
            sums += terms & finite;
            infinite += nonfinite & ONE;
        }
        let [min_sum, max_sum] = sums.to_array();
        let [min_infinite, max_infinite] = infinite.to_array();
        Self {
            min: Extreme {
                sum: min_sum,
                infinite: min_infinite as usize,
            },
            max: Extreme {
                sum: max_sum,
                infinite: max_infinite as usize,
            },
        }
    }
    /// The activity without the term of column `exclude`.
    pub fn compute_excluding(
        row: impl IntoIterator<Item = (usize, f64)>,
        bounds: &[Bounds],
        exclude: usize,
    ) -> Self {
        let mut out = Self::default();
        for (j, a) in row {
            if j == exclude {
                continue;
            }
            let (min, max) = Self::terms(a, bounds[j]);
            out.min.add(min);
            out.max.add(max);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incremental_updates_handle_infinity_and_request_recomputation_after_cancellation() {
        let mut extreme = Extreme {
            sum: 3.0,
            infinite: 1,
        };
        assert!(extreme.replace(f64::INFINITY, 4.0));
        assert_eq!(extreme.value(), Some(7.0));
        assert!(extreme.replace(4.0, 5.0));
        assert_eq!(extreme.value(), Some(8.0));
        let mut lost = Extreme {
            sum: 1e20,
            infinite: 0,
        };
        assert!(!lost.replace(1e20, 1.0));
    }

    #[test]
    fn residual_activities_handle_infinite_terms_and_cancellation() {
        let bounds = [
            Bounds {
                lower: 2.0,
                upper: f64::INFINITY,
            },
            Bounds {
                lower: 3.0,
                upper: 4.0,
            },
        ];
        let act = Activity::compute([(0, 2.0), (1, -1.0)], &bounds);
        assert_eq!(act.min.value(), Some(0.0));
        assert_eq!(act.max.value(), None);
        assert_eq!(act.max.excluding(f64::INFINITY), Some(-3.0));
        let lost = Extreme {
            sum: 1e20,
            infinite: 0,
        };
        assert_eq!(lost.excluding(1e20), None);
        assert_eq!(
            Locks::contribution(
                -2.0,
                RowDomain::Linear(Bounds {
                    lower: 0.0,
                    upper: f64::INFINITY
                })
            ),
            Locks { up: 1, down: 0 }
        );
    }

    #[test]
    fn paired_extremes_match_scalar_bits() {
        let mut seed = 23781882u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for case in 0..3000 {
            let n = case % 2048;
            let mut row = Vec::with_capacity(n);
            let mut bounds = Vec::with_capacity(n);
            let mut scalar = Activity::default();
            for j in 0..n {
                let values = [
                    -f64::INFINITY,
                    -f64::MAX,
                    -1e200,
                    -1e-200,
                    -1.0,
                    -0.0,
                    0.0,
                    1e-200,
                    1.0,
                    1e200,
                    f64::MAX,
                    f64::INFINITY,
                ];
                let x = values[next() as usize % values.len()];
                let y = values[next() as usize % values.len()];
                let b = Bounds {
                    lower: x.min(y),
                    upper: x.max(y),
                };
                let coeffs = [
                    -f64::MAX,
                    -1e200,
                    -1e-200,
                    -1.0,
                    -f64::MIN_POSITIVE,
                    f64::MIN_POSITIVE,
                    1e-200,
                    1.0,
                    1e200,
                    f64::MAX,
                ];
                let a = coeffs[next() as usize % coeffs.len()];
                row.push((j, a));
                bounds.push(b);
                let (min, max) = Activity::terms(a, b);
                scalar.min.add(min);
                scalar.max.add(max);
            }
            let paired = Activity::compute(row, &bounds);
            assert_eq!(
                (
                    paired.min.sum.to_bits(),
                    paired.min.infinite,
                    paired.max.sum.to_bits(),
                    paired.max.infinite
                ),
                (
                    scalar.min.sum.to_bits(),
                    scalar.min.infinite,
                    scalar.max.sum.to_bits(),
                    scalar.max.infinite
                )
            );
        }
    }
    #[test]
    fn single_bound_updates_match_the_two_extreme_reference() {
        let domains = [
            Bounds::FREE,
            Bounds::fixed(0.),
            Bounds::fixed(1.),
            Bounds {
                lower: 0.,
                upper: f64::INFINITY,
            },
            Bounds {
                lower: f64::NEG_INFINITY,
                upper: 0.,
            },
            Bounds {
                lower: -2.,
                upper: 3.,
            },
            Bounds {
                lower: -2.,
                upper: 4.,
            },
            Bounds {
                lower: -1.,
                upper: 3.,
            },
            Bounds {
                lower: -f64::MAX,
                upper: f64::MAX,
            },
        ];
        for a in [
            -f64::MAX,
            -2.,
            -f64::MIN_POSITIVE,
            f64::MIN_POSITIVE,
            2.,
            f64::MAX,
        ] {
            for old in domains {
                for new in domains {
                    let initial =
                        Activity::compute([(0, a), (1, 1.)], &[old, Bounds::fixed(1e300)]);
                    let mut reference = initial;
                    let (old_min, old_max) = Activity::terms(a, old);
                    let (new_min, new_max) = Activity::terms(a, new);
                    let expected = reference.min.replace(old_min, new_min)
                        && reference.max.replace(old_max, new_max);
                    let mut actual = initial;
                    assert_eq!(actual.replace_bound(a, old, new), expected);
                    if expected {
                        assert_eq!(actual.min.sum.to_bits(), reference.min.sum.to_bits());
                        assert_eq!(actual.max.sum.to_bits(), reference.max.sum.to_bits());
                        assert_eq!(actual.min.infinite, reference.min.infinite);
                        assert_eq!(actual.max.infinite, reference.max.infinite);
                    }
                }
            }
        }
    }
}
