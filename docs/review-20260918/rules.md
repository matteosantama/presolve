# Refactoring review: reduction rules and scheduling

Scope: `src/core/rules/*.rs`, `src/core/schedule.rs`, `src/settings.rs` (read in
full); `src/core/model.rs`, `src/core/activity.rs`, `src/core/queues.rs`,
`src/postsolve/tape.rs`, `src/presolve.rs` and the README rule catalog and
scheduling sections (read for context). All line numbers refer to the working
tree at commit `803ecea` on branch `new-presolve-rules`. Every claim below was
checked against the code; nothing is inferred from names alone.

The constraint that shapes every recommendation: reductions must stay
bit-identical on the benchmark corpus. Each proposal therefore states whether
it can change (a) floating-point arithmetic or its order, (b) the order in
which candidates are visited, or (c) the order in which queues are filled or
drained. Proposals that touch none of the three are marked "zero behavioral
risk"; the remaining risk is only that of an editing mistake, which the
benchmark comparison catches.

## Ranking (value / risk)

| # | Proposal | Value | Risk | Effort |
|---|----------|-------|------|--------|
| 1 | One `linear_column` predicate instead of `dual_column` + `linear_column` | small | none | 10 min |
| 2 | One residual-activity kernel behind `propagate_bounds`, `bound_implied`, `open_sides`, `dual_propagation` | medium-high | none (same arithmetic, same fallbacks) | 2-3 h |
| 3 | Certificate constructors replacing nine hand-built `Point`s | medium | none | 1 h |
| 4 | Keep run settings on `Model`; delete `Limits`, the `max_fill` parameters and the field/method name collisions | medium-high | none | 2-3 h |
| 5 | One drain loop and one substitution trio in `schedule.rs` | medium | none (identical call order) | 1 h |
| 6 | One stats sink on `Model`; uniform rule return types | medium | none | 1-2 h |
| 7 | Split the four long functions (`simplify_cones`, `dual_propagation`, `short_equalities`, `sparsify_rows`) | medium | low, except `simplify_cones` (no test coverage in repo) | 4-6 h |
| 8 | Per-rule structs for the two scratch-carrying rules | small-medium | none | 1-2 h |
| 9 | A `Budget` type for count-down work limits | small | low; do not convert `sparsify_rows` | 1 h |

Assessments that are not proposals: the `impl Model`-per-family layout is the
right shape (section "Module layout"); tests are reasonably placed except that
`cones.rs` has no coverage in the repository at all (section "Tests"); two
documentation errors in `DualPropagationSettings` (section "Documentation
errors found").

---

## 1. One `linear_column` predicate

**Files.** `src/core/rules/dual_propagation.rs:99-109` (`dual_column`),
`src/core/rules/dominated_columns.rs:56-66` (`linear_column`).

**Problem.** The two private methods are character-for-character identical:
alive, absent from the Hessian column, and (when cones exist) every incident
row is `RowDomain::Linear`. They only differ in doc comment and name. Both
rules are the ones enabled by the aggressive preset and both will need the same
predicate again if a third dual-side rule is added.

**Proposed shape.** Move it to `model.rs` as
`pub(crate) fn linear_column(&self, j: usize) -> bool` with the
`dominated_columns.rs` doc comment ("absent from the Hessian and from every
conic row"), delete both private copies, and replace the nine call sites
(`dual_column` x5, `linear_column` x4).

**Risk.** Zero behavioral risk; the function is pure and the bodies are equal.

**Effort.** Ten minutes. Do it together with proposal 2.

---

## 2. One residual-activity kernel for implied bounds

**Files and functions.**

- `src/core/rules/bounds.rs:70-125` (`propagate_bounds`, per-column loop),
  `bounds.rs:170-211` (`bound_implied`), `bounds.rs:162-166`
  (`remove_redundant_bound`).
- `src/core/rules/dominated_columns.rs:88-126` (`open_sides`).
- `src/core/rules/dual_propagation.rs:221-252` (per-entry loop inside
  `dual_propagation`).
- `src/core/rules/substitution.rs:151-178` (`singleton_range`,
  `non_implied_bounds`).
- `src/core/activity.rs:90-99` (`Extreme::excluding`).

**Problem.** The same three-step computation is written out five times:

1. take the cached extreme (`act.max` for a lower-side implication, `act.min`
   for an upper-side one) and the column's own term from `Activity::terms`;
2. subtract the term with `Extreme::excluding`, and when that reports
   cancellation, fall back to a direct recomputation excluding the column, but
   only if `extreme.infinite == usize::from(!term.is_finite())`;
3. divide `(rhs - residual) / a`.

`propagate_bounds` writes steps 1-3 twice (lines 74-86 and 87-98) with the
fallback `self.residual_activity(i, j)`. `dual_propagation` writes the same two
blocks (lines 223-252) with fallback `Activity::compute(column, &scratch.y,
Some(i))`. `bound_implied` writes them once with the fallback, but hoists the
infinite-count test to a `continue` before calling `excluding` (line 188).
`open_sides` writes them once with no fallback at all (lines 112-118), which
the comment explains as deliberate conservatism. `singleton_range` skips the
cache and always recomputes.

`bound_implied` and `open_sides` are in addition the same *loop*: for a side of
column `j`, walk the column, skip non-linear rows, compute the implied bound
from each row, and report whether any row implies the existing side. The only
differences are (i) the fallback and (ii) `open_sides` computing both sides in
one call with a bit mask. The comparison `bound >= b.lower` / `bound <= b.upper`
and the `is_finite` filter are identical.

**Proposed shape.** Two small additions, no new arithmetic:

```rust
// activity.rs
impl Extreme {
    /// `excluding`, with a caller-supplied recomputation for the cancellation
    /// case. Recomputation is skipped when another term is infinite, because
    /// it could not produce a finite value either.
    pub fn excluding_or(self, term: f64, recompute: impl FnOnce() -> Option<f64>) -> Option<f64> {
        self.excluding(term).or_else(|| {
            (self.infinite == usize::from(!term.is_finite()))
                .then(recompute)
                .flatten()
        })
    }
}
```

```rust
// bounds.rs
impl Model {
    /// Bound on `side` of column `j` implied by linear row `i` (coefficient `a`)
    /// and the other columns' current bounds. `recompute` allows the direct
    /// residual recomputation on cancellation; `open_sides` passes `false`.
    pub(super) fn row_implied_bound(&mut self, i: usize, j: usize, a: f64, side: Side, recompute: bool) -> Option<f64> {
        let RowDomain::Linear(row) = self.rows[i] else { return None };
        let act = self.activity(i);
        let (min, max) = Activity::terms(a, self.bounds[j]);
        let from_lower = (side == Side::Lower) == (a > 0.0);
        let (rhs, extreme, term) = if from_lower { (row.lower, act.max, max) } else { (row.upper, act.min, min) };
        if !rhs.is_finite() { return None; }
        extreme
            .excluding_or(term, || {
                if !recompute { return None; }
                let r = self.residual_activity(i, j);
                if from_lower { r.max.value() } else { r.min.value() }
            })
            .map(|v| (rhs - v) / a)
            .filter(|v| v.is_finite())
    }

    fn bound_implied(&mut self, j: usize, side: Side, recompute: bool) -> bool {
        let b = self.bounds[j];
        if !side.value(b).is_finite() { return false; }
        let mut cursor = self.a.column(j).cursor();
        while let Some((i, a)) = cursor.next(&self.a) {
            let Some(bound) = self.row_implied_bound(i, j, a, side, recompute) else { continue };
            let implied = match side { Side::Lower => bound >= b.lower, Side::Upper => bound <= b.upper };
            if implied { return true; }
        }
        false
    }
}
```

`open_sides` then becomes a five-line loop over
`[(Side::Upper, UPPER_OPEN), (Side::Lower, LOWER_OPEN)]` calling
`bound_implied(j, side, false)`. In `propagate_bounds`, the two 13-line
blocks at lines 74-98 collapse to two calls of `excluding_or` (the row bounds
there are the *relaxed* `bounds`, not `self.rows[i]`, so use `excluding_or`
directly rather than `row_implied_bound`; that is also why the helper takes
`recompute` rather than being folded into `propagate_bounds` wholesale).
`dual_propagation` lines 223-252 use `excluding_or` with the
`Activity::compute(column, &scratch.y, Some(i))` closure.

**Why it is behavior-preserving.** `Extreme::excluding` already returns `None`
whenever `self.infinite != usize::from(!term.is_finite())` (activity.rs:95 and
97), so hoisting that test into `excluding_or` and evaluating it only on the
`None` path gives the same value as `bound_implied`'s early `continue`. The
fallback closure runs under exactly the same condition as today. Division and
the `is_finite` filter are unchanged and applied in the same order. Column
cursor iteration order is unchanged. `open_sides` evaluates sides in the order
Upper, Lower today; keep that array order.

One subtlety worth keeping in mind: `propagate_bounds` re-reads
`self.activity(i)` on every column entry (line 72) because `tighten_bound` may
have updated the cache mid-row. `row_implied_bound` reads the activity inside,
so that behavior is preserved automatically.

**Risk.** Zero behavioral risk by construction; the change is mechanical.

**Effort.** Two to three hours including a targeted unit test for
`excluding_or` in `activity.rs`, where `excluding` is already tested.

**Not proposed.** Folding `singleton_range`/`non_implied_bounds` into the same
helper. They deliberately bypass the cache (always `residual_activity`) and
produce both sides at once; unifying them would change which path computes the
residual and therefore could change rounding. Leave them.

---

## 3. Certificate constructors

**Files.** `src/core/rules/rows.rs:58-68` (`row_certificate`),
`src/core/rules/variables.rs:67-74` (`recession_certificate`), and ad-hoc
construction at `parallel.rs:264-270` (primal, two `y` entries),
`parallel.rs:349-355` (dual, two `x`), `dominated_columns.rs:244-250` (dual,
two `x`), `dual_propagation.rs:517-524` (dual, `x` from `lambda_touched`),
`cones.rs:80-87` and `cones.rs:103-108` (primal, `y` entries),
`dependencies.rs:174-181` (primal, `y` from `proof`).

**Problem.** Nine sites write `let mut point = Point::zeros(self.bounds.len(),
self.rows.len()); point.<x|y>[..] = ..; Certificate { mode: Recovery::<..>,
point }`. The two named helpers cover only the single-row and single-column
cases, so every rule that needs a two-entry certificate inlines the dance.

**Proposed shape.** Put two constructors on `Model` (they need `n` and `m`) in
`model.rs`, or on `Certificate` taking dimensions:

```rust
impl Model {
    /// Farkas certificate with multipliers `y` on the given rows and `z` on
    /// the given columns; everything else zero.
    pub(super) fn primal_certificate(
        &self,
        y: impl IntoIterator<Item = (usize, f64)>,
        z: impl IntoIterator<Item = (usize, f64)>,
    ) -> Certificate { .. }

    /// Recession direction with the given nonzero components.
    pub(super) fn dual_certificate(&self, x: impl IntoIterator<Item = (usize, f64)>) -> Certificate { .. }
}
```

Then `row_certificate` is
`self.primal_certificate([(row, m)], self.a.row(row).iter().map(|(j, a)| (j, -a * m)))`,
`recession_certificate` is `self.dual_certificate([(column, direction)])`, and
the seven inline sites become one expression each, e.g.
`self.dual_certificate([(dominant, 1.0), (dominated, -1.0)])` and
`self.primal_certificate(rows.iter().zip(dual).map(|(&i, v)| (i, -v)), [])`.

**Risk.** Zero behavioral risk. Every site is a terminal `return Err(..)`;
no arithmetic, ordering, or queue state is involved. The dense `Point::zeros`
allocation is unchanged.

**Effort.** One hour. Keep `row_certificate` and `recession_certificate` as
thin wrappers, since they are the names used by seven callers (five and two).

---

## 4. Run settings on `Model`; delete `Limits` and the `max_fill` parameters

**Files.** `src/core/schedule.rs:13-22` (`Limits`), `schedule.rs:66-76`
(`run`), every `limits.fill` / `limits.propagation` / `limits.sparsification` /
`limits.progress` / `limits.sparsify` use in `schedule.rs`; `src/core/model.rs:51-68`
(the settings fields); `src/presolve.rs:98-137` (the copy-out into both places).

**Problem.** Run configuration is split across two containers that are copies
of the same `Settings`:

- `Model` carries `rules`, `numerics`, `propagation`, `dual_propagation`,
  `dominated_columns`, `equalities`, `dependencies`, `allow_hessian_growth`,
  `fill`, `deadline` (model.rs:51-68).
- `Limits` carries `time`, `fill`, `sparsify`, `progress`, `propagation`,
  `sparsification` (schedule.rs:14-22).

Verified redundancies: `Limits.fill` is copied into `self.fill` at
schedule.rs:67 *before* `run_phases`, and `cleanup` already uses `self.fill`
(line 45), yet the other twelve call sites still thread `limits.fill` through
`singleton_columns`, `doubleton_equalities`, `short_equalities`
(lines 95-102, 123, 225-231, 249-255, 275-281). `Limits.sparsify` is
`settings.rules.sparsification` (presolve.rs:133), which is also
`self.rules.sparsification` (presolve.rs:98). `Limits.propagation` is the raw
`settings.propagation` while `self.propagation` is the gain-validated copy; the
schedule reads only `work_limit` and `additional_rounds` from it
(schedule.rs:202-206), which validation does not touch, so `self.propagation`
gives the same values. `Limits.sparsification` and `Limits.progress` have no
`Model` counterpart only because nobody added one; `equalities` and
`dependencies` are already on `Model` for the analogous rules.

Side effects of the split: `sparsify_rows(deadline, options: SparsificationSettings)`
takes its settings as a parameter while `short_equalities` and
`equality_dependencies` read `self.equalities` / `self.dependencies`; and the
`Model` fields `dominated_columns` and `dual_propagation` collide with the rule
methods of the same name, giving schedule.rs:132-133 the pair
`self.dominated_columns.general_search` / `self.dominated_columns()`.

**Proposed shape.** Store one `settings: Settings` on `Model` (all fields of
`Settings` are `Copy` types; it is `Clone`), applied once in `presolve.rs` with
the gain validation written into `settings.propagation` before storing:

```rust
// model.rs
pub(crate) struct Model {
    ..
    pub settings: crate::settings::Settings,   // replaces rules, numerics, propagation,
                                               // dual_propagation, dominated_columns,
                                               // equalities, dependencies, allow_hessian_growth
    pub fill: usize,                           // keep: it is settings.substitution_fill
    ..
}

// schedule.rs
pub fn run(&mut self, deadline: Instant, executor: &Executor) -> Result<Stats, Certificate>
```

`Limits` disappears; `run` takes the deadline (or reads `self.deadline`, see
note) and the schedule reads `self.settings.progress`,
`self.settings.rules.sparsification`, `self.settings.propagation`,
`self.settings.sparsification`. The rule entry points lose their `max_fill`
parameter: `empty_columns()`, `singleton_columns()`, `doubleton_equalities()`,
`short_equalities(deadline)` read `self.fill`. `sparsify_rows(deadline)` reads
`self.settings.sparsification`.

If storing the whole `Settings` feels too broad, a crate-private
`RunSettings` with the same fields minus `threads`/`time_limit` is equivalent;
the point is one container, not its exact type.

**Risk.** Zero behavioral risk for the settings and fill moves: every value is
a verified copy. The only place to be careful is time. Today `run_phases`
takes `start = Instant::now()` at schedule.rs:79 and uses
`start + limits.time` where `limits.time = time_limit - elapsed_so_far`, so the
schedule's deadline is (up to microseconds) the same instant as
`self.deadline`. Switching the schedule to `self.deadline` changes
`Option<Instant>` handling and the `checked_add`/`saturating_sub` edge cases
(a `Duration::MAX` limit would currently overflow-panic in `start + limits.time`
and would instead become "no deadline"). That cannot change reductions unless
the time limit is hit, but it is a semantic change and should land as its own
commit, after the pure plumbing.

**Effort.** Two to three hours; many lines, all mechanical. Compile errors
enumerate the sites.

---

## 5. One drain loop and one substitution trio in the schedule

**Files.** `src/core/schedule.rs:200-237` (`propagate_rounds`), `239-263`
(`sparsify_cleanup`), `265-289` (`substitution_cleanup`), `92-104` (fast
phase).

**Problem.** `sparsify_cleanup` and `substitution_cleanup` are the same 20-line
loop except that the former calls `propagate_bounds` after the first
`cleanup`. The trio "singleton_columns, doubleton_equalities, short_equalities"
with its three `if self.rules.*` guards appears verbatim four times (lines
94-103 with an extra `cleanup` between the first two, 224-232, 248-256,
274-282).

**Proposed shape.**

```rust
/// Singleton columns, doubleton equalities, short equalities, in that order.
fn substitutions(&mut self, deadline: Instant) {
    if self.rules.singleton_columns { self.singleton_columns(); }
    if self.rules.doubleton_equalities { self.doubleton_equalities(); }
    if self.rules.short_equalities { self.short_equalities(deadline); }
}

/// Drain consequences of a final-pass rule until the model stops changing.
fn drain(&mut self, deadline: Instant, propagate: bool) -> Result<(), Certificate> {
    while Instant::now() < deadline {
        let before = self.revision;
        self.cleanup()?;
        if propagate && self.rules.bound_propagation { self.propagate_bounds()?; }
        self.substitutions(deadline);
        self.cleanup()?;
        if self.revision == before { break; }
    }
    Ok(())
}
```

`propagate_rounds` lines 223-233 become `self.cleanup()?; self.substitutions(deadline); self.cleanup()?;`.
The fast phase keeps its own sequence because it has a `cleanup` between
`singleton_columns` and `doubleton_equalities` (line 97); it must not be
rewritten to call `substitutions` unless that intermediate cleanup is
preserved, since cleanup between them changes which doubletons exist when the
doubleton queue is drained.

**Risk.** Zero behavioral risk if the fast phase is left alone: the three
helpers call exactly the same rules in exactly the same order with the same
guards. The `propagate: bool` flag is the only difference between the two
loops today.

**Effort.** One hour. Do it after proposal 4 so the signatures are already
parameter-free.

---

## 6. One stats sink; uniform rule return types

**Files.** `src/core/model.rs:60` (`equality_stats: EqualityStats`),
`src/core/rules/substitution.rs:37-141` (ten inline increments and a
five-arm `match` on `SubstitutionFailure`), `src/core/schedule.rs:24-30`
(`Stats`), `schedule.rs:126-131` (`parallel_comparisons += ...`),
`src/core/rules/parallel.rs:163,290` and `dominated_columns.rs:277-280`
(return `Result<usize, _>` comparison counts), `src/presolve.rs:119-141,176`
(merging both channels).

**Problem.** Two channels carry statistics out of the run: `equality_stats`
lives on `Model` and is mutated in place; `fast_phases`, `medium_phases`,
`parallel_comparisons`, `time_limit` live in `schedule::Stats` and are returned
by `run`. Comparison counts are threaded through return values
(`parallel_columns` -> `dominated_support_groups` -> `+=`), so the two parallel
rules and the support-group scan return `Result<usize, Certificate>` purely to
feed a counter. The `SubstitutionFailure` to stats mapping (substitution.rs:128-141)
is a `match` that belongs next to `EqualityStats`.

Return types across the rule entry points today:

| Returns | Rules |
|---|---|
| `Result<usize, Certificate>` (count used by the schedule) | `propagate_bounds`, `equality_dependencies`, `dual_propagation` |
| `Result<usize, Certificate>` (count is a stat or ignored) | `parallel_rows`, `parallel_columns`, `dominated_columns` (ignored at schedule.rs:133) |
| `Result<(), Certificate>` | `empty_rows`, `singleton_rows`, `empty_columns`, `simple_dual_fix`, `simplify_cones` |
| `usize` (used by the schedule) | `sparsify_rows` |
| `()` | `fixed_variables`, `coupled_dual_fix`, `singleton_columns`, `doubleton_equalities`, `short_equalities` |

Note that `propagate_bounds`'s count is *not* replaceable by a revision check:
it deletes and relaxes rows (bumping `revision`) without counting them as
tightenings, and `propagate_rounds` stops on `tightened == 0` (schedule.rs:207).
So the counts that the schedule consumes must stay.

**Proposed shape.**

```rust
// model.rs (or a small stats.rs)
#[derive(Debug, Default)]
pub(crate) struct RunStats {
    pub equalities: crate::result::EqualityStats,
    pub parallel_comparisons: usize,
    pub fast_phases: usize,
    pub medium_phases: usize,
    pub time_limit: bool,
}
// Model { pub stats: RunStats, .. }

// result.rs
impl EqualityStats {
    pub(crate) fn record_rejection(&mut self, failure: SubstitutionFailure) { .. }  // the match at substitution.rs:129-141
}
```

`parallel_rows`, `parallel_columns`, `dominated_support_groups` add to
`self.stats.parallel_comparisons` directly and return
`Result<(), Certificate>`; `run` returns `Result<(), Certificate>` and
`presolve.rs` reads `model.stats` once. Convention to document in
`rules/mod.rs`: a rule returns `Result<usize, Certificate>` only when the
count is a scheduling signal (reductions applied), `Result<(), Certificate>`
when it can produce a certificate, `()` otherwise. Under that convention
`dominated_columns` becomes `Result<(), _>` (its count is discarded) and
`sparsify_rows` could stay `usize` since it cannot fail.

**Risk.** Zero behavioral risk; counters and return plumbing only.

**Effort.** One to two hours.

---

## 7. Split the long functions

Lengths, verified: `simplify_cones` cones.rs:51-295 (245 lines),
`dual_propagation` dual_propagation.rs:134-377 (244) plus `apply_direction`
381-535 (155), `sparsify_rows` sparsification.rs:80-241 (162),
`short_equalities` substitution.rs:14-147 (134), `run_phases` schedule.rs:78-195
(118), `propagate_bounds` bounds.rs:18-128 (111, largely fixed by proposal 2).
`parallel_rows` (64 lines) is not long; its only awkward part is the
"unmatched" canonical-order fallback (lines 186-223), which reads fine as-is.

**7a. `simplify_cones`.** The `match cone` at cones.rs:92-292 has arms of 6,
90, 16 and 85 lines using `continue` inside the arm (lines 119, 158, 224, 251).
Extract `fn simplify_soc(&mut self, block, rows: Vec<usize>) -> Result<(), Certificate>`,
`simplify_rsoc`, `simplify_psd`; `continue` becomes `return Ok(())`. Also the
constant-block classification (lines 70-91) is a natural `fn classify_constant_block`.
Risk: mechanical, but see "Tests": no test in this repository builds a
`SecondOrder` or `PositiveSemidefinite` cone, so only the benchmark corpus
would catch a slip. Do this arm-by-arm and rerun the corpus after each.

**7b. `dual_propagation`.** The function has two halves with a clear seam at
line 284: propagation over the dual rows (141-283, producing `consistent`,
`complete`, and the filled `scratch`), then extraction (287-374, looping rows
then columns). Split into `fn propagate_multipliers(&mut self, scratch: &mut DualScratch, budget: usize) -> (bool /*consistent*/, bool /*complete*/)`
and `fn extract_directions(&mut self, scratch: &mut DualScratch, budget: usize, complete: bool) -> Result<usize, Certificate>`.
The `(|| { .. })()` closure wrapping extraction exists only so `scratch` is
put back on the error path; with proposal 8 it goes away, and without it a
plain inner function achieves the same. The per-entry tightening at 221-281
shrinks to a third of its size with proposal 2. Risk: none if loop bodies are
moved verbatim; the round order (`swap_round`, `for at in 0..round.len()`)
must be kept exactly.

**7c. `short_equalities`.** Two extractions: candidate collection (lines
44-87) into `fn pivot_candidates(&mut self, i: usize, options: &EqualitySettings, relative: f64, out: &mut Vec<Candidate>)`
and the attempt loop (89-142) into `fn try_pivots(..) -> bool /*deferred*/`.
Replace the anonymous 5-tuple `(bool, usize, usize, usize, usize)` with

```rust
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Candidate { bounded: bool, score: usize, column: usize, degree: usize, cost: usize }
```

Derived `Ord` compares fields in declaration order, which is the tuple order,
so `select_nth_unstable`, `sort_unstable`, and the `candidate < candidates[0]`
comparison at line 76 produce the same order. Risk: none, provided the field
order is kept as above; the work accounting (`work`, `deferred`) is passed as
`&mut`.

**7d. `sparsify_rows`.** Extract the per-reference scan (129-167) as
`fn scan_cancellation_ratios(&self, reference, base, minimum_ratio, scan: &mut [Scan], candidates: &mut Vec<usize>, work: &mut usize, limit) -> bool /*exhausted*/`
and the evaluation/plan (168-211) as `fn plan_cancellations(..)`. Keep the
count-up `work` accounting untouched (see proposal 9). Risk: low; the
`break 'scan` becomes an early return that the caller checks.

**Effort.** Four to six hours for all four; each is independent.

---

## 8. Per-rule structs for the two scratch-carrying rules

**Files.** `dominated_columns.rs:30-53, 260-272, 281-320, 327-377`;
`dual_propagation.rs:31-57, 141, 287, 375`; `model.rs:56-57`.

**Problem.** `DualScratch` and `DominatedScratch` live as `Model` fields only
so they survive between passes. Each pass does `std::mem::take(&mut self.x_scratch)`,
runs the body inside an immediately-invoked closure so `?` cannot skip the
restore, then writes the scratch back (three occurrences). It is the one place
where the `impl Model` layout works against the code.

**Proposed shape.** Move the scratch into a rule struct owned by the presolver
(or by `Model` as today, but borrowed disjointly):

```rust
pub(crate) struct DualPropagation { scratch: DualScratch }
impl DualPropagation {
    pub fn run(&mut self, model: &mut Model) -> Result<usize, Certificate> { .. }
}
```

The take/closure/restore pattern disappears because `self.scratch` and `model`
are separate borrows. The rule body stays in the same file; only the receiver
changes. For the support-group entry point called from `parallel_columns`,
pass `&mut model.dominated` (or keep the field on `Model` and split the borrow
with `let Model { dominated_scratch, .. } = self` at the top, which also
removes the closure).

**Risk.** Zero behavioral risk; ownership only.

**Effort.** One to two hours. Lower priority than 1-6; worth doing when
touching these files for proposal 7b.

---

## 9. A `Budget` type for count-down work limits

**Files.** `dominated_columns.rs:267-271, 300, 310, 322, 332, 355`;
`dual_propagation.rs:154-158, 187-191, 293, 407-410`; `dependencies.rs:72-73,
83, 93-96, 103-106, 119-122`; `substitution.rs:21, 24, 97-100`;
`schedule.rs:202-222`; `sparsification.rs:95-105, 130-134, 176-179`.

**Problem.** Five rules resolve `WorkLimit` against a rule-specific default
and then decrement a `usize` with three slightly different idioms: `work -= 1`,
`work = work.saturating_sub(k)`, and `if cost > work { stop } else { work -= cost }`.
`sparsify_rows` alone counts *up* toward `limit`. `propagate_rounds` counts up
in a `try_fold` and then subtracts.

**Proposed shape.** A tiny `struct Budget(usize)` in `settings.rs` next to
`WorkLimit` with `exhausted()`, `charge(cost) -> bool` (fails without
charging when `cost > remaining`), and `charge_saturating(cost)`. Apply it to
the count-down users only.

**Risk.** Low for the count-down sites, since each idiom maps one-to-one.
Converting `sparsify_rows` is *not* recommended: its `work >= limit` checks
happen at different points than its increments (line 105 before the reference,
line 133 after the `work += 1` at 130), and while I convinced myself the
count-down equivalent is exact (including the `usize::MAX` unlimited case and
the `work == limit + 1` overshoot), it is the kind of change that buys nothing
and needs a proof. Leave it.

**Effort.** One hour. Small value; include only if the team wants one vocabulary.

---

## Module layout: is `impl Model` per family the right shape?

Yes. The honest assessment is that the current organization is correct for
this library and a trait or per-rule struct would not read better:

- Rules need intimate, mutable access to `Model` internals (the activity
  cache, `equations`, `queues`, `locks`, the tape). Every rule is a sequence
  of "inspect, decide, call a `Model` primitive" and the primitives
  (`fix`, `substitute`, `tighten_bound`, `replace_row`, `restrict_row`,
  `aggregate`) already enforce the invariants. That is the right boundary.
- The schedule is hand-ordered and rules have heterogeneous inputs
  (`executor`, `deadline`, counts). A `trait Rule { fn apply(&mut Model) }`
  would force a lowest-common-denominator signature and gain nothing, because
  nothing iterates over rules generically.
- One file per family matches the README's rule catalog and keeps related
  helpers (`fingerprint`, `proportional`, `exact_product`, `subtract`) private
  to their family.

The two exceptions are covered above: the scratch-carrying rules
(proposal 8) and the one duplicated predicate (proposal 1). A small
improvement to the existing layout: `rules/mod.rs` could carry the
return-type convention from proposal 6 and a one-line index of which schedule
phase each family belongs to, since that is the first question a reader asks.

## Tests

Unit tests exist where pure helpers exist: `dependencies.rs` (exact
arithmetic), `dual_propagation.rs` (`multiplier_range`, `used_bound`),
`parallel.rs` (fingerprint equivalence, executor threading),
`sparsification.rs` (`subtract`, `shifted`). Integration tests under `tests/`
cover aggressive, dependencies, dominated columns, dual propagation, parallel
structure, and quadratic elimination end to end. That split is sensible: rule
behavior is a property of the whole pipeline and belongs in `tests/`.

Gaps that matter for refactoring safety:

- **`cones.rs` has no test in the repository.** No file under `tests/` or
  `src/` constructs a `Cone::SecondOrder`, `Cone::RotatedSecondOrder`, or
  `Cone::PositiveSemidefinite` (checked with grep; the only `#[test]` in
  `presolve.rs` is about auxiliary columns). Proposal 7a should be preceded by
  a handful of small conic problems in `tests/cones.rs`, otherwise the corpus
  is the only safety net.
- `schedule.rs::significant_progress` has no unit test although it is pure
  and has NaN/clamp edge cases; a five-line test belongs beside it.
- `bounds.rs`, `substitution.rs`, `dual_fixing.rs`, `rows.rs`, `variables.rs`
  have no unit tests, but their logic is not separable from `Model`, so
  integration tests are the right home; nothing to move.

## Documentation errors found

Not refactors, but discovered while verifying defaults; both in
`src/settings.rs:302-307`:

- `DualPropagationSettings::work_limit` says "Default: twice the constraint
  nonzeros". The code resolves `self.a.nnz().saturating_mul(4)`
  (`dual_propagation.rs:154-157`) and the README table says "four times". The
  doc comment is wrong.
- The struct doc says "one dual propagation pass per medium phase". Since
  `803ecea` the schedule runs it once after the final passes
  (`schedule.rs:180-193`), which the README describes correctly.

All other `work_limit` doc defaults match their code (`equalities` x2,
`dominated_columns` x2, `dependencies` x4, `propagation` `max(nnz/4, 256)`,
`sparsification` `max(8*(nnz+bound entries), 1024)`).

## The three to do first

1. **Proposal 4 (settings on `Model`, delete `Limits` and `max_fill`).** It
   touches the most lines but every value is a verified copy, so it is the
   safest large win: it removes a second configuration container, four
   parameters from rule entry points, an inconsistent `sparsify_rows`
   signature, and the `dominated_columns`/`dual_propagation` field-versus-
   method collisions. Proposals 5 and 6 become smaller once it is in.

2. **Proposal 2 (residual-activity kernel), with proposal 1 folded in.** This
   is the only real *algorithmic* duplication in the area: five hand-written
   copies of the excluding-with-fallback computation that decides every implied
   bound. `Extreme::excluding_or` and `row_implied_bound` make `bound_implied`
   and `open_sides` one function with a flag, shrink `propagate_bounds` and
   `dual_propagation` by about 40 lines, and the equivalence argument is short
   enough to be checked line by line.

3. **Proposals 5 and 6 together (schedule drain loops and stats sink).**
   Both are an hour each, zero risk, and they are what makes the schedule read
   as the README describes it: one `substitutions` trio, one `drain` loop, and
   rule signatures that return a count only when the schedule uses it.

Proposal 3 (certificates) is a good warm-up commit for whoever starts, and
proposal 7a (`simplify_cones`) should wait for cone tests.
