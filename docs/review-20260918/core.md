# Refactoring review: working model and storage layer

Scope: `src/core/{model,activity,queues,objective,execution}.rs`, `src/matrix/{linked,sparse,csc,quadratic}.rs`, read completely; `src/core/rules/*.rs`, `src/core/schedule.rs`, `src/postsolve/tape.rs`, `src/presolve.rs` read for how the primitives are consumed. All line numbers refer to the working tree at commit `803ecea` (branch `new-presolve-rules`). `cargo test --lib` passes (36 tests) at this commit.

General impression: the layer is in good shape. Comments overwhelmingly explain *why* (e.g. `relax_bound` at `model.rs:379-382`, the lock-merge comment at `model.rs:290-295`, `replace_row_bounds` at `model.rs:682`). There is very little dead code. The problems are concentrated in three places: (1) the implied-bound arithmetic is written out by hand in three rule files instead of once next to `Activity`; (2) `Model`'s mutation primitives each re-spell the same three bookkeeping sequences (queue push, cached-activity shift, lock add/remove); (3) `Model` is configured by poking ten public fields from `presolve.rs`, and two of those settings exist in two copies.

Proposals are ranked by value divided by risk. "Bit-identical" below means: same floating-point operations in the same order, and same push order into every `Worklist`/`ActivityRows`.

---

## 1. One implied-bound helper on `Extreme`/`Activity` (six hand-written copies)

**Files/functions.**
- `src/core/rules/bounds.rs:70-125` (`propagate_bounds` inner loop, two copies: lines 74-86 and 87-98)
- `src/core/rules/bounds.rs:170-211` (`bound_implied`, lines 188-199)
- `src/core/rules/dual_propagation.rs:221-252` (`dual_propagation` inner loop, two copies: lines 223-237 and 238-252)
- `src/core/rules/dominated_columns.rs:88-126` (`open_sides`, lines 112-118)
- Helper would live in `src/core/activity.rs` next to `Extreme::excluding` (line 90).

**The duplication.** Every one of these sites computes "the bound on one variable implied by one row" as

```rust
extreme.excluding(term)
    .or_else(|| (extreme.infinite == usize::from(!term.is_finite()))
        .then(|| /* recompute residual excluding this column */)
        .flatten())
    .map(|v| (rhs - v) / a)
```

`propagate_bounds` (bounds.rs:74-98) and `dual_propagation` (dual_propagation.rs:223-252) are the same 25 lines twice each, differing only in whether the residual comes from `self.residual_activity(i, j)` or `Activity::compute(column, &scratch.y, Some(i))`. `bound_implied` (bounds.rs:188-199) writes the same thing with the `infinite` guard hoisted into a `continue`; that is equivalent, because `Extreme::excluding` (activity.rs:90-99) already returns `None` whenever `infinite != usize::from(!term.is_finite())`, so the guard only ever gates the recompute. `open_sides` (dominated_columns.rs:112-118) is the same expression with no recompute fallback. The comment explaining the guard ("Any other infinite contribution prevents a finite bound", bounds.rs:78) appears once out of six sites.

**Proposed shape.**

```rust
// activity.rs
impl Extreme {
    /// Bound implied for the excluded term's variable by `rhs`:
    /// `(rhs - residual) / a`. When the cached sum cannot exclude the term
    /// (cancellation), `recompute` supplies the residual directly, but only
    /// if this term is the sole infinite one: any other infinite term keeps
    /// the residual infinite.
    pub fn implied(self, term: f64, rhs: f64, a: f64,
                   recompute: impl FnOnce() -> Option<f64>) -> Option<f64> {
        self.excluding(term)
            .or_else(|| (self.infinite == usize::from(!term.is_finite()))
                .then(recompute).flatten())
            .map(|v| (rhs - v) / a)
    }
}

impl Activity {
    /// Both bounds one row implies for the variable with coefficient `a`
    /// and bounds `b`, as (from the row's lower side, from its upper side).
    /// Each side is skipped when its rhs is infinite.
    pub fn implied(self, a: f64, b: Bounds, rhs: Bounds,
                   residual: impl Fn() -> Activity) -> (Option<f64>, Option<f64>) {
        let (min_term, max_term) = Self::terms(a, b);
        let lower = rhs.lower.is_finite()
            .then(|| self.max.implied(max_term, rhs.lower, a, || residual().max.value()))
            .flatten();
        let upper = rhs.upper.is_finite()
            .then(|| self.min.implied(min_term, rhs.upper, a, || residual().min.value()))
            .flatten();
        (lower, upper)
    }
}
```

`propagate_bounds` becomes

```rust
let (lower, upper) = self.activity(i).implied(a, self.bounds[j], bounds,
                                              || self.residual_activity(i, j));
```

and the dual loop becomes `activity.implied(a, scratch.y[i], range, || Activity::compute(column, &scratch.y, Some(i)))`. Note `residual` may be called twice (once per side) in the cancellation case; today each side already recomputes independently (bounds.rs:81 and 93), so this is the same work. `bound_implied` and `open_sides` use `Extreme::implied` directly; they are then ~8 lines each and visibly the same test with/without recompute (see proposal 2b for merging them).

**Risk.** Very low. The arithmetic `(rhs - v) / a` is unchanged and evaluated in the same order; the recompute closure runs under the same condition; no queue is touched. The `is_finite` gate on the rhs is preserved by `.then()`. Watch for one thing: the closures must stay lazy (`then`, not `then_some`).

**Effort.** ~2 hours including reading each site against the helper.

---

## 2. `Bounds`/`Side` vocabulary: move `Side` next to `Bounds` and add the four helpers that are re-spelled everywhere

**Files/functions.**
- `Side` is defined in `src/postsolve/tape.rs:38-57` but imported by seven `core` files (`model.rs:13`, `rules/{bounds,rows,dual_fixing,dual_propagation,dominated_columns,parallel}.rs`); it is a bounds concept, not a tape concept.
- "set one side": `model.rs:385-388` (`relax_bound`), `413-417` (`tighten_bound`), `694-697` (`tighten_row`), `dual_propagation.rs:572-575` (`tighten_multiplier`).
- "opposite side": `bounds.rs:225-232` (`implied_bound`), `dual_propagation.rs:543-550` (`tighten_multiplier`), which are the same nine lines.
- "orient by sign" `if s > 0.0 { (b.lower, b.upper) } else { (b.upper, b.lower) }`: `activity.rs:115-119` (`Activity::terms`), `model.rs:733-737` (`aggregate`), `model.rs:555-559` (`substitute_equation`, with the orientation reversed), `rules/rows.rs:42-46`, `rules/parallel.rs:241-245`, `postsolve/tape.rs:419-423`.
- "insignificant gain" threshold `(factor * feasibility).max(relative * old.abs())`: `bounds.rs:251-255` and `dual_propagation.rs:565-569`, identical.

**Proposed shape.**

```rust
// problem/bounds.rs (Side moves here; tape.rs does `pub use crate::problem::Side`)
impl Side {
    pub fn opposite(self) -> Self { .. }
}
impl Bounds {
    pub fn side(self, side: Side) -> f64 { .. }          // replaces Side::value(bounds)
    pub fn set_side(&mut self, side: Side, value: f64) { .. }
    pub fn with_side(self, side: Side, value: f64) -> Self { .. }
    /// (lower, upper) for a positive multiplier, swapped for a negative one.
    pub fn oriented(self, sign: f64) -> (f64, f64) {
        if sign > 0.0 { (self.lower, self.upper) } else { (self.upper, self.lower) }
    }
}
// model.rs
impl Model {
    pub(super) fn insignificant_gain(&self, old: f64, gain: f64) -> bool {
        gain <= (self.propagation.minimum_gain_factor * self.numerics.feasibility)
            .max(self.propagation.minimum_relative_gain * old.abs())
    }
}
```

`Activity::terms` becomes `let (l, u) = b.oriented(a); (a * l, a * u)`; `aggregate` and `merge_parallel_rows` become `c.oriented(ratio)`; `rows.rs:42-46` becomes `let (l, u) = b.oriented(a); (l / a, u / a)`. `substitute_equation:555-559` wants the swapped orientation; write it as `effective_bounds.oriented(-pivot)` with a one-line comment, or keep it literal. `implied_bound` and `tighten_multiplier` both start with `let (lower, upper) = old.with_side(side, value).into()`.

Keep `Side::value(bounds)` as a thin alias if you prefer not to touch 20 call sites at once; the point is that the concept has one home.

**Risk.** Zero when the helper returns exactly the tuple the site built; `sign > 0.0` is the test used at every site (NaN and zero go to the "else" branch in all of them today). The `oriented(-pivot)` rewrite is the one site to double-check.

**Effort.** ~2 hours; mostly mechanical.

### 2b. `bound_implied` and `open_sides` are one function

`bounds.rs:170-211` and `dominated_columns.rs:88-126` walk a column, pick `(rhs, extreme, term)` by `from_lower = (side == Side::Lower) == (a > 0.0)`, and ask whether the implied bound is at least as tight as the variable's own. The only difference is the recompute fallback (`bound_implied` recomputes on cancellation; `open_sides` conservatively does not, per its doc comment at lines 86-87). After proposal 1 both reduce to a loop over `Extreme::implied`; give `bound_implied` a `recompute: bool` parameter (or a `Fallback` enum) and have `open_sides` call it. The `cursor()` walk (bounds.rs:175-176 and dominated_columns.rs:96-97) is needed in both because `self.activity(i)` takes `&mut self`.

Risk: zero for `bound_implied`; for `open_sides` the merged function must pass `|| None` as the fallback so it stays conservative.

---

## 3. Three bookkeeping helpers inside `Model` (queue push x10, cached-activity shift x3, lock delta x6)

**Files/functions.** All in `src/core/model.rs`.
- `self.queues.column_changed(j, self.a.column(j).len())` at lines 187, 223, 320, 333, 374-375, 448, 652, 678, 761, 800.
- The cached-activity shift
  ```rust
  if activity.min.infinite != STALE && !activity.replace_bound(a, old, new) { *activity = stale(); }
  ```
  at 354-357 (`set_bounds`), 391-394 (`relax_bound`), 461-464 (`fix`), and the `STALE`/`stale()` sentinel handling at 74-96, 255, 274-277.
- `set_bounds` (347-377) and `relax_bound` (383-397) share the same activity loop; `relax_bound` deliberately skips the queues (documented at 379-382).
- Lock add/remove: `from_parts:180-183`, `replace_rows_batch:211-216`, `replace_row:296-317`, `fix:455`, `replace_row_bounds:674-680`, `aggregate:757`.

**Proposed shape.**

```rust
impl Model {
    fn column_changed(&mut self, j: usize) {
        self.queues.column_changed(j, self.a.column(j).len());
    }

    /// Shift every cached activity of `column` from `old` to `new`, marking
    /// rows stale where the incremental update would lose precision.
    fn shift_activities(&mut self, column: usize, old: Bounds, new: Bounds) {
        for (i, a) in self.a.column(column) {
            self.activities[i].shift(a, old, new);
        }
    }
}
```

with the sentinel logic moved onto a small newtype so `STALE` stops leaking into three functions:

```rust
// activity.rs or model.rs
struct CachedActivity(Activity);   // min.infinite == STALE means "recompute"
impl CachedActivity {
    fn shift(&mut self, a: f64, old: Bounds, new: Bounds) {
        if !self.is_stale() && !self.0.replace_bound(a, old, new) { *self = Self::STALE; }
    }
}
```

`set_bounds` then reads: replace bounds; push fixed-column; `for (i, a) in column { activities[i].shift(..); changed_activities.push(i); match row_kinds[i] { .. } }`; `column_changed`. `relax_bound` is `replace; shift_activities; revision += 1`. `fix` keeps its "shift around `changed_row`" trick (lines 461-466), but as `let mut kept = self.activities[i]; kept.shift(a, bounds, Bounds::fixed(0.0)); self.changed_row(i); self.activities[i] = kept;`.

For locks, a `fn move_lock(&mut self, j: usize, a: f64, from: RowDomain, to: RowDomain)` covers `replace_rows_batch` (which does remove-all then add-all, lines 211-216, and must keep doing so; see "not proposed" below) and `replace_row_bounds:675-676`. `replace_row`'s merged walk (296-317) is the optimised path and can keep its shape (or use the merge iterator from proposal 5).

**Risk.** Zero: pure extraction, same loops, same push order. `set_bounds` must keep interleaving the `changed_activities.push(i)` and the kind-based pushes *inside* the same loop as the shift (as today at 354-372), so do not replace its loop with `shift_activities` plus a second loop; only `relax_bound` and `fix` use the standalone helper.

**Effort.** ~1 hour.

---

## 4. `RowKind` enum, computed once, shared by `changed_row` and `Queues::row_changed`

**Files/functions.** `model.rs:79-86` (`mod row_kind` with four `u8` constants), `model.rs:250-268` (`changed_row` classifies the row twice: once into `row_kinds[row]`, once implicitly via `queues.row_changed(row, length, equality)`), `queues.rs:162-171` (`row_changed` re-derives empty/singleton/doubleton/longer from `(size, equality)`), `model.rs:361-372` (`set_bounds` matches on the byte).

**Problem.** The eligibility classification lives in two places with two encodings. `row_kind::OTHER` conflates empty, singleton, and ranged rows, so `set_bounds` cannot distinguish them and `Queues::row_changed` must redo the length test. The `u8` constants plus `Vec<u8>` are an enum written by hand.

**Proposed shape.**

```rust
#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum RowKind { Deleted, Cone, Empty, Singleton, DoubletonEquality, LongerEquality, Other }

impl RowKind {
    fn of(domain: RowDomain, length: usize) -> Self { .. }   // the match at model.rs:259-264, extended
}

// queues.rs
pub fn row_changed(&mut self, row: usize, kind: RowKind) {
    match kind {
        RowKind::Empty => self.empty_rows.push(row),
        RowKind::Singleton => self.singleton_rows.push(row),
        RowKind::DoubletonEquality => self.doubleton_rows.push(row),
        RowKind::LongerEquality => self.short_equalities.push(row),
        _ => {}
    }
    self.changed_activities.push(row);
}
```

`changed_row` computes `let kind = RowKind::of(domain, length); self.row_kinds[row] = kind; if kind != RowKind::Deleted { self.queues.row_changed(row, kind) }`. `Vec<RowKind>` is still one byte per row (the doc comment at 42-45 stays true).

**Risk.** Zero if the mapping is table-checked: today `row_changed` pushes `empty_rows` for length 0 regardless of domain (including cone rows), `singleton_rows` for length 1 regardless of domain, `doubleton_rows` only for equality length 2, `short_equalities` only for equality length >= 3. `RowKind::of` must therefore rank `Empty`/`Singleton` by length *before* checking `Cone`, or keep a separate `is_cone` bit; and `set_bounds`'s `CONE` arm must still fire for a cone row of length 0 or 1. Simplest: keep `Cone` as a top-priority kind and let `row_changed` also take `length` for the 0/1 cases, or add `ConeEmpty`/`ConeSingleton`. Write the truth table in the test.

**Effort.** ~1.5 hours including a unit test of `RowKind::of`.

---

## 5. One sorted-merge iterator, and split `substitute_equation`

### 5a. `merge_sorted`

**Files/functions.** The same "two-pointer walk over two sorted `(index, f64)` sequences" appears in:
- `model.rs:593-634` (`substitute_equation`, row update)
- `objective.rs:120-161` (`try_substitute`, Hessian column vs slopes)
- `rules/sparsification.rs:43-77` (`subtract`)
- `rules/dependencies.rs:40-67` (`subtract`, index-based)
- `rules/dominated_columns.rs:135-158` (`dominates`, `match`-based)
- `model.rs:296-317` and `323-335` (`replace_row`, two index-based walks over `old` vs `entries`)

Five of them use the identical idiom `x.peek().map_or(usize::MAX, |e| e.0)`; two use indices. Each is 15-25 lines of cursor management around 3-5 lines of actual arithmetic.

**Proposed shape.** In `src/matrix/sparse.rs` (where `Entries` lives):

```rust
/// Walk two ascending-index sequences together, yielding every index in
/// either with the value each side holds for it.
pub(crate) fn merge_sorted<L, R>(left: L, right: R)
    -> impl Iterator<Item = (usize, Option<f64>, Option<f64>)>
where L: IntoIterator<Item = (usize, f64)>, R: IntoIterator<Item = (usize, f64)>
{ /* the peekable loop, once */ }
```

`sparsification::subtract` then reads

```rust
for (column, b, a) in merge_sorted(base, other) {
    let (a, b) = (a.unwrap_or(0.0), b.unwrap_or(0.0));
    let v = a - alpha * b;
    ...
}
```

and `substitute_equation`'s loop becomes

```rust
for (j, old, change) in merge_sorted(original, slopes.iter().copied()) {
    let Some(v) = change else { entries.push((j, old.unwrap())); continue };
    let old = old.unwrap_or(0.0);
    let change = a * v;
    let value = old + change;
    ... // unchanged cancellation / fill / push logic
}
```

`replace_row` collapses to one merge for locks and one for the "new-only" queue pass, keeping the two-pass push order exactly (old-support columns first in old order, then columns new to the row in entry order, as the comment at 290-295 requires).

Yield `Option<f64>` rather than `0.0` for the missing side: `objective.rs:138-144` adds `offset * curvature` only when the Hessian entry exists, and `-0.0`/`0.0` differences in `c[j]` could otherwise change `signum()` downstream (`variables.rs:61`, `dual_fixing.rs:50`). With `Option`, every site keeps its exact branch structure.

**Risk.** Low but not zero, and this is the hot path of substitution. The iteration order and the arithmetic per index are unchanged by construction, so reductions stay bit-identical; the risk is performance (the closure-free loops today are easy for LLVM to keep in registers). Benchmark before and after on the corpus; if `substitute_equation` regresses, keep its hand-written loop and apply the helper only to the four cold sites.

**Effort.** Half a day including the benchmark.

### 5b. Split `substitute_equation` (`model.rs:523-668`, 145 lines)

It does five things in sequence: (1) pivot/offset/slopes and the finiteness gate (532-546); (2) the retained-row domain from `effective_bounds` (553-571); (3) the per-row merge with deadline, fill, and cancellation checks (572-636); (4) the objective transaction (637-650); (5) commit: queues, `alive`, `replace_row` calls, tape (651-667). Also note `slopes` (538-543) and `remaining` (547-552) are two passes over the same filtered iterator.

Extract (2) as `fn retained_domain(equation: &Equation, pivot: f64, effective: Bounds) -> Option<RowDomain>` and (3) as `fn substituted_row(&self, i, a, offset, slopes, column, max_fill, fill: &mut usize) -> Result<(Entries, RowDomain), SubstitutionFailure>`. The deadline check (581-587) can stay in the caller's loop. Returning `Result<_, SubstitutionFailure>` also removes the `self.substitution_failure = ...` writes scattered through the body (531, 585, 627, 647): set it once in the caller from the `Err`.

**Risk.** Zero for behaviour; the `SubstitutionFailure::Numerical` reset at line 531 must still happen first, and `Numerical` must remain the reason for the early `return false`s at 535, 545, 563, 589 (today they rely on the reset).

**Effort.** ~2 hours.

---

## 6. Configure `Model` once, and stop carrying two copies of the settings

**Files/functions.**
- `Model` has seven `pub` settings fields plus `allow_hessian_growth`, `fill`, `deadline` (`model.rs:51-68`), all defaulted in `from_parts` (161-175) and then overwritten one by one from `presolve.rs:97-118`.
- `schedule::Limits` (`schedule.rs:14-22`) carries `fill`, `progress`, `propagation`, `sparsification`, `sparsify`; `run` copies `limits.fill` into `self.fill` (67) and then still reads `limits.fill` everywhere (95, 99, 102, 116, ...). `Limits.sparsify` is `settings.rules.sparsification` (`presolve.rs:133`), which `self.rules.sparsification` already holds.
- `Limits.propagation` is the **raw** `settings.propagation` (`presolve.rs:135`), while `model.propagation` is the **sanitised** copy (`presolve.rs:100-111`, `valid_gain`). They agree today only because `propagate_rounds` (`schedule.rs:202-206`) reads just `work_limit` and `additional_rounds`, which are not sanitised. That is a trap for the next person who adds a field.
- `Objective` is constructed as a struct literal at `problem/mod.rs:244-251` and eight times in `objective.rs` tests, each spelling `scratch: Default::default()`; `scratch` is `pub` only for that reason (`objective.rs:63`).

**Proposed shape.**

```rust
// model.rs
pub(crate) struct Config {           // everything a rule reads, sanitised once
    pub rules: Rules, pub numerics: Numerics, pub propagation: PropagationSettings,
    pub dual_propagation: .., pub dominated_columns: .., pub equalities: .., pub dependencies: ..,
    pub allow_hessian_growth: bool, pub fill: usize, pub progress: Progress,
    pub sparsification: SparsificationSettings,
}
impl Config { pub fn from_settings(s: &Settings) -> Self { /* valid_gain lives here */ } }

impl Model {
    pub fn from_parts(a, objective, rows, bounds, config: Config, deadline: Option<Instant>) -> Self
}
// schedule.rs: Limits shrinks to { time: Duration } or disappears; rules read self.config.*
// objective.rs
impl Objective { pub fn new(p: SymmetricMatrix, c: Vec<f64>, constant: f64) -> Self }
```

`presolve.rs:97-118` becomes one call. `Model` keeps `pub a/objective/rows/bounds/alive/locks/queues/postsolve/revision` (all read by `presolve.rs` for packing; encapsulating those is not worth it), but the settings block becomes one field.

Do **not** fold `self.deadline` (`Option<Instant>`, `presolve.rs:117`) together with the `start + limits.time` deadline the scheduler computes (`schedule.rs:102, 116, 164, 171, 184`) in this pass: they are the same instant only up to the microseconds between `start.elapsed()` and the scheduler's own `Instant::now()`, and time-limit behaviour on the corpus is the one thing here that is not deterministic. Note it as a follow-up.

**Risk.** Zero for reductions: every value the rules read is copied verbatim; the sanitised `propagation` replaces the raw copy only in `work_limit`/`additional_rounds`, which are identical in both. Deadline unification is explicitly excluded.

**Effort.** ~3 hours; touches `presolve.rs`, `schedule.rs`, every `limits.` read, and the tests that build an `Objective` literal.

---

## 7. Row deletion has one name and one implementation

**Files/functions.**
- `Model::delete_row` (`model.rs:341-345`) asserts `Linear`, calls `replace_row(row, &[], RowDomain::Deleted)`, pushes `Rule::DeletedRow`.
- `replace_row(i, &[], RowDomain::Deleted)` is spelled at `model.rs:343`, `model.rs:832`, `cones.rs:27`, `cones.rs:111`, `cones.rs:216`, `dependencies.rs:157`, `parallel.rs:281`.
- `cones::remove_cone_row` (`cones.rs:18-29`) is `delete_row` for a cone row, re-implemented because of the `Linear` assert.

**Proposed shape.**

```rust
impl Model {
    /// Remove the row's coefficients and retire its domain; no tape record.
    pub(super) fn clear_row(&mut self, row: usize) {
        debug_assert!(self.rows[row] != RowDomain::Deleted);
        self.replace_row(row, &[], RowDomain::Deleted);
    }
    /// Delete a row whose multiplier is simply zeroed on recovery.
    pub fn delete_row(&mut self, row: usize) {
        self.clear_row(row);
        self.postsolve.rules.push(Rule::DeletedRow(row));
    }
}
```

`remove_cone_row` becomes "push `ConeSlack` if needed; `self.delete_row(i)`". The five other sites that push their own record (`DependentRow`, `MergedRow`, `SocFace`, `PsdZeroFace`, `Unlocked`) call `clear_row`. The `Linear` assert in `delete_row` is dropped (its callers `empty_rows`, `singleton_rows`, `propagate_bounds` all already matched `Linear`).

**Risk.** Zero; same calls in the same order.

**Effort.** 30 minutes.

---

## 8. Tests in the wrong place, and the one module with none

- `src/core/model.rs` has **no** unit tests. `shifted` (`model.rs:100-122`) is tested from `rules/sparsification.rs:254`; that assertion belongs in `model.rs`. The lock-merge invariant of `replace_row` (290-295), the `fix` activity-warming trick (457-466), and the queue order of `replace_rows_batch` are exactly the kind of thing a small in-module test should pin, especially before proposals 3 and 5.
- `src/core/activity.rs:161-192`: `residual_activities_handle_infinite_terms_and_cancellation` ends by asserting `Locks::contribution` (182-191). Split into a `Locks` test.
- `src/matrix/csc.rs:267-303`: `mod arithmetic_tests` sits between `CscMatrix` and `CscMatrixRef` (defined at 307). Move to the end of the file with the other `#[cfg(test)]` helpers (214-265), or into a `tests` module.
- `src/matrix/quadratic.rs` has no tests; `into_csc` (80-104) and `column` (36-52) both implement the same "compact upper triangle of a stable-indexed `SymmetricMatrix`" filter and would benefit from one round-trip test against `working()`.

**Risk.** Zero. **Effort.** ~1 hour.

---

## 9. Smaller items (do opportunistically)

1. **Naming.** `Model::changed_row` (`model.rs:250`) vs `Queues::row_changed`/`column_changed` (`queues.rs:162,173`): pick one word order. `changed_cones` (`model.rs:71`) is a `Worklist` living outside `Queues` and named unlike `changed_activities`; move it into `Queues` as `cones`. `relax_row` (`model.rs:709-711`) is a public alias of the private `replace_row_bounds`; either make `replace_row_bounds` the public name or fold the comment into `relax_row` and delete the indirection.
2. **`Objective::substitute`** (`objective.rs:84-95`) is `try_substitute(..).ok()`; its only non-test caller is `fix` (`model.rs:441-443`). Delete it and write `.ok()` at that call site; tests use `.is_ok()`/`.is_err()`.
3. **`remove_unlocked`** (`model.rs:816-829`): the `.collect::<Vec<_>>().into_iter().map(..)` chain exists only to end the borrow of `self.a` before calling `self.equation(i)`. Two `let` statements say that; better, make `Rule::Unlocked { rows: Vec<Arc<Equation>> }` (`tape.rs:161-165`) and drop the `.as_ref().clone()` deep copy. Tape reads only `equation.row/coefficient/activity`, so recovery is unchanged.
4. **Dead lock bookkeeping.** `fix` (`model.rs:455`) and `aggregate` (`model.rs:757`) decrement `locks[column]` for the column being retired. Nothing reads a dead column's locks (`simple_dual_fix:15-16` and `coupled_dual_fix:78` test `alive` first). Harmless; keep or drop, but say which invariant you intend.
5. **`revision += 1`** appears 15 times. It is only ever compared for equality (`schedule.rs`, `presolve.rs:160,178`), so multiple increments per transaction are harmless; if proposal 3 lands, incrementing once in `changed_row`/`set_bounds` and once per primitive that bypasses them would be enough. Not worth doing alone.
6. **`Worklist::push` grows `queued` lazily** (`queues.rs:21-23`) so that `add_variable` (`model.rs:195-205`) need not touch `Queues`; `ActivityRows::push` does not grow. It works because rows are never added, but a `Queues::add_column()` next to the explicit `locks.push`/`elimination_rejected.push` would make the invariant visible.
7. **Types.** The `u32` boundary in `LinkedMatrix` (ids, `List.len`, `row_singletons`) is well contained behind `usize` accessors and asserted at `zeros:128` and `from_columns:174`. `dual_propagation::Record` (`dual_propagation.rs:25-29`) relies on that assert for its `as u32` casts (`267-269`); a `debug_assert` or `u32::try_from` there would document the dependency. `Locks { up, down }` as `usize` pairs cost 16 bytes per column and `Extreme.infinite` uses `usize` only to host the `STALE` sentinel; shrinking either is a performance experiment, not a readability change, so it is out of scope here.
8. **`csc.rs`** public helpers `is_symmetric`, `into_parts`, `to_owned`, `identity`, `values_mut` have no internal callers but are public API (`lib.rs:36`, `matrix/mod.rs:9`) and `values_mut` is used by the integration tests; not dead. `Quadratic::working` is used by `problem/mod.rs:247`.

---

## Not proposed: unifying `replace_rows_batch` with `replace_row`

They overlap on paper (`model.rs:207-234` vs `286-339`), but they are deliberately different: the batch path lets `LinkedMatrix::replace_rows` sort the updates and exploit column cursors (`linked.rs:536-545`), updates locks by remove-all/add-all instead of merging, and produces a different queue order (touched columns sorted and deduplicated, then `unlocked_columns` for every column of every changed row, then `changed_row` per row). Expressing one in terms of the other would change queue order and therefore reduction order on the corpus. Leave it; proposal 3's `move_lock` helper is enough to make the two visibly parallel.

---

## The three I would do first

1. **Proposal 1 (implied-bound helper).** Highest value per line: it takes the core numerical idea of the library out of three rule files and puts it next to `Extreme::excluding`, where the cancellation guard is explained once. Risk is essentially nil because the expression is moved verbatim. It also makes 2b a five-line change.
2. **Proposal 3 (Model bookkeeping helpers) together with proposal 7.** These are the cheapest changes with the widest readability payoff: every mutation primitive in `model.rs` shrinks, `STALE` stops appearing outside the cache, and "delete a row" has one name. Zero behavioural risk, an hour and a half total, and they make proposals 4 and 5 safer by first pinning the primitives with the in-module tests from proposal 8.
3. **Proposal 6 (single `Config`).** Not the biggest code reduction, but it removes the only latent correctness hazard I found: two copies of `PropagationSettings`, one sanitised and one not, chosen per call site. Doing it now, while `Limits` still has few fields, is much cheaper than after the next setting is added.

Proposals 2, 4, and 5a are worth doing next in that order; 5a is the only one that needs a benchmark run before merging.
