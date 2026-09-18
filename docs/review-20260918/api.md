# Refactoring review: public surface and core plumbing

Scope: `src/lib.rs`, `src/presolve.rs`, `src/result.rs`, `src/settings.rs` (shape only),
`src/problem/{mod,bounds,cone,conic}.rs`, `src/postsolve/{mod,solution,tape}.rs`,
`benchmark/src/*.rs`, `tests/*.rs`, `README.md`. Every file was read in full; every
line reference below was checked against the working tree at commit `803ecea`
(branch `new-presolve-rules`). Nothing under the repository was modified.

Conventions: effort S = under two hours, M = half a day, L = one to two days.
"Behavior risk" is the risk of changing observable results; "API risk" is the risk
of breaking downstream code that uses the public surface.

## Ranking (value divided by risk, best first)

| # | Proposal | Value | Behavior risk | API risk | Effort |
| --- | --- | --- | --- | --- | --- |
| 1 | Delete the duplicate `original_rows`, centralize index translation and the conic-dual sign in one place | High | Very low | None | S |
| 2 | Collapse `Matrix::{Vacant, Moved}` into `Linked { matrix: Option<LinkedMatrix>, .. }` | Medium | Very low | None | S |
| 3 | Split `presolve_owned` and `pack` into model construction, finishing, survivor selection and compaction | High | Low | None | M |
| 4 | Shared `tests/common` module: runner, KKT/Farkas checkers, point constructors | Medium | None | None | M |
| 5 | Benchmark: send `Tuning` to the worker as one JSON payload; derive CLI names from `ValueEnum`; generic `--without` | Medium | None (library) | None | S–M |
| 6 | Recovery tape: merge identical arms, `Recovery::gradient`, `Model::point()`; do not introduce a trait per variant | Low–Medium | Low | None | S |
| 7 | README: remove the settings table and benchmark CLI text that duplicate rustdoc and have already drifted | Medium | None | None | S–M |
| 8 | `into_conic` has no tests or callers; `Solution`/`Certificate` type family small fixes | Medium (coverage) | None | Additive only | S |

Section 9 lists public API inconsistencies that are not refactors on their own.

---

## 1. Delete `original_rows` and centralize index translation

**Files and functions**

- `src/presolve.rs:258-268` `original_rows` — a verbatim copy of
  `src/problem/mod.rs:115-125` `row_indices`, which `presolve.rs` already imports
  (`presolve.rs:13`).
- `src/presolve.rs:146-147` calls `original_rows(&problem.rows)` to obtain the
  linear/conic row positions of the input. But `problem.linear_rows` and
  `problem.conic_rows` (computed by `Problem::from`, `problem/mod.rs:96`) are
  still intact at that point: `presolve_owned` only reads `problem.rows`
  (`presolve.rs:78-85`), and the two vectors are first `take`n inside `pack`
  (`presolve.rs:297-298`). So the call recomputes data the struct already holds.
- `src/postsolve/mod.rs:38-61` `original_point` (stable -> original, allocating),
  `:71-84` `scatter` (compact -> stable), `:85-89` `lift`, `:120-148` `recover_into`
  (compact -> stable, then stable -> original into caller buffers),
  `:177-219` `reduce_warm_start` (original -> stable, then stable -> compact).

**Problem**

Four translation directions (compact<->stable, stable<->original) are written
inline in five places. The convention that a conic dual is stored as `-y[i]` in
the working point appears at `mod.rs:47`, `:82`, `:143`, `:185`, `:211`; a sign
mistake in any one silently breaks one entry point and not the others. The
`x`/`z` prefix truncation (`original_columns`) is repeated at `:45-46`, `:137-138`,
`:179-180`.

**Proposed shape**

```rust
// postsolve/mod.rs
impl Coordinates {
    /// Compact reduced coordinates into the stable working point.
    fn scatter(&self, p: &mut Point, s: SolutionRef<'_>) { /* body of today's scatter */ }
    /// Stable working point back to compact coordinates (used by reduce_warm_start).
    fn gather(&self, p: &Point) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) { ... }
}

/// The caller's original layout: a column prefix plus input row positions.
pub(crate) struct OriginalMap<'a> { pub columns: usize, pub linear: &'a [usize], pub conic: &'a [usize] }
impl OriginalMap<'_> {
    fn scatter(&self, p: &mut Point, s: SolutionRef<'_>)                       // reduce_warm_start input
    fn gather_into(&self, p: &Point, slacks: &[f64], out: SolutionMut<'_>)     // recover_into output
    fn gather(&self, p: Point, slacks: Vec<f64>) -> Solution                    // replaces original_point
}
impl Postsolve {
    fn original(&self) -> OriginalMap<'_> {
        OriginalMap { columns: self.original_columns, linear: &self.input_linear, conic: &self.input_conic }
    }
}
```

`presolve.rs:146-147` becomes
`OriginalMap { columns: n, linear: &problem.linear_rows, conic: &problem.conic_rows }.gather(certificate.point, Vec::new())`
and `original_rows` is deleted. `recover_primal_certificate` (`mod.rs:161-176`) then
gathers straight into a `PrimalCertificate` instead of building a `Solution` and
discarding `x`/`conic_slack`. Keep `original_point`'s allocation-free fast path
(`mod.rs:48-53`, truncate `y` when there are no cones) inside `gather`.

**Risk**: pure code motion; the existing unit test
`recovery_preserves_coordinate_order_and_dual_signs_with_reused_buffers`
(`mod.rs:227-276`) pins every sign and ordering. No public API change.

**Effort**: S.

---

## 2. Collapse `Matrix::{Vacant, Moved}` into an `Option` inside `Linked`

**Files and functions**

- `src/problem/mod.rs:78-93` `enum Matrix { Vacant, Moved{..}, Csc, Linked{..} }`.
- `:254-270` `restore_working_matrix`, `:271-323` `working_matrix`.
- `Vacant` is constructed only as the placeholder in two `std::mem::replace`
  calls (`:260`, `:294`) and matched as unreachable in three places (`:187`, `:218`,
  `:305`). `Moved` exists only between `working_matrix` and
  `restore_working_matrix` and carries the same two maps as `Linked` minus the
  matrix.

**Problem**

Two of four variants encode "a `Linked` whose matrix is temporarily lent to the
`Model`". Each transition needs a `mem::replace` plus a `let ... else { unreachable!() }`
destructure (`:256-263`, `:290-297`), and every reader must list `Moved | Vacant`
as unreachable.

**Proposed shape**

```rust
pub(crate) enum Matrix {
    Csc(CscMatrix),
    Linked {
        /// None while the working model owns the storage (see working_matrix).
        matrix: Option<LinkedMatrix>,
        compact_to_stable_rows: Vec<usize>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    },
}

pub(crate) fn working_matrix(&mut self) -> LinkedMatrix {
    let (m, n) = (self.row_count(), self.variable_count());
    match &mut self.a {
        Matrix::Linked { matrix, compact_to_stable_rows, stable_to_compact_columns }
            if is_identity(compact_to_stable_rows, m) && is_identity(stable_to_compact_columns, n) =>
        {
            matrix.take().expect("working storage is lent at most once")
        }
        Matrix::Csc(a) => LinkedMatrix::from_columns(a.rows(), a.columns(), |j| a.as_ref().column(j)),
        Matrix::Linked { matrix, compact_to_stable_rows, stable_to_compact_columns } => {
            let a = pack_rows(/* as today, lines 314-319 */);
            LinkedMatrix::from_columns(a.rows(), a.columns(), |j| a.as_ref().column(j))
        }
    }
}

pub(crate) fn restore_working_matrix(&mut self, m: LinkedMatrix) {
    if let Matrix::Linked { matrix: slot @ None, .. } = &mut self.a { *slot = Some(m); }
}
```

`matrix_row` and `into_csc` use `matrix.as_ref().expect("editable input belongs to the working model")`
(same message as today's `unreachable!`). `ConstraintMatrix::from_columns` and
`pack` wrap the matrix in `Some`.

**Risk**: none observable. `Outcome::Unchanged` keeps returning the input
allocations exactly as now (the doc at `result.rs:59` promises this).
Lines of code drop by roughly 25 and the `unreachable!` count from 4 to 0 in
this file.

**Effort**: S.

---

## 3. Split `presolve_owned` and `pack`

**Files and functions**

- `src/presolve.rs:71-194` `presolve_owned` (124 lines).
- `src/presolve.rs:270-411` `pack` (142 lines).

**Problems, each verified**

1. Model configuration (`:97-118`) copies ten `Settings` fields onto `Model`
   fields one by one and sanitizes two propagation gains inline (`:101-111`).
   Both `Model::from_parts` (`core/model.rs:139-193`) and this block know the
   full list of configuration fields, so adding a setting means editing both.
   The gain sanitization implements a documented `PropagationSettings` contract
   (`settings.rs:261-266`: "other values use the default") but lives in
   `presolve.rs`.
2. The `Unchanged` path is duplicated verbatim: `:160-166` and `:178-184`
   (restore matrix, restore `c`, restore bounds unless free, return
   `Outcome::Unchanged`). The `stats.after`/`pack` pairing around it is also
   duplicated (`:159,168` vs `:186-187`).
3. `Limits` (`:130-137`) is a second copy of settings fields, assembled here
   rather than by the type that owns it.
4. `pack` interleaves five concerns: surviving column map (`:271-276`), surviving
   row/cone selection with block renumbering (`:277-296` and again `:374-388`),
   `Postsolve` assembly including `direct_slacks` (`:297-330`), the `Solved`
   shortcut (`:332-340`), and objective/bounds/constraint compaction into a
   `Problem` (`:341-410`). `domains` is declared at `:280` and first used at
   `:372`. The `stable_to_compact_columns` inverse map (`:272-276`) is the same
   construction as `ConstraintMatrix::from_columns` (`problem/mod.rs:59`) with a
   non-identity source.

**Proposed shape**

```rust
// core/model.rs — the model owns knowledge of its own configuration fields.
impl Model {
    pub fn configure(&mut self, settings: &Settings, deadline: Option<Instant>) {
        self.rules = settings.rules;
        self.numerics = settings.numerics;
        self.propagation = settings.propagation.sanitized();   // moves valid_gain into settings.rs
        self.dual_propagation = settings.dual_propagation;
        /* dominated_columns, equalities, dependencies, allow_hessian_growth */
        self.deadline = deadline;
    }
}
// core/schedule.rs
impl Limits { pub fn new(settings: &Settings, start: Instant) -> Self { ... } }

// presolve.rs
fn presolve_owned(mut problem: Problem, settings: &Settings, executor: &Executor, start: Instant) -> PresolveResult {
    let (mut model, free_bounds) = working_model(&mut problem, settings, start);
    let before = size(&model);
    let mut stats = Stats::new(before);                         // or keep the literal
    let phases = model.run(Limits::new(settings, start), executor);
    let outcome = match phases {
        Ok(phases) => { /* copy counters */ finish(model, problem, free_bounds, before, &mut stats) }
        Err(certificate) => match certificate.mode {
            Recovery::PrimalInfeasibility => { stats.after = None; infeasible(certificate.point, &problem, n) }
            Recovery::DualInfeasibility => match feasible_point(&model) {
                Some(point) => { stats.after = None; unbounded(point, certificate.point, ...) }
                None => finish(model, problem, free_bounds, before, &mut stats),
            },
            Recovery::Solution => unreachable!(),
        },
    };
    /* elapsed, time_limit */
}

/// Either hand the input back untouched or pack the reduced model.
fn finish(model: Model, problem: Problem, free_bounds: bool, before: Size, stats: &mut Stats) -> Outcome {
    stats.after = Some(size(&model));      // equals `before` when revision == 0, so both paths agree
    if model.revision == 0 { unchanged(model, problem, free_bounds) } else { pack(model, problem, before) }
}
fn unchanged(model: Model, mut problem: Problem, free_bounds: bool) -> Outcome { /* lines 161-166 */ }

// pack, split by concern
struct Survivors { columns: Vec<usize>, linear_rows: Vec<usize>, conic_rows: Vec<usize>, cones: Vec<Cone> }
fn survivors(model: &Model) -> Survivors                              // lines 271, 277-296
fn inverse(map: &[usize], len: usize) -> Vec<usize>                   // lines 272-275; reuse in ConstraintMatrix::from_columns
fn build_postsolve(model_tape, s: Survivors, original: &mut Problem, before: Size) -> Postsolve   // 297-330
fn compact_problem(model: Model, original: Problem, postsolve: &Postsolve, stable_to_compact: Arc<Vec<usize>>, cones: Vec<Cone>, before: Size) -> Problem  // 341-409
```

Note on `finish`: in the current code the Ok/revision==0 path leaves
`stats.after = Some(before)` (set at `:125`) while the DualInfeasibility/revision==0
path sets `Some(size(&model))` (`:159`). With `revision == 0` no edit has
happened, so both are the same value; the unified helper is behavior-preserving.

**Risk**: low. All moves are of contiguous blocks; the only semantic merge is the
`stats.after` observation above. The `#[cfg(test)]` test at `presolve.rs:418-452`
calls `pack` directly and keeps working if `pack` keeps its signature as the
composition of the three helpers. No public API change (`Settings` gains a
`pub(crate) fn sanitized` at most).

**Effort**: M. Do proposal 1 first; it removes `original_rows` and shrinks the
certificate branch.

---

## 4. Shared test-support module

**Files**: `tests/aggressive.rs`, `dependencies.rs`, `dominated_columns.rs`,
`dual_propagation.rs`, `parallel.rs`, `quadratic_elimination.rs`.

**Duplication, counted**

- `Presolver::new(s).unwrap().presolve(presolve::Problem::from(p.clone()))`:
  44 occurrences across the six files (17 in `aggressive.rs` alone).
- `Solution { x, y, z, conic_dual: vec![], conic_slack: vec![] }` literals: 14.
- `fn column(&CscMatrix, j)` is written twice (`aggressive.rs:59-62`,
  `parallel.rs:257-264`) and both re-implement the public
  `CscMatrixRef::column` (`src/matrix/csc.rs:332`) reachable as
  `a.as_ref().column(j)`.
- Three optimality checkers with overlapping content: `aggressive.rs:64-104`
  `check_kkt` (full primal/dual feasibility, complementarity, objective value),
  `dependencies.rs:63-71` `stationarity` (gradient only, hard-codes `P = I`),
  `parallel.rs:117-150` inline Farkas and recession-ray checks. Every future
  rule test will want at least one of these.

**Proposed shape**

```rust
// tests/common/mod.rs   (each test file: `mod common;`)
#![allow(dead_code)]
pub fn run(settings: Settings, data: &ProblemData) -> PresolveResult {
    Presolver::new(settings).unwrap().presolve(Problem::from(data.clone()))
}
pub fn linear_point(x: Vec<f64>, y: Vec<f64>, z: Vec<f64>) -> Solution { Solution { x, y, z, conic_dual: vec![], conic_slack: vec![] } }
pub fn column(a: &CscMatrix, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ { a.as_ref().column(j) }
pub fn check_kkt(data: &ProblemData, s: &Solution) -> f64        // moved from aggressive.rs, tolerance as a parameter
pub fn check_farkas(data: &ProblemData, c: &PrimalCertificate)   // moved from parallel.rs:122-137
pub fn check_ray(data: &ProblemData, ray: &[f64])                 // moved from parallel.rs:142-149
```

`dependencies::stationarity` becomes `check_kkt` with the `P = I` assumption
dropped (it is a strict generalization: `check_kkt` computes `P x` from `data.p`).

**Risk**: none to the library. Because each integration test is its own crate,
unused helpers warn per binary; hence the `allow(dead_code)`.

**Effort**: M (mostly mechanical). This makes every other proposal cheaper to
verify, so it is worth doing early even though it changes no library code.

---

## 5. Benchmark: forward `Tuning` as one payload

**Files and functions**

- `benchmark/src/main.rs:67-93` `struct Tuning` (10 clap flags).
- `benchmark/src/run.rs:210-282` `trial`: rebuilds the argument list by hand
  (`:225-271`), one tuple per flag, formatting each `Option` back to a string.
- `main.rs:40-49` `Preset::as_str`, `:57-65` `SparsificationMode::as_str`,
  `:95-102` `Kind::as_str`, and `run.rs:221-224` a `match` for `PoolMode`: all
  hand-maintained copies of the names clap already derives with `ValueEnum`
  (`to_possible_value().unwrap().get_name()`). `Preset::as_str` and
  `SparsificationMode::as_str` are used only by `trial` (`run.rs:225`, `:261`).
- `run.rs:103-115` `--without` accepts eight hard-coded rule names although the
  `rules!` macro (`run.rs:15-43`) already enumerates all eighteen.

**Problem**: adding a tuning flag requires editing the struct, the `trial`
forwarding table, and the README paragraph at `README.md:564-575`; forgetting the
second silently runs the worker with defaults for that flag, and nothing would
catch it (the worker is spawned from `current_exe`, so no test covers the
round trip).

**Proposed shape**

```rust
// main.rs
#[derive(Clone, Debug, Default, Args, Serialize, Deserialize)]
struct Tuning { ... }                       // Preset, SparsificationMode gain Serialize/Deserialize too
#[derive(Serialize, Deserialize)]
struct WorkerRequest { path: PathBuf, rule: String, threads: usize, pool_mode: PoolMode, tuning: Tuning }
// Command::Worker becomes flagless: the request arrives as JSON on stdin.

// run.rs
fn trial(executable: &Path, request: &WorkerRequest) -> Result<Measurement> {
    let mut child = Command::new(executable).arg("worker").stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped()).spawn()?;
    serde_json::to_writer(child.stdin.take().unwrap(), request)?;
    let output = child.wait_with_output()?;
    /* status check and from_slice as today */
}

// rules! macro also emits the setter used by --without
macro_rules! rules { ($($field:ident),+ $(,)?) => {
    const RULES: [(&str, Rules); 19] = [ ... ];
    fn disable(rules: &mut Rules, name: &str) -> bool {
        match name { $(stringify!($field) => { rules.$field = false; true },)+ _ => false }
    }
}}
```

This deletes `Preset::as_str`, `SparsificationMode::as_str`, the `PoolMode`
match, and the 45-line forwarding table.

**Risk**: none to the library. Behavior notes for the benchmark: the hidden
`worker` subcommand changes shape (it is `#[command(hide = true)]` and only ever
invoked by the parent, `main.rs:131-142`); the generic `--without` would accept
all eighteen rule names instead of eight. If the restriction was intentional,
keep a `const DISABLEABLE: &[&str]` filter, but nothing in the code or README
explains it.

**Effort**: S–M.

---

## 6. Recovery tape: small cleanups, not a trait per variant

**Files**: `src/postsolve/tape.rs:215-438` `recover_with_slacks`, `:443-551`
`reduce_point`.

**Assessment**: the two matches over `Rule` are not duplication. They are the
reverse and forward transforms of each rule; they share structure but almost no
bodies (compare `Substituted` at `:288-317` vs `:478-497`, or `TightenedBound`
at `:328-355` vs `:498-505`). A trait with `recover`/`reduce` per variant struct
would touch the sixteen construction sites in `src/core` (`model.rs` 7,
`rules/cones.rs` 8 including four `ConeSlack`, `parallel.rs`, `dependencies.rs`,
`sparsification.rs` 1 each) and lose the exhaustive match that today forces a
new rule to implement both directions. Keep the enum and the two matches.

**Real, verified repetition worth removing**

1. `Rule::DependentRow { row, .. }` (`:318-322`) and
   `Rule::DeletedRow(row) | Rule::MergedRow { removed: row, .. }` (`:323-327`)
   have identical bodies; merge into one three-pattern arm.
2. `let g = if mode == Recovery::Solution { gradient.evaluate(&point.x) } else { 0.0 }`
   appears in `Fixed` (`:279-283`) and `Substituted` (`:302-306`). Add
   `fn gradient(self, g: &Gradient, x: &[f64]) -> f64` on `Recovery` next to
   `offset`/`bounds` (`:26-35`), which already encode the same "only in
   `Solution` mode" idea.
3. `Point::zeros(self.bounds.len(), self.rows.len())` is written at nine call
   sites in `src/core/rules/*` (`cones.rs:80,103`, `variables.rs:68`,
   `rows.rs:59`, `dependencies.rs:174`, `parallel.rs:264,349`,
   `dominated_columns.rs:244`, `dual_propagation.rs:517` with local `n, m`) plus
   `presolve.rs:199`. A `Model::point(&self) -> Point` (or
   `Certificate::new(mode, &Model)`) removes the dimension pairing from every
   rule. Cross-cutting with the core reviewer's area; flagging it here because
   `Point` is a tape type.
4. `Certificate.mode: Recovery` admits `Recovery::Solution`, which
   `presolve.rs:172` has to mark `unreachable!()`. A two-variant
   `CertificateKind { PrimalInfeasibility, DualInfeasibility }` with
   `impl From<CertificateKind> for Recovery` makes the impossible state
   unrepresentable. Optional; touches the same nine sites as item 3.

An optional larger step, if the per-arm `if !mode.dual() { continue; }` guards
(eight of them, `:219`, `:226`, `:245`, `:254`, `:262`, `:334`, `:362`, and the
`mode.primal()`/`mode.dual()` pairs) are judged noisy: split the reverse pass
into `recover_primal(rule, point, mode)` and `recover_dual(rule, point, mode)`
called in that order per rule (the dual step of `Fixed`/`Substituted` reads the
`x` the primal step just wrote, so per-rule sequencing must be kept). The cost
is either `_ => ()` arms (losing exhaustiveness) or listing every variant twice.
I would not do this now.

`tape.rs` carries the PSLP-derived SPDX header (`:1-2`, see `NOTICE`); any
restructuring should keep it.

**Risk**: low; the five unit tests in `tape.rs:555-732` cover the touched arms.
**Effort**: S.

---

## 7. README: stop duplicating rustdoc, and move benchmark and internals text

**Verified drift (the reason this matters)**: README `:446` says the dual
propagation work limit default "is four times the constraint nonzeros";
`settings.rs:305` says "twice the constraint nonzeros"; the code
(`src/core/rules/dual_propagation.rs:156-157`) uses `nnz * 4`. The README is
right and the rustdoc is wrong, which is exactly what happens when the same
table is kept in two places.

**Sections that duplicate or belong elsewhere**

- `:428-455` settings table: every row restates a field doc in `settings.rs`.
  Replace with a pointer to `Settings` rustdoc (and fix `settings.rs:305`).
- `:564-587` benchmark CLI flags and invocation: belongs in
  `benchmark/README.md` (there is none) or the clap `--help` text, which
  already documents the flags (`main.rs:69-93`, `:161-180`).
- `:538-562` parallel-execution internals (radix sort tie-breaking, the 1,024 /
  32,768 scan thresholds, "no global pool cache"): implementation notes for
  `src/core/execution.rs` module docs, not user guidance. The user-facing
  content (`threads` semantics, `InitError`, `Presolver` reuse, `:544-555`) is
  already in `settings.rs:24-29` and `presolve.rs:25-27,45-46`.
- `:589-624` "Results and postsolve" restates `result.rs` and `lib.rs:12-15`.
- The README has no end-to-end usage example (the only code blocks are
  settings literals); the `lib.rs:18-33` doctest is the obvious candidate.
- Minor: `lib.rs:16-17` has a stray empty `//!` line splitting the crate doc
  paragraph from the example.

**Proposed shape**: README = overview, problem model, quick-start example,
rule catalog with the "Source:" links (this is the part rustdoc cannot express
well), aggressive preset paragraph, pointers to `Settings`, `Outcome`,
`Postsolve`, `into_conic`, and a one-line pointer to `benchmark/README.md`.

**Risk**: none. **Effort**: S–M.

---

## 8. `into_conic` coverage; `Solution` type family

**`into_conic` has no callers and no tests.** `grep` for `into_conic`,
`ConicMap`, `ConicExport`, `ConicData` finds only `src/problem/conic.rs` and the
re-export at `problem/mod.rs:9`; there is no doctest, unit test, integration
test or benchmark use. That is 158 lines (`conic.rs:119-277`) with two distinct
packing paths (CSC scatter `:192-228`, linked `pack_rows` `:229-257`), an
identity fast path (`:166-191`), and three multiplier maps (`:64-110`). Before any
refactor in this file, add a test that exports a small ranged problem through
both storage kinds and checks `map.dual(...)` / `warm_dual` round trips; the
`Coordinates` unit test in `postsolve/mod.rs:227-276` is a good template. Both
paths are needed for performance (`Problem::matrix_row` on CSC input scans every
column per row, `problem/mod.rs:191-195`), so unification is not the goal.

**Solution / SolutionRef / SolutionMut / PrimalCertificate / CertificateRef**
(`postsolve/solution.rs`, 70 lines). The five types are small and each has a
distinct ownership role, so I would not merge them. Two cheap improvements:

- `pack` builds `SolutionRef { x: &[], y: &[], z: &[], conic_dual: &[], conic_slack: &[] }`
  (`presolve.rs:333-339`); a `SolutionRef::EMPTY` const (or `Default`) documents
  the intent.
- `ConicMap::dual` (`conic.rs:64`) returns a `PrimalCertificate` for ordinary
  optimal multipliers, and `recover_primal_certificate` builds a `Solution` only
  to strip it (`postsolve/mod.rs:164-175`). Adding
  `impl From<Solution> for PrimalCertificate` and `impl Solution { fn certificate(&self) -> CertificateRef }`
  is additive. Renaming `PrimalCertificate` to something sign-neutral
  (`Multipliers`) would be clearer but is a breaking change; a doc line on
  `ConicMap::dual` is the non-breaking alternative.

**ProblemData / Problem / ConstraintMatrix / ConicData**: coherent enough to
keep. `ProblemData<M>` is generic only so that `ConstraintMatrix::from_columns`
can feed `Problem::from` without a CSC copy; `ConstraintMatrix` is a one-method
newtype around the private `Matrix`. If proposal 2 lands, `ConstraintMatrix`
becomes `Matrix::Linked { matrix: Some(..), identity maps }` and the `inverse`
helper from proposal 3 replaces `(0..columns).collect()` duplication. Simplifying
the generic away would break `data.a.values_mut()` used by tests
(`aggressive.rs:110`, `:442`) and the benchmark loader, so leave it.

**Risk**: additive only. **Effort**: S for the test, S for the additions.

---

## 9. Public API inconsistencies (no refactor proposed unless noted)

1. **Entry point ergonomics.** `Presolver::presolve(&self, Problem)` forces
   `presolve::Problem::from(data)` at every call (44 times in `tests/`,
   `benchmark/src/run.rs:136,150`). `pub fn presolve(&self, problem: impl Into<Problem>)`
   is source-compatible for all existing callers and removes the boilerplate.
2. **Borrowed vs owned inputs.** `Postsolve::{recover_into, recover_solution, recover_primal_certificate, reduce_warm_start}`
   take `SolutionRef`/`CertificateRef`; `ConicMap::warm_dual` takes `&Solution`
   (`conic.rs:87`) and `ConicMap::dual` takes a bare `&[f64]`. `warm_dual` should
   take `SolutionRef<'_>` (breaking, but the type has no callers in the
   workspace).
3. **Allocation pairs.** `recover_into` (buffers) / `recover_solution` (allocates)
   and `dual_into` / `dual` are paired; `recover_primal_ray`,
   `recover_primal_certificate`, `reduce_warm_start`, `warm_dual` allocate with
   no buffer variant. Either document that only the hot path
   (`recover_into`) is allocation-free, or add `_into` variants uniformly.
4. **Naming.** `Postsolve::solution()` (`mod.rs:103`) allocates a zeroed
   output; the name reads as an accessor. `zeroed_solution()` or
   `Solution::zeros_for(&Postsolve)` would be clearer (breaking).
   `Problem::row(i)` (`mod.rs:198`) means "the i-th *linear* row" while
   `constraints()` and `matrix_row(i)` are indexed by shared row position; and
   `variable_bounds(j)`/`row_bounds(i)`/`conic_rhs(i)` are per-index accessors
   next to slice accessors `c()`, `constraints()`, `cones()`. `variable_bounds`
   cannot be a slice because empty means free, but `row_bounds` could be
   replaced by iterating `constraints()`.
5. **Exhaustiveness.** `Cone` (`cone.rs:7`) and `MatrixError` (`csc.rs:5`) are
   `#[non_exhaustive]`; `Outcome`, `Constraint`, `Progress`, `WorkLimit` are not,
   so adding an outcome or a constraint kind is a breaking change while adding
   a cone is not. Decide one policy; `Outcome` is the one most likely to grow.
6. **Missing docs / derives.** `Size`, `Stats`, `EqualityStats`, `PresolveResult`,
   `ReducedProblem`, `UnboundednessCertificate` (`result.rs`) have field docs
   but no type-level docs; `Stats` has no `Default` although `presolve.rs:120-128`
   builds one by hand.
7. **Error types.** `InitError` and `MatrixError` are both `thiserror` types and
   `Problem::from` cannot fail by design (`problem/mod.rs:22-28`). Consistent.
   `Presolver::new` returns `Result` while `Default` cannot fail; documented at
   `presolve.rs:34`. Fine.
8. **Executor setting placement.** `Settings.threads` (`settings.rs:24-29`) is
   the only field that affects `Presolver::new` rather than a call; benchmark
   metadata has to zero it out to compare algorithm settings
   (`run.rs:341-349`). A `Presolver::new(settings, threads)` split would be
   cleaner but is breaking; note only.

---

## Do these three first

1. **Proposal 1** (delete `original_rows`, centralize translation): it removes a
   verbatim duplicate function and an unnecessary recomputation, puts the
   `-y` conic-dual convention in one place, and is covered by an existing unit
   test. Under two hours, no API change.
2. **Proposal 2** (`Option<LinkedMatrix>` in `Matrix::Linked`): two variants and
   four `unreachable!` arms disappear from a public type's private core, the
   lend/restore protocol becomes a `take`/`Some`, and it shrinks the `Unchanged`
   path that proposal 3 extracts. Under two hours, no observable change.
3. **Proposal 3** (split `presolve_owned` / `pack`): the biggest readability win
   in the file every contributor reads first; it moves model configuration next
   to the model, deletes the duplicated `Unchanged` block, and gives `pack` named
   phases. Do it after 1 and 2 so the blocks being moved are already smaller.

Proposal 4 (shared test support) costs nothing in risk and should be started in
parallel: every proposal above is verified by the integration tests, and a
`run`/`check_kkt` helper makes adding a regression test for each refactor a
one-liner.
