# Data-structure and memory-traffic pass on the presolve

Scope: `src/matrix/linked.rs`, `src/matrix/sparse.rs`, `src/core/model.rs`,
`src/core/queues.rs`, `src/core/execution.rs`, plus minimal call-site changes in
`rules/bounds.rs`, `rules/parallel.rs`, `rules/sparsification.rs`,
`rules/substitution.rs`. `rules/dual_propagation.rs` and `core/schedule.rs` are
untouched. No public API, README or default-setting changes. The full diff is
in `optimization.patch` (`git diff HEAD`, 9 files, +690/-185).

## Where the time went (baseline, samply + atos line attribution)

Profiles were taken with a scratch harness that loops `Presolver::presolve` on
one problem (it reproduces `benchmark run time --pool-mode reused` numbers).

* **BOYD2** (m=186k, n=93k, 78k bound sides removed): `set_bounds` 36% inclusive
  (`tighten_bound` 32%), of which `Worklist::push` 17% self. Every incident row
  of a changed column touched a scattered node (32 B), `Option<Activity>` (40 B),
  `RowDomain` (24 B), the list header (16 B) and two separate `queued` flag
  arrays. Parallel-row/column candidate sorting (`sort_unstable` on fingerprint
  tuples) another ~10%.
* **BOYD1** (18 rows x 93k columns, diagonal P): propagation 30%, sparsification
  scan 10%, parallel rules 8%, `Quadratic::working` 5% (93k one-`Vec`-per-row
  mallocs), `equation()` row snapshots 2% (18 x 31k entries copied per round),
  tape write/drop 4%.
* **CONT-300**: parallel rules 18% (sort ~10%), redundant-bound removal 17%,
  propagation 17%, `Activity::compute` 10%.
* **MAROS-R7**: singleton columns 29% (`Activity::compute` 19%), sparsification
  12%, `replace_row` 7% (full-row lock and queue updates for a one-entry change).
* **Typical netlib LPs** (10-problem profile): linked traversal (`Iter::next` +
  `Cursor::next`) ~31% inclusive, `Extreme::add` 5%, propagation 20%, dual
  propagation 20%, sparsification 16%, parallel rules 12%.
* **EXDATA** (dense 3000x3000 P): building the working Hessian ~55% of the run.

## Changes (all keep the exact sequence of reductions and numerics)

1. **Const-generic axis for `View`/`Iter`/`Cursor`** (`linked.rs`). The axis was
   a runtime field selected on every node; it is now a type parameter, so row and
   column traversal compile to straight-line loads. Same traversal order.
2. **`from_columns` generic visitor** (`linked.rs`). The two construction passes
   called a `&mut dyn FnMut` per nonzero; a nested generic `fn` lets both closures
   inline. Same node layout.
3. **`List.len: u32`** (`linked.rs`). List headers shrink from 16 to 12 bytes;
   `View::len()` still returns `usize`.
4. **Counting sort by column length** (`LinkedMatrix::sort_columns_by_length`,
   used by `remove_redundant_bounds`). Stable counting sort on ascending indices
   gives exactly the `(length, index)` order of the previous comparison sort.
5. **Radix sort for parallel candidate groups** (`execution.rs`, `parallel.rs`).
   `filter_map_sorted` now takes a 64-bit key prefix (the packed support and
   coefficient hashes, which order exactly like the `Ord` prefix), does four
   16-bit counting passes on `(key, position)`, gathers, then sorts each run of
   equal keys by the full `Ord`. The result is identical to `sort_unstable()` on
   the values (values are distinct because the index is part of them); a unit
   test checks this across the threshold. Inputs below 32 768 elements use the
   library sort because the 65536-bucket passes have a fixed cost. The parallel
   (`Executor::Parallel`) branch is unchanged.
6. **`ActivityRows`** (`queues.rs`). `changed_activities` and
   `singleton_activity_rows` received every push together (in `row_changed` and
   `set_bounds`); they are now one structure with two entry lists and one flag
   byte per row (two bits). Each list keeps its own first-push order and is drained
   independently (`take_round`, `take_singleton_round`, `iter`), so the rounds seen
   by propagation and the singleton scan are unchanged; one scattered access per
   incident row instead of two.
7. **`set_bounds` reads a per-row kind byte** (`model.rs`). `row_kinds` is
   refreshed in `changed_row`, which every row edit already ends in (verified for
   `replace_row`, `replace_rows_batch`, `fix`, `aggregate`, `remove_unlocked`,
   `replace_row_bounds`; the cone-to-cone block rewrite in `cones.rs` keeps kind
   `CONE` and the block is still read from `rows[i]`). Queue pushes and their
   order per list are unchanged.
8. **Stale sentinel instead of `Option<Activity>`** (`model.rs`). The cache is
   32 bytes per row instead of 40; `activity()` and `set_bounds` behave exactly
   as before (`replace_bound` is still applied incrementally to cached rows).
9. **`replace_row` merge-walk and scratch reuse** (`model.rs`, `linked.rs`).
   The old row is written into a model-owned scratch buffer instead of a fresh
   `Vec`. Locks of a retained column with an unchanged contribution are left alone
   (the previous code subtracted and re-added the same integers). Queue pushes:
   old-support columns first, in order, then only columns new to the row; a
   retained column was already queued by the first pass with the same column
   length, so the skipped pushes were guaranteed no-ops.
10. **Row snapshot cache for the tape** (`model.rs`). `equation()` caches the
    `Arc<Equation>` per row and `changed_row` invalidates it, so repeated
    propagation rounds over an unchanged row share one copy. `substitute_on_side`
    builds its own snapshot (it needs modified bounds) so the shared `Arc` is never
    cloned by `make_mut`. `implied_bound`'s proof closure takes `&mut Self`.
11. **`take_round` pre-sizes the next buffer** (`queues.rs`) instead of regrowing
    from empty (measured neutral within noise; kept as it removes reallocations).
12. **Arena-backed `SymmetricMatrix`** (`sparse.rs`). One `Vec<(usize, f64)>`
    with a `(start, len, capacity)` slot per variable replaces one `Vec` per row.
    Construction is two counting passes and one allocation; row contents and
    order (lower entries, diagonal, upper entries) are identical to before. Inserts
    shift within a slot, or relocate the row to the end with doubled capacity (as a
    `Vec` would reallocate); compaction runs when dead space covers half the arena
    and only changes slot positions. A randomized test checks every row against a
    dense reference through growth, removal and compaction.
13. **AoS scan record in `sparsify_rows`** (`sparsification.rs`): `seen`,
    `ratio`, `count` for a candidate row sit on one cache line. Pure layout change.

## Measurements

All timings are single-threaded (`threads: 1`), reused presolver, min of
repeated in-process runs, base and candidate binaries alternated (A/B/A/B).
Two other engineers were benchmarking on the machine at the same time, so
expect a few percent of noise on individual problems; a bisect of batch one
showed build-to-build layout noise of about +-2% on problems under 25 ms.

Interleaved A/B, 13 slowest problems (15 iterations x 2 rounds, min):

| problem | base ms | new ms | change |
|---|---|---|---|
| maros-meszaros/BOYD1 | 91.22 | 77.68 | -14.8% |
| maros-meszaros/BOYD2 | 65.58 | 53.13 | -19.0% |
| maros-meszaros/CONT-300 | 47.66 | 42.26 | -11.3% |
| maros-meszaros/CONT-200 | 19.53 | 16.70 | -14.5% |
| maros-meszaros/EXDATA | 9.58 | 7.63 | -20.4% |
| netlib/MAROS-R7 | 24.12 | 21.76 | -9.8% |
| netlib/PILOT87 | 19.15 | 18.05 | -5.7% |
| netlib/FIT2P | 18.32 | 17.15 | -6.4% |
| netlib/FIT2D | 16.60 | 14.88 | -10.4% |
| netlib/PILOT | 13.70 | 12.37 | -9.7% |
| netlib/STOCFOR3 | 13.36 | 12.37 | -7.4% |
| netlib/QAP15 | 13.08 | 11.43 | -12.6% |
| netlib/DFL001 | 13.26 | 12.38 | -6.6% |
| total | 365.16 | 317.79 | -13.0% |

Interleaved A/B over the whole corpus (236 problems, 5 iterations x 2 rounds,
min per problem, measured before the radix threshold was raised to 32 768):

* total 680.6 ms -> 611.9 ms (**-10.1%**); netlib -8.8%, Maros-Meszaros -11.1%.
* 172 problems faster by more than 3%, 48 within +-3%, 16 slower by more than
  3%. The slower ones are small (0.2-4.8 ms): the LISWET family (+4-7% at
  ~1.4 ms), CVXQP1_L/CVXQP3_L, QSHELL, STADAT3. Their sizes (10k rows/columns)
  fell into the radix path whose four 65536-bucket passes cost more than
  `sort_unstable` at that size; the threshold was raised from 4096 to 32 768
  afterwards (the sort result is identical by construction and unit-tested), but
  there was no time to re-measure the corpus after that change.
* The 24 problems above 5 ms: 466.5 ms -> 410.9 ms (-11.9%).

`benchmark run time --trials 5 --pool-mode reused` was also run for both
trees (`benchmark/results/agent/time/base-time.jsonl`, `cand-time.jsonl`), but
the candidate run coincided with the other engineers' benchmarks (for example
BOYD2 median 83 ms with p10-p90 of 60-182 ms against a quiet baseline of 67.6
ms), so `compare time base-time cand-time` is not meaningful; the interleaved
numbers above are the honest measurement.

Sort microbenchmark (186k random `((u32,u32),usize)` keys): `sort_unstable`
5.8 ms, packed-u64 `sort_unstable` 4.1 ms, stable radix 2.2 ms.

## Verification

* `cargo run --release -p benchmark -- run size --name base ...` was taken before
  any edit; `compare size base <candidate>` was run after each batch
  (`cand1`, `cand2`, `cand3`) and on the final tree (`final`):
  **236 matched cases; 0 changed** (variables, rows, nonzeros, bound sides and
  outcome identical for every problem).
* `cargo test --release`: all passing (35 lib tests including 4 new ones, plus
  the integration tests). `cargo clippy --all-targets --release`: no warnings.
  `cargo fmt --all --check`: clean.

## Ideas tried or considered and rejected

* **Pooling the per-row `Entries` in `substitute_equation`/`fix`**: allocation
  is a few percent of a substitution's cost (the sorted merge dominates); not
  worth the plumbing.
* **A shared push log for the two activity consumers** (instead of two entry
  lists with shared flags): a column pushed after one consumer drained would
  have to appear in the log twice, which breaks the first-push order the other
  consumer must see; the two-list/one-flag design was used instead.
* **Software prefetch in column traversal**: needs inline asm or unstable
  intrinsics on aarch64; not acceptable for a portable stable-Rust crate.
* **Uninitialised node arena in `from_columns`** (skipping the `vec![Node; nnz]`
  fill, ~2% on the largest problems): requires `unsafe`, which the crate does
  not use.
* **SoA node layout**: would cut row-traversal bytes but turn each scattered
  column step into three cache misses; column traversal already dominates on
  the LP set.
* **Boxing the large `Rule` variants** (`Rule` is 80 bytes; tape write + drop
  is ~4% on BOYD1): estimated under 2% and touches every tape match; skipped.
* **Merging the column worklists' flags** like `ActivityRows`: estimated ~1%
  on substitution-heavy problems; skipped for time.
* **Pre-sizing `take_round`**: measured neutral (within +-2% noise); kept only
  because it removes reallocations.
* The bisect variants (reverting the merge-walk, the `set_bounds`
  restructure or the `take_round` change one at a time) all landed within the
  +-2% build noise of each other on MAROS-R7/FIT2P/PILOT87, so none of the
  batch-one changes regress anything measurable; the large effects came from
  the const-generic axis, the shared activity flags, the kind byte, the radix
  sort and the Hessian arena.
