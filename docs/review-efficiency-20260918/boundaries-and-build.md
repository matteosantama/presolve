# Efficiency review: entry/exit boundaries, postsolve tape, build configuration

Branch `new-presolve-rules`, commit 2528782. Scope: `src/presolve.rs`,
`LinkedMatrix::from_columns`/`pack`, `SymmetricMatrix::from_upper_columns`/
`pack_upper`, `CscMatrix`, `Problem::take_objective`, `Model::from_parts`,
the recovery tape and `Postsolve`, and `Cargo.toml`.

No file in the main tree was modified. All measurements come from an
instrumented copy under the session scratchpad
(`scratchpad/presolve`, phase timers in `presolve.rs`/`model/mod.rs`, a
`benchmark/src/bin/phases.rs` driver, a counting global allocator test, and
five release builds of the benchmark binary under different profiles). The
in-process timings below were taken while other agents were using the
machine; they are minimums over repeated calls, so the *shares* are
trustworthy to a few percent but absolute numbers are not corpus numbers.

## Summary

| Rank | Finding | Est. corpus impact | Risk to reductions | Effort |
| --- | --- | --- | --- | --- |
| 1 | `Model::from_parts` walks the arena through the linked cursor for locks, and initialises queues one row/column at a time (`src/model/mod.rs:179-191`) | -1.5% to -2% of corpus sum; -2% to -4% on LISWET/AUG2D/CONT/QSHIP families; geomean ~0.98 | none (integer sums, identical push order) | low |
| 2 | Tape records are 80 B and each `TightenedBound` clones an `Arc<Equation>` (`src/model/tape.rs:84-176`, `src/model/mod.rs:482`) | -0.5% to -1% of corpus sum; -2% to -3% on BOYD1, ~-1% on BOYD2/QAP15/CONT-200; also -1.5 ms of `drop(result)` on BOYD1 (outside the timed region) | none (bookkeeping only) | medium |
| 3 | Build profile: keep the reference benchmark on the default profile; add `#[inline]` to the hot leaf functions that the 16-CGU build leaves out of line; treat `codegen-units = 1` as a measurement-stability option, not a user-facing gain | 0 to -3% (needs the A/B described in section "Compiler settings") | none | low |
| 4 | `LinkedMatrix::from_columns` fills 32 B per nonzero with a non-zero `Node` before overwriting every slot (`src/matrix/linked.rs:175-185`) | -0.3% to -0.5%; ~-0.7% on BOYD1/2, CONT-300 | none | low, but needs `bytemuck` or `unsafe` |
| 5 | `size()` is an O(m) walk run twice per call; `before` is computable in O(1) without cones and `after` can ride on the `survivors()` walk (`src/presolve.rs:223-243`) | -0.2% (up to -0.65% on BOYD2/CONT-300) | none, but the reported `Size` must stay byte-identical (`run size` checks it) | low |
| 6 | Exit path does seven separate O(n)/O(m) passes and keeps an identity `input_linear` map (`src/presolve.rs:253-399`, `src/problem/mod.rs:51-61`) | -0.2% to -0.3% | none | low-medium |
| 7 | `a.pack` scatter uses `usize` maps; `next`/`stable_to_compact` could be `u32` (`src/matrix/linked.rs:234-270`) | -0.1% to -0.3% (pack is 2.3% of corpus) | none | low |
| 8 | Hessian copy is bandwidth/page-fault bound (`src/matrix/sparse.rs:52-95`); EXDATA spends 67% of its time here | only via a SoA arena (out of scope); otherwise nothing safe to gain | n/a | high |
| 9 | Nine separate `queued` flag arrays and a non-zero `stale()` fill in `Queues::new`/`from_parts` (`src/model/queues.rs:147-159`, `src/model/mod.rs:166`) | -0.1% to -0.2%, mostly page-fault/first-touch cost on m > 50k | none | low-medium |
| 10 | `Coordinates.compact_to_stable_columns` is an `Arc<Vec<usize>>` that is never cloned (`src/postsolve/mod.rs:15`) | none measurable; cleanup | none | trivial |

The boundary code (working model in, packed problem out) is **13% of the
corpus time** (78 ms of 600 ms in-process; median share per problem 12%,
mean share 16% for the 141 problems under 1 ms, which drive the geometric
mean). Cutting it by 30% would be about -4% on the sum and 0.95 on the
geomean; findings 1, 2, 4, 5, 6 together are a realistic 25-35% cut.

## Where the boundary time goes

In-process, min of 3 calls per problem, all 236 problems, default settings,
`Presolver` reused (so no thread-pool cost; `threads = 1` never builds a
pool anyway, `src/executor.rs:12-14`):

| Phase | ms | % of corpus |
| --- | ---: | ---: |
| `working_model` total | 55.5 | 9.2 |
| &nbsp;&nbsp;`LinkedMatrix::from_columns` | 21.6 | 3.6 |
| &nbsp;&nbsp;`take_objective` (`from_upper_columns`) | 10.7 | 1.8 |
| &nbsp;&nbsp;`Model::from_parts` | 21.2 | 3.5 |
| &nbsp;&nbsp;&nbsp;&nbsp;lock sweep over `a.row(i)` | 10.1 | 1.7 |
| &nbsp;&nbsp;&nbsp;&nbsp;`changed_row(i)` loop | 6.3 | 1.0 |
| &nbsp;&nbsp;&nbsp;&nbsp;`column_changed(j)` loop | 3.3 | 0.6 |
| &nbsp;&nbsp;&nbsp;&nbsp;allocations | 1.5 | 0.2 |
| `size()` (both calls) | 1.0 | 0.2 |
| `survivors` + `build_postsolve` | 2.1 | 0.4 |
| `compact_problem` | 19.5 | 3.3 |
| &nbsp;&nbsp;`LinkedMatrix::pack` | 13.9 | 2.3 |
| &nbsp;&nbsp;`pack_upper` | 1.3 | 0.2 |
| `model.run` | 522.2 | 87.0 |

Per problem (quiet-ish machine, min of 9; the totals match the corpus
medians within a few percent):

| Problem | total ms | from_columns | take_objective | from_parts | a.pack | boundary share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BOYD1 (18 x 93k, nnz 559k) | 76.0 | 2.5 | 0.6 | 1.7 | 2.4 | 9.5% |
| BOYD2 (186k x 93k, nnz 424k) | 45.0 | 2.1 | 0.4 | 2.7 | 2.0 | 17% |
| CONT-300 (90k x 90k, nnz 449k) | 33.3 | 2.1 | 0.4 | 2.0 | 1.3 | 16% |
| CONT-200 | 15.8 | 0.9 | 0.2 | 0.9 | 0.6 | 18% |
| EXDATA (dense 3000 x 3000 P) | 7.2 | 0.04 | 5.0 | 0.07 | - | 67% |
| MAROS-R7 | 15.1 | 0.6 | 0.0 | 0.4 | 0.3 | 9% |
| FIT2D | 12.0 | 0.5 | 0.0 | 0.3 | 0.5 | 11% |
| QAP15 | 8.6 | 0.5 | 0.0 | 0.3 | 0.3 | 12% |
| STOCFOR3 | 8.2 | 0.3 | 0.0 | 0.3 | 0.2 | 10% |
| LISWET1 (10k x 10k) | 0.96 | 0.13 | 0.05 | 0.17 | - | 21% |
| AUG2DC | 1.98 | 0.18 | 0.11 | 0.23 | - | 16% |

Fixed overhead on a 2 x 1 problem (HS21-like): 51 allocations, ~1.1 KB, 4.0 µs
warm in-process. The corpus number for HS21 is 15 µs per fresh process, so
for problems under ~50 µs the timing is cold-start (page faults, cold
caches), not the presolve. Nothing in this scope changes that except smaller
code (see the build section).

Type sizes (release, `size_of`): `Rule` 80 B, `Equation` 48 B, `Node` 32 B,
`List` 12 B, `Activity` 32 B, `Locks` 16 B, `RowDomain` 24 B,
`Constraint` 24 B, `Bounds` 16 B, `Slot` 24 B, `Model` 1992 B (moved by
value four times per call; ~100 ns, irrelevant), `Settings` 296 B with no
heap (the `clone` in `configure` is free).

## Findings

### 1. `Model::from_parts`: lock sweep through the linked cursor, per-row queue initialisation

`src/model/mod.rs:179-191`:

```rust
for i in 0..m {
    for (j, a) in model.a.row(i) {
        model.locks[j].add(Locks::contribution(a, model.rows[i]));
    }
    model.changed_row(i);
}
for j in 0..n {
    model.column_changed(j);
    ...
}
```

What happens now. Right after `from_columns` the arena is laid out row by
row, so `a.row(i)` reads consecutive memory, but it does so through
`Cursor::next` (`src/matrix/linked.rs:60-70`): each step loads
`nodes[self.next]` and the next id comes from that load, a dependent-load
chain at about 2.5 ns per nonzero (1.4 ms on BOYD1, 1.1 ms on BOYD2,
0.85 ms on CONT-300). The `changed_row` loop then costs 7-10 ns per row
(1.37 ms for BOYD2's 186k rows, 0.88 ms for CONT-300): it rewrites
`activities[i] = stale()` and `equations[i] = None`, which `from_parts`
just initialised to exactly those values (lines 166, 169), computes the
kind byte, and pushes into `changed_activities` (two `Vec::push`es with
doubling growth from empty, plus a flag byte) and into one size-class
worklist. The column loop pushes every column into `unlocked_columns`,
again growing from empty. Corpus: lock sweep 1.7%, row loop 1.0%, column
loop 0.6%.

Change.

* Locks: iterate the arena in storage order instead of following links.
  Add `LinkedMatrix::for_each_entry(&self, f: impl FnMut(row, col, value))`
  that walks `self.nodes` by index (assert `self.free == NONE`, which is
  true right after construction, or skip nodes whose `row == NONE` if the
  free list is ever non-empty). Nodes are in row order, so `rows[row]`
  reads are monotone and `locks[col]` is the only scattered access. Lock
  contributions are integer counts, so any order gives identical values.
  Expected 2.5 ns/nnz -> about 1 ns/nnz: -0.8 ms BOYD1, -0.6 ms BOYD2,
  -0.5 ms CONT-300, about -1% corpus.
* Row initialisation: replace the `changed_row(i)` loop with a dedicated
  initialisation that (a) computes `row_kinds[i]` in a tight loop,
  (b) fills `changed_activities` in bulk (`propagation = (0..m).collect()`,
  `singleton = (0..m).collect()`, `flags = vec![PROPAGATION | SINGLETON;
  m]`), (c) pushes the size-class queues in the same `i` order as today
  with `Vec::with_capacity` sized from a first counting pass or simply
  reserved to `m`. Push order per list is unchanged, so every rule sees the
  same rounds. Expected 7 ns/row -> 2-3 ns/row: -0.8 ms BOYD2, -0.5 ms
  CONT-300, -4% on the twelve LISWET problems, -2% on AUG2D*.
* Column loop: `unlocked_columns.entries = (0..n).collect()`,
  `queued = vec![true; n]`, then the singleton/empty/fixed pushes in `j`
  order. Small (-0.3%).

Bit-identity: locks identical; queue contents and order identical; the
`changed_cones` pushes for cone rows must be kept (they are, since
`set_cones` runs afterwards and re-pushes anyway, but keep the order).
Effort: low; one afternoon including a unit test that compares the queue
state after `from_parts` against the old loop on a random model.

### 2. Tape records: 80 B each, one `Arc` clone and one `Arc` drop per `TightenedBound`

`src/model/tape.rs:84-176` (`Rule`), `src/model/mod.rs:257-269`
(`equation()` cache), `src/model/mod.rs:462-490` (`tighten_bound`).

What happens now. `Rule` is 80 B because `Substituted` carries
`{column, Arc<Equation>, Gradient (32 B), Entries (24 B), bool}`. The
dominant record on the corpus is `TightenedBound {column, Arc<Equation>,
Side, old}` which needs 32 B of payload but occupies 80. Counts from the
instrumented run:

| Problem | rules | of which TightenedBound | distinct equations | tape bytes |
| --- | ---: | ---: | ---: | ---: |
| BOYD1 | 244 376 | 244 369 | 17 | 19.5 MB |
| BOYD2 | 43 982 | 43 982 | 55 | 3.5 MB |
| CONT-200 | 39 601 | 39 601 | 39 601 | 3.2 MB + 5.7 MB of equations |
| QAP15 | 22 275 | 22 275 | 1 590 | 1.8 MB |
| STOCFOR3 | 23 858 | 19 447 | 14 197 | 1.9 MB + 1.8 MB |
| CONT-300 | 14 950 | 13 156 | 5 382 | 1.2 MB |
| PILOT87 | 12 548 | 12 202 | 1 161 | 1.0 MB |

On BOYD1 the tape is written through `Vec::push` with doubling (about
19.5 MB of realloc copies plus 39 MB of fresh pages), each push does one
relaxed atomic increment on one of 17 `Arc`s, and `drop(result)` does
244k release-ordered decrements plus the enum drop loop: 1.7 ms measured,
outside the benchmark's timed region but inside a user's. The structures
report estimated tape write plus drop at ~4% of BOYD1; boxing was skipped
then because it "touches every tape match".

Change. Two independent steps:

* Replace `Arc<Equation>` by an index into a tape-owned
  `equations: Vec<Equation>` (`RecoveryTape { rules, equations }`). The
  model's cache becomes `equations: Vec<u32>` (NONE when stale) pointing
  into that table; `equation(row)` pushes a snapshot when the cache is
  stale and returns the index; `substitute_on_side` and the four
  `Arc::new(Equation {..})` sites in `tape.rs` tests push their own
  snapshot. `Recovery` code looks up `self.equations[idx]`. The proof
  closure in `rules/bounds.rs:185` returns the index. No atomics anywhere,
  `Postsolve` stays `Send + Sync` (an `Rc` would not).
* Shrink `Rule` to 32 B: `u32` for column/row indices (dimensions are
  already required to be below `u32::MAX`, `src/problem/mod.rs:38-39`),
  `Side` as a byte, and `Box` the payload of `Substituted`, `Fixed`,
  `Unlocked`, `RowCombination`, `Eliminated`, `DependentRow`, `ConeSlack`
  (each of those already allocates a `Vec`, so the extra `Box` is one more
  small allocation per record on rules that are 1-3 orders of magnitude
  rarer than `TightenedBound`). `TightenedBound {column: u32, equation:
  u32, side: u8, old: f64}` then fits in 24 B including the tag.

Expected effect: 2.5x less tape traffic and no atomics. BOYD1 -2% to -3%
of the timed region and -1.5 ms of drop; BOYD2, QAP15, CONT-200, PILOT87,
STOCFOR3 about -1%; corpus sum -0.5% to -1%. Replay loops
(`recover_with_slacks`, `reduce_point`) become denser too, though postsolve
is not timed by the benchmark.

Risk: none to reductions (the tape is write-only during presolve); the
tape tests in `tape.rs` cover replay. Effort: medium, half a day; the
`Rule` match arms change shape but not logic.

### 3. Build configuration (see the dedicated section below)

### 4. `from_columns` fills every `Node` before overwriting it

`src/matrix/linked.rs:175-185`:

```rust
out.nodes = vec![Node { value: 0.0, row: 0, col: 0, prev: [NONE; 2], next: [NONE; 2] }; nnz];
```

Pass 2 (lines 186-211) writes every slot exactly once with all five fields
set (`offsets[row]` advances through each row's contiguous range and the
ranges partition `0..nnz`), and the only other write, `nodes[col_list.tail].next[1] = id`,
targets a slot already written. So the initial contents are never read.
Because `NONE` is `u32::MAX`, `vec![...]` cannot use `alloc_zeroed` and
emits a 32 B-per-nonzero fill: 18 MB on BOYD1, 13.5 MB on BOYD2/CONT-300,
about 0.4-0.6 ms each (the first-touch page faults are unavoidable and
would simply move to pass 2). Corpus: -0.3% to -0.5%.

Change: obtain zeroed, lazily-mapped memory without a fill. Without
`unsafe` in this crate that means `bytemuck::zeroed_vec::<Node>(nnz)` with
`#[derive(Pod, Zeroable)]` on `Node` (all-`u32`/`f64` fields, no padding:
32 B, `repr(C)` needed). With `unsafe`, `Vec::with_capacity(nnz)` plus
`set_len` after pass 2 and a debug assertion that `offsets[row]` reached
each row's tail. The structures report rejected the unsafe route; the
`bytemuck` route adds a small, `no_std`, widely audited dependency. Same
trick applies to `SymmetricMatrix::from_upper_columns` only in principle:
`vec![(0, 0.0); nnz]` already goes through `alloc_zeroed` (std's `IsZero`
covers tuples), so nothing to gain there.

Also in `from_columns`: `offsets: Vec<usize>` (line 164) could be `Vec<u32>`
(halves 1.5 MB on BOYD2; negligible), and the two visitor passes are already
inlined (the `visit` generic is instantiated once per closure; confirmed by
the symbol list). Pass 1 must stay: CSC carries column counts, and the
row-major layout needs row counts.

Note on node order: nodes are inserted row-contiguous with column links
threaded through them (line 202-208), which is what makes the arena-order
sweep in finding 1 valid and keeps row walks sequential until rules start
recycling nodes from the free list. Nothing to change.

### 5. `size()` walks every row, twice per call

`src/presolve.rs:223-243`, called at line 77 (`before`) and 165 (`after`).

Each call reads `rows[i]` (24 B) and `lists[0][i]` (12 B) for all m rows:
0.11-0.18 ms per call on BOYD2, 0.06-0.09 ms on CONT-300, i.e. 0.3-0.65% of
those problems, 0.2% of the corpus. `alive.iter().filter().count()` is a
third O(n) pass.

Change:

* `before`: when `problem.cones.is_empty()`, `variables = n`,
  `linear_rows = m`, `a_nonzeros = model.a.nnz()` (explicit zeros are
  dropped by `from_columns`, so the arena count equals the walked sum),
  `p_nonzeros = p.nnz()`, `conic_rows = g_nonzeros = 0`. With cones, keep
  the walk (it is one pass; `input_rows` already tells how many rows are
  conic, so only `g_nonzeros`/`a_nonzeros` need the walk).
* `after` on the unchanged path equals `before` (the comment at line 164
  already says so): skip the call when `model.revision == 0`.
* `after` on the reduced path: `survivors()` (line 253) already visits
  every row and every `alive` flag; accumulate `a_nonzeros`/`g_nonzeros`
  and the counts there and return them alongside `Survivors`.

Risk: the reported `Size` must not change by a single count, since
`compare size` is the regression gate; the identities above hold as long
as the linked matrix never stores an explicit zero (`insert` is only
reached with non-zero values; worth a `debug_assert!`). Effort: low.

### 6. Exit path: seven separate linear passes and an identity row map

`src/presolve.rs:253-399` and `src/problem/mod.rs:51-61`.

Passes today, all O(n) or O(m): `survivors` columns filter, `survivors`
rows filter, `inverse` for columns (line 335), the `.eq(0..before.variables)`
identity check (338-341), the `c`/`bounds` gather (356-361), `packed_rows`
(362-367), `domains` (368-389), plus `row_indices` at entry (line 75:
two allocations, 1.5 MB on BOYD2, identity map when there are no cones).
Total `survivors + build_postsolve` is 0.4% of the corpus and the
`compact_problem` scalar loops another ~0.5% (pack is the rest). Each pass
alone is tiny; the point is that they are separate allocations and
separate sweeps over arrays that just left cache.

Change (all trivially order-preserving):

* Build `columns`, `stable_to_compact_columns` and the identity check in
  one loop over `alive`; build `linear_rows`, `packed_rows`, and `domains`
  for linear rows in one loop over `rows`.
* Keep `row_indices` lazy when `cones.is_empty()`: `OriginalMap` gets an
  `Identity(m)` variant (or `linear: Option<&[usize]>`), `gather`/`scatter`
  use `copy_from_slice`, and `Postsolve::original_linear_rows()`
  materialises on first request (`OnceLock<Vec<usize>>`). This also
  removes 8 B x m from every `Postsolve` on LP inputs.
* The `rows: Vec<RowDomain>` copy at `working_model` (lines 133-140) can
  reuse the caller's allocation because `Constraint` and `RowDomain` are
  both 24 B/8-aligned: `std::mem::take(&mut problem.rows).into_iter().map(..).collect()`
  is done in place by std. The unchanged path then has to convert back
  (one O(m) pass, same cost as the copy it replaces), so this is a win only
  on reduced outcomes (207 of 236). Marginal; list it for completeness.

Expected: -0.2% to -0.3% corpus, mostly on m > 20k. Effort: low-medium
(the `OriginalMap` change touches `postsolve/mod.rs` and the certificate
path at `presolve.rs:96-100`).

### 7. `LinkedMatrix::pack`: index widths

`src/matrix/linked.rs:234-270`. The pack is 2.3% of the corpus (2.4 ms on
BOYD1, 2.0 ms on BOYD2). It reads rows sequentially (good: the arena is
still mostly row-ordered after presolve) and scatters into `ri`, `values`,
and `next[j]`, plus a `stable_to_compact[j]` lookup: four scattered
touches per nonzero, about 4 ns each on BOYD1. `next` and
`stable_to_compact` are `usize` (0.75 MB each on BOYD); as `u32` they are
half that and stay in L2 more reliably. The compaction branch
(`filled != nnz`, lines 256-268) copies everything a second time when any
row outside `rows` still holds nodes; it would be worth counting how often
it triggers on the corpus (a `Deleted` row whose entries were not removed)
before investing further. The output format (separate `Vec<usize>` row
indices and `Vec<f64>` values) fixes the two write streams. Expected -0.1%
to -0.3%; effort low.

### 8. Hessian copy (`from_upper_columns`)

`src/matrix/sparse.rs:52-95`. Two passes over the input upper triangle, a
zeroed arena (`alloc_zeroed`, no fill cost), and a transpose scatter.
The arena is 16 B per stored entry with full symmetric storage, so EXDATA
(1.125M input entries) writes a 36 MB arena: 5.0 ms of its 7.2 ms, and
since EXDATA comes back unchanged the copy is pure overhead. The cost is
bandwidth and first-touch page faults, not instruction count; the scatter's
3000 concurrent write streams fit in L2. The only real lever is layout:
a `(u32, f64)` tuple is still 16 B because of alignment, so halving the
traffic requires SoA (`Vec<u32>` + `Vec<f64>`, 12 B per entry, 27 MB) and
that changes `SymmetricMatrix::row()`'s return type, which every rule
consumes. Out of scope here; note that `Slot` (`start, len, capacity` as
`usize`, 24 B) could be three `u32`s without touching the rules. Elsewhere
the Hessian copy is 1.8% of the corpus, diagonal P mostly (LISWET, AUG2D,
CONT: 0.05-0.4 ms), where it is dominated by the per-column closure and
the slot arithmetic; not worth separate work.

### 9. Queue flag arrays and the `stale()` fill

`src/model/queues.rs:147-159` allocates nine `queued: Vec<bool>` arrays
(four of size m, five of size n) plus `ActivityRows.flags`; `from_parts`
adds `alive`, `locks` (16 B x n), `activities` (32 B x m, non-zero
`stale()` so it is a real fill: 6 MB on BOYD2, 0.2 ms), `row_kinds`,
`equations`, `elimination_rejected`. Allocation itself is 0.2% of the
corpus; the first-touch cost lands in the loops that use them. Two row-side
and two column-side bit-mask arrays (`u8` per index) would replace the
nine `bool` arrays and put the flags that `row_changed` and
`column_changed` update together on one cache line, which is what the
`ActivityRows` merge did for two of them (that merge was one of the larger
wins in the structures report). A zero-representable stale marker (for
example storing `infinite + 1`) would make `activities` `alloc_zeroed`.
Both are small (-0.1% to -0.2%) and the flag merge belongs with the queue
owner; listed so it is not lost.

### 10. Unneeded `Arc` in `Coordinates`

`src/postsolve/mod.rs:15` and `src/presolve.rs:296`:
`compact_to_stable_columns: Arc<Vec<usize>>`. Nothing clones it (grep shows
no `Arc::clone` of it anywhere); a plain `Vec<usize>` removes one
allocation and one indirection per `Postsolve`. No measurable effect;
cleanup.

## Things checked and found fine

* `CscMatrix::from_parts` (`src/matrix/csc.rs:39-53`) adopts arrays with
  no scan, as documented; validation lives only in `CscMatrixRef::new`
  (lines 342-380), which the presolve entry never calls. `from_triplets`
  sorts and canonicalises, but it is a caller-side constructor. The
  `column(j)` iterator (line 336) slices twice per column; inlined, one
  bounds check per column, no per-nonzero cost beyond the zip.
* `Problem::take_objective` moves `c` and copies only the Hessian; the
  bounds vector is moved, not copied; `cones.clone()` is a few words.
* `configure` clones `Settings` (296 B, no heap); free.
* `Presolver::new(1)` builds no Rayon pool; `Executor::Serial` has no
  runtime cost. Rayon is linked into the binary regardless (symbols
  present in every build variant) but never called on the serial path.
* `Model` (1992 B) is moved by value into `finish`, `pack`,
  `compact_problem`; four 2 KB memcpys per call, irrelevant.
* Tape replay loops (`recover_with_slacks`, `reduce_point`) touch memory
  in tape order and per-record work is proportional to the record's own
  entries; nothing quadratic. Not in the timed region.
* `feasible_point` runs only on a dual-infeasibility certificate; rare.
* Per-call allocation count grows slowly with size (51 on a 2 x 1 problem,
  446 at n = 200, 1 867 at n = 2 000, 2 605 at n = 20 000) so no
  per-element allocation pattern remains on the boundary path.

## Compiler settings

`Cargo.toml` has no `[profile.*]` section. Release therefore uses
`opt-level = 3`, `codegen-units = 16`, `lto = false` (which in practice
means thin *local* LTO across the crate's own 16 codegen units),
`panic = "unwind"`, `debug-assertions = false`, no `target-cpu`. The
`benchmark` crate is a workspace member with no profile of its own, and
Cargo only honours the workspace root's profiles, so the benchmark binary
and the library are already built identically: **timings reflect what a
user with a default release profile gets.**

That last point is the important one for deciding what to change:
profile settings in a library's `Cargo.toml` are ignored by downstream
builds. Whatever is put under `[profile.release]` here changes only this
workspace's binaries. So:

* The fix for "small helpers not inlined across codegen units" (the +5%
  regression in `docs/benchmarks/refactor-20260918/REPORT.md`) is
  `#[inline]` in source, which reaches users' builds. It was applied to
  the six activity helpers; the symbol comparison below shows more
  candidates.
* `codegen-units = 1` and `lto` change the benchmark's numbers but not a
  user's, unless the README tells users to set them in their own profile.

### Measured build variants

Five builds of the `benchmark` binary from the scratch copy (`cargo build
--release -p benchmark --config ...`), Rust 1.98.1, arm64:

| Profile | build s | binary bytes | `__text` bytes | text symbols | distinct `presolve` text symbols |
| --- | ---: | ---: | ---: | ---: | ---: |
| default (16 CGU, thin-local LTO) | 16 | 2 706 464 | 1 494 556 | 5 281 | 456 |
| `codegen-units = 1` | 24 | 2 316 528 | 1 350 548 | 3 749 | 337 |
| `lto = "thin"` | 17 | 2 673 792 | 1 512 716 | 5 110 | 437 |
| `lto = "fat"`, `codegen-units = 1` | 48 | 2 085 392 | 1 316 276 | 2 886 | 283 |
| fat + cgu 1 + `panic = "abort"` | 33 | 1 761 360 | 1 200 028 | 2 165 | 237 |

Functions that exist as out-of-line symbols in the default build but are
fully inlined under `codegen-units = 1` (crate hashes stripped; excluding
Rayon/std plumbing):

* `matrix::linked::Iter<0>::next` (row iteration through the arena; the
  hottest loop primitive in the crate)
* `model::activity::Activity::replace_bound`
* `model::activity::Activity::compute::<View<0>>` and `::<View<1>>`
* `model::queues::Queues::row_changed`, `Queues::column_changed`,
  `Worklist::reset`
* `rules::sparsification::subtract::<Iter<0>>`
* `rules::parallel::proportional::<0>` and `::<1>`, `groups_from`,
  `candidate_groups`
* `matrix::sparse::SymmetricMatrix::set`, `add_variable`, `pack_upper`
* `matrix::linked::LinkedMatrix::pack`, `replace_rows`, `replace_row_into`,
  `add_column`, `sort_columns_by_length`
* `Model::from_parts`, `set_cones`, `singleton_range`, `equality_fingerprint`

A symbol being present does not prove every call site went through a
call, but for `Iter::next`, `replace_bound`, `row_changed` and
`column_changed` it means at least one hot caller in some codegen unit did
not inline it, and that which caller that is depends on the hash-based
CGU partition, which shifts whenever a function is added or renamed. That
is the mechanism behind the +/-2% build-to-build noise on sub-25 ms
problems noted in the structures report.

### Recommendations

1. **Keep the reference benchmark on the default profile.** It measures
   what users get.
2. **Add `#[inline]`** (not `always`) to: `Cursor::next`, `Iter::next`,
   `View::len/iter/cursor` (`src/matrix/linked.rs:56-105`),
   `Activity::replace_bound` and `compute` (`src/model/activity.rs`),
   `Queues::row_changed`, `column_changed`, `Worklist::push`, `pop`,
   `ActivityRows::push` (`src/model/queues.rs`), `Locks::add/remove/contribution`,
   `sparsification::subtract`, `parallel::proportional`. These are 5-20
   line leaf functions; the hint costs nothing when LLVM already inlined
   them and makes the decision independent of the CGU partition. Portable
   to users. Expected 0 to -2% depending on which sites are currently
   out-of-line; also reduces the layout noise.
3. **`codegen-units = 1` in `[profile.release]`: adopt for measurement
   stability, not as a speed claim.** It removes the partition noise and
   shrinks text by 10%, at +50% build time. If adopted, note in the README
   that the reported numbers use it and that users should set it in their
   own release profile. If not adopted, run the A/B once anyway to learn
   how far the default build is from the single-CGU ceiling; if the gap is
   more than 2%, more `#[inline]` hints are needed.
4. **`lto = "thin"`: no effect** (thin-local LTO is already on; the
   binary is within 1% of default in every metric) and it would only
   matter for cross-crate inlining into std/rayon, which the hot loops do
   not do. Skip.
5. **`lto = "fat"` + `codegen-units = 1`**: text -12% vs default. Fat LTO
   can help by inlining `alloc`/`dealloc` shims and by whole-program
   devirtualisation, but this crate has no trait objects on hot paths and
   its hot loops are already monomorphic. Expect 0 to -2% beyond
   `codegen-units = 1` alone; measure, do not assume. Build time 48 s vs
   16 s.
6. **`panic = "abort"`**: safe for this workspace. No `catch_unwind`,
   `resume_unwind` or `#[should_panic]` anywhere in `src/`, `tests/`,
   `benchmark/`; Cargo ignores `panic` for `cargo test` builds (verified:
   `cargo test --release --config 'profile.release.panic="abort"'` runs);
   the benchmark parent process checks `status.success()` of the worker
   (`benchmark/src/run.rs:272-279`), so an aborted worker (SIGABRT after
   the panic message) is still reported as a failure with its stderr.
   `Executor::Parallel` would abort the process instead of propagating a
   worker panic through `rayon::join`, which is acceptable for a
   benchmark binary. Text -9%, fewer landing pads, typically <=1% runtime.
   But it is a final-binary decision: a user's `panic = "unwind"` build
   compiles presolve with unwinding regardless. Only worth it if the
   owner accepts that the benchmark then measures a configuration users
   cannot get from the crate alone.
7. **`-C target-cpu=native`: do not use, and it would not help here.**
   The `aarch64-apple-darwin` target already sets `cpu = "apple-m1"` with
   `neon`, `lse`, `fp16`, `dotprod` etc. enabled (`rustc --print
   target-spec-json`), so `native` on an M-series machine changes at most
   the scheduling model. `wide` on aarch64 implements `f64x4` as two NEON
   `f64x2` halves and `u32x4` as one NEON register with or without
   `native` (`wide-1.7.0/src/f64x4_.rs:4-24`); the `avx` branch only
   exists for x86_64. On a generic x86_64 build (`x86-64` baseline, no
   AVX) `f64x4` becomes two SSE2 ops, which is still correct and the hash
   in `rules/parallel.rs:44-60` is not bandwidth-bound. For distribution,
   leave the target CPU to the user; `.cargo/config.toml` rustflags do not
   propagate through crates.io either.

### Timing experiment to run (when the machine is quiet)

Build three `benchmark` binaries from the same commit into separate target
directories:

```
cargo build --release -p benchmark --target-dir target-default
cargo build --release -p benchmark --target-dir target-cgu1 --config 'profile.release.codegen-units=1'
cargo build --release -p benchmark --target-dir target-fat  --config 'profile.release.codegen-units=1' --config 'profile.release.lto="fat"'
```

Then alternate them (`A B C A B C ...`) with `run time --trials 5 --rule all`
fresh-process per trial, three rounds each, and compare sums of medians
and geomeans exactly as in `docs/benchmarks/plain-csc-20260918/REPORT.md`.
Decision rule: if `cgu1` is within +/-1% of default, keep the default and
rely on `#[inline]` hints; if it is 2% or more faster, either adopt it and
document it, or treat the gap as a list of missing `#[inline]`s (compare
the two symbol lists with the `nm | sed` recipe above; the functions that
disappear under `cgu1` are the candidates). Add a fourth binary with
`panic = "abort"` only if the owner is open to the caveat in item 6.

## Scratch artifacts

Everything used for this review is under the session scratchpad:
`presolve/` (instrumented copy: `src/presolve.rs` and `src/model/mod.rs`
phase timers, `src/postsolve/mod.rs::tape_summary`, `tests/alloc_count.rs`,
`src/lib.rs::scratch_sizes`, `benchmark/src/bin/phases.rs`),
`profiles.sh` and `profiles.out` (the five builds), `target-*/` (their
outputs), `phases-corpus.txt` (per-problem phase breakdown for all 236
problems). The scratch build reads the corpus from the main tree's
`benchmark/data` and never writes to the main tree.
