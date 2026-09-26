# benchmark

Regression tooling for the presolve library. It records what presolve does to
every instance in a corpus of MPS files and compares those records between two
versions of the code, so a pull request shows exactly which reductions it
changes.

## Corpus

Each subdirectory of `benchmark/data` is a family, and each `.mps` or `.mps.gz` file in
it is an instance, identified as `family/name`:

| Family | Instances | Contents |
| --- | --: | --- |
| `kennington` | 16 | Kennington LPs |
| `maros-meszaros` | 138 | Maros–Mészáros convex QPs |
| `miplib` | 236 | MIPLIB 2017 benchmark LP relaxations |
| `netlib` | 98 | Netlib LPs |
| `qplib` | 19 | QPLIB quadratic programs |

Every data file stays in the repository. A snapshot skips MIPLIB files above
4 MiB compressed by default to keep runs short, which leaves 471 of the 507
instances. `--max-size-mib FAMILY=MIB` replaces that default with per-family
limits, `--all-sizes` runs the whole corpus, and `--family` restricts a run to
some families. `--data DIR` points the tool at another directory laid out the
same way, such as a collection of your own models.

The MPS reader in [`src/mps.rs`](src/mps.rs) documents the format conventions
it follows. It reads free format and falls back to fixed format for files
whose names contain spaces.

## Snapshots

A snapshot presolves every selected instance and writes one JSON record per
line, sorted by instance. A record holds the status (presolved, load error, or
panic), the outcome, the sizes before and after, the reductions applied by rule
and kind, and work counters such as phase counts and equality pivot attempts.
It holds no timings.

Profiles run each presolve call on one thread without a time limit, so a
snapshot depends only on the code, the corpus, and the platform. Instances run
in parallel, one per core by default.

```sh
cargo run --release -p benchmark -- snapshot --profile default --out benchmark/results/base.jsonl
```

`--profile aggressive` uses `Settings::aggressive`. The command reports wall
time, the slowest presolve calls, and any instances that failed to load or
panicked. On 8 cores, the default selection takes about 8 s under the default
profile, and the whole corpus takes about 12 s under the default profile and
26 s under the aggressive one.

## Comparing snapshots

```sh
cargo run --release -p benchmark -- compare benchmark/results/base.jsonl benchmark/results/head.jsonl
```

The comparison lists outcome counts, sizes after presolve, reductions by rule
and kind summed over the corpus, and every changed field of every instance. It
exits with status 1 when statuses, outcomes, sizes, reductions, or the set of
instances differ, 0 when only work counters or nothing differ, and 2 on error.
Work counters are reported but never cause status 1. A field that only one
snapshot has, such as a counter added by the newer binary, is listed as not
compared. `--format markdown` writes the report for a CI job summary, and
`--max-rows N` limits the per-instance table.

Snapshots can differ between platforms, so compare only snapshots taken on the
same machine. They are build outputs; write them under `benchmark/results`,
which git ignores.

## Comparing two source trees

[`compare-trees.sh`](compare-trees.sh) builds the benchmark binary from two
source trees, snapshots both profiles with each binary against the head tree's
corpus, and compares them with the head tree's binary:

```sh
git worktree add ../presolve-main main
benchmark/compare-trees.sh ../presolve-main . benchmark/results/compare
```

It writes the four snapshots and a markdown summary to the output directory.
It exits with 0 for no changes, 1 for changed reductions, 2 on error, and 3
when the base tree has no snapshot command to compare against. Extra
arguments after the output directory select profiles.

## Continuous integration

The Reductions workflow runs `compare-trees.sh` on every pull request, on
Linux x86, Linux arm, and macOS arm runners. The base is the first parent of
the pull request's merge commit, and both sides build and run on the same
runner. The comparison goes to the job summary and the snapshots are uploaded
as artifacts. The job fails when reductions change, unless the pull request
carries the `reductions-change` label that accepts the change.

## Timing

Time presolve only with the workspace release profile, which builds one
codegen unit with full link-time optimization. With the default sixteen
units, code placement alone moves corpus timings by about 1% between builds,
which is as large as the effects worth measuring. A standalone timing crate
must set the same profile, and both sides of an A/B comparison must share one
`Cargo.lock`.
