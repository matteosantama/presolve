# benchmark

Records presolve's reductions on the MPS corpus in `data`, one subdirectory
per family, and compares them between versions.

```sh
cargo run --release -p benchmark -- snapshot --profile default --out benchmark/results/base.jsonl
cargo run --release -p benchmark -- compare benchmark/results/base.jsonl benchmark/results/head.jsonl
```

`compare` exits 1 when outcomes, sizes, or reductions differ. By default MIPLIB
files above 4 MiB are skipped; `--all-sizes` runs everything.

To compare two source trees, as CI does on every pull request:

```sh
benchmark/compare-trees.sh ../presolve-main . benchmark/results/compare
```

CI fails on changed reductions unless the pull request has the
`reductions-change` label.

Time presolve only with the workspace release profile (one codegen unit, fat
LTO).
