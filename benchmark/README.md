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

CI fails when reductions change. `benchmark/compare-instructions.sh` does the
same for the Callgrind instruction counts of `cargo bench -p benchmark --bench
presolve` (Linux with Valgrind); CI fails when one rises by more than 1%.

`verify --out FILE` takes the same selection flags as `snapshot`, checks each
recovered solution against Clarabel on the original, and exits 1 on a failure.

Time presolve only with the workspace release profile (one codegen unit, fat
LTO). To profile, build with the `profiling` profile, which adds debug info,
and run one instance repeatedly with the `profile` subcommand:

```sh
cargo build --profile profiling -p benchmark
```

On macOS, profilers read the debug info from the object files under
`target/profiling`, so leave that directory in place. Pick `--repeat` so the
run lasts a few seconds. With [samply](https://github.com/mstange/samply)
(`cargo install samply`), the profile opens in the Firefox Profiler:

```sh
samply record target/profiling/benchmark profile --instance miplib/rail507 --repeat 50
```

Instruments needs a full Xcode install, not only the Command Line Tools. Use
the Time Profiler template for where time goes, or CPU Counters for hardware
events such as cache misses, chosen in the Instruments app:

```sh
xcrun xctrace record --template 'Time Profiler' --output rail507.trace --launch -- target/profiling/benchmark profile --instance miplib/rail507 --repeat 50
open rail507.trace
```
