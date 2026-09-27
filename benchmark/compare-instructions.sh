#!/usr/bin/env bash
# Run the Callgrind benches built from two source trees and compare the
# instruction counts of each benchmark.
#
# Usage: benchmark/compare-instructions.sh BASE_TREE HEAD_TREE OUT_DIR
#
# Needs Valgrind and jq. Each tree builds into its own target directory with
# its own Cargo.lock and runs with the gungraun-runner version that Cargo.lock
# names, installed under HEAD_TREE/target. A benchmark fails when its count
# rises by more than MAX_INCREASE_PCT percent (default 1).
#
# Writes OUT_DIR/{base,head}.jsonl (gungraun summaries), OUT_DIR/{base,head}.tsv
# (benchmark and instruction count), and OUT_DIR/summary.md (markdown).
# Exit status: 0 no regressions, 1 a count rose above the threshold, 2 error,
# 3 BASE_TREE has no benches to compare against.
set -euo pipefail

max_increase=${MAX_INCREASE_PCT:-1}
bench=presolve

if [ $# -ne 3 ]; then
    sed -n '2,16p' "$0" >&2
    exit 2
fi
base=$(cd "$1" && pwd)
head=$(cd "$2" && pwd)
mkdir -p "$3"
out=$(cd "$3" && pwd)
summary="$out/summary.md"
: > "$summary"
export LC_ALL=C

# Install the gungraun-runner that a tree's Cargo.lock asks for and print its
# directory.
runner() {
    local version root
    version=$(awk '$0 == "name = \"gungraun\"" { getline; gsub(/version = |"/, ""); print }' \
        "$1/Cargo.lock")
    [ -n "$version" ] || return 1
    root="$head/target/gungraun-runner-$version"
    [ -x "$root/bin/gungraun-runner" ] ||
        cargo install --locked --root "$root" gungraun-runner --version "$version" >&2 ||
        return 1
    echo "$root/bin"
}

# Run a tree's benches, writing their summaries to OUT_DIR/SIDE.jsonl and
# benchmark ids with instruction counts to OUT_DIR/SIDE.tsv.
run() {
    local tree=$1 side=$2 bin
    bin=$(runner "$tree") || return 1
    echo "== $side" >&2
    PATH="$bin:$PATH" CARGO_TARGET_DIR="$tree/target" cargo bench --locked -p benchmark \
        --bench "$bench" --manifest-path "$tree/Cargo.toml" \
        -- --output-format=json --parallel=auto > "$out/$side.jsonl" || return 1
    jq -r '[.module_path + (if .id then "." + .id else "" end),
            ([.profiles[] | select(.tool == "Callgrind") | .data.total.metrics.Ir.values.new][0]
             // 0)]
           | @tsv' "$out/$side.jsonl" | sort > "$out/$side.tsv"
}

if ! run "$head" head; then
    echo "The head tree's benches failed to build or run." | tee -a "$summary" >&2
    exit 2
fi
if [ ! -s "$out/head.tsv" ]; then
    echo "The head tree's benches reported no results." | tee -a "$summary" >&2
    exit 2
fi
# Without benches in the base, list the head's counts and exit 3.
missing=0
if [ ! -f "$base/benchmark/benches/$bench.rs" ]; then
    missing=1
    : > "$out/base.tsv"
elif ! run "$base" base; then
    echo "The base tree's benches failed to build or run." | tee -a "$summary" >&2
    exit 2
fi

# Rows of benchmark, base, head, and change, then a verdict. A benchmark on one
# side only is listed but never fails.
set +e
join -t $'\t' -a 1 -a 2 -e - -o 0,1.2,2.2 "$out/base.tsv" "$out/head.tsv" | awk -F '\t' \
    -v max="$max_increase" -v summary="$summary" '
    BEGIN { status = 0 }
    {
        if ($3 == "-") change = "removed"
        else if ($3 == 0 || $2 == 0) { change = "no count"; status = 2 }
        else if ($2 == "-") change = "new"
        else {
            pct = ($3 - $2) * 100 / $2
            change = sprintf("%+.3f%%", pct)
            if (pct > max) { change = change " regression"; if (status == 0) status = 1 }
        }
        rows = rows sprintf("| %s | %s | %s | %s |\n", $1, $2, $3, change)
        text = text sprintf("%-55s %14s %14s  %s\n", $1, $2, $3, change)
    }
    END {
        verdict = status == 1 ? "Instruction counts rose by more than " max "%." \
            : status == 2 ? "A benchmark reported no instruction count." \
            : "No instruction count rose by more than " max "%."
        printf "**%s**\n\n| benchmark | base | head | change |\n| :-- | --: | --: | --: |\n%s", \
            verdict, rows > summary
        printf "%s\n%s", verdict, text
        exit status
    }'
status=$?
set -e
if [ "$missing" -eq 1 ]; then
    printf '\nThe base tree has no instruction count benches, so there is nothing to compare.\n' \
        | tee -a "$summary" >&2
    [ "$status" -eq 2 ] && exit 2
    exit 3
fi
exit "$status"
