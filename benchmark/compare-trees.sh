#!/usr/bin/env bash
# Snapshot the corpus with the benchmark binaries built from two source trees,
# then compare the snapshots with the head tree's binary.
#
# Usage: benchmark/compare-trees.sh BASE_TREE HEAD_TREE OUT_DIR [PROFILE...]
#
# Both binaries read HEAD_TREE/benchmark/data, so only code differs. Profiles
# default to "default aggressive". Each tree builds into its own target
# directory with its own Cargo.lock and the workspace release profile.
#
# Writes OUT_DIR/{base,head}-PROFILE.jsonl and OUT_DIR/summary.md (markdown).
# Exit status: 0 no changes, 1 statuses, outcomes, sizes, or reductions
# changed, 2 error, 3 BASE_TREE has no snapshot command to compare against.
set -euo pipefail

if [ $# -lt 3 ]; then
    sed -n '2,15p' "$0" >&2
    exit 2
fi
base=$(cd "$1" && pwd)
head=$(cd "$2" && pwd)
mkdir -p "$3"
out=$(cd "$3" && pwd)
shift 3
if [ $# -gt 0 ]; then profiles=("$@"); else profiles=(default aggressive); fi
data="$head/benchmark/data"
summary="$out/summary.md"
: > "$summary"

# Build a tree's benchmark binary and print its path.
build() {
    CARGO_TARGET_DIR="$1/target" cargo build --release --locked -p benchmark \
        --manifest-path "$1/Cargo.toml" >&2
    echo "$1/target/release/benchmark"
}

head_bin=$(build "$head") || exit 2
if ! base_bin=$(build "$base"); then
    echo "The base tree's benchmark crate failed to build." | tee -a "$summary" >&2
    exit 2
fi
if ! "$base_bin" snapshot --help 2>/dev/null | grep -q -- '--out'; then
    echo "The base tree has no snapshot command, so there is nothing to compare." \
        | tee -a "$summary" >&2
    exit 3
fi

status=0
for profile in "${profiles[@]}"; do
    for side in base head; do
        bin=${side}_bin
        echo "== $side, $profile profile" >&2
        "${!bin}" snapshot --profile "$profile" --data "$data" --slowest 5 \
            --out "$out/$side-$profile.jsonl" || exit 2
    done
    printf '### %s profile\n\n' "$profile" >> "$summary"
    set +e
    "$head_bin" compare "$out/base-$profile.jsonl" "$out/head-$profile.jsonl" \
        --format markdown >> "$summary"
    code=$?
    "$head_bin" compare "$out/base-$profile.jsonl" "$out/head-$profile.jsonl" --max-rows 20
    set -e
    echo >> "$summary"
    case $code in
        0) ;;
        1) [ "$status" -eq 0 ] && status=1 ;;
        *) exit 2 ;;
    esac
done
exit "$status"
