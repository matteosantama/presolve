#!/usr/bin/env python3
"""Compare benchmark runs whose settings differ (the built-in compare refuses).

Usage:
  analyze.py size  <old> <new>            # reduced-size totals and per-problem changes
  analyze.py time  <old> <new> [<new2>..] # sum of medians and geometric mean of ratios
Runs are read from ./size or ./time next to this script.
"""
import json, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(kind, name):
    cases = {}
    with open(os.path.join(HERE, kind, f"{name}.jsonl")) as f:
        for line in f:
            rec = json.loads(line)
            if rec["type"] == "case":
                cases[rec["id"]] = rec["result"]
    return cases


def sizes(old, new):
    keys = ["variables", "linear_rows", "a_nonzeros", "p_nonzeros"]
    tot = {k: [0, 0] for k in keys + ["bound_sides"]}
    changed = []
    outcomes = {}
    for pid in sorted(set(old) | set(new)):
        a, b = old.get(pid), new.get(pid)
        if not a or not b or a["measurements"][0]["outcome"] != "reduced" or b["measurements"][0]["outcome"] != "reduced":
            outcomes[pid] = (a["measurements"][0]["outcome"] if a else None, b["measurements"][0]["outcome"] if b else None)
            continue
        ma, mb = a["measurements"][0], b["measurements"][0]
        for k in keys:
            tot[k][0] += ma["after"][k]
            tot[k][1] += mb["after"][k]
        tot["bound_sides"][0] += ma["after_bound_sides"]
        tot["bound_sides"][1] += mb["after_bound_sides"]
        d = tuple(mb["after"][k] - ma["after"][k] for k in keys)
        if any(d) or mb["after_bound_sides"] != ma["after_bound_sides"]:
            changed.append((pid, ma["after"]["variables"], mb["after"]["variables"], ma["after"]["linear_rows"], mb["after"]["linear_rows"], ma["after"]["a_nonzeros"], mb["after"]["a_nonzeros"], mb["after_bound_sides"] - ma["after_bound_sides"]))
    print("| Measure | old | new | change |")
    print("| --- | ---: | ---: | ---: |")
    for k, (x, y) in tot.items():
        print(f"| {k} | {x} | {y} | {y - x:+d} ({100 * (y - x) / max(x, 1):+.2f}%) |")
    diff = {p: o for p, o in outcomes.items() if o[0] != o[1]}
    print(f"\nproblems with a changed outcome: {diff}")
    print(f"\nproblems with changed sizes: {len(changed)}")
    print("| Problem | vars old | vars new | rows old | rows new | nnz old | nnz new | bound sides |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in sorted(changed, key=lambda r: r[1] - r[2], reverse=True):
        print("| " + " | ".join(str(x) for x in row[:7]) + f" | {row[7]:+d} |")


def medians(run):
    out = {}
    for pid, r in run.items():
        ts = [m["elapsed_ns"] for m in r["measurements"] if m.get("elapsed_ns") is not None]
        if ts and not r.get("errors"):
            out[pid] = statistics.median(ts) / 1e6
    return out


def times(names):
    runs = {n: medians(load("time", n)) for n in names}
    common = set.intersection(*(set(r) for r in runs.values()))
    base = runs[names[0]]
    print(f"{len(common)} common problems; base = {names[0]}")
    print("| Run | Sum of medians | vs base | Geomean ratio | Faster (>5%) | Slower (>5%) | Slowest problem change |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for n, r in runs.items():
        s = sum(r[p] for p in common)
        sb = sum(base[p] for p in common)
        g = math.exp(sum(math.log(r[p] / base[p]) for p in common) / len(common))
        faster = sum(1 for p in common if r[p] < 0.95 * base[p])
        slower = sum(1 for p in common if r[p] > 1.05 * base[p])
        worst = max(common, key=lambda p: r[p] / base[p])
        print(f"| {n} | {s:.1f} ms | {100 * (s - sb) / sb:+.2f}% | {g:.4f} | {faster} | {slower} | {worst} {100 * (r[worst] / base[worst] - 1):+.0f}% |")
    big = sorted(common, key=lambda p: -base[p])[:20]
    print("\nLargest problems (median ms):")
    print("| Problem | " + " | ".join(names) + " |")
    print("| --- |" + " ---: |" * len(names))
    for p in big:
        print(f"| {p} | " + " | ".join(f"{runs[n][p]:.2f}" for n in names) + " |")


if __name__ == "__main__":
    kind = sys.argv[1]
    if kind == "size":
        sizes(load("size", sys.argv[2]), load("size", sys.argv[3]))
    else:
        times(sys.argv[2:])
