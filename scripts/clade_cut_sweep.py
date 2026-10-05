#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sweep the clade-inclusion threshold of the mean-by-phylum estimator and report
how the estimate moves with it.

The mean-by-phylum estimator averages the per-phylum proportions without
weighting, over the phyla holding at least min-n genomes. That min-n is a
choice: the published analysis uses 10, the estimate moves by several points
across the range swept here, and the published value sits high in that range.
The sweep is published so that two things can be checked. First, the raw
proportion, which is the headline figure, does not depend on min-n at all,
because it is not a per-clade estimator. Second, where the mean-by-phylum
estimator sits relative to the raw proportion: it is above it at every cut
except the smallest, min-n 1, at which every phylum enters however few genomes
it holds and the unweighted mean falls below the raw proportion.

Usage:
  python3 scripts/clade_cut_sweep.py --analysis results/gtdb/baseline_archaea.json
"""
import argparse, json, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--cuts", default="1,2,3,5,8,10,12,15,20,30")
    ap.add_argument("--published-cut", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    R = json.load(open(a.analysis))
    ph = R["by_phylum"]
    raw = R["raw_proportion"]["prop"]
    cuts = [int(c) for c in a.cuts.split(",")]

    print(f"raw proportion (independent of min-n): {100*raw:.2f} %\n")
    print(f"{'min-n':>6} {'phyla':>6} {'zero-positive':>14} {'mean-by-phylum':>16} {'> raw?':>7}")
    rows = []
    for c in cuts:
        sel = [p for p in ph if p["N"] >= c]
        if not sel:
            continue
        m = sum(p["prop"] for p in sel) / len(sel)
        z = sum(1 for p in sel if p["k"] == 0)
        rows.append(dict(min_n=c, n_phyla=len(sel), n_zero=z, mean_by_phylum=m))
        star = "*" if c == a.published_cut else " "
        print(f"{c:>5}{star} {len(sel):>6} {z:>14} {100*m:>15.2f}% {'yes' if m > raw else 'NO':>7}")

    vals = [r["mean_by_phylum"] for r in rows if r["min_n"] >= 3]
    pub = [r["mean_by_phylum"] for r in rows if r["min_n"] == a.published_cut]
    print()
    print(f"  range over min-n 3..{max(cuts)} : {100*min(vals):.1f} - {100*max(vals):.1f} %"
          f"  (amplitude {100*(max(vals)-min(vals)):.1f} points)")
    if pub:
        print(f"  published value (min-n {a.published_cut}) : {100*pub[0]:.1f} %")
    print(f"  estimator exceeds the raw proportion at every cut: "
          f"{'YES' if all(r['mean_by_phylum'] > raw for r in rows) else 'NO'}")
    print()
    print("  Note: phyla with zero positives are retained. Dropping them would")
    print("  bias the mean upward; they are real clades with a measured zero.")

    if a.out:
        json.dump(dict(raw_proportion=raw, published_cut=a.published_cut,
                       sweep=rows), open(a.out, "w"), indent=2)
        print(f"\n[ok] {a.out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
