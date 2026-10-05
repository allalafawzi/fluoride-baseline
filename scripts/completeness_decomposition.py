#!/usr/bin/env python3
"""
Decompose the 9.10-point gap the 95 % completeness filter opens in the Fluc
carrier proportion.

Applying a 95 % completeness filter after the fact to the panel already drawn
moves the carrier proportion from 13.03 % to 22.13 %, a gap of 9.10 points
(section 2.5). The raw gap mixes two effects the rest of the analysis keeps
apart: detection, because an incomplete genome has lost genes and the channel
can be absent from the assembly without being absent from the organism; and
composition, because the filter removes whole lineages, and not the lineages
where the channel is common.

The script separates them by holding phylum, then family, fixed:

  1. nothing held fixed  -> the raw gap.
  2. phylum held fixed   -> the whole panel is restricted to the phyla that
                            survive the filter, so what is left of the gap is
                            no longer attributable to losing phyla.
  3. family held fixed   -> only families holding both genomes below 95 % and
                            genomes above it are compared. This is the closest
                            comparison to a pure completeness effect.

Measured here, the four lost phyla account for 0.15 of the 9.10 points.

Genome size is not held fixed at any level. It is the strongest predictor of
the channel (section 2.6) and it is correlated with completeness, so the
residual within-family gap still contains that part.
"""
import argparse
import collections
import csv
import json
import sys


def load(plan_path, scan_path, prefix, marker):
    plan = list(csv.DictReader(open(plan_path), delimiter="\t"))[:prefix]
    pos = set()
    for r in csv.DictReader(open(scan_path), delimiter="\t"):
        if r["marker"] == marker and int(r["n_copies"] or 0) > 0:
            pos.add(r["genome_id"])
    for r in plan:
        r["carrier"] = r["accession"] in pos
        r["comp"] = float(r["completeness"])
    return plan


def prop(rows):
    return (sum(r["carrier"] for r in rows), len(rows),
            sum(r["carrier"] for r in rows) / len(rows) if rows else float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="data/gtdb/ar53_plan.tsv")
    ap.add_argument("--scan", default="results/gtdb/scan_archaea.tsv")
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--threshold", type=float, default=95.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    plan = load(a.plan, a.scan, a.prefix, a.marker)
    hi = [r for r in plan if r["comp"] >= a.threshold]
    if not hi:
        print("no genome above the threshold", file=sys.stderr)
        return 2

    k0, n0, p0 = prop(plan)
    k1, n1, p1 = prop(hi)

    # level 2: phylum held fixed
    survivors = {r["phylum"] for r in hi}
    sub = [r for r in plan if r["phylum"] in survivors]
    k2, n2, p2 = prop(sub)
    lost = sorted({r["phylum"] for r in plan} - survivors)

    # level 3: family held fixed
    byfam = collections.defaultdict(list)
    for r in plan:
        byfam[r["family"]].append(r)
    lo_rows, hi_rows, nfam = [], [], 0
    for rs in byfam.values():
        lo = [r for r in rs if r["comp"] < a.threshold]
        hh = [r for r in rs if r["comp"] >= a.threshold]
        if lo and hh:
            nfam += 1
            lo_rows += lo
            hi_rows += hh
    kl, nl, pl = prop(lo_rows)
    kh, nh, pheq = prop(hi_rows)

    comp_effect = 100 * (p2 - p0)
    det_effect = 100 * (p1 - p2)
    within = 100 * (pheq - pl)

    print(f"marker {a.marker}   completeness threshold {a.threshold} %\n")
    print(f"  whole panel                     {k0:>4}/{n0:<5} = {100*p0:6.2f} %")
    print(f"  completeness >= {a.threshold:g} %           {k1:>4}/{n1:<5} = {100*p1:6.2f} %")
    print(f"  raw gap                                          {100*(p1-p0):+6.2f} points\n")
    print(f"  --- phylum held fixed ({len(survivors)} surviving phyla) ---")
    print(f"  panel restricted to those phyla {k2:>4}/{n2:<5} = {100*p2:6.2f} %")
    print(f"  composition effect (phyla)                       {comp_effect:+6.2f} points")
    print(f"  residual at fixed phylum                         {det_effect:+6.2f} points")
    print(f"  lost phyla : {', '.join(lost) if lost else 'none'}\n")
    print(f"  --- family held fixed ({nfam} mixed families) ---")
    print(f"  genomes < {a.threshold:g} %                   {kl:>4}/{nl:<5} = {100*pl:6.2f} %")
    print(f"  genomes >= {a.threshold:g} %                  {kh:>4}/{nh:<5} = {100*pheq:6.2f} %")
    print(f"  within-family gap                                {within:+6.2f} points")
    print("\nGenome size is not held fixed at any level and is correlated with "
          "completeness, so the within-family gap still contains it.")

    rec = dict(marker=a.marker, threshold=a.threshold,
               panel=dict(k=k0, n=n0, prop=p0),
               filtered=dict(k=k1, n=n1, prop=p1),
               raw_gap_points=100 * (p1 - p0),
               phylum_fixed=dict(k=k2, n=n2, prop=p2,
                                 n_surviving_phyla=len(survivors),
                                 lost_phyla=lost),
               composition_effect_points=comp_effect,
               residual_at_fixed_phylum_points=det_effect,
               family_fixed=dict(n_mixed_families=nfam,
                                 below=dict(k=kl, n=nl, prop=pl),
                                 above=dict(k=kh, n=nh, prop=pheq),
                                 gap_points=within),
               caveat="Genome size is not held fixed at any level; it is the "
                      "strongest predictor of the channel and is correlated with "
                      "completeness, so the within-family gap still contains it.")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(rec, fh, indent=2)
        print(f"\n[ok] {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
