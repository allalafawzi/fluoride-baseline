#!/usr/bin/env python3
"""
Reach of the Fluc/CrcB detector into the phyla absent from its seed.

The headline prevalence is measured with a model whose seed holds 25 bacterial
and 5 archaeal entries covering 4 of the panel's 21 phyla, and the phyla with
the lowest prevalence are the ones the seed does not cover. Absence of the
channel and blindness of the model are therefore not distinguishable from the
census table alone.

This script measures what the deposited files allow: for each phylum, the
distribution of the best unthresholded Fluc score among the undetected genomes.
A model that half-saw a divergent subfamily would leave a mass of intermediate
scores, above background and below the threshold.

A sufficiently divergent homolog can score at background level, so the absence
of intermediate scores rules out partial detection, not total blindness. The
decisive test -- removing the five archaeal seed entries, recalibrating and
rescanning -- needs the proteomes (about 700 GB), which are not deposited.
"""
import argparse
import csv
import json
import statistics
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="data/gtdb/ar53_plan.tsv")
    ap.add_argument("--scan", default="results/gtdb/scan_archaea_all_scores.tsv")
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--ga", type=float, default=40.95)
    ap.add_argument("--nc", type=float, default=13.8,
                    help="background threshold: score of the best confirmed negative")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    plan = list(csv.DictReader(open(a.plan), delimiter="\t"))[:a.prefix]
    phylum = {r["accession"]: r["phylum"] for r in plan}
    score = {}
    for r in csv.DictReader(open(a.scan), delimiter="\t"):
        if r["marker"] == a.marker and r["genome_id"] in phylum:
            try:
                score[r["genome_id"]] = float(r["best_score_unthresholded"] or 0)
            except ValueError:
                score[r["genome_id"]] = 0.0
    if not score:
        print("no score read: wrong marker or wrong table?", file=sys.stderr)
        return 2

    rows = []
    for p in sorted(set(phylum.values())):
        g = [x for x in phylum if phylum[x] == p and x in score]
        if not g:
            continue
        neg = [score[x] for x in g if score[x] < a.ga]
        k = sum(1 for x in g if score[x] >= a.ga)
        rows.append(dict(
            phylum=p, n=len(g), k=k, prevalence=k / len(g),
            undetected_median=statistics.median(neg) if neg else None,
            undetected_max=max(neg) if neg else None,
            undetected_above_nc=sum(1 for x in neg if x >= a.nc)))
    rows.sort(key=lambda r: (r["prevalence"], -r["n"]))

    allv = sorted(score.values())
    neg = [x for x in allv if x < a.ga]
    pos = [x for x in allv if x >= a.ga]

    print(f"marker {a.marker}   N = {len(allv)}   GA = {a.ga}   NC = {a.nc}\n")
    print(f"{'phylum':28} {'N':>4} {'k':>4} {'prev':>6} | "
          f"undetected :   {'med':>6} {'max':>7} {'>=NC':>5}")
    for r in rows:
        print(f"{r['phylum']:28} {r['n']:>4} {r['k']:>4} {100*r['prevalence']:>5.1f}% |"
              f" {r['undetected_median']:>6.1f} {r['undetected_max']:>7.1f}"
              f" {r['undetected_above_nc']:>5}")
    gap = min(pos) - max(neg) if pos and neg else None
    print(f"\nundetected : median {statistics.median(neg):.1f} bits, "
          f"max {max(neg):.1f} ; {sum(1 for x in neg if x >= a.nc)} of {len(neg)} "
          f"reach the background threshold")
    print(f"detected   : min {min(pos):.1f} bits, median {statistics.median(pos):.1f}")
    print(f"empty from {max(neg):.1f} to {min(pos):.1f} bits : {gap:.1f} bits, 0 genome")
    print("\nA sufficiently divergent homolog can score at background level, so "
          "\nthe absence of intermediate scores rules out partial detection, not "
          "total blindness.")

    rec = dict(marker=a.marker, n_genomes=len(allv), ga=a.ga, nc=a.nc,
               by_phylum=rows,
               undetected_median=statistics.median(neg), undetected_max=max(neg),
               undetected_above_nc=sum(1 for x in neg if x >= a.nc),
               detected_min=min(pos), detected_median=statistics.median(pos),
               empty_band_bits=gap,
               caveat="A divergent homolog can score at background level; the "
                      "absence of intermediate scores rules out partial detection, "
                      "not total blindness. The decisive test (leave-one-out on the "
                      "five archaeal seed entries + rescan) needs the proteomes, "
                      "which are not deposited.")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(rec, fh, indent=2)
        print(f"\n[ok] {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
