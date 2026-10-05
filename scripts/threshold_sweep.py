#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sweep the gathering threshold and report how the prevalence moves.

The threshold rule is the methodological claim of this work, so the first
question a reader has is how much the published number depends on where the
threshold sits. The sweep answers it from the deposited sub-threshold scores.

Positivity follows the census rule by default: `hmmsearch --cut_ga`, which
produced the census, requires the sequence score and at least one domain score
to clear the threshold. Counting on the sequence score alone is a different
rule; the two agree at the published threshold (164 either way) and diverge at
the edges of the window, where at NC they give 209 against 195. The rule in
force is therefore written into the output, and check_deposit.py re-derives a
point of the sweep from the raw scores.

Usage:
  python3 scripts/threshold_sweep.py \
      --scan results/gtdb/scan_archaea_all_scores.tsv \
      --plan data/gtdb/ar53_plan.tsv --prefix 1259 --marker Fluc_CrcB \
      --out results/gtdb/threshold_sensitivity.json
"""
import argparse, csv, json, math, re, sys
from pathlib import Path

import os as _os_sc
import sys as _sys_sc
_sys_sc.path.insert(0, _os_sc.path.dirname(_os_sc.path.abspath(__file__)))
from stats_common import wilson as _wilson, Z95


def hmm_ga(path):
    """The GA line of a model header, or None."""
    for line in open(path, errors="ignore"):
        m = re.match(r"^GA\s+([-\d.]+)", line)
        if m:
            return float(m.group(1))
        if line.startswith("HMM "):
            break
    return None


def wilson(k, n, z=Z95):
    """Delegate to scripts/stats_common.wilson: one definition, one z value."""
    return _wilson(k, n, z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", required=True,
                    help="uncensored scan table (needs both score columns)")
    ap.add_argument("--plan", required=True)
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--hmm", default=None,
                    help="model, to read GA and the window; defaults to results/hmm/<marker>.hmm")
    ap.add_argument("--thresholds", default=None,
                    help="comma-separated; defaults to a grid spanning the window")
    ap.add_argument("--rule", choices=["cut_ga", "sequence_only"], default="cut_ga",
                    help="'cut_ga' (default, and what the census used): sequence and "
                         "domain must clear the threshold. 'sequence_only' is provided "
                         "for comparison; it is not the census rule.")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.scan), delimiter="\t"))
    if not rows:
        print("empty scan table", file=sys.stderr); return 1
    need = ("best_score_unthresholded", "best_domain_score_unthresholded")
    missing = [c for c in need if c not in rows[0]]
    if missing and a.rule == "cut_ga":
        print(f"ERROR: {a.scan} lacks {missing}.", file=sys.stderr)
        print("  The census rule needs both scores. Regenerate the table with", file=sys.stderr)
        print("  scripts/uncensor_scan.sh, or pass --rule sequence_only (not the", file=sys.stderr)
        print("  census rule, and it will be recorded as such).", file=sys.stderr)
        return 1

    plan = [r["accession"] for r in
            list(csv.DictReader(open(a.plan), delimiter="\t"))[:a.prefix]]
    keep = set(plan)

    seq, dom = {}, {}
    for r in rows:
        if r["marker"] != a.marker or r["genome_id"] not in keep:
            continue
        s = (r.get("best_score_unthresholded") or "").strip()
        d = (r.get("best_domain_score_unthresholded") or "").strip()
        if s:
            seq[r["genome_id"]] = float(s)
        if d:
            dom[r["genome_id"]] = float(d)
    if not seq:
        print(f"no scores for marker {a.marker}", file=sys.stderr); return 1

    hp = Path(a.hmm) if a.hmm else Path("results/hmm") / f"{a.marker}.hmm"
    ga = hmm_ga(hp) if hp.exists() else None
    nc = tc = None
    if hp.exists():
        for line in open(hp, errors="ignore"):
            m = re.match(r"^NC\s+([-\d.]+)", line)
            if m: nc = float(m.group(1))
            m = re.match(r"^TC\s+([-\d.]+)", line)
            if m: tc = float(m.group(1))
            if line.startswith("HMM "): break

    if a.thresholds:
        grid = [float(x) for x in a.thresholds.split(",")]
    else:
        grid = sorted({v for v in (nc, 20.0, 27.2, 35.0, ga, 50.0, 61.29, tc, 80.0,
                                   100.0, 120.0) if v is not None})

    def positive(g, t):
        if a.rule == "sequence_only":
            return seq.get(g, -1e9) >= t
        return seq.get(g, -1e9) >= t and dom.get(g, -1e9) >= t

    N = len(plan)
    out = []
    print(f"marker {a.marker}   N = {N}   rule = {a.rule}")
    if ga is not None:
        print(f"window [{nc} ; {tc}]   published GA {ga}")
    print(f"\n{'threshold':>10} {'positives':>10} {'prevalence':>11}  {'95% Wilson':>16}")
    for t in grid:
        k = sum(1 for g in plan if positive(g, t))
        _, lo, hi = wilson(k, N)
        tag = "  <- published GA" if ga is not None and abs(t - ga) < 1e-9 else ""
        print(f"{t:>10.2f} {k:>10} {100*k/N:>10.2f} %  [{100*lo:>5.2f} ; {100*hi:>5.2f}]{tag}")
        out.append(dict(threshold=t, positives=k, prevalence=k / N,
                        ci_lo=lo, ci_hi=hi))

    inside = [r for r in out if nc is not None and tc is not None
              and nc - 1e-9 <= r["threshold"] <= tc + 1e-9]
    amp = (max(r["prevalence"] for r in inside) -
           min(r["prevalence"] for r in inside)) if inside else None
    if amp is not None:
        print(f"\namplitude inside the window [{nc} ; {tc}] : "
              f"{100*amp:.2f} points")

    rec = dict(marker=a.marker, n_genomes=N, rule=a.rule,
               rule_description=("sequence score AND domain score must clear the "
                                 "threshold, as hmmsearch --cut_ga does"
                                 if a.rule == "cut_ga" else
                                 "sequence score only -- NOT the rule the census used"),
               window=[nc, tc], published_GA=ga,
               amplitude_in_window=amp, sweep=out,
               note=("Recomputable from the deposited uncensored scan table. "
                     "check_deposit.py re-derives a point of this sweep."))
    if a.out:
        Path(a.out).write_text(json.dumps(rec, indent=2))
        print(f"\n[ok] {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
