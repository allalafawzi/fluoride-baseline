#!/usr/bin/env python3
"""
Measure the detection floor of a model: how much recall it loses on homologs
more divergent than its seed.

A gathering threshold calibrated on characterised enzymes is calibrated on
sequences that are close to one another, so it under-detects divergent
environmental homologs by an amount that has to be measured rather than
assumed.

The model is re-thresholded at a series of margins of its lowest leave-one-out
score. At each margin the recall on a validation set of true homologs absent
from the seed is scored against the false positives among the negatives, giving
a trade-off table in <outdir>/<name>.sensitivity.tsv and the highest-recall
threshold that still admits no false positive.

With --confirmed-prefix, only negatives whose header begins with that prefix
(for example 'sp|') count as false positives; the others are reported
separately as unlabelled, since a sequence with no established function is
evidence neither way.

USAGE
  python3 scripts/sensitivity_audit.py --positives seeds/fluc_pos.faa \
      --negatives seeds/negatives.faa \
      --truth-positives seeds/Fluc_CrcB_truth.faa \
      --name Fluc_CrcB --outdir results/hmm --confirmed-prefix 'sp|'
"""
import argparse, subprocess, sys, json, tempfile, os, csv
from pathlib import Path

def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode: sys.exit(f"[FAIL] {c}\n{r.stderr[:1200]}")
    return r.stdout

def scored(hmm, fa, ga):
    with tempfile.NamedTemporaryFile(suffix=".tbl", delete=False) as t: tbl = t.name
    sh(f"hmmsearch --max -T {ga} --tblout {tbl} -o /dev/null {hmm} {fa}")
    s = set()
    for l in open(tbl):
        if not l.startswith("#"):
            f = l.split()
            if len(f) > 5: s.add(f[0])
    os.unlink(tbl); return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positives", required=True)
    ap.add_argument("--negatives", required=True)
    ap.add_argument("--truth-positives", required=True,
                    help="FASTA of true homologs not in the seed (validation)")
    ap.add_argument("--name", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--margins", default="1.00,0.90,0.80,0.70,0.60,0.50,0.40,0.30")
    ap.add_argument("--confirmed-prefix", default=None,
                    help="count as false positives only those negatives whose "
                         "header begins with this prefix (e.g. 'sp|'). The "
                         "others are reported separately, as unlabelled.")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    calib = json.loads((out / f"{a.name}.calib.json").read_text())
    loo_min = calib["loo_min"]
    hmm = out / f"{a.name}.hmm"
    n_truth = sum(1 for l in open(a.truth_positives) if l.startswith(">"))
    n_neg = sum(1 for l in open(a.negatives) if l.startswith(">"))

    rows = []
    for m in [float(x) for x in a.margins.split(",")]:
        ga = m * loo_min
        tp = len(scored(hmm, a.truth_positives, ga))
        all_fp = scored(hmm, a.negatives, ga)
        if a.confirmed_prefix:
            confirmed = [x for x in all_fp if str(x).startswith(a.confirmed_prefix)]
            fp, unlab = len(confirmed), len(all_fp) - len(confirmed)
        else:
            fp, unlab = len(all_fp), 0
        rows.append(dict(margin=m, GA=round(ga, 1),
                         recall=round(tp / n_truth, 4), TP=tp, N_truth=n_truth,
                         false_positives=fp, unlabelled_retained=unlab,
                         N_neg=n_neg, FP_rate=round(fp / n_neg, 4)))
    with open(out / f"{a.name}.sensitivity.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"\n{a.name} - recall / false positives trade-off (loo_min={loo_min:.1f})")
    print(f"{'margin':>6} {'GA':>8} {'recall':>8} {'TP':>6} {'FP':>5} "
          f"{'unlabelled':>10} {'FP rate':>9}")
    for r in rows:
        print(f"{r['margin']:>6.2f} {r['GA']:>8.1f} {100*r['recall']:>7.1f}% "
              f"{r['TP']:>6} {r['false_positives']:>5} "
              f"{r['unlabelled_retained']:>10} {100*r['FP_rate']:>8.2f}%")
    best = [r for r in rows if r["false_positives"] == 0]
    if best:
        b = max(best, key=lambda r: r["recall"])
        print(f"\n  -> optimal GA at 0 false positives: {b['GA']:.1f} (margin {b['margin']:.2f}), "
              f"recall {100*b['recall']:.1f} %")
        print(f"  -> detection floor at that threshold: {100*(1-b['recall']):.1f} % "
              "of the divergent homologs are missed.")

if __name__ == "__main__":
    main()
