#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Audit of a completeness filter: how much it would remove, and from where.

  python3 scripts/completeness_audit.py \
      --meta data/gtdb/ar53_metadata_r226.tsv.gz --domain Archaea \
      --out results/gtdb/completeness_archaea.json

sample_design.py can drop genomes whose CheckM completeness falls below a
threshold, on the grounds that an absent gene should be a biological absence
rather than an assembly gap. Any such threshold is nevertheless a choice, and
in archaea an expensive one: many uncultivated lineages are known only from
metagenome-assembled genomes (MAG) of middling completeness, so a strict cut
removes the uncultivated extremophiles the work is about.

This script measures the two quantities that decide the threshold:
  1. how much is lost at each candidate threshold;
  2. where the loss falls. A loss concentrated on a few phyla is not noise
     being cleaned up but whole lineages being removed.

It downloads nothing and writes nothing but its own output: it reads the GTDB
metadata file already on disk.
"""

import argparse, glob, gzip, json, os, sys
from collections import Counter, defaultdict

THRESHOLDS = [0, 50, 70, 80, 85, 90, 95, 99]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", help="*_metadata_r*.tsv.gz ; otherwise auto-detected")
    ap.add_argument("--meta-dir", default="data/gtdb")
    ap.add_argument("--domain", default="Archaea", choices=["Archaea", "Bacteria"])
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    meta = a.meta
    if not meta:
        pre = "ar53" if a.domain == "Archaea" else "bac120"
        c = sorted(glob.glob(os.path.join(a.meta_dir, f"{pre}_metadata*.tsv.gz")))
        if not c:
            sys.exit(f"[ERROR] no {pre}_metadata*.tsv.gz file in {a.meta_dir}")
        meta = c[-1]
    print(f"[reading] {meta}", file=sys.stderr)

    reps = []                       # (completeness, contamination, phylum, class)
    with gzip.open(meta, "rt", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(head)}
        comp_c = next((c for c in ("checkm2_completeness", "checkm_completeness")
                       if c in idx), None)
        cont_c = next((c for c in ("checkm2_contamination", "checkm_contamination")
                       if c in idx), None)
        if not comp_c or not cont_c:
            sys.exit("[ERROR] CheckM columns not found: " + ", ".join(head[:40]))
        print(f"[columns] completeness={comp_c}  contamination={cont_c}", file=sys.stderr)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(head):
                continue
            if p[idx["gtdb_representative"]].strip().lower() not in ("t", "true"):
                continue
            tax = p[idx["gtdb_taxonomy"]]
            if not tax.startswith(f"d__{a.domain}"):
                continue
            r = {t[:3]: t[3:] for t in tax.split(";") if len(t) > 3}
            try:
                comp = float(p[idx[comp_c]]); cont = float(p[idx[cont_c]])
            except ValueError:
                continue
            reps.append((comp, cont, r.get("p__", "?"), r.get("c__", "?")))

    n = len(reps)
    if n == 0:
        sys.exit("[ERROR] 0 representative read")
    print(f"[total] {n} {a.domain.lower()} species representatives\n", file=sys.stderr)

    # 1. How much is lost, threshold by threshold.
    glob_rows = []
    print(f"{'thresh':>7} {'retained':>9} {'%':>7}   (contamination <= "
          f"{a.max_contamination} % in every case)", file=sys.stderr)
    for s in THRESHOLDS:
        k = sum(1 for c, ct, _, _ in reps if c >= s and ct <= a.max_contamination)
        glob_rows.append(dict(threshold=s, retained=k, share=round(100 * k / n, 1)))
        print(f"{s:>6}% {k:>9} {100*k/n:>6.1f}%", file=sys.stderr)

    # 2. Where the loss falls: retention by phylum.
    tot_p = Counter(p for _, _, p, _ in reps)
    by_phylum = {}
    for s in (80, 90, 95):
        keep = Counter(p for c, ct, p, _ in reps
                       if c >= s and ct <= a.max_contamination)
        by_phylum[s] = {p: dict(total=tot_p[p], retained=keep.get(p, 0),
                                share=round(100 * keep.get(p, 0) / tot_p[p], 1))
                        for p in tot_p}

    s0 = 95
    lost = sorted(((p, d) for p, d in by_phylum[s0].items() if d["total"] >= 5),
                  key=lambda t: t[1]["share"])
    print(f"\n[phyla most amputated at the {s0} % threshold] "
          f"(>= 5 representatives only)", file=sys.stderr)
    print(f"  {'phylum':<34} {'total':>6} {'retained':>8} {'%':>7}", file=sys.stderr)
    for p, d in lost[:18]:
        print(f"  {p[:34]:<34} {d['total']:>6} {d['retained']:>8} {d['share']:>6.1f}%",
              file=sys.stderr)

    removed = [p for p, d in by_phylum[s0].items() if d["retained"] == 0]
    halved = [p for p, d in by_phylum[s0].items()
              if d["total"] >= 5 and 0 < d["share"] < 50]
    print(f"\n[verdict] at the {s0} % threshold: {len(removed)} phyla entirely "
          f"removed, {len(halved)} amputated by more than half", file=sys.stderr)

    # 3. Detection bias: is a gene that is present actually observed?
    # In a genome that is c % complete, a gene that is really present is seen
    # with probability about c/100, so the observed prevalence is
    # underestimated by roughly the mean completeness.
    det = {}
    for s in (0, 50, 80, 90, 95):
        sel = [c for c, ct, _, _ in reps if c >= s and ct <= a.max_contamination]
        if sel:
            m = sum(sel) / len(sel)
            det[s] = dict(n=len(sel), mean_completeness=round(m, 2),
                          relative_underestimation_pct=round(100 * (1 - m / 100), 2))
    print(f"\n[detection bias] mean completeness and expected "
          f"underestimation of the prevalence", file=sys.stderr)
    for s, d in det.items():
        print(f"  threshold {s:>3}% : n={d['n']:>5}  mean completeness "
              f"{d['mean_completeness']:>6.2f}%  -> prevalence underestimated "
              f"by about {d['relative_underestimation_pct']:.2f}%", file=sys.stderr)

    out = dict(domain=a.domain, file=os.path.basename(meta),
               total_representatives=n, max_contamination=a.max_contamination,
               by_threshold=glob_rows, by_phylum=by_phylum,
               phyla_entirely_removed_at_95=sorted(removed),
               phyla_more_than_halved_at_95=sorted(halved),
               detection_bias=det)
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        json.dump(out, open(a.out, "w"), indent=2, ensure_ascii=False)
        print(f"\n[ok] written: {a.out}", file=sys.stderr)
        print("     this is the measurement the threshold choice rests on.", file=sys.stderr)


if __name__ == "__main__":
    main()
