#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
How many genomes have to be downloaded, at each completeness threshold.

The power calculation showed that the information comes from the number of
families sampled, not from the number of genomes: beyond one genome per family,
each additional genome is largely redundant.

This script crosses the two. The completeness threshold sets the number of
visible families, and a high threshold does not only remove genomes, it removes
whole families, so the information they carry cannot be recovered by
downloading more. The number of layers sets the number of genomes to download.
The script also reports how many are already on disk, hence the real remaining
volume.

No network. Reads the GTDB metadata and, if present,
data/gtdb/log/done_protein.txt.

Usage:
  python3 scripts/plan_sizing.py --domain Archaea
"""

import argparse, glob, gzip, json, os, sys
from collections import Counter, defaultdict

THRESHOLDS = [50, 70, 80, 85, 90, 95]
MB_PER_GENOME = 0.436          # measured: 658 MB for 1510 proteomes


def kish(counts):
    c = [x for x in counts if x > 0]
    return sum(x * x for x in c) / sum(c) if c else 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-dir", default="data/gtdb")
    ap.add_argument("--domain", default="Archaea", choices=["Archaea", "Bacteria"])
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--done", default="data/gtdb/log/done_protein.txt")
    ap.add_argument("--icc", type=float, default=0.15,
                    help="ICC assumed at planning time. It is not a "
                         "measurement: see --measured-icc-from.")
    ap.add_argument("--measured-icc-from",
                    default="results/gtdb/design_effect.json",
                    help="JSON carrying icc_family, the ICC measured on the "
                         "drawn panel. If it exists and differs from --icc, the "
                         "script says so and records both.")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    # The --icc default of 0.15 is a planning assumption, fixed before any
    # draw; it is not a measurement. The ICC measured on the drawn panel is
    # 0.730: with 1 259 genomes and m_eff 2.53 the assumption gives N_eff 1 024
    # and the measurement gives 595, so the planned information yield was 72 %
    # too high. The script therefore compares its assumption against the
    # measurement whenever the measurement is on disk, warns when the two
    # differ, and records both values in its output.
    icc_measured = None
    icc_cli = any(x.startswith("--icc") for x in sys.argv[1:])
    if a.measured_icc_from and os.path.exists(a.measured_icc_from):
        try:
            icc_measured = json.load(open(a.measured_icc_from))["icc_family"]
        except (KeyError, ValueError) as e:
            print(f"[warn] {a.measured_icc_from} unreadable ({e})", file=sys.stderr)
    if icc_measured is not None and abs(icc_measured - a.icc) > 0.05:
        print(f"[ATTENTION] assumed ICC {a.icc:.3f}, measured ICC {icc_measured:.3f} "
              f"({os.path.basename(a.measured_icc_from)}).\n"
              f"            The deff and n_eff columns below use the "
              f"assumption, not the measurement.\n"
              f"            The JSON carries both, and n_eff_at_measured_icc "
              f"gives the corrected value.", file=sys.stderr)

    pre = "ar53" if a.domain == "Archaea" else "bac120"
    c = sorted(glob.glob(os.path.join(a.meta_dir, f"{pre}_metadata*.tsv.gz")))
    if not c:
        sys.exit(f"[ERROR] no {pre}_metadata*.tsv.gz in {a.meta_dir}")
    meta = c[-1]

    already = set()
    if os.path.exists(a.done):
        already = {l.strip() for l in open(a.done) if l.strip()}
    print(f"[reading] {meta}", file=sys.stderr)
    print(f"[disk]    {len(already)} genomes already retrieved\n", file=sys.stderr)

    # accession -> (completeness, family)
    recs = []
    with gzip.open(meta, "rt", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(head)}
        comp_c = next((x for x in ("checkm2_completeness", "checkm_completeness")
                       if x in idx), None)
        cont_c = next((x for x in ("checkm2_contamination", "checkm_contamination")
                       if x in idx), None)
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
            fam = r.get("f__")
            if not fam:
                continue
            try:
                comp = float(p[idx[comp_c]]); cont = float(p[idx[cont_c]])
            except (ValueError, TypeError):
                continue
            if cont > a.max_contamination:
                continue
            acc = p[idx["accession"]]
            for pfx in ("RS_", "GB_"):
                if acc.startswith(pfx):
                    acc = acc[len(pfx):]
            recs.append((acc, comp, fam))

    rows = []
    print(f"{'thresh':>6} {'genomes':>8} {'families':>9} | "
          f"{'layers':>7} {'to draw':>8} {'have':>6} {'LEFT':>7} {'MB':>7} "
          f"{'DEFF':>6} {'compl.':>7} {'<80%':>7}", file=sys.stderr)
    print("  (compl. = mean completeness of the genomes actually chosen ;",
          file=sys.stderr)
    print("   <80% = how many of them are below 80 %, that is to say the "
          "families where nothing better exists)", file=sys.stderr)
    print("-" * 100, file=sys.stderr)
    for s in THRESHOLDS:
        sel = [(acc, comp, fam) for acc, comp, fam in recs if comp >= s]
        byfam = defaultdict(list)
        for acc, comp, fam in sel:
            byfam[fam].append((acc, comp))
        for v in byfam.values():
            v.sort(key=lambda t: (-t[1], t[0]))     # the most complete first
        n_fam = len(byfam)
        for k in (1, 2, 3):
            chosen = [(acc, comp) for v in byfam.values() for acc, comp in v[:k]]
            pick = [acc for acc, _ in chosen]
            comps = sorted(c for _, c in chosen)
            mean = sum(comps) / len(comps)
            med = comps[len(comps) // 2]
            below80 = sum(1 for c in comps if c < 80)
            below90 = sum(1 for c in comps if c < 90)
            counts = [min(len(v), k) for v in byfam.values()]
            m = kish(counts)
            d = 1 + (m - 1) * a.icc
            have = len(set(pick) & already)
            remaining = len(pick) - have
            rows.append(dict(threshold=s, visible_genomes=len(sel), families=n_fam,
                             layers=k, to_draw=len(pick), already=have,
                             remaining=remaining, mb=round(remaining * MB_PER_GENOME),
                             m_eff=round(m, 2), deff=round(d, 2),
                             n_eff=round(len(pick) / d),
                             n_eff_at_measured_icc=(
                                 round(len(pick) / (1 + (m - 1) * icc_measured))
                                 if icc_measured is not None else None),
                             mean_completeness=round(mean, 1),
                             median_completeness=round(med, 1),
                             n_below_80=below80, n_below_90=below90))
            print(f"{s:>5}% {len(sel):>8} {n_fam:>9} | {k:>7} {len(pick):>8} "
                  f"{have:>6} {remaining:>7} {remaining*MB_PER_GENOME:>7.0f} {d:>6.2f} "
                  f"{mean:>7.1f} {below80:>7}", file=sys.stderr)
        print("-" * 100, file=sys.stderr)

    # what a high threshold costs, in families: a lost family is information
    # that no further download makes up for
    fam_by_threshold = {}
    for s in THRESHOLDS:
        fam_by_threshold[s] = len({fam for _, comp, fam in recs if comp >= s})
    base = fam_by_threshold[50]
    print(f"\n[families visible by threshold]  "
          f"a family that is lost is information that no download "
          f"makes up for", file=sys.stderr)
    for s in THRESHOLDS:
        n = fam_by_threshold[s]
        print(f"  threshold {s:>3}% : {n:>4} families  ({100*n/base:>5.1f} % of what "
              f"exists at 50 %)", file=sys.stderr)

    out = dict(domain=a.domain, file=os.path.basename(meta),
               already_on_disk=len(already), icc=a.icc,
               icc_role="planning assumption, made before any draw",
               icc_source=("command line" if icc_cli
                           else "script default (0.15)"),
               icc_measured=icc_measured,
               icc_measured_from=(a.measured_icc_from
                                  if icc_measured is not None else None),
               mb_per_genome=MB_PER_GENOME,
               families_by_threshold=fam_by_threshold, sizing=rows)
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        json.dump(out, open(a.out, "w"), indent=2, ensure_ascii=False)
        print(f"\n[ok] written: {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
