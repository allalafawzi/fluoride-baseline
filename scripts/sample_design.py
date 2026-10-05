#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the genome sampling order from the GTDB metadata.

  python3 scripts/sample_design.py --meta data/gtdb/ar53_metadata_r226.tsv.gz \
      --domain Archaea --stratum family --out data/gtdb/ar53_plan.tsv

scripts/power.py shows that under proportional sampling (drawn at random among
the species representatives) the effective sample size saturates around
N_eff ~ 400-470 however many genomes are downloaded, because the largest
families absorb the extra draws. The same budget spent on one genome per
family gives N_eff ~ 2400 at N = 4000.

The output is therefore an ordered list in which every prefix is itself a
valid balanced sample:
  - layer 1: one genome per family, families in a reproducible pseudo-random
             order (FNV-1a hash of the family name, plus the seed)
  - layer 2: a second genome per family, same order
  - and so on.
In practice: download 500 genomes, estimate pi (the prevalence of FAcD), then
extend the same list up to N. The first 500 remain part of the final sample.

Completeness is a covariate, contamination a filter. An absent gene has to be
a biological absence rather than an assembly gap, and the obvious response is
to require completeness >= 95 %. In archaea that is the wrong move: many
lineages are known only from metagenome-assembled genomes (MAG) of middling
completeness, and a 95 % cut removes whole phyla, which are exactly the
uncultivated extremophiles of interest. scripts/completeness_audit.py
measures that loss.

So the two quality measures are treated differently:
  - contamination <= 5 % is a hard filter. Contamination produces false
    positives, a gene credited to a genome that does not carry it, and here a
    false positive costs more than a false negative.
  - completeness keeps a low floor (50 %, GTDB's own floor) and the value is
    written into the plan to be used as a covariate. The analysis stratifies
    by completeness and tests for a trend, and a sensitivity analysis
    restricted to genomes >= 90 % accompanies the main result. A gene present
    in a genome that is c % complete is seen with probability about c/100, so
    the raw prevalence is underestimated by about that factor; the correction
    belongs in the analysis, not in discarding genomes.

Within a family the most complete genome comes first, so layer 1 holds the
best representative of each family.
"""

import argparse, gzip, os, sys
from collections import defaultdict


def fnv1a(s, seed=0):
    """Deterministic hash: same seed -> same order, on any machine."""
    h = 2166136261 ^ (seed * 2654435761 & 0xFFFFFFFF)
    for ch in s.encode("utf-8"):
        h ^= ch
        h = (h * 16777619) & 0xFFFFFFFF
    return h


RANKS = {"d__": "domain", "p__": "phylum", "c__": "class",
         "o__": "order", "f__": "family", "g__": "genus", "s__": "species"}


def parse(meta, min_comp, max_cont, want_domain, log=sys.stderr):
    rows, seen_cols = [], {}
    kept = dropped_qual = 0
    with gzip.open(meta, "rt", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(head)}
        for need in ("accession", "gtdb_taxonomy", "gtdb_representative"):
            if need not in idx:
                print(f"[ERROR] missing column: {need}", file=log)
                return None, None
        # The CheckM column names change from one GTDB release to the next.
        size_c = next((c for c in ("genome_size", "total_gap_length")
                       if c in idx), None)
        prot_c = next((c for c in ("protein_count",) if c in idx), None)
        comp_c = next((c for c in ("checkm2_completeness", "checkm_completeness")
                       if c in idx), None)
        cont_c = next((c for c in ("checkm2_contamination", "checkm_contamination")
                       if c in idx), None)
        seen_cols = dict(completeness=comp_c, contamination=cont_c)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(head):
                continue
            if p[idx["gtdb_representative"]].strip().lower() not in ("t", "true"):
                continue
            tax = p[idx["gtdb_taxonomy"]]
            r = {}
            for tok in tax.split(";"):
                if tok[:3] in RANKS and len(tok) > 3:
                    r[RANKS[tok[:3]]] = tok[3:]
            if r.get("domain") != want_domain:
                continue
            if not r.get("family"):
                continue
            comp = cont = ""
            if comp_c and cont_c:
                try:
                    comp = float(p[idx[comp_c]]); cont = float(p[idx[cont_c]])
                    if comp < min_comp or cont > max_cont:
                        dropped_qual += 1
                        continue
                except ValueError:
                    comp = cont = ""
            r["_comp"], r["_cont"] = comp, cont
            # Genome size and protein count. DPANN genomes (Nanobdellota,
            # Micrarchaeota, Aenigmatarchaeota and others) are strongly
            # reduced, about 0.6-1.2 Mb. An absent gene there can mean
            # "symbiont with a reduced genome" rather than "no fluoride
            # pressure", and without these columns genome reduction would be
            # confounded with the signal being looked for.
            r["_size"] = p[idx[size_c]] if size_c else ""
            r["_nprot"] = p[idx[prot_c]] if prot_c else ""
            acc = p[idx["accession"]]
            for pre in ("RS_", "GB_"):
                if acc.startswith(pre):
                    acc = acc[len(pre):]
            rows.append((acc, r))
            kept += 1
    print(f"[quality] {kept} representatives kept, {dropped_qual} discarded "
          f"(completeness < {min_comp}% or contamination > {max_cont}%)", file=log)
    if not (seen_cols["completeness"] and seen_cols["contamination"]):
        print("[quality] WARNING: CheckM columns not found, filter not applied",
              file=log)
    return rows, seen_cols


def order(rows, seed=20260826, stratum="family"):
    """Order the genomes so that every prefix of the list is balanced."""
    by = defaultdict(list)
    for acc, r in rows:
        by[r[stratum]].append((acc, r))
    prio = {k: fnv1a(k, seed) for k in by}
    out = []
    for k, v in by.items():
        # Most complete first, so layer 1 holds the best representative of
        # the family; then RefSeq, for homogeneous annotation; then accession.
        v.sort(key=lambda t: (-(t[1].get("_comp") or 0),
                              0 if t[0].startswith("GCF_") else 1, t[0]))
        for j, (acc, r) in enumerate(v):
            out.append((j, prio[k], acc, r))
    out.sort(key=lambda t: (t[0], t[1], t[2]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True, help="*_metadata_r*.tsv.gz from GTDB")
    ap.add_argument("--domain", default="Bacteria", choices=["Bacteria", "Archaea"])
    ap.add_argument("--stratum", default="family",
                    choices=["family", "genus", "order"])
    ap.add_argument("--min-completeness", type=float, default=50.0,
                    help="low floor: completeness is a covariate, not a hard filter")
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows, _ = parse(a.meta, a.min_completeness, a.max_contamination, a.domain)
    if not rows:
        sys.exit(1)
    ordered = order(rows, a.seed, a.stratum)

    # Column names are fixed here and read by name everywhere else, so a
    # reordering or a rename would surface immediately rather than silently.
    with open(a.out, "w") as fh:
        fh.write("accession\tlayer\tphylum\tclass\torder\tfamily\tgenus"
                 "\tcompleteness\tcontamination\tgenome_size\tn_proteins\n")
        for j, _, acc, r in ordered:
            fh.write("\t".join([acc, str(j + 1), r.get("phylum", ""),
                                r.get("class", ""), r.get("order", ""),
                                r.get("family", ""), r.get("genus", ""),
                                str(r.get("_comp", "")), str(r.get("_cont", "")),
                                str(r.get("_size", "")), str(r.get("_nprot", ""))]) + "\n")

    n_strat = len({r[a.stratum] for _, r in rows})
    print(f"[plan] {len(ordered)} genomes ordered over {n_strat} {a.stratum} strata",
          file=sys.stderr)
    print(f"[plan] capacity at 1/{a.stratum} = {n_strat} genomes "
          f"(= sample size at design effect ~ 1)", file=sys.stderr)
    # Kish effective cluster size for a few prefix lengths.
    import math
    fam = [r[a.stratum] for _, _, _, r in ordered]
    for N in (250, 500, 1000, 2000, 3000, 4000, 6000, 8000):
        if N > len(fam):
            break
        c = defaultdict(int)
        for x in fam[:N]:
            c[x] += 1
        v = list(c.values())
        m = sum(x * x for x in v) / sum(v)
        print(f"[plan]   N={N:<6} {len(v):>5} {a.stratum} strata touched  "
              f"m_eff={m:5.2f}  DEFF(ICC=.15)={1+(m-1)*.15:5.2f}", file=sys.stderr)
    print(f"[plan] written: {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
