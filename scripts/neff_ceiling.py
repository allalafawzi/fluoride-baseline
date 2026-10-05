#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ceiling on N_eff: how far a proportional draw keeps adding information, and on
which taxonomy the answer is computed.

For each draw size the script reports N_eff under a proportional draw, as the
median over R replicates with its range, and N_eff under the balanced draw,
until the capacity of the balanced plan is exhausted. The proportional column is
a random multinomial draw, so a single draw is not reproducible: the replicate
seeds are derived deterministically from the draw size and the replicate index,
and two runs give the same table.

Which taxonomy the numbers describe has to be read off the output. power.py
calls load_taxo(), which looks for bac120_metadata*.tsv.gz, the bacterial
metadata; without that file it falls back to modelled_taxo(), a synthetic
taxonomy of 4 000 families and 113 073 species with a power law on family
sizes. A ceiling computed that way describes the synthetic taxonomy and not the
archaeal domain, which holds 6 657 species representatives in 591 families. When
ar53_metadata*.tsv.gz is available, this script computes the archaeal domain as
well, so the two can be compared.
"""
import argparse
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import power as P                                  # noqa: E402

DRAWS = [1000, 4000, 20000, 40000]


def archaeal_taxo(meta_dir, min_comp=50.0, max_cont=5.0):
    """Same logic as power.load_taxo, but on ar53_metadata: the archaeal
    domain, for contrast with the synthetic taxonomy.

    The two quality filters are those of the plan, deliberately. power.load_taxo
    does not filter on quality; applied as is to the archaea it gives 603
    families and 6 968 representatives, whereas the rest of the analysis works
    on the 591 families and 6 657 representatives that pass completeness >= 50 %
    and contamination <= 5 % (sample_design.py, plan_sizing.py). The filters of
    the plan are therefore applied here, and the JSON records them."""
    import glob
    import gzip
    from collections import Counter
    cand = sorted(glob.glob(os.path.join(meta_dir, "ar53_metadata*.tsv.gz")))
    if not cand:
        return None
    fam, fo = Counter(), {}
    with gzip.open(cand[-1], "rt", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_tax, i_rep = head.index("gtdb_taxonomy"), head.index("gtdb_representative")
        i_comp = next((head.index(c) for c in ("checkm2_completeness",
                                               "checkm_completeness") if c in head), None)
        i_cont = next((head.index(c) for c in ("checkm2_contamination",
                                               "checkm_contamination") if c in head), None)
        if i_comp is None or i_cont is None:
            sys.exit("[ERROR] CheckM columns not found in " + cand[-1])
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(i_tax, i_rep, i_comp, i_cont):
                continue
            if p[i_rep].strip().lower() not in ("t", "true"):
                continue
            try:
                if float(p[i_comp]) < min_comp or float(p[i_cont]) > max_cont:
                    continue
            except ValueError:
                continue
            r = p[i_tax].split(";")
            f = next((x for x in r if x.startswith("f__")), None)
            o = next((x for x in r if x.startswith("o__")), None)
            if f and f != "f__" and o and o != "o__":
                fam[f] += 1
                fo[f] = o
    if not fam:
        return None
    keys = sorted(fam)
    ords_ = sorted({fo[k] for k in keys})
    oi = {o: i for i, o in enumerate(ords_)}
    return P.Taxo([fam[k] for k in keys], [oi[fo[k]] for k in keys],
                  "ARCHAEA %s (completeness >= %g %%, contamination <= %g %%)"
                  % (os.path.basename(cand[-1]), min_comp, max_cont))


def table(tx, icc, split, cap, reps, base_seed):
    i_f, i_o = icc * split, icc * (1 - split)
    rows = []
    for N in DRAWS:
        vals = []
        for s in range(reps):
            # seed derived from N and the replicate index: deterministic, not shared
            rng = np.random.default_rng(base_seed + 1_000_003 * s + N)
            D, _ = P.deff_at(tx, N, "proportional", None, i_f, i_o, rng)
            vals.append(N / D)
        prop = dict(median=round(statistics.median(vals)),
                    low=round(min(vals)), high=round(max(vals)))
        capacity = tx.capacity(cap)
        if N > capacity:
            bal = dict(verdict="capacity exceeded", capacity=capacity)
        else:
            rng = np.random.default_rng(base_seed + N)
            Db, cb = P.deff_at(tx, N, "balanced", cap, i_f, i_o, rng)
            bal = dict(n_eff=round(N / Db), deff=round(Db, 3),
                       m_fam=round(P.kish(cb), 2))
        rows.append(dict(n_drawn=N, proportional=prop, balanced=bal))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-dir", default="data/gtdb",
                    help="if ar53_metadata*.tsv.gz is there, the archaeal domain "
                         "is computed alongside the synthetic taxonomy")
    ap.add_argument("--icc", type=float, default=0.15,
                    help="assumed total ICC, split between family and order")
    ap.add_argument("--icc-split", type=float, default=0.5)
    ap.add_argument("--cap", type=int, default=3, help="layers of the balanced plan")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261001)
    ap.add_argument("--min-completeness", type=float, default=50.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing output even if it carries an "
                         "archaeal row this run cannot recompute")
    a = ap.parse_args()

    out = dict(icc_total=a.icc, icc_split=a.icc_split, cap=a.cap,
               replicates=a.reps, seed=a.seed, draws=DRAWS,
               archaea_filters=dict(min_completeness=a.min_completeness,
                                    max_contamination=a.max_contamination),
               taxonomies=[])

    todo = [("model", P.modelled_taxo())]
    arch = (archaeal_taxo(a.meta_dir, a.min_completeness, a.max_contamination)
            if os.path.isdir(a.meta_dir) else None)
    if arch is not None:
        todo.append(("archaea", arch))
    else:
        print("[note] no ar53_metadata*.tsv.gz in %s : the archaeal contrast "
              "is not computed (see data/gtdb/README.md)" % a.meta_dir,
              file=sys.stderr)

    for name, tx in todo:
        rows = table(tx, a.icc, a.icc_split, a.cap, a.reps, a.seed)
        rec = dict(key=name, source=tx.source, n_families=tx.n_fam,
                   n_orders=tx.n_ord, n_representatives=int(tx.fam_size.sum()),
                   capacity_at_cap=tx.capacity(a.cap), rows=rows)
        out["taxonomies"].append(rec)

        print(f"\n{tx.source}")
        print(f"  {tx.n_fam} families, {tx.n_ord} orders, "
              f"{tx.fam_size.sum()} representatives ; "
              f"capacity at {a.cap}/family = {tx.capacity(a.cap)}")
        print(f"  {'drawn':>7} | {'proportional : N_eff median (range)':<38} | "
              f"balanced {a.cap} layers")
        for r in rows:
            p = r["proportional"]
            b = r["balanced"]
            bal = (b["verdict"] if "verdict" in b
                   else f"N_eff {b['n_eff']}  (DEFF {b['deff']})")
            print(f"  {r['n_drawn']:>7} | {p['median']:>6} "
                  f"({p['low']} to {p['high']}){'':<20} | {bal}")

    print("\nUnder a proportional draw N_eff plateaus ; under a balanced draw it "
          "\nfollows N until the capacity of the plan is exhausted.")
    if len(todo) == 1:
        print("\nThis table is computed on the synthetic taxonomy, not on the "
              "\narchaeal domain.")

    if a.out:
        # The archaeal row needs the GTDB metadata, which the deposit does not
        # redistribute. Running this command without them would otherwise replace
        # a deposited output that carries the archaeal row with a model-only one,
        # with nothing in the output to say so. Such a write is refused, the file
        # is kept, and the refusal is reported.
        if (os.path.exists(a.out) and arch is None and not a.force):
            try:
                prev = json.load(open(a.out))
                had = any(t.get("key") == "archaea"
                          for t in prev.get("taxonomies", []))
            except (OSError, ValueError):
                had = False
            if had:
                print("\n[kept] %s is not rewritten : it carries an archaeal "
                      "row\n       this run cannot recompute (GTDB metadata "
                      "absent).\n       Add ar53_metadata*.tsv.gz to --meta-dir, "
                      "or pass --force\n       to write the model-only version."
                      % a.out, file=sys.stderr)
                return 0
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        json.dump(out, open(a.out, "w"), indent=2, ensure_ascii=False)
        print(f"\n[ok] {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
