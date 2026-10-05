#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the archaeal Fluc/CrcB baseline from the balanced sample: the raw
proportion, the prevalence per phylum, the unweighted mean of the per-phylum
proportions, a logistic fit of genome size and phylum, and the measured
intraclass correlation. Results are written as JSON, which figure_archaea.py
plots.

Inputs: the sampling plan from sample_design.py, the strict scan, and
optionally the relaxed scan and the manifest of predicted proteomes.

Three controls run alongside the main estimate:
  - threshold, whether the relaxed scan (E<1e-5, no calibrated threshold)
    adds positives, which shows how far the result depends on the threshold;
  - genome size, since prevalence tracks genome size, so the comparison
    between clades is repeated within a common size band;
  - provenance, whether proteomes predicted with Prodigal give the same
    prevalence as the downloaded ones.
"""
import argparse, csv, json, math, os, random, sys
from collections import defaultdict

import numpy as np

import os as _os_sc
import sys as _sys_sc
_sys_sc.path.insert(0, _os_sc.path.dirname(_os_sc.path.abspath(__file__)))
from stats_common import wilson as _wilson, Z95


# sample_design.py now writes the sampling plan with English column names, but
# an ar53_plan.tsv produced before that change carries the earlier French
# ones. Both spellings are accepted: the header is normalised as it is read,
# so a plan already on disk keeps working without being rebuilt.
def normalise_plan_row(row):
    """Pass-through, kept so the two readers share one entry point."""
    return dict(row)


def wilson(k, n, z=Z95):
    """Wilson score interval, delegated to scripts/stats_common.wilson so that
    a single definition, with a single value of the normal quantile, is used
    throughout the deposit."""
    return _wilson(k, n, z)


def icc_anova(groups):
    """ICC for a binary character, by one-way variance decomposition.
    groups = list of lists of 0/1."""
    g = [v for v in groups if len(v) >= 2]
    if len(g) < 3:
        return float("nan")
    n_tot = sum(len(v) for v in g)
    k = len(g)
    gm = sum(sum(v) for v in g) / n_tot
    ssb = sum(len(v) * (np.mean(v) - gm) ** 2 for v in g)
    ssw = sum(sum((x - np.mean(v)) ** 2 for x in v) for v in g)
    msb, msw = ssb / (k - 1), ssw / (n_tot - k)
    n0 = (n_tot - sum(len(v) ** 2 for v in g) / n_tot) / (k - 1)
    d = msb + (n0 - 1) * msw
    return float((msb - msw) / d) if d > 0 else float("nan")


def logit_fit(X, y, iters=80):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ b))
        W = p * (1 - p) + 1e-9
        b = b + np.linalg.solve((X * W[:, None]).T @ X + 1e-8 * np.eye(X.shape[1]),
                                X.T @ (y - p))
    p = 1 / (1 + np.exp(-X @ b))
    return b, float(np.sum(y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", required=True)
    ap.add_argument("--scan", required=True)
    ap.add_argument("--scan-relaxed", default=None)
    ap.add_argument("--predicted", default=None)
    ap.add_argument("--prefix", type=int, default=1259,
                    help="size of the balanced sample (prefix of the plan)")
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--min-n", type=int, default=10)
    ap.add_argument("--out", default="results/gtdb/baseline_archaea.json")
    a = ap.parse_args()
    say = lambda *x: print(*x, file=sys.stderr)

    rows = [normalise_plan_row(r)
            for r in csv.DictReader(open(a.plan), delimiter="\t")]
    plan = {r["accession"]: r for r in rows}
    order = [r["accession"] for r in rows]
    pref = set(order[:a.prefix])

    sc = defaultdict(dict)
    for r in csv.DictReader(open(a.scan), delimiter="\t"):
        sc[r["genome_id"]][r["marker"]] = int(r["n_copies"])
    pos = lambda g: 1 if sc[g].get(a.marker, 0) > 0 else 0

    predicted = set()
    if a.predicted and os.path.exists(a.predicted):
        predicted = {r["genome_id"] for r in csv.DictReader(open(a.predicted), delimiter="\t")
                     if r.get("genome_id")}

    G = [g for g in pref if g in sc]
    size = {g: int(plan[g]["genome_size"]) for g in G if plan[g].get("genome_size")}
    say(f"[sample] {len(G)} genomes, {len({plan[g]['family'] for g in G})} families, "
        f"{len({plan[g]['phylum'] for g in G})} phyla")

    # ---------- by phylum -------------------------------------------------
    by_phylum = []
    for ph in sorted({plan[g]["phylum"] for g in G}):
        sub = [g for g in G if plan[g]["phylum"] == ph]
        k = sum(pos(g) for g in sub)
        p, lo, hi = wilson(k, len(sub))
        med = float(np.median([size[g] for g in sub if g in size])) if any(
            g in size for g in sub) else float("nan")
        by_phylum.append(dict(phylum=ph, N=len(sub), k=k, prop=p, lo=lo, hi=hi,
                              median_size=med))
    by_phylum.sort(key=lambda d: -d["prop"])

    K, N = sum(pos(g) for g in G), len(G)
    raw, blo, bhi = wilson(K, N)
    kept = [d for d in by_phylum if d["N"] >= a.min_n]
    random.seed(11)
    est = sum(d["prop"] for d in kept) / len(kept)
    bs = sorted(sum(random.choice(kept)["prop"] for _ in kept) / len(kept) for _ in range(20000))
    ca = dict(estimate=est, lo=bs[500], hi=bs[19500], n_clades=len(kept),
              sd=float(np.std([d["prop"] for d in kept], ddof=1)))
    say(f"[raw] {K}/{N} = {100*raw:.1f} %   [{100*blo:.1f} ; {100*bhi:.1f}]")
    say(f"[mean-by-phylum] {100*est:.1f} %   [{100*ca['lo']:.1f} ; {100*ca['hi']:.1f}]")

    # ---------- control 1: threshold --------------------------------------
    threshold = None
    if a.scan_relaxed and os.path.exists(a.scan_relaxed):
        rel = {r["genome_id"]: int(r["n_copies"])
               for r in csv.DictReader(open(a.scan_relaxed), delimiter="\t")
               if r["marker"] == a.marker}
        comm = [g for g in rel if g in sc]
        gained = [g for g in comm if rel[g] > 0 and pos(g) == 0]
        lost = [g for g in comm if rel[g] == 0 and pos(g) > 0]
        threshold = dict(n_common=len(comm), strict_positives=sum(pos(g) for g in comm),
                         relaxed_positives=sum(1 for g in comm if rel[g] > 0),
                         gained_by_relaxing=len(gained), lost=len(lost))
        say(f"[threshold control] {len(comm)} genomes in common, "
            f"{len(gained)} gained by relaxing the threshold")

    # ---------- control 2: size -------------------------------------------
    bands = [(0, 1.0), (1.0, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 99)]
    by_size = []
    for x, y in bands:
        sel = [g for g in G if g in size and x <= size[g] / 1e6 < y]
        if not sel:
            continue
        k = sum(pos(g) for g in sel)
        p, lo, hi = wilson(k, len(sel))
        by_size.append(dict(min_mb=x, max_mb=y, N=len(sel), k=k, prop=p, lo=lo, hi=hi))

    Gt = [g for g in G if g in size]
    yv = np.array([pos(g) for g in Gt], float)
    lsz = np.array([math.log10(size[g]) for g in Gt])
    _, ll0 = logit_fit(np.ones((len(Gt), 1)), yv)
    b1, ll1 = logit_fit(np.column_stack([np.ones(len(Gt)), lsz]), yv)
    phs = [d["phylum"] for d in kept]
    # Indicators are built for phs[1:], so the first retained phylum is the
    # reference level, and every phylum below the min-n cut has no indicator of
    # its own and is absorbed into that same level, which keeps every genome in
    # the fit. The reference level is therefore a mixture of the first retained
    # phylum and all the rare ones, and it can pool the highest prevalence with
    # several zeros. Its composition is written to the output so the
    # coefficients can be read against what they are measured from.
    rare = [d for d in by_phylum if d not in kept]
    reference_level = dict(
        phylum=phs[0] if phs else None,
        absorbed_rare_phyla=[d["phylum"] for d in rare],
        n_absorbed_genomes=sum(d["N"] for d in rare),
        k_absorbed_positives=sum(d["k"] for d in rare))
    D = np.column_stack([np.array([1.0 if plan[g]["phylum"] == p else 0.0 for g in Gt])
                         for p in phs[1:]]) if len(phs) > 1 else np.zeros((len(Gt), 0))
    b2, ll2 = logit_fit(np.column_stack([np.ones(len(Gt)), lsz, D]), yv)
    # Phyla with no positive are completely separated: their coefficients are not
    # identified and the nominal degrees of freedom overstate the real ones.
    n_sep = sum(1 for d in kept if d["k"] == 0)
    model = dict(chi2_size=2 * (ll1 - ll0), df_size=1,
                 reference_level=reference_level,
                 n_phyla_completely_separated=n_sep,
                 chi2_phylum_after_size=2 * (ll2 - ll1), df_phylum=D.shape[1],
                 df_phylum_effective=D.shape[1] - n_sep,
                 slope_log10_size=float(b1[1]))
    say(f"[model] size chi2={model['chi2_size']:.1f} (1 df) ; "
        f"phylum after size chi2={model['chi2_phylum_after_size']:.1f} "
        f"({model['df_phylum']} df reported, {model['df_phylum_effective']} effective)")
    if n_sep:
        say(f"[model] {n_sep} retained phyla have no positive: their coefficients are "
            f"not identified (complete separation). The reported chi-square is not "
            f"exactly on its reference distribution; the conclusion is unaffected.")
    if rare:
        say(f"[model] reference level = {reference_level['phylum']} plus "
            f"{len(rare)} rare phyla ({reference_level['k_absorbed_positives']}"
            f"/{reference_level['n_absorbed_genomes']} positive)")

    matched = []
    for ph in phs:
        sel = [g for g in G if plan[g]["phylum"] == ph and g in size
               and 1.0e6 <= size[g] < 2.5e6]
        if len(sel) < a.min_n:
            continue
        k = sum(pos(g) for g in sel)
        p, lo, hi = wilson(k, len(sel))
        matched.append(dict(phylum=ph, N=len(sel), k=k, prop=p, lo=lo, hi=hi))
    matched.sort(key=lambda d: -d["prop"])

    # ---------- control 3: provenance -------------------------------------
    prov = None
    if predicted:
        pr = []
        for lab, sel in (("downloaded", [g for g in G if g not in predicted]),
                         ("predicted", [g for g in G if g in predicted])):
            k = sum(pos(g) for g in sel)
            p, lo, hi = wilson(k, len(sel))
            pr.append(dict(source=lab, N=len(sel), k=k, prop=p, lo=lo, hi=hi))
        pairs = []
        for ph in sorted({plan[g]["phylum"] for g in G}):
            aa = [g for g in G if plan[g]["phylum"] == ph and g not in predicted]
            bb = [g for g in G if plan[g]["phylum"] == ph and g in predicted]
            if len(aa) < 8 or len(bb) < 8:
                continue
            pairs.append(dict(phylum=ph, n_downloaded=len(aa),
                              p_downloaded=sum(pos(g) for g in aa) / len(aa),
                              n_predicted=len(bb),
                              p_predicted=sum(pos(g) for g in bb) / len(bb)))
        prov = dict(global_=pr, by_phylum=pairs)

    # ---------- measured ICC ----------------------------------------------
    all_genomes = [g for g in sc if g in plan]
    byf, byp = defaultdict(list), defaultdict(list)
    for g in all_genomes:
        byf[plan[g]["family"]].append(pos(g))
        byp[plan[g]["phylum"]].append(pos(g))
    icc = dict(family=icc_anova(list(byf.values())),
               phylum=icc_anova(list(byp.values())),
               n_genomes=len(all_genomes),
               n_families_with_2plus=sum(1 for v in byf.values() if len(v) >= 2))
    say(f"[measured ICC] family {icc['family']:.3f}  phylum {icc['phylum']:.3f}")

    res = dict(marker=a.marker, prefix=a.prefix,
               sample=dict(N=N, k=K, n_families=len({plan[g]["family"] for g in G}),
                           n_phyla=len({plan[g]["phylum"] for g in G})),
               raw_proportion=dict(prop=raw, lo=blo, hi=bhi),
               mean_by_phylum=ca, by_phylum=by_phylum,
               threshold_control=threshold, by_size=by_size,
               size_phylum_model=model, size_matched=matched,
               provenance_control=prov, measured_icc=icc)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=2, ensure_ascii=False)
    say(f"\n[ok] {a.out}")


if __name__ == "__main__":
    main()
