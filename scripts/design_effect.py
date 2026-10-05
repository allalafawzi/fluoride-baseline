#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Measure the design effect of the balanced plan, and report cluster-corrected
confidence intervals.

Genomes drawn one per family are not independent observations: members of a
family share a recent ancestor and therefore most of their gene content. A
Wilson interval computed on N = 1259 pretends otherwise. The design effect
(DEFF) is how much information that costs:

    DEFF = 1 + (m_eff - 1) * ICC        N_eff = N / DEFF

where m_eff is Kish's effective cluster size, sum(n_i^2)/sum(n_i), and ICC is
the intraclass correlation of the trait at family rank.

m_eff is counted in the plan file, not assumed. The most even allocation of
1259 genomes over 591 families (514 families of 2 plus 77 of 3) gives
m_eff = 2.18, which is a lower bound on the real value; treating that bound as
a measurement overstates the information available. The plan file gives 293
families of 3, 82 of 2 and 216 of 1, so m_eff = 2.53 and the corrected
interval is wider.

Usage:
  python3 scripts/design_effect.py --plan data/gtdb/ar53_plan.tsv \\
      --analysis results/gtdb/baseline_archaea.json --prefix 1259
"""
import argparse, collections, csv, json, math, sys

import os as _os_sc
import sys as _sys_sc
_sys_sc.path.insert(0, _os_sc.path.dirname(_os_sc.path.abspath(__file__)))
from stats_common import wilson as _wilson, Z95

def wilson(k, n, z=Z95):
    """Delegate to scripts/stats_common.wilson: one definition, one z value."""
    return _wilson(k, n, z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", required=True)
    ap.add_argument("--analysis", required=True,
                    help="baseline_archaea*.json, for k, N and the measured ICC")
    ap.add_argument("--prefix", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    R = json.load(open(a.analysis))
    N = R["sample"]["N"]
    k = R["sample"]["k"]
    icc = R["measured_icc"]["family"]
    prefix = a.prefix or R.get("prefix") or N

    rows = list(csv.DictReader(open(a.plan), delimiter="\t"))
    fam = collections.Counter(r["family"] for r in rows[:prefix])

    n_tot = sum(fam.values())
    if n_tot != N:
        print(f"[warn] plan prefix has {n_tot} genomes, analysis says N={N}",
              file=sys.stderr)

    sizes = collections.Counter(fam.values())
    m_eff = sum(v * v for v in fam.values()) / n_tot
    deff = 1 + (m_eff - 1) * icc
    n_eff = n_tot / deff

    # the lower bound an even allocation would have given, for comparison
    F = len(fam)
    q, r = divmod(n_tot, F)
    even = [q + 1] * r + [q] * (F - r)
    m_even = sum(v * v for v in even) / n_tot
    d_even = 1 + (m_even - 1) * icc

    print("=== plan composition (counted, not assumed) ===")
    print(f"  genomes {n_tot}   families {F}")
    for s in sorted(sizes):
        print(f"    {sizes[s]:>4} families contribute {s} genome(s)")
    print()
    print(f"  m_eff (Kish)  = {m_eff:.4f}     [even allocation would give {m_even:.4f}]")
    print(f"  ICC (family)  = {icc:.4f}")
    print(f"  DEFF          = {deff:.4f}     [even allocation would give {d_even:.4f}]")
    print(f"  N_eff         = {n_eff:.1f}       [even allocation would give {n_tot/d_even:.1f}]")
    print()
    print("=== intervals ===")
    _, lo, hi = wilson(k, n_tot)
    _, clo, chi = wilson(k * n_eff / n_tot, n_eff)
    print(f"  marker positive {k}/{n_tot} = {100*k/n_tot:.2f} %")
    print(f"    naive    [{100*lo:.2f} ; {100*hi:.2f}]")
    print(f"    corrected[{100*clo:.2f} ; {100*chi:.2f}]")
    _, zlo, zhi = wilson(0, n_tot)
    _, czlo, czhi = wilson(0, n_eff)
    print(f"  zero count 0/{n_tot}")
    print(f"    naive    [{100*zlo:.2f} ; {100*zhi:.2f}]")
    print(f"    corrected[{100*czlo:.2f} ; {100*czhi:.2f}]")

    out = dict(n_genomes=n_tot, n_families=F,
               family_size_distribution={str(s): sizes[s] for s in sorted(sizes)},
               m_eff_kish=m_eff, m_eff_even_allocation_lower_bound=m_even,
               icc_family=icc, deff=deff, n_eff=n_eff,
               marker_k=k,
               interval_naive=[lo, hi], interval_cluster_corrected=[clo, chi],
               zero_interval_naive=[zlo, zhi],
               zero_interval_cluster_corrected=[czlo, czhi])
    if a.out:
        json.dump(out, open(a.out, "w"), indent=2)
        print(f"\n[ok] {a.out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
