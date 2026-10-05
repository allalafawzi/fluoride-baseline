#!/usr/bin/env python3
"""
Comparison of the Fluc/CrcB channel between archaeal lifestyles: halophiles,
thermophiles(-acidophiles), and the rest.

GTDB carries no lifestyle annotation. Lifestyle is therefore assigned here, by
hand, from the taxonomic rank, using only groups whose physiology is established
in the literature and homogeneous at that rank. It is an assignment by taxonomic
proxy, not a measurement: it applies to a minority of the genomes, and any
genome whose group is not physiologically homogeneous is left unclassified
rather than guessed.

Three rules were fixed before looking at the results:
  1. A group is kept only if the literature gives it one lifestyle throughout
     (for example Halobacteria: obligate halophiles with no known exception).
  2. No genome belongs to two categories. The Sulfolobales and the
     Thermoplasmatales are both thermophilic and acidophilic, so they form a
     category of their own, 'thermoacidophile', instead of being counted twice.
  3. Unclassified groups are not redistributed: they stay unclassified and their
     share is reported, so that the real coverage is visible.

The analysis is confounded with clade by construction: the halophiles are the
Halobacteria, so comparing halophiles with thermophiles compares two clades. The
script therefore also reports genome size by group, which is the most obvious
competing explanation. No p-value is computed: at this degree of confounding a
test would give a number with no interpretation.

Output: results/gtdb/lifestyle_fluc.json
"""
import argparse, csv, json, math, os
from collections import defaultdict

import os as _os_sc
import sys as _sys_sc
_sys_sc.path.insert(0, _os_sc.path.dirname(_os_sc.path.abspath(__file__)))
from stats_common import wilson as _wilson, Z95

# --- assignment by taxonomic proxy, declared in full --------------------------
# GTDB class -> lifestyle ; otherwise GTDB order -> lifestyle.
CLASS_LIFESTYLE = {
    # Obligate halophiles: growth at molar NaCl.
    "Halobacteria":  "halophile",
    "Nanosalinia":   "halophile",   # nanohaloarchaea, symbionts of haloarchaea
    # Thermophiles that are not acidophiles
    "Thermococci":   "thermophile",
    "Methanobacteria": None,        # too heterogeneous: mesophiles and thermophiles
}
ORDER_LIFESTYLE = {
    "Halobacteriales":    "halophile",
    "Halorutilales":      "halophile",
    "Thermococcales":     "thermophile",
    "Archaeoglobales":    "thermophile",
    "Methanopyrales":     "thermophile",
    "Thermoproteales":    "thermophile",
    "Thermofilales":      "thermophile",
    # thermoacidophiles: hot and acid, a category of their own (rule 2)
    "Sulfolobales":       "thermoacidophile",
    "Thermoplasmatales":  "thermoacidophile",
    "Aciduliprofundales": "thermoacidophile",
}


def lifestyle(row):
    o = ORDER_LIFESTYLE.get(row["order"])
    if o:
        return o
    c = CLASS_LIFESTYLE.get(row["class"])
    return c if c else "unclassified"


def wilson(k, n, z=Z95):
    """Delegates to scripts/stats_common.wilson so the deposit holds a single
    definition of the interval, and no second copy can drift from it.
    """
    return _wilson(k, n, z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="data/gtdb/ar53_plan.tsv")
    ap.add_argument("--scan", default="results/gtdb/scan_archaea.tsv")
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--out", default="results/gtdb/lifestyle_fluc.json")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.plan), delimiter="\t"))[:a.prefix]
    pos, score = set(), {}
    for r in csv.DictReader(open(a.scan), delimiter="\t"):
        if r["marker"] == a.marker and int(r["n_copies"]) > 0:
            pos.add(r["genome_id"])
            try:
                score[r["genome_id"]] = float(r["best_score"])
            except ValueError:
                pass

    grp = defaultdict(lambda: {"N": 0, "k": 0, "sizes": [], "scores": [],
                               "phyla": set(), "orders": set()})
    for r in rows:
        g = grp[lifestyle(r)]
        g["N"] += 1
        g["sizes"].append(float(r["genome_size"]))
        g["phyla"].add(r["phylum"])
        g["orders"].add(r["order"])
        if r["accession"] in pos:
            g["k"] += 1
            if r["accession"] in score:
                g["scores"].append(score[r["accession"]])

    def med(v):
        v = sorted(v)
        return None if not v else (v[len(v) // 2] if len(v) % 2
                                   else (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2)

    out = {
        "marker": a.marker,
        "N_total": len(rows),
        "assignment": "taxonomic proxy, hand-curated; see module docstring",
        "caveat": "lifestyle is confounded with clade by construction; no p-value "
                  "is computed and none should be",
        "class_map": {k: v for k, v in CLASS_LIFESTYLE.items() if v},
        "order_map": ORDER_LIFESTYLE,
        "groups": [],
    }
    for name in ("halophile", "thermophile", "thermoacidophile", "unclassified"):
        g = grp.get(name)
        if not g:
            continue
        p, lo, hi = wilson(g["k"], g["N"])
        out["groups"].append({
            "lifestyle": name, "N": g["N"], "k": g["k"],
            "prop": p, "lo": lo, "hi": hi,
            "median_genome_size": med(g["sizes"]),
            "median_best_score": med(g["scores"]),
            "n_phyla": len(g["phyla"]), "n_orders": len(g["orders"]),
            "orders": sorted(g["orders"])[:12],
        })
    classified = sum(x["N"] for x in out["groups"] if x["lifestyle"] != "unclassified")
    out["coverage"] = {"classified": classified,
                       "fraction": classified / len(rows)}

    # --- matching on genome size ---------------------------------------------
    # Genome size explains much of Fluc status (a factor of 69 from the smallest
    # genome to the largest). Comparing lifestyles without matching on size
    # would measure size, so the comparison is stratified into size bands.
    # Within a band size is only approximately matched, not held fixed.
    def band(sz):
        mb = sz / 1e6
        for lo, hi in ((0, 1.0), (1.0, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 99)):
            if lo <= mb < hi:
                return f"{lo}-{hi}"
        return "3.5-99"

    strat = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in rows:
        cell = strat[band(float(r["genome_size"]))][lifestyle(r)]
        cell[0] += 1
        if r["accession"] in pos:
            cell[1] += 1
    out["size_matched"] = {
        b: {g: {"N": v[0], "k": v[1], "prop": (v[1] / v[0]) if v[0] else None}
            for g, v in d.items()}
        for b, d in strat.items()
    }

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write("\n")

    print(f"{'lifestyle':<18}{'N':>5}{'k':>5}{'%':>8}   CI95        "
          f"{'med. size':>12}{'med. score':>11}")
    for x in out["groups"]:
        print(f"{x['lifestyle']:<18}{x['N']:>5}{x['k']:>5}{100*x['prop']:>7.1f}%"
              f"  [{100*x['lo']:4.1f};{100*x['hi']:5.1f}]"
              f"{(x['median_genome_size'] or 0)/1e6:>11.2f} Mb"
              f"{x['median_best_score'] if x['median_best_score'] else 0:>10.1f}")
    print(f"\ncoverage : {classified}/{len(rows)} "
          f"({100*classified/len(rows):.1f} %) of the genomes are classified")

    print("\nMatched on genome size (size explains much of the status) :")
    cols = ("halophile", "thermophile", "thermoacidophile", "unclassified")
    print(f"{'band (Mb)':<12}" + "".join(f"{c:>20}" for c in cols))
    for b in ("0-1.0", "1.0-1.5", "1.5-2.5", "2.5-3.5", "3.5-99"):
        line = f"{b:<12}"
        for c in cols:
            v = strat[b].get(c)
            line += f"{(f'{v[1]}/{v[0]} = {100*v[1]/v[0]:.0f}%' if v and v[0] else '-'):>20}"
        print(line)


if __name__ == "__main__":
    main()
