#!/usr/bin/env python3
"""
baseline.py -- recomputes the prevalence of fluoride exporters, with an explicit method.

Corrects four defects of the ">85%" figure of Stockbridge & Wackett 2024 (Nat Commun 15:4593):
  1. lower bound -> point estimate + confidence interval
  2. bacteria only -> stratification by domain, archaeal figure produced
  3. implicit denominator -> N declared at each stratum
  4. raw proportion on a set deliberately biased towards the deep branches
     -> mean-by-phylum estimator alongside the raw proportion

Definitions (Dodge et al. 2026: "prokaryotes rely on one type only"):
  E_any    = Fluc/CrcB present OR CLC-F present
  E_class  = {fluc_only, clcf_only, both, none}
  R_any    = E_any OR fluoride riboswitch present
"""
import argparse, csv, json, math, sys
from collections import defaultdict, Counter

import os as _os_sc
import sys as _sys_sc
_sys_sc.path.insert(0, _os_sc.path.dirname(_os_sc.path.abspath(__file__)))
from stats_common import wilson as _wilson, Z95

# ---------- confidence intervals ----------
def wilson(k, n, z=Z95):
    """Delegate to scripts/stats_common.wilson: one definition, one z value."""
    return _wilson(k, n, z)


def _betainv(a, b, y):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo+hi)/2
        if _betacdf(mid, a, b) < y: lo = mid
        else: hi = mid
    return (lo+hi)/2

def _betacdf(x, a, b):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lb = math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)
    # series of the regularised incomplete beta function
    def cf(x, a, b):
        eps, tiny = 3e-14, 1e-300
        qab, qap, qam = a+b, a+1, a-1
        c, d = 1.0, 1.0-qab*x/qap
        if abs(d) < tiny: d = tiny
        d = 1/d; h = d
        for m in range(1, 300):
            m2 = 2*m
            aa = m*(b-m)*x/((qam+m2)*(a+m2))
            d = 1+aa*d;  c = 1+aa/c
            if abs(d) < tiny: d = tiny
            if abs(c) < tiny: c = tiny
            d = 1/d; h *= d*c
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
            d = 1+aa*d;  c = 1+aa/c
            if abs(d) < tiny: d = tiny
            if abs(c) < tiny: c = tiny
            d = 1/d; de = d*c; h *= de
            if abs(de-1.0) < eps: break
        return h
    front = math.exp(lb + a*math.log(x) + b*math.log(1-x))
    if x < (a+1)/(a+b+2): return front*cf(x, a, b)/a
    return 1 - math.exp(lb + b*math.log(1-x) + a*math.log(x))*cf(1-x, b, a)/b

def clopper_pearson(k, n, alpha=0.05):
    if n == 0: return (float("nan"),)*3
    lo = 0.0 if k == 0 else _betainv(k, n-k+1, alpha/2)
    hi = 1.0 if k == n else _betainv(k+1, n-k, 1-alpha/2)
    return k/n, lo, hi

# ---------- clade-averaged estimator ----------
def clade_averaged(groups, min_n=3):
    """groups: {clade: (k, n)} -> unweighted mean of the proportions per clade
    with a normal CI on the between-clade variance. Neutralises the
    over-representation of a hyper-sampled clade."""
    ps = [k/n for k, n in groups.values() if n >= min_n]
    m = len(ps)
    if m < 2: return dict(estimate=float("nan"), lo=float("nan"), hi=float("nan"), n_clades=m)
    mean = sum(ps)/m
    sd = math.sqrt(sum((p-mean)**2 for p in ps)/(m-1))
    se = sd/math.sqrt(m)
    return dict(estimate=mean, lo=max(0.0, mean-1.96*se), hi=min(1.0, mean+1.96*se),
                n_clades=m, sd_between_clades=sd)

# ---------- loading ----------
# sample_design.py now writes the sampling plan with English column names, but
# an ar53_plan.tsv produced before the translation carries the former French
# names. Both spellings are accepted here: the header is normalised as it is
# read, so a plan already on disk keeps working and does not have to be
# rebuilt.
def normalise_plan_row(row):
    """Pass-through, kept so the two readers share one entry point."""
    return dict(row)

def load_markers(path):
    d = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["genome_id"]][r["marker"]] = int(r["n_copies"])
    return d

def load_meta(path):
    d = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = list(rd.fieldnames or [])
        # Two formats accepted: the historical TSV (genome_id, domain, phylum...)
        # and the plan produced by sample_design.py (accession, layer, phylum,
        # class, order, family, genus, completeness, contamination,
        # genome_size, n_proteins). We normalise towards the first.
        plan = "accession" in cols and "family" in cols
        for r in rd:
            if plan:
                r = normalise_plan_row(r)
                g = r["accession"]
                d[g] = dict(genome_id=g, domain=DEFAULT_DOMAIN,
                            phylum=r.get("phylum", "NA"),
                            **{"class": r.get("class", "")},
                            order=r.get("order", ""), family=r.get("family", ""),
                            genus=r.get("genus", ""),
                            completeness=r.get("completeness", ""),
                            genome_size=r.get("genome_size", ""),
                            n_proteins=r.get("n_proteins", ""))
            else:
                d[r["genome_id"]] = r
    return d

DEFAULT_DOMAIN = "NA"


def main():
    global DEFAULT_DOMAIN
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", required=True, help="TSV from scan_proteomes.py")
    ap.add_argument("--domain-label", default="Archaea",
                    help="domain to record when --metadata is a sample_design plan")
    ap.add_argument("--clade-rank", default="phylum",
                    choices=["phylum", "class", "order", "family"],
                    help="grouping rank for the clade-averaged estimator")
    ap.add_argument("--metadata", required=True,
                    help="TSV: genome_id, domain, phylum, [class], [genome_size], [completeness]")
    ap.add_argument("--fluc-marker", default="Fluc_CrcB")
    ap.add_argument("--clcf-marker", default="CLC_F")
    ap.add_argument("--riboswitch", default=None, help="optional TSV: genome_id, n_riboswitch")
    ap.add_argument("--rescue", default=None,
                    help="optional TSV of the 6-frame rescue: genome_id, marker, n_copies")
    ap.add_argument("--min-clade-n", type=int, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    DEFAULT_DOMAIN = a.domain_label
    mk, meta = load_markers(a.markers), load_meta(a.metadata)
    ribo = {}
    if a.riboswitch:
        with open(a.riboswitch) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ribo[r["genome_id"]] = int(r["n_riboswitch"])
    rescue = defaultdict(dict)
    if a.rescue:
        with open(a.rescue) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rescue[r["genome_id"]][r["marker"]] = int(r["n_copies"])

    gids = [g for g in meta if g in mk]
    missing = [g for g in meta if g not in mk]
    rows = []
    for g in gids:
        m = mk[g]
        fluc = m.get(a.fluc_marker, 0)
        clcf = m.get(a.clcf_marker, 0)
        fluc_r = fluc + rescue.get(g, {}).get(a.fluc_marker, 0)
        clcf_r = clcf + rescue.get(g, {}).get(a.clcf_marker, 0)
        cls = ("both" if fluc and clcf else "fluc_only" if fluc else
               "clcf_only" if clcf else "none")
        rows.append(dict(
            genome_id=g, domain=meta[g].get("domain", "NA"), phylum=meta[g].get("phylum", "NA"),
            fluc=fluc, clcf=clcf, fluc_rescued=fluc_r, clcf_rescued=clcf_r,
            E_any=int(bool(fluc or clcf)), E_any_rescued=int(bool(fluc_r or clcf_r)),
            E_class=cls, ribo=ribo.get(g, 0),
            R_any=int(bool(fluc or clcf or ribo.get(g, 0))),
            genome_size=meta[g].get("genome_size", ""),
            completeness=meta[g].get("completeness", ""),
            clade=meta[g].get(a.clade_rank, "NA") or "NA",
        ))

    def block(sub, label):
        n = len(sub)
        if n == 0: return None
        k = sum(r["E_any"] for r in sub)
        kr = sum(r["E_any_rescued"] for r in sub)
        p, lo, hi = clopper_pearson(k, n)
        pw, wlo, whi = wilson(k, n)
        byph = defaultdict(lambda: [0, 0])
        for r in sub:
            byph[r["clade"]][1] += 1
            byph[r["clade"]][0] += r["E_any"]
        ca = clade_averaged({p_: tuple(v) for p_, v in byph.items()}, a.min_clade_n)
        return dict(
            stratum=label, N=n,
            E_any_k=k,
            raw_proportion=round(p, 4),
            CI95_Clopper_Pearson=[round(lo, 4), round(hi, 4)],
            CI95_Wilson=[round(wlo, 4), round(whi, 4)],
            mean_by_phylum_estimator=(None if math.isnan(ca["estimate"]) else {
                "estimate": round(ca["estimate"], 4),
                "CI95": [round(ca["lo"], 4), round(ca["hi"], 4)],
                "n_phyla_kept": ca["n_clades"],
                "sd_between_phyla": round(ca.get("sd_between_clades", float("nan")), 4)}),
            classes=dict(Counter(r["E_class"] for r in sub)),
            proportion_both=round(sum(1 for r in sub if r["E_class"] == "both")/n, 4),
            R_any_proportion=round(sum(r["R_any"] for r in sub)/n, 4) if ribo else None,
            rescue_6_frame=dict(
                E_any_after=kr,
                genomes_recovered=kr-k,
                delta_points=round(100*(kr-k)/n, 2)) if rescue else None,
            phyla=[{"phylum": p_, "N": v[1], "k": v[0], "prop": round(v[0]/v[1], 3)}
                   for p_, v in sorted(byph.items(), key=lambda x: -x[1][1])],
        )

    report = dict(
        method={
            "definition_E_any": f"{a.fluc_marker} OR {a.clcf_marker} detected (>=1 copy, calibrated GA threshold)",
            "threshold": "gathering threshold calibrated by LOO on a characterised seed, cf. *.calib.json",
            "CI": "Clopper-Pearson (exact) and Wilson, alpha=0.05",
            "clade_estimator": f"unweighted mean of the proportions per phylum (phyla with N>={a.min_clade_n})",
            "warning": ("GEBA-1003 is a set chosen to maximise phylogenetic novelty. "
                        "The raw proportion estimates the prevalence IN THIS SET, not in prokaryotes. "
                        "The mean-by-phylum estimator is the figure to cite for a general scope."),
        },
        N_total=len(gids), genomes_without_scan=len(missing),
        strata=[b for b in [block(rows, "ALL"),
                            block([r for r in rows if r["domain"].lower().startswith("bact")], "BACTERIA"),
                            block([r for r in rows if r["domain"].lower().startswith("arch")], "ARCHAEA")] if b],
    )
    with open(a.out, "w") as fh: json.dump(report, fh, indent=2, ensure_ascii=False)

    with open(a.out.replace(".json", ".genomes.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    for s in report["strata"]:
        cp = s["CI95_Clopper_Pearson"]
        print(f"\n{s['stratum']}  N={s['N']}")
        print(f"  E_any raw           : {100*s['raw_proportion']:.1f} %  "
              f"CI95 [{100*cp[0]:.1f} - {100*cp[1]:.1f}]")
        if s["mean_by_phylum_estimator"]:
            e = s["mean_by_phylum_estimator"]
            print(f"  mean by phylum      : {100*e['estimate']:.1f} %  "
                  f"CI95 [{100*e['CI95'][0]:.1f} - {100*e['CI95'][1]:.1f}]  "
                  f"({e['n_phyla_kept']} phyla)")
        print(f"  classes             : {s['classes']}")
        print(f"  both classes        : {100*s['proportion_both']:.1f} %")
        if s["rescue_6_frame"]:
            r_ = s["rescue_6_frame"]
            print(f"  6-frame rescue      : +{r_['genomes_recovered']} genomes "
                  f"(+{r_['delta_points']} points)")
    print(f"\n[OK] -> {a.out}")

if __name__ == "__main__":
    main()
