#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Power calculation for the FAcD x fluoride exporter association test.

  python3 scripts/power.py --out results/power

N bacterial genomes are drawn and scored on two binary characters:
  FAcD+/-    does the genome carry the defluorinase (equivalog PF00561)?
  EXPORT+/-  does it carry a fluoride export system
             (Fluc/CrcB PF02537, CLC-F, riboswitch RF01734)?

The test has to separate two hypotheses of opposite sign:
  H1 co-retention : P(EXPORT+|FAcD+) > P(EXPORT+|FAcD-)
  H2 purge        : P(EXPORT+|FAcD+) < P(EXPORT+|FAcD-)
                    (O'Connor et al. 2024: adaptive evolution deletes the
                     defluorinase once export is sufficient)

The script answers three separate questions, labelled E1 to E3 in the output
files and in the "exp" field of validation_simulation.json.

E1, how many genomes: variance inflation from relatedness. Clade effects are
drawn independently for the two characters, so there is pseudo-replication but
no confounding. The design-effect formula is validated by simulation, the type
I error rate is checked, and N is tabulated.

E2, which test: confounding by clade. The two characters share their clade
effects, with beta = 0, so no real association exists. The false positive rate
of a co-occurrence test on the terminal states measures the cost of ignoring
the phylogeny, and decides whether a contingency table suffices or a
branch-event method (Pagel, EvolCCM) is required.

E3, which sampling design: the currency is events, not genomes. A branch-based
method draws its information from the number of independent origins of FAcD
sampled. The expected number of positive clades detected is
(N/k)*phi*(1-(1-psi)^k), strictly decreasing in k, so the best design is one
genome per family spread as widely as possible.

No input data are required. load_taxo looks for the bacterial GTDB metadata in
--meta-dir (data/gtdb/bac120_metadata_r*.tsv.gz) and uses the real family and
order sizes when it finds the file. When it does not, the run falls back to a
synthetic taxonomy of 4000 families and 113 073 species, reported in the
output as source = MODEL. Numbers produced without the metadata file describe
that model and not any real domain, so they are not interchangeable with the
ones in the paper.
"""

import argparse, gzip, json, math, os, sys
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import norm, chi2

Z = norm.ppf


# Taxonomy and sampling designs.

class Taxo:
    """Hierarchy order > family > species. `fam_size[i]` species in the
    family i, `fam_order[i]` index of its order."""

    def __init__(self, fam_size, fam_order, source):
        self.fam_size = np.asarray(fam_size, dtype=int)
        self.fam_order = np.asarray(fam_order, dtype=int)
        self.n_ord = int(self.fam_order.max()) + 1
        self.source = source

    @property
    def n_fam(self):
        return self.fam_size.size

    def capacity(self, cap):
        return int(np.minimum(self.fam_size, cap).sum()) if cap else int(self.fam_size.sum())

    def alloc(self, N, design="balanced", cap=1, rng=None):
        """Number of genomes drawn per family."""
        if design == "proportional":
            p = self.fam_size / self.fam_size.sum()
            idx = rng.choice(self.n_fam, size=int(N), replace=True, p=p)
            return np.bincount(idx, minlength=self.n_fam)
        limit = np.minimum(self.fam_size, int(cap))
        N = int(min(N, limit.sum()))
        lo, hi = 0, int(limit.max())
        while lo < hi:                       # largest uniform layer
            mid = (lo + hi + 1) // 2
            if np.minimum(limit, mid).sum() <= N:
                lo = mid
            else:
                hi = mid - 1
        out = np.minimum(limit, lo)
        rest = N - int(out.sum())
        if rest > 0:
            out[np.flatnonzero(out < limit)[:rest]] += 1
        return out

    def order_counts(self, fam_counts):
        return np.bincount(self.fam_order, weights=fam_counts,
                           minlength=self.n_ord)


def kish(counts):
    """Effective cluster size, sum(n^2)/sum(n). With clusters of unequal size
    it is this quantity, not the mean size, that drives the inflation."""
    c = np.asarray(counts, dtype=float)
    c = c[c > 0]
    return float((c ** 2).sum() / c.sum()) if c.size else 1.0


def deff_multi(m_fam, m_ord, icc_fam, icc_ord):
    """Design effect with two nested levels (generalised Kish)."""
    return 1.0 + max(0.0, m_fam - 1) * icc_fam + max(0.0, m_ord - 1) * icc_ord


def load_taxo(meta_dir="data/gtdb"):
    if not os.path.isdir(meta_dir):
        return None
    cand = [f for f in os.listdir(meta_dir)
            if f.startswith("bac120_metadata") and f.endswith(".tsv.gz")]
    if not cand:
        return None
    path = os.path.join(meta_dir, sorted(cand)[-1])
    fam = Counter(); fo = {}
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            try:
                i_tax = head.index("gtdb_taxonomy")
                i_rep = head.index("gtdb_representative")
            except ValueError:
                return None
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(i_tax, i_rep):
                    continue
                if p[i_rep].strip().lower() not in ("t", "true"):
                    continue
                r = p[i_tax].split(";")
                f = next((x for x in r if x.startswith("f__")), None)
                o = next((x for x in r if x.startswith("o__")), None)
                if f and f != "f__" and o and o != "o__":
                    fam[f] += 1
                    fo[f] = o
    except OSError:
        return None
    if not fam:
        return None
    keys = sorted(fam)
    ords_ = sorted({fo[k] for k in keys})
    oi = {o: i for i, o in enumerate(ords_)}
    return Taxo([fam[k] for k in keys], [oi[fo[k]] for k in keys],
                os.path.basename(path))


def modelled_taxo(n_fam=4000, n_ord=1500, n_species=113000, seed=0):
    """Tracing of GTDB: power law on the family sizes (a few enormous
    families, long tail), families distributed among orders."""
    rng = np.random.default_rng(seed)
    raw = rng.pareto(1.15, n_fam) + 1.0
    sizes = np.maximum(1, np.round(raw / raw.sum() * n_species)).astype(int)
    w = rng.pareto(1.5, n_ord) + 1.0
    fam_ord = rng.choice(n_ord, size=n_fam, p=w / w.sum())
    fam_ord = np.unique(fam_ord, return_inverse=True)[1]
    return Taxo(sizes, fam_ord, "MODEL (4000 families / 1500 orders / 113000 spp.)")


# E1: how many genomes, given the variance inflation from relatedness.

def power_two_prop(N, pi, p0, p1, alpha, D=1.0):
    """Two-sided power for two proportions under unbalanced allocation.
    Checked against Fleiss: 0.50/0.60 at power .80 gives 388 per group,
    0.30/0.50 at power .90 gives 124 per group (the formula without the
    continuity correction)."""
    n1, n0 = N * pi / D, N * (1 - pi) / D
    if n1 < 2 or n0 < 2:
        return 0.0
    pb = (n1 * p1 + n0 * p0) / (n1 + n0)
    se0 = math.sqrt(pb * (1 - pb) * (1 / n1 + 1 / n0))
    se1 = math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    return float(norm.cdf((abs(p1 - p0) - Z(1 - alpha / 2) * se0) / se1)) if se1 > 0 else 0.0


def deff_at(tx, N, design, cap, icc_fam, icc_ord, rng):
    c = tx.alloc(N, design, cap, rng)
    return deff_multi(kish(c), kish(tx.order_counts(c)), icc_fam, icc_ord), c


def n_required(tx, pi, p0, p1, alpha, target, icc_fam, icc_ord,
               design="balanced", cap=1, rng=None):
    if not (0 < p1 < 1) or abs(p1 - p0) < 1e-9:
        return None, None
    hi = tx.capacity(cap) if design == "balanced" else 300_000
    pw = lambda N: power_two_prop(N, pi, p0, p1, alpha,
                                  deff_at(tx, N, design, cap, icc_fam, icc_ord, rng)[0])
    best = pw(hi)
    if best < target:
        return None, round(best, 4)
    lo = 20
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if pw(mid) >= target:
            hi = mid
        else:
            lo = mid
    return int(hi), None


def cluster_robust_chi2(cl, x, y):
    """Difference of proportions with a cluster-linearised (Taylor) variance.
    Every genome is kept, singleton clusters included, which is how a GEE with
    exchangeable correlation behaves; a stratified test would discard them."""
    C = int(cl.max()) + 1
    n1 = np.bincount(cl, weights=x, minlength=C)
    a = np.bincount(cl, weights=x * y, minlength=C)
    n0 = np.bincount(cl, weights=1 - x, minlength=C)
    b = np.bincount(cl, weights=(1 - x) * y, minlength=C)
    N1, N0 = n1.sum(), n0.sum()
    if N1 < 5 or N0 < 5:
        return 0.0
    p1, p0 = a.sum() / N1, b.sum() / N0
    u = (a - p1 * n1) / N1 - (b - p0 * n0) / N0
    occ = int((n1 + n0 > 0).sum())
    if occ < 2:
        return 0.0
    v = occ / (occ - 1.0) * float((u ** 2).sum())
    return float((p1 - p0) ** 2 / v) if v > 0 else 0.0


def _sim_draw(tx, N, design, cap, rng):
    c = tx.alloc(N, design, cap, rng)
    fam_id = np.repeat(np.arange(tx.n_fam), c)
    ord_id = tx.fam_order[fam_id]
    fam_id = np.unique(fam_id, return_inverse=True)[1]
    ord_id = np.unique(ord_id, return_inverse=True)[1]
    return fam_id, ord_id


def _sigma(icc):
    """Logistic threshold model: ICC = s^2/(s^2 + pi^2/3)."""
    return math.sqrt(icc / (1 - icc)) * (math.pi / math.sqrt(3)) if icc > 0 else 0.0


def simulate(tx, N, pi, p0, p1, icc_fam, icc_ord, alpha, n_rep,
             design="balanced", cap=1, seed=1, shared_clade=False):
    """With shared_clade=False (E1) each character gets its own clade effects,
    giving pseudo-replication without confounding. With shared_clade=True (E2)
    the two characters share them, giving phylogenetic confounding."""
    rng = np.random.default_rng(seed)
    fam_id, ord_id = _sim_draw(tx, N, design, cap, rng)
    n = fam_id.size
    Cf, Co = int(fam_id.max()) + 1, int(ord_id.max()) + 1
    crit = chi2.ppf(1 - alpha, 1)
    sf, so = _sigma(icc_fam), _sigma(icc_ord)
    l0, l1 = math.log(p0 / (1 - p0)), math.log(p1 / (1 - p1))
    beta, lpi = l1 - l0, math.log(pi / (1 - pi))

    hits = 0
    for _ in range(n_rep):
        ux = rng.normal(0, sf, Cf)[fam_id] + rng.normal(0, so, Co)[ord_id]
        uy = ux if shared_clade else (rng.normal(0, sf, Cf)[fam_id]
                                      + rng.normal(0, so, Co)[ord_id])
        x = (rng.random(n) < 1 / (1 + np.exp(-(lpi + ux)))).astype(float)
        y = (rng.random(n) < 1 / (1 + np.exp(-(l0 + beta * x + uy)))).astype(float)
        if cluster_robust_chi2(fam_id, x, y) >= crit:
            hits += 1
    return hits / n_rep


# E3: the currency is events, so count the positive clades detected.

def detected_clades(N, k, phi, psi, n_fam_max):
    """Expected number of FAcD+ families detected by drawing k genomes in each
    of N/k families.
      phi = fraction of families carrying FAcD (a proxy for the number of
            independent origins)
      psi = prevalence of FAcD within a positive family
    E = (N/k) * phi * (1 - (1-psi)^k), strictly decreasing in k."""
    F = min(N / k, n_fam_max)
    return F * phi * (1 - (1 - psi) ** k)


# Main program.

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/power")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--power", type=float, default=0.80)
    ap.add_argument("--sim", type=int, default=2000)
    ap.add_argument("--meta-dir", default="data/gtdb")
    ap.add_argument("--icc-split", type=float, default=0.5,
                    help="share of the total ICC assigned to the family level")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(7)
    say = lambda *x: print(*x, file=sys.stderr)

    tx = load_taxo(a.meta_dir) or modelled_taxo()
    real = not tx.source.startswith("MODEL")
    say(f"[taxo] {tx.source}")
    say(f"[taxo] {tx.n_fam} families, {tx.n_ord} orders, {tx.fam_size.sum()} representatives")
    say(f"[taxo] capacity 1/family={tx.capacity(1)}  2/family={tx.capacity(2)}  "
        f"5/family={tx.capacity(5)}")

    split = a.icc_split
    P0, DELTA = [0.75, 0.85, 0.92], [0.05, 0.10, 0.15]
    PI, ICC = [0.02, 0.05, 0.10, 0.20], [0.0, 0.15, 0.30]
    DESIGNS = [("balanced", 1), ("balanced", 3), ("proportional", None)]

    # E1: the table of required N.
    rows = []
    for p0 in P0:
        for d in DELTA:
            for sgn, lab in ((+1, "H1"), (-1, "H2")):
                p1 = round(p0 + sgn * d, 4)
                for pi in PI:
                    for icc in ICC:
                        i_f, i_o = icc * split, icc * (1 - split)
                        for design, cap in DESIGNS:
                            base = dict(p0=p0, delta=d, direction=lab, p1=p1,
                                        pi=pi, icc=icc, design=design,
                                        cap=cap or "")
                            if not (0 < p1 < 1):
                                rows.append({**base, "verdict": "arithmetic ceiling"})
                                continue
                            N, best = n_required(tx, pi, p0, p1, a.alpha, a.power,
                                                 i_f, i_o, design, cap, rng)
                            if N is None:
                                rows.append({**base, "verdict":
                                    "saturated" if design == "proportional"
                                    else "insufficient capacity",
                                    "max_power": best})
                                continue
                            D, c = deff_at(tx, N, design, cap, i_f, i_o, rng)
                            rows.append({**base, "N": N, "verdict": "ok",
                                         "deff": round(D, 2),
                                         "m_fam": round(kish(c), 2),
                                         "m_ord": round(kish(tx.order_counts(c)), 2),
                                         "n_facd_pos": int(round(N * pi)),
                                         "residual_H0": int(round(N * pi * (1 - p0))),
                                         "residual_Ha": int(round(N * pi * (1 - p1)))})
    cols = ["p0", "delta", "direction", "p1", "pi", "icc", "design", "cap", "N",
            "deff", "m_fam", "m_ord", "n_facd_pos", "residual_H0", "residual_Ha",
            "max_power", "verdict"]
    with open(f"{a.out}/power_grid.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    json.dump(dict(source=tx.source, real_taxonomy=real, alpha=a.alpha,
                   target_power=a.power, icc_split=split, rows=rows),
              open(f"{a.out}/power_grid.json", "w"), indent=2)

    # E1b: validation of the design-effect formula by simulation.
    val = []
    if a.sim > 0:
        say("\n=== E1  independent clade effects (pseudo-replication only) ===")
        for icc in (0.0, 0.15, 0.30):
            i_f, i_o = icc * split, icc * (1 - split)
            for design, cap in DESIGNS:
                t1 = simulate(tx, 3000, .05, .85, .85, i_f, i_o, a.alpha, a.sim,
                              design, cap, seed=101, shared_clade=False)
                val.append(dict(exp="E1", type="type_I_error", icc=icc,
                                design=design, cap=cap, N=3000,
                                nominal=a.alpha, observed=round(t1, 4)))
                say(f"[H0] ICC={icc:<4} {design}/{cap}: alpha observed = {t1:.4f}")
        for p0, d, pi, icc in ((.85, .10, .05, .15), (.85, .10, .05, .30),
                               (.85, .15, .05, .15), (.85, .10, .10, .15),
                               (.92, .10, .05, .15)):
            i_f, i_o = icc * split, icc * (1 - split)
            N, _ = n_required(tx, pi, p0, p0 - d, a.alpha, a.power, i_f, i_o,
                              "balanced", 1, rng)
            if not N:
                continue
            emp = simulate(tx, N, pi, p0, p0 - d, i_f, i_o, a.alpha, a.sim,
                           "balanced", 1, seed=42, shared_clade=False)
            val.append(dict(exp="E1", type="power", p0=p0, delta=d, pi=pi,
                            icc=icc, N=N, analytic=a.power, simulated=round(emp, 3)))
            say(f"[H2] p0={p0} d={d} pi={pi} ICC={icc} N={N} -> power {emp:.3f}")

        # E2: confounding by clade.
        say("\n=== E2  shared clade effects, beta=0 (no real association) ===")
        for icc in (0.05, 0.15, 0.30):
            i_f, i_o = icc * split, icc * (1 - split)
            for design, cap in (("balanced", 1), ("proportional", None)):
                fp = simulate(tx, 3000, .05, .85, .85, i_f, i_o, a.alpha, a.sim,
                              design, cap, seed=202, shared_clade=True)
                val.append(dict(exp="E2", type="false_positives_confounding", icc=icc,
                                design=design, cap=cap, N=3000,
                                nominal=a.alpha, observed=round(fp, 4)))
                say(f"[H0 confounded] ICC={icc:<4} {design}/{cap}: "
                    f"false positives = {fp:.4f}  (nominal {a.alpha})")
    json.dump(val, open(f"{a.out}/validation_simulation.json", "w"), indent=2)

    # E3: expected number of clades detected, per design.
    ev = []
    for N in (1000, 2000, 3000, 4000):
        for phi in (0.02, 0.05, 0.10):
            for psi in (0.10, 0.30, 0.70):
                for k in (1, 2, 3, 5, 10, 20):
                    ev.append(dict(N=N, phi=phi, psi=psi, k=k,
                                   clades_detected=round(
                                       detected_clades(N, k, phi, psi, tx.n_fam), 2)))
    json.dump(ev, open(f"{a.out}/events.json", "w"), indent=2)
    with open(f"{a.out}/events.tsv", "w") as fh:
        fh.write("N\tphi\tpsi\tk\tclades_detected\n")
        for r in ev:
            fh.write(f"{r['N']}\t{r['phi']}\t{r['psi']}\t{r['k']}\t{r['clades_detected']}\n")

    # The arithmetic ceiling on the detectable margin, given p0.
    json.dump([dict(p0=p, max_margin_H1=round(1 - p, 3), max_margin_H2=round(p, 3))
               for p in (.70, .75, .80, .85, .88, .90, .92, .95)],
              open(f"{a.out}/margin_ceiling.json", "w"), indent=2)

    say(f"\n[ok] {a.out}/power_grid.tsv ({len(rows)} rows), events.tsv, "
        f"validation_simulation.json")


if __name__ == "__main__":
    main()
