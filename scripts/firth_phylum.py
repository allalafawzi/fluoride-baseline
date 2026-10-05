#!/usr/bin/env python3
"""
Firth penalised logistic regression: Fluc status ~ log10(genome size) + phylum.

The unpenalised model suffers complete separation: three of the retained phyla
(B1Sed10-29, Korarchaeota, Nanohalarchaeota) hold no positive genome. For those
phyla the maximum likelihood estimate does not exist -- the coefficient runs to
minus infinity -- and the likelihood-ratio statistic is not on its reference
chi-square distribution. No adjustment of the degrees of freedom repairs that.

Firth (1993) penalisation adds the term (1/2) * log det I(beta) to the
log-likelihood, where I is the Fisher information. That term is the Jeffreys
prior; it makes the maximum finite even under complete separation, and it
removes the O(1/n) bias of the estimates. The estimates reported here are
therefore finite under separation, but the test that goes with them is the
penalised likelihood-ratio test, not the ordinary one.

References: Firth D. (1993) Bias reduction of maximum likelihood estimates,
Biometrika 80:27-38. Penalised likelihood-ratio test: Heinze & Schemper (2002)
Stat Med 21:2409-2419.

Output: results/gtdb/firth_phylum.json
"""
import argparse, csv, json, math, os, sys
import numpy as np

# Firth core, written out rather than imported: statsmodels carries no Firth
# implementation, and third-party packages (firthlogist) are not guaranteed to
# be installed. The lines below are short enough to check by reading.

def _fit_firth(X, y, max_iter=500, tol=1e-9):
    """Return (beta, penalised log-likelihood). Damped Newton-Raphson."""
    n, p = X.shape
    beta = np.zeros(p)
    pll_old = -np.inf
    for _ in range(max_iter):
        eta = X @ beta
        # numerical stability: clip eta so that exp() cannot overflow
        eta = np.clip(eta, -30.0, 30.0)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = mu * (1.0 - mu)
        w = np.maximum(w, 1e-12)
        XW = X * w[:, None]
        I = X.T @ XW                       # Fisher information
        try:
            Iinv = np.linalg.pinv(I)
        except np.linalg.LinAlgError:
            break
        # h = leverages (diagonal of the weighted hat matrix)
        h = np.einsum('ij,jk,ik->i', XW, Iinv, X)
        # Firth penalised score: U* = X'(y - mu + h*(1/2 - mu))
        U = X.T @ (y - mu + h * (0.5 - mu))
        step = Iinv @ U
        # damping: bound the step so the penalised likelihood increases
        ss = np.max(np.abs(step))
        if ss > 5.0:
            step = step * (5.0 / ss)
        for _halve in range(30):
            cand = beta + step
            pll = _penalised_loglik(X, y, cand)
            if pll >= pll_old - 1e-12:
                break
            step = step / 2.0
        beta = cand
        pll = _penalised_loglik(X, y, beta)
        if abs(pll - pll_old) < tol:
            pll_old = pll
            break
        pll_old = pll
    return beta, pll_old


def _penalised_loglik(X, y, beta):
    eta = np.clip(X @ beta, -30.0, 30.0)
    mu = 1.0 / (1.0 + np.exp(-eta))
    eps = 1e-12
    ll = float(np.sum(y * np.log(mu + eps) + (1 - y) * np.log(1 - mu + eps)))
    w = np.maximum(mu * (1.0 - mu), 1e-12)
    I = X.T @ (X * w[:, None])
    sign, logdet = np.linalg.slogdet(I)
    if sign <= 0:
        return -np.inf
    return ll + 0.5 * logdet


def _fit_firth_constrained(X, y, free_idx, max_iter=500, tol=1e-9):
    """Maximise the penalised log-likelihood with the coefficients outside
    free_idx constrained to zero, the penalty still computed on the full X."""
    n, p = X.shape
    beta = np.zeros(p)
    free = np.array(free_idx, dtype=int)
    pll_old = -np.inf
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -30.0, 30.0)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(mu * (1.0 - mu), 1e-12)
        XW = X * w[:, None]
        I = X.T @ XW
        Iinv = np.linalg.pinv(I)
        h = np.einsum('ij,jk,ik->i', XW, Iinv, X)
        U = X.T @ (y - mu + h * (0.5 - mu))
        Isub = I[np.ix_(free, free)]
        step_sub = np.linalg.pinv(Isub) @ U[free]
        step = np.zeros(p)
        step[free] = step_sub
        ss = np.max(np.abs(step))
        if ss > 5.0:
            step = step * (5.0 / ss)
        cand = beta + step
        for _halve in range(30):
            if _penalised_loglik(X, y, cand) >= pll_old - 1e-12:
                break
            step = step / 2.0
            cand = beta + step
        beta = cand
        pll = _penalised_loglik(X, y, beta)
        if abs(pll - pll_old) < tol:
            pll_old = pll
            break
        pll_old = pll
    return beta, pll_old


def firth_lrt(X_full, X_red, y, red_idx=None):
    """Penalised likelihood-ratio test between two nested models.

    The penalty term (1/2) log det I(beta) depends on the design matrix. Fitting
    the full model and the reduced model separately and subtracting their
    penalised log-likelihoods compares two quantities that do not carry the same
    penalty: the difference in penalty, which measures no signal, enters the
    statistic in full and inflates it. Simulated under the null, that version
    rejects at 100 % for a nominal 5 % level.

    The test used here (Heinze & Schemper 2002) evaluates both models under the
    same penalty, that of the full design: the reduced model is fitted by
    constraining the tested coefficients to zero, without changing the design.
    """
    if red_idx is None:
        red_idx = list(range(X_red.shape[1]))
    _, pll_f = _fit_firth(X_full, y)
    _, pll_r = _fit_firth_constrained(X_full, y, red_idx)
    stat = 2.0 * (pll_f - pll_r)
    df = X_full.shape[1] - len(red_idx)
    return stat, df, pll_f, pll_r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="data/gtdb/ar53_plan.tsv")
    ap.add_argument("--scan", default="results/gtdb/scan_archaea.tsv")
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--min-n", type=int, default=10,
                    help="phyla below this count are pooled into 'other'")
    ap.add_argument("--out", default="results/gtdb/firth_phylum.json")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.plan), delimiter="\t"))[:a.prefix]
    pos = set()
    for r in csv.DictReader(open(a.scan), delimiter="\t"):
        if r["marker"] == a.marker and int(r["n_copies"]) > 0:
            pos.add(r["genome_id"])

    from collections import Counter
    cnt = Counter(r["phylum"] for r in rows)
    kept = sorted(p for p, n in cnt.items() if n >= a.min_n)
    lvl = kept + ["other"]

    y = np.array([1.0 if r["accession"] in pos else 0.0 for r in rows])
    size = np.array([math.log10(max(float(r["genome_size"]), 1.0)) for r in rows])
    ph = [r["phylum"] if r["phylum"] in kept else "other" for r in rows]

    # reference level = the most populated phylum, so that the contrasts are
    # estimated against a well populated category
    ref = max(kept, key=lambda p: cnt[p])
    others = [p for p in lvl if p != ref]

    n = len(rows)
    X_red = np.column_stack([np.ones(n), size])                 # size alone
    D = np.column_stack([[1.0 if p == q else 0.0 for p in ph] for q in others])
    X_full = np.column_stack([X_red, D])                        # size + phylum

    stat_ph, df_ph, pll_f, pll_r = firth_lrt(X_full, X_red, y, red_idx=[0, 1])

    # size effect within the full model: column 1 (log10_size) is constrained
    # to zero while the design stays complete
    free_nosize = [0] + list(range(2, X_full.shape[1]))
    stat_sz, df_sz, _, _ = firth_lrt(X_full, None, y, red_idx=free_nosize)

    from scipy.stats import chi2
    beta, _ = _fit_firth(X_full, y)

    sep = [p for p in lvl
           if sum(1 for r, q in zip(rows, ph) if q == p and r["accession"] in pos) == 0]

    out = {
        "method": "Firth penalised logistic regression (Jeffreys prior), "
                  "penalised likelihood-ratio test",
        "implementation": "scripts/firth_phylum.py, Newton-Raphson with step halving",
        "marker": a.marker,
        "N": n,
        "k": int(y.sum()),
        "min_n_phylum": a.min_n,
        "reference_level": ref,
        "n_levels": len(lvl),
        "levels_absorbed_into_other": sorted(p for p, c in cnt.items() if c < a.min_n),
        "levels_with_zero_positives": sep,
        "phylum_after_size": {
            "lrt_stat": stat_ph, "df": df_ph,
            "p_value": float(chi2.sf(stat_ph, df_ph)),
        },
        "size_after_phylum": {
            "lrt_stat": stat_sz, "df": df_sz,
            "p_value": float(chi2.sf(stat_sz, df_sz)),
        },
        "slope_log10_size": float(beta[1]),
        "odds_ratio_per_decade_size": float(math.exp(beta[1])),
        "penalised_loglik_full": pll_f,
        "penalised_loglik_size_only": pll_r,
        "coefficients": {
            lab: float(v) for lab, v in
            zip(["intercept", "log10_size"] + [f"phylum={p}" for p in others], beta)
        },
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, sort_keys=False)
        fh.write("\n")
    print(json.dumps({k: out[k] for k in
          ("N","k","reference_level","n_levels","levels_with_zero_positives",
           "phylum_after_size","size_after_phylum","odds_ratio_per_decade_size")},
          indent=1))


if __name__ == "__main__":
    main()
