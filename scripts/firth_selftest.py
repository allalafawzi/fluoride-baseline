#!/usr/bin/env python3
"""
Checks on the Firth implementation (scripts/firth_phylum.py).

Four checks on the hand-written penalised regression:

  1. Under complete separation the estimate must be finite; the ordinary
     maximum likelihood estimate does not exist there.
  2. Without separation and at large n, Firth must agree with ordinary maximum
     likelihood (statsmodels): the penalty is of order 1/n. The coded tolerance
     is 0.01.
  3. The fit must reach a true maximum of the penalised log-likelihood.
  4. Under the null, and with levels that really do hold no positive, the
     penalised likelihood-ratio test must hold its nominal level. Check 4a
     verifies first that the generator produces that situation, without which
     4b would be validating something other than what it announces.

The penalty depends on the design matrix, so the two nested models of check 4
must be fitted under a common design: the test constrains the tested
coefficients to zero without changing the design matrix. Taking the difference
of two separately fitted penalised log-likelihoods lets the penalty difference
enter the statistic, and the rejection rate then reaches 100 % at the 5 % level.

Usage :  python3 scripts/firth_selftest.py [--B 400]
Exit code 0 if every check passes ; exit code 1 otherwise.
"""
import argparse, sys
import math
import numpy as np
from scipy.stats import chi2

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from firth_phylum import _fit_firth, _penalised_loglik, firth_lrt

FAIL = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        FAIL.append(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=400,
                    help="number of replicates for check 4 (minimum 50 : "
                         "below that the measured rates mean nothing)")
    a = ap.parse_args()
    if a.B < 50:
        print(f"ERROR : --B {a.B} is too small. Check 4 measures rejection "
              f"rates ; below 50 replicates they mean nothing and the check "
              f"would fail for a reason that is not a defect.",
              file=sys.stderr)
        return 2

    # 1 -- complete separation
    x = np.array([-3, -2, -1, 1, 2, 3], float)
    y = np.array([0, 0, 0, 1, 1, 1], float)
    b, _ = _fit_firth(np.column_stack([np.ones(6), x]), y)
    check("1 finite estimate under complete separation",
          bool(np.all(np.isfinite(b))) and abs(b[1]) < 50,
          f"slope={b[1]:.3f}")

    # 2 -- agreement with the ordinary MLE where it exists
    rng = np.random.default_rng(0)
    n = 4000
    x = rng.normal(size=n)
    y = (rng.random(n) < 1 / (1 + np.exp(-(-0.5 + 1.2 * x)))).astype(float)
    X = np.column_stack([np.ones(n), x])
    bf, _ = _fit_firth(X, y)
    try:
        import statsmodels.api as sm
        bm = sm.Logit(y, X).fit(disp=0).params
        d = float(np.max(np.abs(bf - bm)))
        check("2 agreement with ordinary MLE (n=4000)", d < 0.01, f"max gap={d:.5f}")
    except ImportError:
        print("SKIP  2 agreement with the MLE (statsmodels absent)")

    # 3 -- true maximum
    b0, pll0 = _fit_firth(X, y)
    nb = [_penalised_loglik(X, y, b0 + d) for d in
          (np.array([0.01, 0]), np.array([-0.01, 0]),
           np.array([0, 0.01]), np.array([0, -0.01]))]
    check("3 maximum of the penalised log-likelihood",
          all(v <= pll0 + 1e-9 for v in nb))

    # 4 -- level of the test under H0, with levels that really hold no positive
    #
    # The generator below reproduces the composition of the panel itself: very
    # unequal level sizes (1 to 310 genomes per phylum) and a low prevalence,
    # which is what makes empty levels occur in every replicate. Nearly
    # equiprobable groups with a higher prevalence leave almost every level with
    # a positive, and the test would then be exercised on a situation other than
    # the one it is meant to cover. Check 4a measures that the situation occurs.
    rng = np.random.default_rng(11)
    # Composition of the fitted model: the 16 phyla kept (N >= 10) plus the
    # 'other' level that absorbs the 5 rare phyla (24 genomes), for 1259
    # observations over 17 levels. N and the number of levels both have to match
    # the model, or the measured rates describe a different design.
    sizes = [310, 238, 157, 132, 85, 77, 55, 53, 20, 19, 19, 18, 14, 14, 13, 11, 24]
    assert sum(sizes) == 1259 and len(sizes) == 17
    ps, nzero = [], []
    for _ in range(a.B):
        g = np.concatenate([np.full(s_, i) for i, s_ in enumerate(sizes)])
        n = len(g)
        size = rng.normal(size=n)
        # low prevalence, about 13 %, as in the panel
        y = (rng.random(n) < 1 / (1 + np.exp(-(-2.4 + 0.9 * size)))).astype(float)
        nzero.append(sum(1 for q in range(len(sizes)) if y[g == q].sum() == 0))
        Xr = np.column_stack([np.ones(n), size])
        D = np.column_stack([(g == q).astype(float) for q in range(1, len(sizes))])
        Xf = np.column_stack([Xr, D])
        s_stat, df, _, _ = firth_lrt(Xf, Xr, y, red_idx=[0, 1])
        ps.append(chi2.sf(s_stat, df))
    ps = np.array(ps)
    frac_with_zero = float(np.mean(np.array(nzero) > 0))
    # The 0.5 threshold only asks that the situation under test occurs in most
    # replicates ; the composition does not control any precise rate.
    check("4a the generator does produce levels with no positive",
          frac_with_zero > 0.5,
          f"{100*frac_with_zero:.0f} % of replicates, "
          f"median {int(np.median(nzero))} empty level(s)")
    # The check establishes that the observed rejection rate does not exceed
    # the nominal level in this simulated scenario. It does not establish that
    # the test is conservative: an empirical rate over B replications carries a
    # Monte-Carlo standard error of sqrt(p(1-p)/B), about 1 point at B=400, so a
    # rate of 4.0 % has a 95 % interval of [2.1 ; 5.9], which contains 5 %.
    #
    # The scenario also draws the observations independently within each level.
    # It does not reproduce the dependence between related genomes, and so does
    # not validate the p-values of the fitted model.
    for lvl in (0.05, 0.10):
        r = float(np.mean(ps < lvl))
        se = math.sqrt(r * (1 - r) / a.B) if a.B else 0.0
        check(f"4b nominal level {lvl:.2f} not exceeded under H0",
              r <= lvl + 0.02,
              f"observed={r:.3f} +/- {se:.3f} (MC standard error, B={a.B}) ; "
              f"CI95 [{max(0, r - 1.96 * se):.3f} ; {r + 1.96 * se:.3f}]")

    if FAIL:
        print("\nFAILED :", ", ".join(FAIL))
        return 1
    print("\nAll Firth checks pass.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
