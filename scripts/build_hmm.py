#!/usr/bin/env python3
"""
Build an equivalog HMM profile and calibrate its gathering threshold (GA).

Steps: align the positive seed (biochemically characterised sequences), run
hmmbuild, then calibrate GA by leave-one-out cross-validation of the positives
against a set of related but functionally distinct families.

With the default rule, GA is the midpoint in bits between the highest-scoring
confirmed negative and the lowest leave-one-out positive. Only the calibration
data enter that choice, so genera held out of the seed stay a test rather than
a tuning criterion. The alternative rule puts GA at a fixed fraction of the
lowest leave-one-out score.

"Separation" is judged against the confirmed negatives alone, those that are
both manually reviewed and carry an established function name. The remaining
negatives are reported as unlabelled rather than counted as false positives,
and the number of them scoring above GA is reported too, since a vague
annotation is not evidence of a different function.

The threshold, the minimum leave-one-out score and the maximum negative score
are written to a .calib.json file next to the model, which is what lets the
published figure be recomputed.

  python3 scripts/build_hmm.py --name CLC_F --positives seeds/clcf_pos.faa \
      --negatives seeds/negatives_clc.faa --outdir results/hmm --rule midpoint
"""
import argparse
import re, json, subprocess, sys, tempfile, os, shutil
from pathlib import Path

def sh(cmd, **kw):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"[FAIL] {cmd}\n{r.stderr[:2000]}")
    return r.stdout

def read_fasta(p):
    name, seq, out = None, [], []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if name: out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name: out.append((name, "".join(seq)))
    return out

def write_fasta(recs, p):
    with open(p, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")

def align(fa, out, threads=4):
    # --thread 1: MAFFT in --localpair mode is not deterministic when
    # multi-threaded. Three runs of the same command gave loo_min = 173.1 /
    # 174.3 / 173.8 and 2,776 / 2,758 / 2,810 scored negatives. One bit changes
    # no conclusion, but it does stop the deposit from reproducing its own
    # figures. Single-threaded costs little on seeds of 22 to 214 sequences.
    sh(f"mafft --maxiterate 1000 --localpair --anysymbol --quiet --thread 1 {fa} > {out}")

def hmmbuild(aln, hmm, name):
    sh(f"hmmbuild --amino -n {name} -o /dev/null {hmm} {aln}")

def best_scores(hmm, fa):
    """Maximum bit score per target sequence."""
    with tempfile.NamedTemporaryFile(suffix=".tbl", delete=False) as t:
        tbl = t.name
    sh(f"hmmsearch --max -E 1000 --tblout {tbl} -o /dev/null {hmm} {fa}")
    best = {}
    for line in open(tbl):
        if line.startswith("#"): continue
        f = line.split()
        if len(f) < 6: continue
        tgt, score = f[0], float(f[5])
        if score > best.get(tgt, -1e9): best[tgt] = score
    os.unlink(tbl)
    return best

def read_alignment(p):
    """(name, aligned sequence); the name is the first word of the header."""
    return read_fasta(p)

def sub_alignment(recs, i, out):
    """Alignment of every record except the i-th, with the columns that have
    become entirely empty removed. Deterministic and immediate: the same
    alignment minus one row, rather than a fresh alignment. Realigning instead
    costs O(n^2) per removal, so O(n^3) for the whole cross-validation, which
    214 sequences make impractical."""
    kept = [r for j, r in enumerate(recs) if j != i]
    L = len(kept[0][1])
    useful = [c for c in range(L) if any(s[c] != "-" for _, s in kept)]
    write_fasta([(n, "".join(s[c] for c in useful)) for n, s in kept], out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positives", required=True, help="seed FASTA, characterised sequences")
    ap.add_argument("--negatives", required=True, help="FASTA of related non-functional families")
    ap.add_argument("--name", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rule", choices=["midpoint", "margin"], default="midpoint",
                    help="'midpoint' (default): GA = the midpoint, in bits, of the window "
                         "[best confirmed negative ; lowest LOO positive]. "
                         "This rule uses ONLY the calibration data, so "
                         "the set of unknown genera remains a test and not a "
                         "tuning criterion. 'margin': former rule, "
                         "GA = margin * min(LOO).")
    ap.add_argument("--margin", type=float, default=0.90,
                    help="used only with --rule margin ; default 0.90")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--confirmed-prefix", default=None,
                    help="header prefix designating the CONFIRMED negatives "
                         "inside --negatives (e.g. 'sp|' for the only "
                         "manually reviewed part of UniProt). The threshold is "
                         "constrained by those alone ; the others are reported "
                         "as unlabelled, not as false positives.")
    ap.add_argument("--loo-mode", choices=["realign", "prune"], default="realign",
                    help="'realign' (default, exact): each removal triggers a "
                         "full MAFFT realignment. Costs O(n^3) and becomes "
                         "impractical beyond a hundred or so sequences. "
                         "'prune': the sub-alignment is obtained by removing "
                         "one row from the full alignment and deleting the "
                         "columns that have become empty. Instantaneous and deterministic. "
                         "to be validated against 'realign' on a seed small enough "
                         "to support both before being used on its own.")
    ap.add_argument("--loo-cache", default=None,
                    help="JSON file in which to store/reload the cross-validation ; "
                         "it costs 1 alignment per seed sequence, so it is worth "
                         "computing it only once")
    ap.add_argument("--confirmed-pattern", default=None,
                    help="REGULAR EXPRESSION applied to the FULL HEADER of the "
                         "negatives to designate those whose DIFFERENT function "
                         "is established. To be preferred over --confirmed-prefix when the "
                         "database section is not enough: Swiss-Prot confirms "
                         "that a protein exists and which family it belongs to, "
                         "not which ion it transports. CLC-F case: sp|Q4VFY6 is a "
                         "Swiss-Prot ClC homolog of UNKNOWN specificity, which scores "
                         "524 bits ; it is not a confirmed negative.")
    ap.add_argument("--confirmed-negatives", default=None,
                    help="optional FASTA of CONFIRMED negatives. When it is "
                         "supplied, it is THAT set which constrains the threshold ; the other "
                         "negatives are only reported. This is the lesson of "
                         "CLC-F: a vague annotation is not a negative.")
    a = ap.parse_args()

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    pos = read_fasta(a.positives)
    if len(pos) < 4:
        sys.exit(f"[FAIL] seed too small ({len(pos)} sequences). Minimum 4, aim for >=8.")

    # Final model, built on the whole seed.
    aln_all = out / f"{a.name}.seed.aln"
    hmm_all = out / f"{a.name}.hmm"
    align(a.positives, aln_all, a.threads)
    hmmbuild(aln_all, hmm_all, a.name)

    # Leave-one-out cross-validation.
    loo, loo_provenance = {}, a.loo_mode
    if a.loo_cache and Path(a.loo_cache).exists():
        cache = json.loads(Path(a.loo_cache).read_text())
        # Current format: {"mode": ..., "scores": {...}}; older caches are flat.
        c_mode = cache.get("mode", "unknown") if "scores" in cache else "unknown"
        c_sc = cache.get("scores", cache) if "scores" in cache else cache
        if set(c_sc) == {n for n, _ in pos}:
            loo, loo_provenance = c_sc, c_mode
            print(f"[LOO] resuming from cache: {len(loo)} scores, mode {c_mode}",
                  file=sys.stderr)
            if c_mode != a.loo_mode:
                print(f"[LOO] WARNING: cache produced in mode '{c_mode}', "
                      f"whereas '{a.loo_mode}' is requested. The published threshold will not be "
                      f"reproducible with the command displayed. Delete "
                      f"{a.loo_cache} to recompute.", file=sys.stderr)
    tmp = Path(tempfile.mkdtemp())
    recs_aln = read_alignment(aln_all) if a.loo_mode == "prune" else None
    to_do = pos if not loo else []
    if to_do:
        print(f"[LOO] mode {a.loo_mode}, {len(to_do)} leave-one-out steps", file=sys.stderr)
    for i, (n, s) in enumerate(to_do):
        f_sub, f_held = tmp / "sub.fa", tmp / "held.fa"
        write_fasta([(n, s)], f_held)
        if a.loo_mode == "prune":
            sub_alignment(recs_aln, i, tmp / "sub.aln")
        else:
            write_fasta([r for j, r in enumerate(pos) if j != i], f_sub)
            align(f_sub, tmp / "sub.aln", a.threads)
        hmmbuild(tmp / "sub.aln", tmp / "sub.hmm", "loo")
        sc = best_scores(tmp / "sub.hmm", f_held)
        loo[n] = sc.get(n, float("-inf"))
        if (i + 1) % 25 == 0:
            print(f"   {i+1}/{len(to_do)}", file=sys.stderr)
    shutil.rmtree(tmp, ignore_errors=True)
    if a.loo_cache and loo:
        Path(a.loo_cache).write_text(json.dumps(
            {"mode": loo_provenance, "scores": loo}, indent=1))

    neg = best_scores(hmm_all, a.negatives)
    n_neg_total = len(read_fasta(a.negatives))
    # First word of the header -> full header. hmmsearch keeps only the first
    # word in column 0 of --tblout, whereas the pattern applies to the
    # description, so both forms are needed.
    headers = {}
    for _l in open(a.negatives):
        if _l.startswith(">"):
            _h = _l[1:].rstrip()
            headers[_h.split()[0]] = _h

    # Where the best-ranked negatives come from. A UniProt header starts with
    # sp| (Swiss-Prot, manually reviewed) or tr| (TrEMBL, automatic
    # annotation). A well-ranked TrEMBL negative is often a mislabelled true
    # positive rather than a false positive of the model.
    def source(name):
        return "SwissProt" if name.startswith("sp|") else (
            "TrEMBL" if name.startswith("tr|") else "other")
    by_source = {}
    for name, sc in neg.items():
        by_source.setdefault(source(name), []).append((sc, name))
    diag = {}
    for src, lst in by_source.items():
        lst.sort(reverse=True)
        diag[src] = dict(n=len(lst), max=lst[0][0],
                         median=lst[len(lst) // 2][0],
                         top10=[[round(x, 1), n] for x, n in lst[:10]])
    print("[negatives] maximum score by annotation source:", file=sys.stderr)
    for src, d in sorted(diag.items()):
        print(f"   {src:<10} n={d['n']:<6} max={d['max']:8.1f}  "
              f"median={d['median']:8.1f}", file=sys.stderr)

    confirmed = None
    if a.confirmed_pattern:
        # A confirmed negative must satisfy two independent conditions:
        #   (a) manually reviewed              -> --confirmed-prefix (e.g. 'sp|')
        #   (b) an established function name   -> --confirmed-pattern
        # Neither condition suffices alone, and CLC-F fails each of them in turn:
        # sp|Q4VFY6 is manually reviewed but of unknown ionic specificity
        # (526 bits); tr|U5MTF6 carries a precise name that was propagated
        # automatically (524 bits). Both are excluded from the negatives.
        rx = re.compile(a.confirmed_pattern, re.I)
        pref = a.confirmed_prefix or ""
        confirmed = {n: v for n, v in neg.items()
                     if n.startswith(pref) and rx.search(headers.get(n, n))}
        others = {n: v for n, v in neg.items() if n not in confirmed}
        omax = max(others.values()) if others else float("-inf")
        print(f"[confirmed negatives] prefix '{pref}' AND pattern: {len(confirmed)} confirmed "
              f"(max {max(confirmed.values()) if confirmed else float('-inf'):.1f}) ; "
              f"{len(others)} others treated as unlabelled (max {omax:.1f})",
              file=sys.stderr)
        rest = sorted(((v, n) for n, v in others.items()), reverse=True)[:10]
        print("   ten best-ranked unlabelled entries (mis-annotated candidates):",
              file=sys.stderr)
        for v, n in rest:
            print(f"     {v:8.1f}  {headers.get(n, n)[:88]}", file=sys.stderr)
    elif a.confirmed_prefix:
        confirmed = {n: v for n, v in neg.items() if n.startswith(a.confirmed_prefix)}
        others = {n: v for n, v in neg.items() if not n.startswith(a.confirmed_prefix)}
        print(f"[confirmed negatives] prefix '{a.confirmed_prefix}': {len(confirmed)} "
              f"confirmed (max {max(confirmed.values()):.1f}) ; {len(others)} others "
              f"treated as unlabelled (max "
              f"{max(others.values()) if others else float('-inf'):.1f})",
              file=sys.stderr)
    elif a.confirmed_negatives:
        confirmed = best_scores(hmm_all, a.confirmed_negatives)
        print(f"[confirmed negatives] {len(confirmed)} scores, max "
              f"{max(confirmed.values()) if confirmed else float('-inf'):.1f}", file=sys.stderr)

    loo_min = min(loo.values())
    constraining = confirmed if confirmed is not None else neg
    neg_max = max(constraining.values()) if constraining else float("-inf")

    # negative_max is computed on the confirmed negatives only, which for CLC-F
    # is 86 sequences out of 3351. Recomputing the maximum over the whole
    # negatives file gives 605.5 instead of 160.4, so the record has to say
    # which subset the number came from: a field that cannot be re-derived from
    # the record is not traceable. The filter itself and the number of
    # sequences it kept are therefore stored alongside the number.
    negative_filter = {
        "confirmed_prefix": a.confirmed_prefix,
        "confirmed_pattern": a.confirmed_pattern,
        "criterion": ("manually reviewed (prefix) AND established function name (pattern)"
                      if a.confirmed_pattern else
                      ("manually reviewed (prefix) only" if a.confirmed_prefix else
                       "none: every negative is constraining")),
        "n_confirmed": (len(constraining) if constraining is not None else None),
        # Named for the cutoff it was measured at: this script reports at
        # -E 1000, while the provenance pass uses -E 1e9 and so sees slightly
        # more sequences. The two are correct counts of different populations,
        # and one name for both would look like a contradiction.
        "n_scored_at_E1000": len(neg),
        "reporting_cutoff": "-E 1000 (build pass)",
        "fraction_kept": (round(len(constraining) / len(neg), 6)
                          if constraining and neg else None),
    }
    # How many negatives of any status score above the shipped threshold, and
    # how high the best one reaches. Reported, not filtered out.
    unfiltered_max = max(neg.values()) if neg else float("-inf")
    if a.rule == "midpoint" and neg_max != float("-inf") and neg_max < loo_min:
        ga = 0.5 * (neg_max + loo_min)
    else:
        if a.rule == "midpoint":
            print("[rule] empty window or no negatives -> falling back on the margin",
                  file=sys.stderr)
        ga = a.margin * loo_min

    verdict = "OK" if ga > neg_max else "SEPARATION IMPOSSIBLE"
    above = {n: v for n, v in neg.items() if v > ga}
    n_above = len(above)
    n_above_reviewed = sum(1 for n in above if n.startswith("sp|"))
    if n_above:
        print(f"[separation] {n_above} negative(s) of any status score above GA "
              f"({n_above_reviewed} manually reviewed); best {max(above.values()):.1f}",
              file=sys.stderr)
    calib = {
        "model": a.name,
        "n_positives": len(pos),
        "n_negatives_tested": n_neg_total,
        "n_negatives_scored": len(neg),
        "negatives_by_source": diag,
        "constraining_negatives": negative_filter["criterion"],
        "negative_filter": negative_filter,
        "threshold_window": (None if confirmed is None else
                             [max(confirmed.values()), min(loo.values())]),
        "loo_scores": loo,
        "loo_min": loo_min,
        "negative_max": None if neg_max == float("-inf") else neg_max,
        "unfiltered_negative_max": (None if unfiltered_max == float("-inf")
                                    else unfiltered_max),
        "gathering_threshold": ga,
        "rule": a.rule,
        "loo_mode": loo_provenance,
        "loo_mode_requested": a.loo_mode,
        "margin_factor": a.margin,
        "separation": verdict,
        "separation_scope": ("no negative of any status scores above GA"
                             if n_above == 0 else
                             "separation is against the CONFIRMED negatives only"),
        "n_negatives_above_GA": n_above,
        "n_negatives_above_GA_reviewed": n_above_reviewed,
        "loo_sensitivity_at_GA": sum(1 for v in loo.values() if v >= ga) / len(loo),
    }
    (out / f"{a.name}.calib.json").write_text(json.dumps(calib, indent=2))

    # Inject GA into the HMM header so hmmsearch --cut_ga can use it.
    txt = hmm_all.read_text().splitlines(True)
    with open(hmm_all, "w") as fh:
        for line in txt:
            fh.write(line)
            if line.startswith("CKSUM"):
                fh.write(f"GA    {ga:.2f} {ga:.2f}\n")
                fh.write(f"TC    {loo_min:.2f} {loo_min:.2f}\n")
                if neg_max != float("-inf"):
                    fh.write(f"NC    {neg_max:.2f} {neg_max:.2f}\n")

    print(json.dumps({k: v for k, v in calib.items() if k != "loo_scores"}, indent=2))
    if verdict != "OK":
        print(f"\n[ALERT] GA={ga:.1f} <= max negative score={neg_max:.1f}. "
              "The clade is not separable by HMM alone: add a residue filter.",
              file=sys.stderr)

if __name__ == "__main__":
    main()
