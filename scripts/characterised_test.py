#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Score a model against enzymes whose activity was measured at the bench:
enzymes with a published assay, listed in seeds/characterised_facd.tsv
together with the reference that measured each one.

Recall reported on a held-out set drawn from UniProt is weaker evidence than it
looks. Those entries are labelled "fluoroacetate dehalogenase" mostly by
automatic pipelines that are themselves profile-HMM based, so a profile HMM
validated against them is close to being validated against itself. This script
scores the model against a small, hard set instead.

The set contains a negative on purpose. Q01399 (DehH2) carries EC 3.8.1.3 in
the manually reviewed part of UniProt, yet its own entry states it does not act
on fluoroacetate. A model that accepts it is following the label; a model that
rejects it is following the biology.

Sequences are fetched from UniProt and cached, so the test is rerunnable
offline once the cache exists.
"""
import argparse, csv, json, os, subprocess, sys, tempfile, urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"

def read_table(p):
    with open(p) as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))
    return [r for r in rows if r.get("accession")]

_AA = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

def valid_fasta(txt):
    """Raise unless txt is a real protein FASTA.

    Without this check, an empty body from UniProt would raise nothing: the
    accession would count as fetched, an empty string would go into the query
    file, hmmsearch would see no sequence, and the entry would be scored as
    wrong, which is indistinguishable from the model rejecting it. An obsolete
    or demerged accession, an error page and a truncated download over a slow
    link all produce such a body. A failed download has to surface as untested,
    never as a verdict about the model."""
    if not txt or not txt.strip():
        raise ValueError("empty body (obsolete, demerged or deleted accession)")
    if not txt.lstrip().startswith(">"):
        raise ValueError("not FASTA (error page or redirect?)")
    lines = txt.strip().splitlines()
    seq = "".join(lines[1:]).replace(" ", "").upper()
    if len(seq) < 50:
        raise ValueError(f"only {len(seq)} residues (truncated download?)")
    bad = set(seq) - _AA
    if bad:
        raise ValueError(f"non-amino-acid characters {sorted(bad)[:5]}")
    return txt

def fetch(acc, cache):
    f = Path(cache) / f"{acc}.fasta"
    if f.exists() and f.stat().st_size:
        try:
            return valid_fasta(f.read_text())
        except ValueError as e:
            # A corrupt cache entry must be discarded, not trusted. Caching a
            # bad download once would otherwise make the error permanent.
            print(f"  [cache] discarding {acc}: {e}", file=sys.stderr)
            f.unlink()
    txt = valid_fasta(
        urllib.request.urlopen(UNIPROT.format(acc), timeout=90).read().decode())
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(txt)
    return txt

def scores(hmm, fasta):
    """Best bit score per target. -E 1000 so that a sequence the model does
    not recognise at all is reported as 'no hit' rather than silently missing
    because the threshold was too strict to say anything."""
    with tempfile.TemporaryDirectory() as d:
        tbl = Path(d) / "t.txt"
        subprocess.run(["hmmsearch", "--max", "-E", "1000", "--tblout", str(tbl),
                        "-o", os.devnull, hmm, fasta], check=True)
        best = {}
        for line in open(tbl):
            if line.startswith("#"):
                continue
            f = line.split()
            tgt, sc = f[0], float(f[5])
            if sc > best.get(tgt, -1e9):
                best[tgt] = sc
    return best

def ga_of(hmm):
    for line in open(hmm):
        if line.startswith("GA"):
            return float(line.split()[1])
        if line.startswith("HMM "):
            break
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="seeds/characterised_facd.tsv")
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--seed", default=None,
                    help="seed FASTA of the model. When given, membership is "
                         "determined from the file rather than from the "
                         "in_seed column of the table, which is hand-entered "
                         "and therefore the weakest link in this test. A "
                         "sequence present in the seed cannot test the model.")
    ap.add_argument("--cache", default="seeds/characterised_cache")
    ap.add_argument("--out", default="results/hmm/characterised_test.json")
    a = ap.parse_args()

    rows = read_table(a.table)
    seed_accs = set()
    if a.seed:
        for line in open(a.seed):
            if line.startswith(">"):
                tok = line[1:].split()[0].split("|")
                seed_accs.add(tok[1] if len(tok) > 2 else tok[0])
        n_declared = sum(1 for r in rows
                         if r.get("in_seed", "no").lower() == "yes")
        n_actual = sum(1 for r in rows if r["accession"] in seed_accs)
        print(f"[seed] {len(seed_accs)} accessions in {a.seed}; "
              f"{n_actual} of the {len(rows)} table entries are in it "
              f"(table declared {n_declared})", file=sys.stderr)
        for r in rows:
            r["in_seed"] = "yes" if r["accession"] in seed_accs else "no"
    ga = ga_of(a.hmm)
    if ga is None:
        sys.exit("[FAIL] no GA line in the HMM header")
    print(f"[test] {len(rows)} characterised entries, GA = {ga:.2f} bits",
          file=sys.stderr)

    # A sequence that could not be downloaded is reported as untested, never as
    # a failure of the model. Without that distinction a network timeout gives
    # a "no hit" line indistinguishable from a genuine rejection, that is a
    # false negative in a validation table.
    fetched, unavailable = {}, []
    for r in rows:
        try:
            fetched[r["accession"]] = fetch(r["accession"], a.cache)
        except Exception as e:
            unavailable.append(r["accession"])
            print(f"  [network] {r['accession']} could not be fetched: {e}",
                  file=sys.stderr)
    if unavailable:
        print(f"  [network] {len(unavailable)} sequence(s) UNTESTED this run; "
              "rerun when the network is back. Successful downloads are cached.",
              file=sys.stderr)
    with tempfile.TemporaryDirectory() as d:
        fa = Path(d) / "characterised.faa"
        with open(fa, "w") as fh:
            for txt in fetched.values():
                fh.write(txt)
        best = scores(a.hmm, str(fa)) if fetched else {}

    # UniProt FASTA headers are "sp|ACC|NAME ..." or "tr|ACC|NAME ...", and
    # hmmsearch reports the first word, so match on the middle field.
    by_acc = {}
    for tgt, sc in best.items():
        parts = tgt.split("|")
        by_acc[parts[1] if len(parts) > 2 else tgt] = sc

    out, n_ok, n_test, n_cand = [], 0, 0, 0
    print(f"\n  {'accession':<10} {'name':<14} {'expected':<9} {'score':>8}  "
          f"{'verdict':<12} note", file=sys.stderr)
    for r in rows:
        acc = r["accession"]
        got = acc in fetched
        sc = by_acc.get(acc)
        called = "positive" if (sc is not None and sc >= ga) else "negative"
        held_out = r.get("in_seed", "no").lower() != "yes"
        exp = r["expected"]

        if not got:
            verdict, ok, note = "UNTESTED", None, "download failed"
        elif exp == "out_of_scope":
            # A real defluorinase of a different fold. An equivalog model is
            # deliberately narrow: detecting one of these would mean the model
            # had left its own family. "not detected" is the correct answer,
            # and it measures specificity of scope, not sensitivity.
            ok = called == "negative"
            verdict = "SCOPE OK" if ok else "OUT OF FAMILY"
            note = "other fold, must not be detected"
            n_test += 1
            n_ok += ok
        elif exp == "candidate":
            # Structure solved, activity not measured on fluoroacetate. The
            # model's answer is a prediction here, not a test of the model.
            verdict, ok = ("predicted +" if called == "positive"
                           else "predicted -"), None
            note = "candidate, activity not measured"
            n_cand += 1
        else:
            ok = called == exp
            verdict = "OK" if ok else "WRONG"
            note = "" if held_out else "seed anchor"
            if held_out:
                n_test += 1
                n_ok += ok
        out.append(dict(accession=acc, name=r["name"],
                        organism=r.get("organism", ""), expected=exp,
                        score=sc, called=called, correct=ok, tested=got,
                        held_out=held_out, evidence=r.get("evidence", "")))
        s = f"{sc:8.1f}" if sc is not None else ("       -" if not got
                                                 else "  no hit")
        print(f"  {acc:<10} {r['name']:<14} {exp:<9} {s}  "
              f"{verdict:<12} {note}", file=sys.stderr)

    print(f"\n  held-out characterised entries: {n_ok}/{n_test} correct",
          file=sys.stderr)
    if n_cand:
        print(f"  {n_cand} candidate(s) reported as PREDICTIONS, not scored: "
              "their activity on fluoroacetate has never been measured.",
              file=sys.stderr)
    if unavailable:
        print(f"  {len(unavailable)} entry(ies) UNTESTED -- network, not model.",
              file=sys.stderr)
    print("  (entries marked 'seed' are anchors of the model and are reported"
          "\n   for completeness only; they cannot test anything.)",
          file=sys.stderr)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(
        dict(hmm=a.hmm, gathering_threshold=ga, n_held_out=n_test,
             n_correct=n_ok, n_untested=len(unavailable),
             entries=out), indent=2))
    print(f"[ok] {a.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
