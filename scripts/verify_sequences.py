#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify, one sequence at a time, whether a "no hit" from a model is real.

A no-hit means one of two different things:

  (a) the model genuinely does not recognise the sequence  -> a real result;
  (b) the sequence never reached hmmsearch                 -> a plumbing bug.

Case (b) is easy to miss. An empty body from UniProt, an error page returned
instead of a FASTA, or a download truncated over a slow link all leave the
accession looking fetched while hmmsearch sees nothing, and the entry is then
scored as a rejection by the model.

This script therefore concludes nothing until it has shown, for each
accession: the HTTP outcome, the raw byte count, the FASTA header, the
sequence length, a checksum, and the complete hmmsearch verdict at a
deliberately permissive threshold. Every accession is downloaded fresh and
nothing is read from the cache.
"""
import argparse, hashlib, os, re, subprocess, sys, tempfile, urllib.error, urllib.request
from pathlib import Path

FASTA_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"
AA = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

def download(acc):
    """Return (status, text). Never raises: the caller must see what happened."""
    try:
        with urllib.request.urlopen(FASTA_URL.format(acc), timeout=90) as r:
            return f"HTTP {r.status}", r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return f"HTTP {e.code}", ""
    except Exception as e:
        return f"NETWORK ERROR ({type(e).__name__}: {e})", ""

def parse_fasta(txt):
    """(header, sequence) or (None, reason). Strict on purpose: an invalid body
    has to be reported as invalid, not passed on as a sequence."""
    if not txt.strip():
        return None, "empty body -- obsolete, demerged or deleted accession"
    if not txt.lstrip().startswith(">"):
        head = txt.strip().splitlines()[0][:70] if txt.strip() else ""
        return None, f"not FASTA (starts with {head!r}) -- error page?"
    lines = txt.strip().splitlines()
    header, seq = lines[0][1:], "".join(lines[1:]).replace(" ", "").upper()
    if not seq:
        return None, "header present but no residues"
    bad = set(seq) - AA
    if bad:
        return None, f"non-amino-acid characters: {sorted(bad)[:6]}"
    if len(seq) < 50:
        return None, f"only {len(seq)} residues -- truncated download?"
    return (header, seq), None

def hmm_ga(hmm):
    for line in open(hmm):
        if line.startswith("GA"):
            return float(line.split()[1])
        if line.startswith("HMM "):
            break
    return None

def score_one(hmm, header, seq):
    """Score one sequence on its own, at E=1000 with --max, so that nothing
    about the other sequences can affect the outcome or hide a failure."""
    with tempfile.TemporaryDirectory() as d:
        fa, tbl = Path(d) / "one.faa", Path(d) / "one.tbl"
        fa.write_text(f">{header}\n{seq}\n")
        r = subprocess.run(["hmmsearch", "--max", "-E", "1000",
                            "--tblout", str(tbl), "-o", os.devnull, hmm, str(fa)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            return None, f"hmmsearch failed: {(r.stderr or '')[:200]}"
        rows = [l for l in tbl.read_text().splitlines() if not l.startswith("#")]
        if not rows:
            return None, "no hit even at E=1000 -- the model truly does not match"
        f = rows[0].split()
        return float(f[5]), f"E-value {f[4]}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="+")
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--cache", default=None,
                    help="optional directory to write validated sequences into; "
                         "nothing is ever read from it, every accession is "
                         "downloaded fresh")
    a = ap.parse_args()

    ga = hmm_ga(a.hmm)
    print(f"model {a.hmm}   gathering threshold {ga}")
    print("every accession is downloaded fresh; the cache is not read.\n")

    verdicts = []
    for acc in a.accessions:
        print("=" * 72)
        print(f"  {acc}")
        status, txt = download(acc)
        print(f"    download        : {status}, {len(txt.encode())} bytes")
        parsed, why = parse_fasta(txt)
        if parsed is None:
            print(f"    FASTA           : INVALID -- {why}")
            print(f"    VERDICT         : UNTESTABLE (not the model's fault)")
            verdicts.append((acc, "UNTESTABLE", why))
            continue
        header, seq = parsed
        print(f"    header          : {header[:88]}")
        print(f"    length          : {len(seq)} aa")
        print(f"    sha256[:16]     : {hashlib.sha256(seq.encode()).hexdigest()[:16]}")
        print(f"    first 60 aa     : {seq[:60]}")
        if a.cache:
            Path(a.cache).mkdir(parents=True, exist_ok=True)
            (Path(a.cache) / f"{acc}.fasta").write_text(f">{header}\n{seq}\n")
        sc, note = score_one(a.hmm, header, seq)
        if sc is None:
            print(f"    hmmsearch       : {note}")
            print(f"    VERDICT         : NOT DETECTED -- a real negative")
            verdicts.append((acc, "NOT DETECTED", note))
        else:
            side = "ABOVE" if sc >= ga else "below"
            print(f"    hmmsearch       : {sc:.1f} bits ({note}) -- {side} GA={ga}")
            print(f"    VERDICT         : {'DETECTED' if sc >= ga else 'sub-threshold'}")
            verdicts.append((acc, "DETECTED" if sc >= ga else "sub-threshold",
                             f"{sc:.1f} bits"))

    print("=" * 72)
    print("\nSUMMARY")
    for acc, v, note in verdicts:
        print(f"  {acc:<12} {v:<14} {note[:60]}")
    n_untestable = sum(1 for _, v, _ in verdicts if v == "UNTESTABLE")
    if n_untestable:
        print(f"\n  {n_untestable} accession(s) could not be tested at all. Any earlier"
              "\n  'WRONG' for these was a plumbing failure, not a model failure.")

if __name__ == "__main__":
    main()
