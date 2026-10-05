#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Describe a seed: how many of its sequences carry experimental evidence, and how
many were added by homology. Reports how many entries are manually reviewed
(sp|) rather than automatic (tr|), and how many carry each level of UniProt
protein existence evidence.

A seed built by extending a small characterised core through homology search is
ordinary practice, but only the core is biochemically characterised, and the
Methods section has to state which part that is. This script measures the split
so the Methods can state it rather than imply it.

PE levels, from UniProt:
  1 evidence at protein level   2 evidence at transcript level
  3 inferred from homology      4 predicted            5 uncertain
Only PE 1 means the protein itself was observed. Nothing below PE 1 supports
a claim of biochemical characterisation.
"""
import argparse, collections, json, sys, time, urllib.request
from pathlib import Path

FIELDS = "accession,protein_existence,reviewed,protein_name,organism_name"
URL = ("https://rest.uniprot.org/uniprotkb/search?query=accession:({})"
       "&fields=" + FIELDS + "&format=tsv&size=500")

def accessions(fasta):
    """UniProt FASTA headers read sp|ACC|NAME or tr|ACC|NAME."""
    out = []
    for line in open(fasta):
        if line.startswith(">"):
            tok = line[1:].split()[0]
            p = tok.split("|")
            out.append((tok, p[1] if len(p) > 2 else tok,
                        p[0] if len(p) > 2 else "other"))
    return out

def fetch_pe(accs, cache, batch=100):
    cache = Path(cache)
    known = json.loads(cache.read_text()) if cache.exists() else {}
    todo = [a for a in accs if a not in known]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        q = "%20OR%20".join(chunk)
        try:
            raw = urllib.request.urlopen(URL.format(q), timeout=120).read().decode()
        except Exception as e:
            print(f"  [error] batch {i//batch}: {e}", file=sys.stderr)
            continue
        for line in raw.splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 3:
                known[f[0]] = dict(pe=f[1], reviewed=f[2],
                                   name=f[3] if len(f) > 3 else "",
                                   organism=f[4] if len(f) > 4 else "")
        print(f"  fetched {min(i+batch, len(todo))}/{len(todo)}", file=sys.stderr)
        time.sleep(0.4)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(known, indent=1))
    return known

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", required=True, help="seed FASTA to describe")
    ap.add_argument("--cache", default="seeds/pe_cache.json")
    ap.add_argument("--out", default="results/hmm/seed_audit.json")
    a = ap.parse_args()

    entries = accessions(a.seed)
    accs = [acc for _, acc, _ in entries]
    print(f"[seed] {len(entries)} sequences in {a.seed}", file=sys.stderr)
    src = collections.Counter(s for _, _, s in entries)
    for k, v in sorted(src.items()):
        print(f"   {k:<8} {v}", file=sys.stderr)

    pe = fetch_pe(accs, a.cache)
    by_pe = collections.Counter(pe.get(x, {}).get("pe", "unknown") for x in accs)
    print("\n[protein existence]", file=sys.stderr)
    for k, v in sorted(by_pe.items(), key=lambda t: -t[1]):
        print(f"   {k:<32} {v:4d}  ({100*v/len(accs):.1f} %)", file=sys.stderr)

    n_pe1 = by_pe.get("Evidence at protein level", 0)
    print(f"\n  Characterised core (PE 1): {n_pe1} of {len(accs)} "
          f"({100*n_pe1/len(accs):.1f} %)", file=sys.stderr)
    print("  Every other entry was added by homology rather than by a"
          "\n  published assay, and the Methods state the split.", file=sys.stderr)

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(
        dict(seed=a.seed, n_sequences=len(entries),
             by_source=dict(src), by_protein_existence=dict(by_pe),
             n_evidence_at_protein_level=n_pe1,
             pe1_accessions=[x for x in accs
                             if pe.get(x, {}).get("pe") == "Evidence at protein level"]),
        indent=2))
    print(f"[ok] {a.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
