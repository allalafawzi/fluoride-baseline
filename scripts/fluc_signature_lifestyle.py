#!/usr/bin/env python3
"""
Compares the sequence signatures of the Fluc/CrcB channel between archaeal
lifestyles: halophiles, thermophiles, thermoacidophiles, and the rest.

The prevalence analysis (scripts/lifestyle_fluc.py) asks how many genomes carry
the channel in each lifestyle, a question of presence and absence. This script
asks a different one: among the genomes that carry it, is the protein the same?
A halophile lives in a molar brine and a thermoacidophile at 80 degrees and
pH 2; if those constraints act on the channel itself, they should leave a trace
in the sequence, and the first place to look is the surface of the protein
rather than its active site.

Method:
1. The diagnostic positions are not chosen by hand. They are read from the model
   itself: each match state of the profile HMM carries an amino-acid
   distribution whose entropy measures directly how constrained the position is.
2. The detected sequences are aligned to the model with `hmmalign`, which places
   each residue against the model state it corresponds to. Without that
   alignment, comparing "position 88" between two proteins of different lengths
   means nothing.
3. Three quantities are reported per lifestyle: the residue composition at the
   most constrained positions, the copy number per genome, and the charge
   composition.

The charge composition has to be read together with the hydrophobicity. The
acidic enrichment of obligate halophiles is an adaptation of the
solvent-exposed surface, where carboxylates hold a hydration layer in molar
brine. Fluc is an integral membrane protein of about 126 residues, almost all of
it buried in the lipid bilayer, so it has very little solvent-exposed surface
and very little room for that adaptation. Measured on the seed, which is
entirely bacterial and so independent of any halophile question: GRAVY
(Kyte-Doolittle hydrophobicity index) = +0.96, a composition dominated by L
(16.9 %), G (13.7 %), A (10.4 %), F (8.0 %), V (7.1 %), and only 2.9 % acidic
residues, which is the profile of a transmembrane helix. A whole bacterial
proteome runs near 12 % acidic residues and an obligate halophile proteome near
15-20 %, so at 2.9 % in bacteria the floor is structural and leaves no margin.
An absence of acidic enrichment in halophilic Fluc is therefore the expected
consequence of Fluc being transmembrane, not a negative biological result. The
script reports hydrophobicity next to charge for that reason, and
--proteome-background compares Fluc with the rest of the proteome of the same
genome, which holds the organism fixed.

Lifestyle is confounded with clade (the archaeal halophiles are the
Halobacteria). Any difference observed can be an adaptation or an inheritance,
and this script does not separate the two. It describes; it does not explain.

Usage
  python3 scripts/fluc_signature_lifestyle.py --proteomes /path/to/proteomes
"""
import argparse
import csv
import gzip
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lifestyle_fluc import lifestyle  # the same assignment as the prevalence analysis

AA = "ACDEFGHIKLMNPQRSTVWY"
ACIDIC = set("DE")
BASIC = set("KR")

# Kyte-Doolittle hydrophobicity scale. The GRAVY index (the mean of these
# values over the sequence) measures how hydrophobic a protein is: positive for
# membrane proteins, negative for soluble ones. It is reported here because it
# bounds the interpretation of the charge composition.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}


def gravy(s):
    v = [KD[c] for c in s if c in KD]
    return sum(v) / len(v) if v else 0.0


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def hmm_match_entropy(hmm_path):
    """Entropy of each match state of the profile, read from the HMM file.

    HMMER stores scores as -log(probability) in nats. Low entropy means a
    constrained position: that is the operational definition of a diagnostic
    position used here, so the positions come from the model and not from a
    hand-made choice.
    """
    ent, pos = [], 0
    with open(hmm_path) as fh:
        started = False
        for line in fh:
            if line.startswith("HMM "):
                started = True
                continue
            if not started:
                continue
            f = line.split()
            # match-state line: index then 20 scores
            if len(f) >= 21 and f[0].isdigit():
                try:
                    p = [math.exp(-float(x)) if x != "*" else 0.0
                         for x in f[1:21]]
                except ValueError:
                    continue
                tot = sum(p) or 1.0
                p = [x / tot for x in p]
                h = -sum(x * math.log(x) for x in p if x > 0)
                pos = int(f[0])
                ent.append((pos, h, AA[p.index(max(p))]))
    return ent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--plan", default="data/gtdb/ar53_plan.tsv")
    ap.add_argument("--scan", default="results/gtdb/scan_archaea.tsv")
    ap.add_argument("--hmm", default="results/hmm/Fluc_CrcB.hmm")
    ap.add_argument("--seed", default="seeds/fluc_pos.faa")
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--top", type=int, default=12,
                    help="number of the most constrained positions to compare")
    ap.add_argument("--proteome-background", action="store_true",
                    help="compare Fluc with the rest of the proteome of the "
                         "same genome, which holds the organism fixed and is "
                         "what makes the halophile control testable")
    ap.add_argument("--out", default="results/gtdb/fluc_signature_lifestyle.json")
    a = ap.parse_args()

    for exe in ("hmmalign",):
        if subprocess.call(["which", exe], stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL) != 0:
            sys.exit(f"ERROR : {exe} not found. Install HMMER 3.4 "
                     f"(conda install -c bioconda hmmer) and rerun.")

    rows = list(csv.DictReader(open(a.plan), delimiter="\t"))[:a.prefix]
    life = {r["accession"]: lifestyle(r) for r in rows}

    wanted = defaultdict(set)
    for r in csv.DictReader(open(a.scan), delimiter="\t"):
        if (r["marker"] == a.marker and r.get("hit_ids")
                and r["genome_id"] in life and int(r["n_copies"]) > 0):
            for h in r["hit_ids"].split(";"):
                if h:
                    wanted[r["genome_id"]].add(h)
    print(f"[plan] {len(wanted)} carrying genomes in the balanced sample")

    seqs = {}          # id -> (sequence, lifestyle)
    missing = 0
    for gid, ids in wanted.items():
        hit = None
        for cand in (f"{gid}.faa.gz", f"{gid}.faa"):
            p = os.path.join(a.proteomes, cand)
            if os.path.exists(p):
                hit = p
                break
        if hit is None:
            missing += 1
            continue
        for name, s in read_fasta(hit):
            if name in ids:
                seqs[f"{gid}|{name}"] = (s, life[gid])
    print(f"[extraction] {len(seqs)} sequences retrieved, "
          f"{missing} proteomes missing from disk")
    if len(seqs) < 30:
        sys.exit("Too few sequences to compare : check --proteomes")

    # --- the most constrained positions, read from the model -----------------
    ent = hmm_match_entropy(a.hmm)
    if not ent:
        sys.exit("ERROR : the HMM match states cannot be read")
    top = sorted(ent, key=lambda t: t[1])[:a.top]
    top_pos = sorted(p for p, _, _ in top)
    cons = {p: c for p, _, c in ent}
    print(f"[model] {len(ent)} match states ; the {a.top} most "
          f"constrained : {top_pos}")

    # --- alignment to the model ---------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "hits.faa")
        with open(fa, "w") as fh:
            for k, (s, _) in seqs.items():
                fh.write(f">{k}\n{s}\n")
        out = subprocess.run(
            ["hmmalign", "--outformat", "A2M", "--amino", a.hmm, fa],
            capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit("hmmalign failed :\n" + out.stderr[:500])
        aln = {}
        name, buf = None, []
        for line in out.stdout.splitlines():
            if line.startswith(">"):
                if name:
                    aln[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if name:
            aln[name] = "".join(buf)

    # In A2M the match columns are upper case or '-', and the insertions are
    # lower case or '.'. Keeping only the match columns gives back exactly the
    # states of the model.
    def match_cols(s):
        return [c for c in s if c.isupper() or c == "-"]

    bygrp = defaultdict(list)
    for k, s in aln.items():
        bygrp[seqs[k][1]].append(match_cols(s))

    groups = ["halophile", "thermophile", "thermoacidophile", "unclassified"]
    result = {"marker": a.marker, "n_sequences": len(seqs),
              "top_positions": top_pos,
              "consensus": {str(p): cons.get(p) for p in top_pos},
              "by_lifestyle": {}}

    print("\n" + "=" * 74)
    print("Composition at the most constrained positions of the model")
    print("=" * 74)
    hdr = f"{'position':>9}{'consensus':>11}"
    for g in groups:
        hdr += f"{g[:12]:>14}"
    print(hdr)
    for p in top_pos:
        line = f"{p:>9}{cons.get(p, '?'):>11}"
        for g in groups:
            col = [s[p - 1] for s in bygrp.get(g, []) if len(s) >= p]
            if not col:
                line += f"{'-':>14}"
                continue
            c = Counter(col).most_common(1)[0]
            line += f"{c[0]} {100*c[1]/len(col):4.0f}%{'':>6}"
        print(line)

    # --- copy number ---------------------------------------------------------
    ncop = defaultdict(list)
    percontig = defaultdict(int)
    for k in seqs:
        percontig[k.split("|")[0]] += 1
    for gid, n in percontig.items():
        ncop[life[gid]].append(n)

    print("\n" + "=" * 74)
    print("Copies of the channel per genome carrying it")
    print("=" * 74)
    print(f"{'lifestyle':<20}{'genomes':>9}{'copies':>9}{'copies/genome':>16}")
    for g in groups:
        v = ncop.get(g, [])
        if not v:
            continue
        print(f"{g:<20}{len(v):>9}{sum(v):>9}{sum(v)/len(v):>16.2f}")
        result.setdefault("by_lifestyle", {}).setdefault(g, {}).update(
            {"n_genomes": len(v), "n_copies": sum(v),
             "copies_per_genome": sum(v) / len(v)})

    # --- charge and hydrophobicity, read together ----------------------------
    print("\n" + "=" * 74)
    print("Charge composition, with the hydrophobicity that bounds it")
    print("Fluc is an integral membrane protein. The acidic enrichment of")
    print("halophiles is an adaptation of the solvent-exposed surface, which")
    print("this protein has almost none of, and a positive GRAVY is the sign")
    print("of that structural constraint.")
    print("=" * 74)
    print(f"{'lifestyle':<20}{'n':>5}{'% acidic':>11}{'% basic':>13}"
          f"{'net charge':>14}{'GRAVY':>9}")
    for g in groups:
        ss = [seqs[k][0] for k in seqs if seqs[k][1] == g]
        if not ss:
            continue
        ac = sum(sum(1 for c in s_ if c in ACIDIC) / len(s_) for s_ in ss) / len(ss)
        ba = sum(sum(1 for c in s_ if c in BASIC) / len(s_) for s_ in ss) / len(ss)
        net = sum((sum(1 for c in s_ if c in BASIC)
                   - sum(1 for c in s_ if c in ACIDIC)) for s_ in ss) / len(ss)
        gv = sum(gravy(s_) for s_ in ss) / len(ss)
        print(f"{g:<20}{len(ss):>5}{100*ac:>10.1f}%{100*ba:>12.1f}%"
              f"{net:>14.1f}{gv:>9.2f}")
        result.setdefault("by_lifestyle", {}).setdefault(g, {}).update(
            {"n": len(ss), "frac_acidic": ac, "frac_basic": ba,
             "net_charge": net, "gravy": gv})

    # --- Fluc against the proteome of its own genome --------------------------
    if a.proteome_background:
        print("\n" + "=" * 74)
        print("Fluc compared with the rest of the proteome of the same genome")
        print("If the halophilic acidic enrichment is present in these organisms")
        print("but not in Fluc, it shows in the 'proteome' column and not in the")
        print("'Fluc' column. That is what separates a structural constraint")
        print("from an absence of adaptation.")
        print("=" * 74)
        bg = defaultdict(list)
        for gid in percontig:
            hit = None
            for cand in (f"{gid}.faa.gz", f"{gid}.faa"):
                pth = os.path.join(a.proteomes, cand)
                if os.path.exists(pth):
                    hit = pth
                    break
            if hit is None:
                continue
            tot_ac = tot_len = 0
            for _n, s_ in read_fasta(hit):
                tot_ac += sum(1 for c in s_ if c in ACIDIC)
                tot_len += len(s_)
            if tot_len:
                bg[life[gid]].append(tot_ac / tot_len)
        print(f"{'lifestyle':<20}{'genomes':>9}"
              f"{'% acidic proteome':>20}{'% acidic Fluc':>16}")
        for g in groups:
            v = bg.get(g, [])
            if not v:
                continue
            ss = [seqs[k][0] for k in seqs if seqs[k][1] == g]
            fl = sum(sum(1 for c in s_ if c in ACIDIC) / len(s_)
                     for s_ in ss) / len(ss)
            print(f"{g:<20}{len(v):>9}{100*sum(v)/len(v):>19.1f}%"
                  f"{100*fl:>15.1f}%")
            result["by_lifestyle"].setdefault(g, {})["proteome_frac_acidic"] = \
                sum(v) / len(v)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="\n") as fh:
        json.dump(result, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\n-> {a.out}")
    print("\nLifestyle is confounded with clade. Any difference observed "
          "can be an\nadaptation or an inheritance ; this script does not "
          "separate the two.")


if __name__ == "__main__":
    main()
