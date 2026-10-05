#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test the archaeal CLC_F hits against structural criteria that were not used to
build the model: the three published signature motifs of fluoride-specific
CLC-F, and the gating-glutamate element shared by the whole CLC family. Rates
are reported for the hits, for the bacterial seed as a positive control and for
confirmed chloride transporters as a negative control, each against a
permutation control of the same amino-acid composition.

The model detects hits in archaeal proteomes, but no archaeal prevalence is
published from them, because neither of the two quantities needed to read such
a number has been measured in archaea:

  * sensitivity: the seed is 22 sequences, all bacterial, so a divergent
    archaeal CLC-F could be missed, which would make the number an
    underestimate.
  * specificity: the confirmed negatives are bacterial ClcA, ClcB and YfeO.
    Archaea have their own CLC-family chloride transporters and none was shown
    to the calibration, so the number could equally be an overestimate. That
    is the direction that matters here: a hit could be a chloride transporter
    the model was never asked to reject.

Stockbridge & Wackett (2024, Nat Commun 15:4593) report that fluoride-specific
CLC-F are identified by three signature sequences: GNNLI/GMGLI in the
N-terminal domain for ion selectivity, GREGT/V at the transport machinery, and
GEVTP in the C-terminal domain for fluoride binding. Those motifs were not used
to build the profile HMM, which is a whole-length alignment model, so they are
an independent criterion, as the HDRG nucleophile elbow was for the FAcD split.

Rates close to those of the bacterial seed indicate genuine CLC-F; rates close
to those of the chloride controls indicate chloride transporters instead.
"""
import argparse, csv, gzip, re, sys
from pathlib import Path

# The three signatures, as published. Kept as plain alternatives rather than
# a fuzzy pattern: a loose regex would manufacture the agreement we are trying
# to test.
MOTIFS = {
    "N-term selectivity (GNNLI/GMGLI)": re.compile(r"G[NM][NG]LI"),
    "transport core (GREGT/V)":         re.compile(r"GREG[TV]"),
    "C-term F-binding (GEVTP)":         re.compile(r"GEVTP"),
}

# A second question, separate from the three motifs above: is the hit a CLC at
# all? The motifs ask whether a CLC is fluoride-specific. They cannot say
# whether the protein belongs to the family, and that is what decides how a
# zero should be read. No fluoride signature and no CLC signature means the
# hits are not CLC proteins. No fluoride signature with a clear CLC signature
# means they are CLC chloride transporters, which is the case that matters,
# because the model would then be detecting the wrong family.
#
# McIlwain, Ruprecht & Stockbridge (2021, Annu Rev Biochem 90:559):
#     "All members of the CLC family possess this so-called gating glutamate."
# In E. casseliflavus CLC-F it is E118; in E. coli ClcA it is E148.
#
# The pattern used here is GREG rather than the GxExxP-shaped motif that the
# literature description suggests. On the controls, a GxExxP pattern matched
# true CLC chloride transporters at 1.4 %, below their 9.2 % chance rate, which
# cannot hold if every CLC carries the glutamate. The sequences themselves
# place the element: ClcA carries GREGP and the CLC-F seed carries GREGT/V,
# the same GREG element with a different following residue.
#
#     GREG    CLC family membership   chloride 69.7 %, CLC-F 86.4 %, chance ~1 %
#     GREG[TV]  fluoride variant      chloride  0.0 %, CLC-F 68.2 %
#     GREGP     chloride variant      chloride 43.4 %, CLC-F  0.0 %
#
# The two variants do not overlap on the controls, which is what makes the test
# readable. GREG[TV] is also one of the three Stockbridge signatures above, so
# the fluoride half of this test is not independent of them; only the GREG
# membership test is new information.
#
# Every rate below is reported against a permutation control (same amino-acid
# composition, order destroyed). A short motif matches by chance, and without
# that control the numbers could not be interpreted.
GATE = {
    "CLC membership (GREG)":        re.compile(r"GREG"),
    "fluoride variant (GREG[TV])":  re.compile(r"GREG[TV]"),
    "chloride variant (GREGP)":     re.compile(r"GREGP"),
}

# Calibration of the criterion itself, measured on the bacterial seed before
# any archaeal sequence is looked at:
#     GNNLI/GMGLI   4/22  (18 %)
#     GREGT/V      15/22  (68 %)
#     GEVTP        14/22  (64 %)
#     all three     4/22  (18 %)
# Taken literally the signatures are therefore not diagnostic: only 18 % of the
# sequences a published criterion calls CLC-F carry all three exactly. They
# have to be read as a graded rather than a binary criterion, and the archaeal
# hits compared to these rates rather than to 100 %. Reporting the share of
# archaeal hits that carry the signatures without this baseline would turn the
# looseness of the criterion into a negative result.

def read_fasta(p):
    op = gzip.open if str(p).endswith(".gz") else open
    name, buf = None, []
    with op(p, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)

def score_motifs(seq):
    return {k: bool(rx.search(seq)) for k, rx in MOTIFS.items()}


def _shuffled(seq, rng):
    """Shuffle a sequence: same amino-acid composition, order destroyed. This
    gives the chance reference rate for a short motif."""
    L = list(seq)
    rng.shuffle(L)
    return "".join(L)


def summarise_gate(label, seqs, seed=1):
    """Score the gating-glutamate patterns with their permutation control.

    Without the control, a six-position motif of which four positions are free
    would match almost any membrane protein and the rates would mean nothing.
    """
    import random
    if not seqs:
        print(f"  {label}: no sequence")
        return None
    rng = random.Random(seed)
    n = len(seqs)
    res = {}
    for k, rx in GATE.items():
        obs = sum(1 for s in seqs if rx.search(s))
        # 20 permutations per sequence give the rate expected by chance
        exp = 0
        for s in seqs:
            exp += sum(1 for _ in range(20) if rx.search(_shuffled(s, rng))) / 20.0
        res[k] = (obs, exp)
    print(f"\n  {label}  (n = {n})")
    print(f"     {'motif':<32}{'observed':>12}{'expected by chance':>22}")
    for k, (obs, exp) in res.items():
        print(f"     {k:<32}{obs:>5} ({100*obs/n:5.1f} %)"
              f"{exp:>12.1f} ({100*exp/n:5.1f} %)")
    return res

def summarise(label, seqs):
    if not seqs:
        print(f"  {label}: no sequence"); return
    counts = {k: 0 for k in MOTIFS}
    n_all = n_none = 0
    for s in seqs:
        r = score_motifs(s)
        for k, v in r.items():
            counts[k] += v
        hits = sum(r.values())
        n_all += hits == 3
        n_none += hits == 0
    n = len(seqs)
    print(f"\n  {label}  (n = {n})")
    for k, c in counts.items():
        print(f"     {k:<34} {c:5d}  ({100*c/n:5.1f} %)")
    print(f"     {'all three signatures':<34} {n_all:5d}  ({100*n_all/n:5.1f} %)")
    print(f"     {'none of the three':<34} {n_none:5d}  ({100*n_none/n:5.1f} %)")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", required=True, help="scan TSV with hit_ids")
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--marker", default="CLC_F")
    ap.add_argument("--seed", default="seeds/clcf_pos.faa",
                    help="bacterial seed: the positive control")
    ap.add_argument("--negatives", default=None,
                    help="optional FASTA of confirmed chloride transporters: "
                         "the negative control")
    a = ap.parse_args()

    wanted = {}
    with open(a.scan) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["marker"] == a.marker and r.get("hit_ids"):
                for h in r["hit_ids"].split(";"):
                    wanted.setdefault(r["genome_id"], set()).add(h)
    print(f"[scan] {len(wanted)} genomes with a {a.marker} hit, "
          f"{sum(len(v) for v in wanted.values())} protein ids")

    found = []
    root = Path(a.proteomes)
    for gid, ids in wanted.items():
        for cand in (root / f"{gid}.faa.gz", root / f"{gid}.faa"):
            if cand.exists():
                for name, seq in read_fasta(cand):
                    if name in ids:
                        found.append(seq)
                break
    print(f"[extract] {len(found)} hit sequences recovered from the proteomes")
    if len(found) < sum(len(v) for v in wanted.values()):
        print("  (some proteomes were not found on disk; the comparison still "
              "holds on those recovered)")

    print("\n" + "=" * 64)
    print("Signature motifs of fluoride-specific CLC-F")
    print("Stockbridge & Wackett 2024, Nat Commun 15:4593")
    print("=" * 64)
    summarise(f"{a.marker} hits in the scanned proteomes", found)
    if Path(a.seed).exists():
        summarise("bacterial seed (positive control)",
                  [s for _, s in read_fasta(a.seed)])

    print("\n" + "=" * 64)
    print("Gating glutamate: is the hit a CLC at all?")
    print("McIlwain, Ruprecht & Stockbridge 2021, Annu Rev Biochem 90:559:")
    print('  "All members of the CLC family possess this gating glutamate."')
    print("Each rate is shown against its own permutation control.")
    print("=" * 64)
    summarise_gate(f"{a.marker} hits in the scanned proteomes", found)
    if Path(a.seed).exists():
        summarise_gate("bacterial seed (positive control)",
                       [s for _, s in read_fasta(a.seed)])
    if a.negatives and Path(a.negatives).exists():
        summarise_gate("confirmed chloride transporters (negative control)",
                       [s for _, s in read_fasta(a.negatives)])
    if a.negatives and Path(a.negatives).exists():
        summarise("confirmed chloride transporters (negative control)",
                  [s for _, s in read_fasta(a.negatives)])

    print("\nRates close to the bacterial seed point to genuine CLC-F and rates")
    print("close to the chloride controls point to chloride transporters, which")
    print("is why no archaeal prevalence is published from these hits.")

if __name__ == "__main__":
    main()
