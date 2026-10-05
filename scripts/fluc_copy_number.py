#!/usr/bin/env python3
"""
Copy number of the Fluc/CrcB channel per genome, and adjacency of the pairs, by
lifestyle.

Fluc channels work as antiparallel dimers, and two architectures exist: a single
gene whose product assembles with itself, or a pair of adjacent genes whose
products assemble with each other. Counting copies per genome, then testing
whether the two copies are neighbours on the contig, separates the two from the
detection data alone.

No coordinate file is consulted, so adjacency can only be decided when the
identifiers themselves carry a position (contig plus gene index). On this panel
that covers 7 of the 46 halophile pairs; the other 39 are protein accessions and
stay undecided. NCBI assigns those accession numbers at submission and not by
gene order on the contig, so consecutive accessions say nothing about position.
The script therefore reports three categories and computes no adjacency
frequency while the undecided pairs dominate.

Output: results/gtdb/fluc_copy_number.json
"""
import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lifestyle_fluc import lifestyle

GROUPS = ("halophile", "thermophile", "thermoacidophile", "unclassified")


def id_format(a):
    """Is the trailing number of this identifier a position?

    Three forms occur in the scan tables:
      CP040679.1_157   -> full contig (with version) + gene index
      PSP77129.1       -> protein accession + version suffix
      WP_008417722.1   -> protein accession, with an underscore in the prefix
    Only the first carries a position. Consecutive protein accessions say
    nothing about where the genes sit: NCBI assigns those numbers at submission,
    not by order on the contig.

    The test therefore requires the full ".<version>_<index>" pattern. Matching
    a trailing "_<number>" alone would read WP_008417722 as a contig index and
    inflate the number of pairs counted as interpretable.
    """
    return "contig_index" if re.search(r"\.\d+_\d+$", a) else "protein_accession"


def adjacent(a, b):
    """True / False / None. None means undecided, for want of coordinates.

    No coordinate file (GFF3, GenBank) is consulted, so the test can only
    conclude when the identifiers themselves carry a position, that is in the
    contig+index form. For a pair of protein accessions the result is None:
    neither adjacent nor non-adjacent, but undecided. Turning those cases into
    an adjacency frequency would build a result out of a naming convention.
    """
    if id_format(a) != "contig_index" or id_format(b) != "contig_index":
        return None
    a = re.sub(r"\.\d+$", "", a)
    b = re.sub(r"\.\d+$", "", b)
    ma = re.match(r"^(.*?)(\d+)$", a)
    mb = re.match(r"^(.*?)(\d+)$", b)
    if not (ma and mb):
        return None
    if ma.group(1) != mb.group(1):
        return False
    return abs(int(ma.group(2)) - int(mb.group(2))) == 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="data/gtdb/ar53_plan.tsv")
    ap.add_argument("--scan", default="results/gtdb/scan_archaea.tsv")
    ap.add_argument("--marker", default="Fluc_CrcB")
    ap.add_argument("--prefix", type=int, default=1259)
    ap.add_argument("--out", default="results/gtdb/fluc_copy_number.json")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.plan), delimiter="\t"))[:a.prefix]
    life = {r["accession"]: lifestyle(r) for r in rows}

    dist = defaultdict(Counter)
    adj = defaultdict(lambda: [0, 0, 0])   # adjacent, not adjacent, undecidable
    ncop = defaultdict(lambda: [0, 0])     # genomes, copies

    for r in csv.DictReader(open(a.scan), delimiter="\t"):
        if r["marker"] != a.marker or r["genome_id"] not in life:
            continue
        k = int(r["n_copies"])
        if k == 0:
            continue
        g = life[r["genome_id"]]
        dist[g][k] += 1
        ncop[g][0] += 1
        ncop[g][1] += k
        if k == 2:
            ids = [x for x in r.get("hit_ids", "").split(";") if x]
            if len(ids) == 2:
                v = adjacent(*ids)
                adj[g][0 if v is True else (1 if v is False else 2)] += 1

    out = {"marker": a.marker, "by_lifestyle": {}}

    print("=" * 72)
    print("Copies per genome carrying the channel")
    print("=" * 72)
    print(f"{'lifestyle':<20}{'genomes':>9}{'copies':>9}"
          f"{'copies/genome':>15}{'distribution':>20}")
    for g in GROUPS:
        if not ncop[g][0]:
            continue
        n, c = ncop[g]
        print(f"{g:<20}{n:>9}{c:>9}{c/n:>15.2f}"
              f"{str(dict(sorted(dist[g].items()))):>20}")
        out["by_lifestyle"][g] = {
            "n_genomes": n, "n_copies": c, "copies_per_genome": c / n,
            "distribution": {str(k): v for k, v in sorted(dist[g].items())}}

    print("\n" + "=" * 72)
    print("Adjacency of the pairs (genomes with exactly 2 copies)")
    print("Only pairs whose identifiers carry a position (contig+index) are")
    print("decidable. A pair of protein accessions stays undecided, and no")
    print("adjacency frequency is derived from it.")
    print("=" * 72)
    print(f"{'lifestyle':<20}{'pairs':>8}{'adjacent':>12}"
          f"{'not adj.':>10}{'undecid.':>10}{'% adj. (decidable)':>18}")
    for g in GROUPS:
        aa, bb, cc = adj[g]
        t = aa + bb + cc
        if not t:
            continue
        dec = aa + bb
        rate = f"{100*aa/dec:.0f}% of {dec}" if dec >= 10 else "not computed"
        print(f"{g:<20}{t:>8}{aa:>12}{bb:>10}{cc:>10}{rate:>18}")
        out["by_lifestyle"].setdefault(g, {}).update(
            {"pairs": t, "pairs_adjacent": aa, "pairs_not_adjacent": bb,
             "pairs_undecidable": cc,
             "frac_adjacent_among_decidable": (aa / (aa + bb)) if (aa + bb) >= 10 else None,
             "note": "adjacency is decidable only for contig+index identifiers"})

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
