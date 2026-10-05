#!/usr/bin/env bash
# Recompute every figure quoted in the Results section from the v2 census.
#
#   bash scripts/rerun_analysis.sh
set -u
cd "$(dirname "$0")/.."
S=results/gtdb/scan_archaea_v2.tsv
P=data/gtdb/ar53_plan.tsv

echo "=== 1. Is the lost genome part of the balanced sample? ==="
python3 - "$P" GCA_003230355.1 << 'PY'
import csv, sys
target = sys.argv[2]
with open(sys.argv[1]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["accession"] == target:
            print(f"  present in the plan -- layer {r.get('layer')}, "
                  f"phylum {r.get('phylum')}, family {r.get('family')}")
            break
    else:
        print("  ABSENT from the plan: the balanced prevalence is unchanged")
PY

echo
echo "=== 2. Archaeal baseline, marker Fluc_CrcB ==="
python3 scripts/analyse_archaea.py --plan "$P" --scan "$S" \
  --marker Fluc_CrcB --prefix 1259 --min-n 10 \
  --out results/gtdb/baseline_archaea_v2.json

echo
echo "=== 3. FAcD: positives and co-occurrence ==="
python3 - "$S" "$P" << 'PY'
import csv, sys, collections
scan, plan = sys.argv[1], sys.argv[2]
best = collections.defaultdict(dict)
with open(scan) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("best_score"):
            best[r["marker"]][r["genome_id"]] = (float(r["best_score"]),
                                                 int(r["n_copies"]))
meta = {}
with open(plan) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        meta[r["accession"]] = r
facd, fluc = best.get("FAcD", {}), best.get("Fluc_CrcB", {})
# The denominator is the set of genomes actually scanned, that is the
# genome_id values present in the scan file, and not its union with the plan.
# The plan lists every genome the sampling design would eventually draw,
# thousands more than have been downloaded, so the union would inflate the
# denominator and deflate every prevalence.
scanned = set()
with open(scan) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        scanned.add(r["genome_id"])
allg = sorted(scanned)
print(f"  genomes scanned          : {len(allg)}")
print(f"  (genomes listed in the plan, scanned or not: {len(meta)})")
print(f"  Fluc+                    : {len(fluc)}  ({100*len(fluc)/len(allg):.1f} %)")
print(f"  FAcD+                    : {len(facd)}  ({100*len(facd)/len(allg):.2f} %)")
print(f"  FAcD+ also carrying Fluc : {sum(1 for g in facd if g in fluc)} / {len(facd)}")
print("\n  The FAcD positives:")
for g, (s, n) in sorted(facd.items(), key=lambda x: -x[1][0]):
    m = meta.get(g, {}); fs = fluc.get(g)
    print(f"    {g}  {s:6.1f}  copies {n}  {m.get('family','?'):<18} "
          f"{float(m.get('genome_size', 0) or 0)/1e6:5.2f} Mb  "
          f"Fluc {'yes x'+str(fs[1]) if fs else 'no'}")
PY
