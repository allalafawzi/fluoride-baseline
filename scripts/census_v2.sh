#!/usr/bin/env bash
# Rerun the archaeal census with the rebuilt gathering thresholds and compare
# the previous result with the new one, genome by genome.
#   Fluc_CrcB   27.20 -> 40.95
#   CLC_F      160.00 -> 167.15
#   FAcD       226.98 -> 227.15
# The comparison names the genomes gained and lost rather than counting them:
# a list can be checked against the score distribution, a count cannot.
#
#   bash scripts/census_v2.sh        # PROTEOMES, JOBS and CPU can be set
set -u
cd "$(dirname "$0")/.."
P=${PROTEOMES:-data/gtdb/proteomes}
OLD=results/gtdb/scan_archaea.tsv
NEW=results/gtdb/scan_archaea_v2.tsv

[ -f "$OLD" ] || { echo "not found: $OLD"; exit 1; }
cp "$OLD" results/gtdb/scan_archaea_v1_backup.tsv
echo "[backup] results/gtdb/scan_archaea_v1_backup.tsv"

rm -f "$NEW"
python3 scripts/scan_proteomes.py --proteomes "$P" \
  --hmm results/hmm/Fluc_CrcB.hmm --hmm results/hmm/CLC_F.hmm \
  --hmm results/hmm/FAcD.hmm \
  --out "$NEW" --jobs "${JOBS:-6}" --cpu "${CPU:-2}"

echo
echo "=== Old / new comparison, per marker ==="
python3 - "$OLD" "$NEW" << 'PY'
import csv, sys, collections
def read(p):
    d = collections.defaultdict(dict)
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("best_score"):
                d[r["marker"]][r["genome_id"]] = float(r["best_score"])
    return d
a, b = read(sys.argv[1]), read(sys.argv[2])
for m in sorted(set(a) | set(b)):
    va, vb = set(a.get(m, {})), set(b.get(m, {}))
    lost, gained = va - vb, vb - va
    print(f"\n  {m}")
    if not va:
        print("    marker absent from the old scan -- 'gained' is an artefact,"
              " not a real gain")
    print(f"    old {len(va):5d}   new {len(vb):5d}   "
          f"lost {len(lost):4d}   gained {len(gained):4d}")
    for g in sorted(lost)[:12]:
        print(f"      lost    {g}  old score {a[m][g]:.1f}")
    if len(lost) > 12:
        print(f"      ... and {len(lost)-12} more")
    for g in sorted(gained)[:12]:
        print(f"      gained  {g}  new score {b[m][g]:.1f}")
PY
echo
echo "Expected for Fluc_CrcB: about 1 genome lost out of 832 -- the only one"
echo "between 27.2 and 41.0 bits. A larger difference needs an explanation."
