#!/usr/bin/env bash
# Does the CLC_F seed itself split, and does a stricter model change the census?
#
#   bash scripts/clcf_strict.sh      # PROTEOMES, JOBS and CPU can be set
#
# The published fluoride signatures of CLC-F (GNNLI/GMGLI, GREGT/V, GEVTP;
# Stockbridge & Wackett 2024) are carried in full by only 4 of the 22 seed
# sequences, and not at all by 7 of them. The seed was drawn from TrEMBL
# entries named "Chloride/fluoride channel protein", which is an automatic
# annotation: the negatives were audited against an independent criterion, the
# positives were not.
#
# The test: build a second model from the seed members carrying at least two of
# the three signatures, calibrate it by the same rule, and rescan. If the
# archaeal hits survive, seed composition was not the explanation. If they
# vanish, the archaeal CLC-F figure is an artefact of seed composition, and
# part of the bacterial figure may be as well.
set -uo pipefail
cd "$(dirname "$0")/.."
S=seeds ; O=results/hmm ; P=${PROTEOMES:-data/gtdb/proteomes}

echo "=== 1. splitting the seed by signature content ==="
python3 - "$S/clcf_pos.faa" "$S/clcf_pos_strict.faa" << 'PY'
import re, sys
M = [re.compile(r"G[NM][NG]LI"), re.compile(r"GREG[TV]"), re.compile(r"GEVTP")]
src, dst = sys.argv[1], sys.argv[2]
recs, h, b = [], None, []
for l in open(src):
    if l.startswith('>'):
        if h: recs.append((h, ''.join(b)))
        h, b = l[1:].rstrip(), []
    else:
        b.append(l.strip())
if h: recs.append((h, ''.join(b)))
keep = [(h, s) for h, s in recs if sum(1 for rx in M if rx.search(s)) >= 2]
with open(dst, 'w') as fh:
    for h, s in keep:
        fh.write(f">{h}\n{s}\n")
print(f"  {len(recs)} sequences in, {len(keep)} carry >= 2 of the 3 signatures")
print(f"  -> {dst}")
PY

echo
echo "=== 2. rebuild and recalibrate the strict model ==="
PATTERN='H\(\+\)/Cl\(-\) exchange transporter|Chloride channel protein ClcB|putative ion channel protein YfeO|Chloride channel protein EriC'
LOG="$O/CLC_F_strict.build.log"
if python3 scripts/build_hmm.py --name CLC_F_strict \
    --positives "$S/clcf_pos_strict.faa" --negatives "$S/negatives_clc.faa" \
    --outdir "$O" --rule midpoint --confirmed-prefix 'sp|' \
    --confirmed-pattern "$PATTERN" --threads 1 > "$LOG" 2>&1; then
  grep -E "separation|gathering_threshold|loo_min|negative_max|threshold_window" "$LOG" | head -10
else
  echo "  BUILD FAILED -- last 15 lines of $LOG:"; tail -15 "$LOG" | sed 's/^/    /'; exit 1
fi
[ -f "$O/CLC_F_strict.hmm" ] || { echo "  no model produced."; exit 1; }

echo
echo "=== 3. rescan with the strict model ==="
python3 scripts/scan_proteomes.py --proteomes "$P" --hmm "$O/CLC_F_strict.hmm" \
  --out results/gtdb/scan_clcf_strict.tsv --jobs "${JOBS:-6}" --cpu "${CPU:-2}"

echo
echo "=== 4. verdict ==="
python3 - << 'PY'
import csv
def pos(path, marker):
    out = set()
    for r in csv.DictReader(open(path), delimiter='\t'):
        if r['marker'] == marker and r.get('best_score'):
            out.add(r['genome_id'])
    return out
loose  = pos('results/gtdb/scan_archaea_v3.tsv', 'CLC_F')
strict = pos('results/gtdb/scan_clcf_strict.tsv', 'CLC_F_strict')
print(f"  loose model  : {len(loose)} archaeal genomes")
print(f"  strict model : {len(strict)} archaeal genomes")
print(f"  kept {len(loose & strict)}   lost {len(loose - strict)}   new {len(strict - loose)}")
print()
if not strict:
    print("  the strict model detects nothing in archaea: the archaeal CLC-F")
    print("  signal came from seed members carrying no fluoride signature.")
elif len(strict) < 0.3 * max(1, len(loose)):
    print("  most of the archaeal signal does not survive a seed restricted to")
    print("  signature-carrying members.")
else:
    print("  the archaeal signal survives, so seed composition is not the")
    print("  explanation and the missing signatures need another one.")
PY
