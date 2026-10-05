#!/usr/bin/env bash
# Does the archaeal FAcD result depend on the archaeal member of the seed?
#
#   bash scripts/seed_composition_control.sh [seeds] [results/hmm] [proteomes]
#
# The 8 FAcD+ archaea found are all haloarchaea, and the model seed itself
# contains an archaeon, Haloarcula, plus one Euryarchaeota entry. The model may
# therefore be finding genuine archaeal FAcD, or only what resembles its single
# archaeal member.
#
# The test: remove the archaea from the seed, rebuild and recalibrate the
# model, and rescan. If the same 8 genomes are still found, detection does not
# depend on that member. If they disappear, the archaeal figure is an artefact
# of seed composition and has to be withdrawn. Same approach as the relaxed
# scan used for the absence of Fluc in DPANN: remove the support and see
# whether the result stands on its own.
set -uo pipefail
S="${1:-seeds}"; R="${2:-results/hmm}"; P="${3:-data/gtdb/proteomes}"
mkdir -p "$R" results/gtdb

echo "=== 1. removing archaea from the seed ==="
python3 - "$S/facd_pos.faa" "$S/facd_pos_no_archaea.faa" <<'PY'
import re, sys
# Archaeal genera and lineages to remove. The pattern is deliberately wider
# than Haloarcula alone: dropping one bacterium too many costs less than
# leaving one archaeon in.
ARCH = re.compile(r"OS=(Halo|Natr|Methano|Sulfolob|Thermococc|Pyrococc|Archaeo|"
                  r"Ferroplasma|Picrophilus|Thermoplasma|Nitroso|Candidatus\s+Nitroso)"
                  r"|_9EURY|_9ARCH|_9CREN|ARCHAEA", re.I)
src, dst = sys.argv[1], sys.argv[2]
kept, removed = [], []
hdr = seq = None
def push():
    if hdr is None: return
    (removed if ARCH.search(hdr) else kept).append((hdr, seq))
for line in open(src):
    if line.startswith(">"):
        push(); hdr, seq = line[1:].rstrip(), ""
    else:
        seq += line
push()
with open(dst, "w") as fh:
    for h, s in kept:
        fh.write(f">{h}\n{s}")
print(f"  original seed : {len(kept)+len(removed)} sequences")
print(f"  removed       : {len(removed)}")
for h, _ in removed:
    print(f"     {h[:100]}")
print(f"  kept          : {len(kept)}  -> {dst}")
PY

echo
echo "=== 2. rebuild and recalibrate WITHOUT the archaea ==="
# The seed has 211 sequences after removal, so exact leave-one-out realignment
# is impractical here (cubic time). This uses the pruning mode, validated
# against realignment on three families by scripts/validate_loo.sh, plus a
# cache so that a rerun costs nothing.
#
# The build output goes to a log file and is only filtered on success.
# Filtering it in the pipeline would swallow a build failure and leave the next
# step to fail on a missing file instead.
LOG="$R/FAcD_no_archaea.build.log"
if python3 scripts/build_hmm.py \
    --positives "$S/facd_pos_no_archaea.faa" --negatives "$S/facd_neg.faa" \
    --confirmed-prefix 'sp|' --name FAcD_no_archaea --outdir "$R" \
    --rule midpoint --loo-mode prune --loo-cache "$R/FAcD_no_archaea.loo.json" \
    --threads 1 > "$LOG" 2>&1; then
  grep -E "separation|gathering_threshold|loo_min|negative_max|threshold_window" "$LOG" \
    | head -12
else
  echo "  BUILD FAILED -- last 15 lines of $LOG:"
  tail -15 "$LOG" | sed 's/^/    /'
  exit 1
fi
[ -f "$R/FAcD_no_archaea.hmm" ] || { echo "  no model produced, stopping."; exit 1; }

echo
echo "=== 3. rescan the archaea with the archaea-free model ==="
python3 scripts/scan_proteomes.py --proteomes "$P" \
  --hmm "$R/FAcD_no_archaea.hmm" \
  --out results/gtdb/scan_facd_no_archaea.tsv --jobs "${JOBS:-6}" --cpu "${CPU:-2}"

echo
echo "=== 4. verdict ==="
python3 - <<'PY'
import csv
a={r['genome_id']:int(r['n_copies']) for r in
   csv.DictReader(open("results/gtdb/scan_facd_archaea.tsv"),delimiter='\t')}
b={r['genome_id']:int(r['n_copies']) for r in
   csv.DictReader(open("results/gtdb/scan_facd_no_archaea.tsv"),delimiter='\t')}
A={g for g,v in a.items() if v>0}; B={g for g,v in b.items() if v>0}
print(f"  with archaea in the seed    : {len(A)} positives")
print(f"  without archaea in the seed : {len(B)} positives")
print(f"  kept: {len(A&B)}   lost: {len(A-B)}   new: {len(B-A)}")
if A-B: print("  lost:", ", ".join(sorted(A-B)))
if B-A: print("  new:", ", ".join(sorted(B-A)))
print()
print("  " + ("the archaeal result does not depend on the archaeal member "
               "of the seed"
               if len(A&B)==len(A) else
               "part of the archaeal result depends on the seed composition"))
PY
