#!/usr/bin/env bash
# Rebuild Fluc_CrcB, CLC_F and FAcD with the current build_hmm.py.
#
#   bash scripts/regenerate_models.sh      # THREADS and LOO_MODE can be set
#
# Two of the three models were first built before the rule that an automatic
# annotation is not a confirmed negative, and their gathering thresholds had
# been written into the .hmm headers by hand, so neither could be reproduced
# from the deposited calibration files. Rebuilding all three from one
# documented command fixes three things at once:
#   - the thresholds become reproducible from the calibration files;
#   - the three models come from the same code;
#   - the noise cutoff NC no longer exceeds the gathering threshold GA. Without
#     the notion of a confirmed negative, NC was max over all negatives, which
#     is how CLC_F came to carry NC 605.2 against GA 160.0.
set -euo pipefail
cd "$(dirname "$0")/.."
S=seeds ; O=results/hmm
mkdir -p "$O/avant"
cp "$O"/Fluc_CrcB.hmm "$O"/Fluc_CrcB.calib.json \
   "$O"/CLC_F.hmm "$O"/CLC_F.calib.json "$O/avant/" 2>/dev/null || true
echo "[backup] previous files in $O/avant/"

echo; echo "=== Fluc_CrcB (all negatives already verified: Swiss-Prot) ==="
python3 scripts/build_hmm.py --name Fluc_CrcB \
  --positives "$S/fluc_pos.faa" --negatives "$S/negatives.faa" \
  --outdir "$O" --rule midpoint --confirmed-prefix 'sp|' --threads "${THREADS:-6}"

echo; echo "=== CLC_F ==="
# A confirmed negative must satisfy two independent conditions: manually
# reviewed (--confirmed-prefix) and carrying the name of an established
# function (--confirmed-pattern). CLC-F needs both:
#   sp|Q4VFY6|SYCA_RHITR  "Symbiosis-assisting ClC homolog", Rhizobium tropici,
#       455 aa, Pfam Voltage_CLC: manually reviewed, but of unknown ionic
#       specificity, scoring 526 bits. Swiss-Prot curation establishes that a
#       protein exists and which family it belongs to, not which ion it moves.
#   tr|U5MTF6|U5MTF6_CLOSA  a precise name, but propagated automatically by
#       UniRule, scoring 524 bits.
# Either criterion alone admits one of the two; together they exclude both.
PATTERN='H\(\+\)/Cl\(-\) exchange transporter|Chloride channel protein ClcB|putative ion channel protein YfeO|Chloride channel protein EriC'
python3 scripts/build_hmm.py --name CLC_F \
  --positives "$S/clcf_pos.faa" --negatives "$S/negatives_clc.faa" \
  --outdir "$O" --rule midpoint --confirmed-prefix 'sp|' \
  --confirmed-pattern "$PATTERN" --threads "${THREADS:-6}"

# FAcD: same code, same rule, same determinism. Its published threshold
# (226.98) came from the earlier margin rule and a multi-threaded MAFFT run, so
# it is rebuilt here with the other two from a single command. With 214 seed
# sequences a full realignment per held-out sequence takes a night, so the
# leave-one-out step uses the pruning approximation, validated against exact
# realignment by scripts/validate_loo.sh.
if [ -f "$S/facd_pos.faa" ] && [ -f "$S/facd_neg.faa" ]; then
  cp "$O"/FAcD.hmm "$O"/FAcD.calib.json "$O/avant/" 2>/dev/null || true
  echo; echo "=== FAcD ==="
  echo "  (check below that the top-scoring Swiss-Prot negatives are"
  echo "   characterised enzymes: named epoxide hydrolases, not 'putative' or"
  echo "   'uncharacterized' entries. If one has no established function, the"
  echo "   Q4VFY6 case is repeating on a third family.)"
  python3 scripts/build_hmm.py --name FAcD \
    --positives "$S/facd_pos.faa" --negatives "$S/facd_neg.faa" \
    --outdir "$O" --rule midpoint --confirmed-prefix 'sp|' \
    --loo-mode "${LOO_MODE:-prune}" --loo-cache "$O/FAcD.loo.json" \
    --threads 1
  MODELS="Fluc_CrcB CLC_F FAcD"
else
  echo; echo "[FAcD] seeds absent from $S/ -- rebuild skipped."
  echo "       Run facd_seed.py first, then this script."
  MODELS="Fluc_CrcB CLC_F"
fi

echo; echo "=== Before / after comparison ==="
for m in $MODELS; do
  echo "-- $m"
  for f in "$O/avant/$m.hmm" "$O/$m.hmm"; do
    [ -f "$f" ] || continue
    printf "   %-28s " "$(basename "$(dirname "$f")")/$(basename "$f")"
    grep -E '^(GA|TC|NC)' "$f" | awk '{printf "%s=%s  ", $1, $2}'; echo
  done
done
echo
echo "If CLC_F still reports SEPARATION IMPOSSIBLE, the ten top-scoring"
echo "unlabelled sequences printed above show where the pattern needs widening."
echo "If GA changed, rerun the census, or justify keeping the previous"
echo "threshold against the score distribution of the census itself."

echo
echo "=== HMMER convention check: NC < GA < TC ==="
for m in $MODELS; do
  awk -v m="$m" '/^NC/{nc=$2} /^GA/{ga=$2} /^TC/{tc=$2} END{
    ok = (nc=="" || (nc+0 < ga+0)) && (ga+0 < tc+0) ? "OK" : "INCOHERENT"
    printf "  %-12s NC=%-8s GA=%-8s TC=%-8s %s\n", m, (nc==""?"-":nc), ga, tc, ok
  }' "results/hmm/$m.hmm"
done
