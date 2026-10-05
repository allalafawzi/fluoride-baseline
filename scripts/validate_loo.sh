#!/usr/bin/env bash
# Validate the fast leave-one-out mode ('prune') against the exact one
# ('realign') on the two seeds small enough to afford both.
#
#   bash scripts/validate_loo.sh
#
# Exact leave-one-out realigns the whole seed once per held-out sequence, which
# costs O(n^3) overall and becomes impractical past about a hundred sequences;
# the FAcD seed has 214. The pruning mode instead removes one row from the
# alignment already computed and drops the columns that become empty. It is
# immediate and deterministic, but it remains an approximation, so it has to be
# shown equivalent on the statistic that sets the threshold: the minimum.
set -u
cd "$(dirname "$0")/.."
S=seeds ; T=results/hmm/loo_test ; mkdir -p "$T"
PATTERN='H\(\+\)/Cl\(-\) exchange transporter|Chloride channel protein ClcB|putative ion channel protein YfeO|Chloride channel protein EriC'

compare () {                 # $1 name  $2 positives  $3 negatives  then args...
  local name=$1 pos=$2 neg=$3; shift 3
  for mode in realign prune; do
    if ! python3 scripts/build_hmm.py --name "${name}_${mode}" \
        --positives "$pos" --negatives "$neg" --outdir "$T" \
        --rule midpoint --loo-mode "$mode" --threads 1 "$@" \
        > "$T/${name}_${mode}.log" 2>&1; then
      echo "  [$name/$mode] FAILED -- last lines:"
      tail -4 "$T/${name}_${mode}.log" | sed 's/^/      /'
      return
    fi
  done
  python3 - "$T" "$name" << 'PY'
import json, sys, pathlib
d, name = pathlib.Path(sys.argv[1]), sys.argv[2]
r = json.loads((d/f"{name}_realign.calib.json").read_text())
e = json.loads((d/f"{name}_prune.calib.json").read_text())
diffs = sorted(((abs(r["loo_scores"][k]-e["loo_scores"][k]), k)
                for k in r["loo_scores"]), reverse=True)
v = [x for x, _ in diffs]
print(f"  {name}")
print(f"    loo_min    realign {r['loo_min']:8.1f}   prune {e['loo_min']:8.1f}"
      f"   diff {e['loo_min']-r['loo_min']:+.1f}")
print(f"    neg_max    realign {r['negative_max']:8.1f}   prune {e['negative_max']:8.1f}")
print(f"    GA         realign {r['gathering_threshold']:8.2f}   "
      f"prune {e['gathering_threshold']:8.2f}"
      f"   diff {e['gathering_threshold']-r['gathering_threshold']:+.2f}")
print(f"    separation {r['separation']} / {e['separation']}")
print(f"    per-sequence difference: median {sorted(v)[len(v)//2]:.1f}  "
      f"max {diffs[0][0]:.1f}  ({diffs[0][1][:44]})")
PY
}

echo "=== Fluc_CrcB (30 sequences) ==="
compare Fluc_CrcB "$S/fluc_pos.faa" "$S/negatives.faa" --confirmed-prefix 'sp|'
echo
echo "=== CLC_F (22 sequences) ==="
compare CLC_F "$S/clcf_pos.faa" "$S/negatives_clc.faa" \
  --confirmed-prefix 'sp|' --confirmed-pattern "$PATTERN"
echo
echo "If the difference on GA is of the order of a bit, 'prune' is usable for"
echo "FAcD (214 sequences), where exact realignment takes a night."
