#!/usr/bin/env bash
# Regenerate the census table without the sub-threshold censoring.
#
#   bash scripts/uncensor_scan.sh <proteome_dir> [out.tsv]
#
# With --cut_ga, hmmsearch reports nothing below the gathering threshold, so in
# the resulting table a genome scoring just under the threshold is
# indistinguishable from a genome with no hit at all: both have an empty score
# cell. That covers 80.9 % of the rows of the shipped scan_archaea.tsv
# (Fluc 65 %, CLC_F 78 %, FAcD 100 %).
#
# The threshold-sensitivity analysis cannot be redone from such a table: there
# is no way to see what the census would look like at 27.2 or 61.3 bits instead
# of 40.95. This script produces the uncensored table instead. It needs the
# proteomes, which are too large to deposit, so run it where they are.
#
# Cost: one extra hmmsearch pass per genome, roughly doubling the scan time.
set -euo pipefail
cd "$(dirname "$0")/.."

# Check for hmmsearch once, here. Without this check a missing hmmsearch is
# reported once per genome, 2408 identical error lines instead of one message.
if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "ERROR: hmmsearch not found in PATH." >&2
  echo "  Install HMMER 3.4, then re-run:" >&2
  echo "    conda env create -f envs/environment.lock.yml   # or:" >&2
  echo "    brew install hmmer        (macOS)" >&2
  echo "    apt-get install hmmer     (Debian/Ubuntu)" >&2
  exit 1
fi
HV=$(hmmsearch -h 2>/dev/null | sed -n '2p')
case "$HV" in
  *"HMMER 3.4"*) echo "[uncensor] $HV" ;;
  *) echo "WARNING: this is not HMMER 3.4." >&2
     echo "  $HV" >&2
     echo "  Bit scores depend on the HMMER version through its background model," >&2
     echo "  so the scores produced here may not match the deposited census." >&2
     echo "  Set ALLOW_HMMER_MISMATCH=1 to proceed anyway." >&2
     [ "${ALLOW_HMMER_MISMATCH:-0}" = "1" ] || exit 1 ;;
esac

PROT="${1:-data/gtdb/proteomes}"
OUT="${2:-results/gtdb/scan_archaea_all_scores.tsv}"

if [ ! -d "$PROT" ]; then
  echo "Proteome directory not found: $PROT" >&2
  echo "Pass it as the first argument." >&2
  exit 1
fi
n=$(find "$PROT" -name '*.faa*' | wc -l | tr -d ' ')
echo "[uncensor] $n proteomes in $PROT"
echo "[uncensor] output -> $OUT"

python3 scripts/scan_proteomes.py \
  --hmm results/hmm/Fluc_CrcB.hmm \
  --hmm results/hmm/CLC_F.hmm \
  --hmm results/hmm/FAcD.hmm \
  --proteomes "$PROT" \
  --out "$OUT" \
  --all-scores \
  --jobs "${JOBS:-4}" --cpu "${CPU:-2}"

echo
echo "=== verification ==="
python3 - "$OUT" <<'PY'
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
if "best_score_unthresholded" not in rows[0]:
    print("FAIL: the uncensored column is absent"); sys.exit(1)
if "best_domain_score_unthresholded" not in rows[0]:
    print("FAIL: the domain-score column is absent. --cut_ga requires both the")
    print("      sequence score and a domain score to clear GA, so without it the")
    print("      census count cannot be re-derived from this table.")
    sys.exit(1)
empty_old = sum(1 for r in rows if not r["best_score"].strip())
empty_new = sum(1 for r in rows if not r["best_score_unthresholded"].strip())
print(f"rows                              : {len(rows)}")
print(f"empty in best_score (censored)    : {empty_old} ({100*empty_old/len(rows):.1f} %)")
print(f"empty in best_score_unthresholded : {empty_new} ({100*empty_new/len(rows):.1f} %)")
# where a hit exists the two columns must agree exactly
bad = 0
for r in rows:
    a, b = r["best_score"].strip(), r["best_score_unthresholded"].strip()
    if a and b and abs(float(a) - float(b)) > 0.05:
        bad += 1
print(f"disagreements where both present  : {bad}")
if bad:
    print("FAIL: the two passes disagree on a scored genome"); sys.exit(1)
# Re-derive the census count from the raw scores, the way --cut_ga does it.
import re as _re
ga = {}
for m in ("Fluc_CrcB", "CLC_F", "FAcD"):
    try:
        for line in open(f"results/hmm/{m}.hmm"):
            g = _re.match(r"^GA\s+([-\d.]+)", line)
            if g:
                ga[m] = float(g.group(1)); break
            if line.startswith("HMM "):
                break
    except OSError:
        pass
import collections as _c
cnt, der = _c.Counter(), _c.Counter()
for r in rows:
    m = r["marker"]
    if m not in ga:
        continue
    if (r["best_score"] or "").strip():
        cnt[m] += 1
    sq = (r["best_score_unthresholded"] or "").strip()
    dm = (r["best_domain_score_unthresholded"] or "").strip()
    if sq and dm and float(sq) >= ga[m] and float(dm) >= ga[m]:
        der[m] += 1
bad_m = [m for m in ga if cnt[m] != der[m]]
for m in sorted(ga):
    print(f"  {m:<12} census {cnt[m]:>5}   re-derived {der[m]:>5}   "
          f"{'OK' if cnt[m] == der[m] else 'MISMATCH'}")
if bad_m:
    print("FAIL: the census cannot be re-derived from the raw scores"); sys.exit(1)
print("OK: uncensored columns complete, consistent, and the census re-derives exactly")
PY
