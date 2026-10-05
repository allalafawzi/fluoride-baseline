#!/usr/bin/env bash
# Report the FAcD figure on the balanced sample rather than on every proteome
# that happens to be on disk, and check two things that go with it:
#   - the seed-composition control quoted thresholds from the earlier margin
#     rule (226.6 / 227.0), so it is rerun here under --rule midpoint;
#   - the effect of hmmsearch --max on CLC_F, whose threshold window is only
#     13.5 bits wide.
#
#   bash scripts/facd_balanced.sh    # N can be set (default: 1259)
set -u
cd "$(dirname "$0")/.."
S=results/gtdb/scan_archaea_v3.tsv
P=data/gtdb/ar53_plan.tsv
N=${N:-1259}

echo "=== 1. FAcD on the balanced sample (first $N rows of the plan) ==="
python3 - "$S" "$P" "$N" << 'PY'
import csv, sys, math
scan, plan, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
pos = {m: set() for m in ("FAcD","Fluc_CrcB","CLC_F")}
for r in csv.DictReader(open(scan), delimiter="\t"):
    if r.get("best_score") and r["marker"] in pos:
        pos[r["marker"]].add(r["genome_id"])
rows = list(csv.DictReader(open(plan), delimiter="\t"))
prefix = [r["accession"] for r in rows[:n]]
inpref = set(prefix)
def wilson(k, N, z=1.959963985):
    if N == 0: return (0.0, 0.0)
    p = k/N; d = 1+z*z/N; c = (p+z*z/(2*N))/d
    h = z*math.sqrt(p*(1-p)/N + z*z/(4*N*N))/d
    return (100*(c-h), 100*(c+h))
print(f"  balanced sample : {len(prefix)} genomes")
for m in ("Fluc_CrcB","CLC_F","FAcD"):
    k = len(pos[m] & inpref)
    lo, hi = wilson(k, len(prefix))
    print(f"    {m:<10} {k:4d}/{len(prefix)} = {100*k/len(prefix):6.2f} %"
          f"   CI95 [{lo:.2f} ; {hi:.2f}]")
print(f"\n  FAcD positives INSIDE the balanced sample:")
for g in sorted(pos['FAcD'] & inpref): print(f"    {g}")
outside = sorted(pos['FAcD'] - inpref)
print(f"  FAcD positives OUTSIDE it ({len(outside)}):")
for g in outside: print(f"    {g}")
PY

echo
echo "=== 2. Seed-composition control, rerun under the midpoint rule ==="
if [ -f seeds/facd_pos.faa ]; then
  bash scripts/seed_composition_control.sh
else
  echo "  seeds/facd_pos.faa absent -- skipped"
fi

echo
echo "=== 3. Effect of --max on CLC_F (the 13.5-bit window) ==="
PROT=$(ls data/gtdb/proteomes/*.faa.gz 2>/dev/null | head -200)
if [ -z "$PROT" ]; then echo "  no proteomes found"; exit 0; fi
tmp=$(mktemp -d); zcat $PROT > "$tmp/lot.faa"
for opt in "--max" ""; do
  hmmsearch $opt -T 100 --tblout "$tmp/t$(echo "$opt"|tr -d ' -').txt" \
            -o /dev/null results/hmm/CLC_F.hmm "$tmp/lot.faa"
done
python3 - "$tmp" << 'PY'
import sys, pathlib
d = pathlib.Path(sys.argv[1])
def read(p):
    s = {}
    for l in open(p):
        if l.startswith("#"): continue
        c = l.split(); s[c[0]] = max(s.get(c[0], -1e9), float(c[5]))
    return s
a, b = read(d/"tmax.txt"), read(d/"t.txt")
com = sorted(set(a) & set(b))
print(f"  above 100 bits: {len(a)} with --max, {len(b)} without")
if com:
    ec = sorted(a[k]-b[k] for k in com)
    print(f"  score difference (--max minus default) over {len(com)} in common:")
    print(f"    median {ec[len(ec)//2]:+.2f}   max {ec[-1]:+.2f}   min {ec[0]:+.2f}")
    print(f"  window width is 13.5 bits: a shift of {max(abs(ec[0]),abs(ec[-1])):.2f}"
          f" bits is {100*max(abs(ec[0]),abs(ec[-1]))/13.5:.1f} % of it")
lost = [k for k in a if k not in b]
print(f"  seen with --max and LOST without: {len(lost)}")
for k in lost[:10]: print(f"    {k}  {a[k]:.1f} bits")
PY
rm -rf "$tmp"
