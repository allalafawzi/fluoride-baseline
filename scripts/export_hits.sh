#!/usr/bin/env bash
# Extract the detected protein sequences and their alignment to the model.
#
# Table 10 of the article is computed from the 239 detected Fluc sequences. The
# proteomes they come from are about 700 GB and cannot be deposited, but the 239
# sequences themselves can, which is enough to redo the table. Run this where the
# proteomes are; it writes the two files the deposit carries:
#
#     results/gtdb/fluc_hits.faa       the 239 sequences
#     results/gtdb/fluc_hits.a2m       their alignment to the profile
#
# USAGE
#   bash scripts/export_hits.sh /path/to/proteomes
set -euo pipefail

PROT="${1:?usage: bash scripts/export_hits.sh /chemin/vers/proteomes}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if ! command -v hmmalign >/dev/null 2>&1; then
    echo "ERROR: hmmalign not found. Install HMMER 3.4 and run again." >&2
    exit 1
fi
if [ ! -d "$PROT" ]; then
    echo "ERROR: $PROT is not a directory." >&2
    exit 1
fi

# Written through temporary files, so a failure cannot replace a complete
# export with a partial one.
OUT_FAA="results/gtdb/fluc_hits.faa"
OUT_ALN="results/gtdb/fluc_hits.a2m"
TMP_FAA="$(mktemp)"
TMP_ALN="$(mktemp)"
trap 'rm -f "$TMP_FAA" "$TMP_ALN"' EXIT

python3 - "$PROT" "$TMP_FAA" <<'PYEOF'
import csv, gzip, os, sys
prot, out = sys.argv[1], sys.argv[2]

plan = list(csv.DictReader(open("data/gtdb/ar53_plan.tsv"), delimiter="\t"))[:1259]
keep = {r["accession"] for r in plan}

wanted = {}
for r in csv.DictReader(open("results/gtdb/scan_archaea.tsv"), delimiter="\t"):
    if (r["marker"] == "Fluc_CrcB" and r["genome_id"] in keep
            and int(r["n_copies"]) > 0 and r.get("hit_ids")):
        wanted.setdefault(r["genome_id"], set()).update(
            x for x in r["hit_ids"].split(";") if x)

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

exp = {f"{g}|{h}" for g, hs in wanted.items() for h in hs}
n, missing = 0, 0
with open(out, "w", newline="\n") as fh:
    for gid in sorted(wanted):
        path = None
        for cand in (f"{gid}.faa.gz", f"{gid}.faa"):
            p = os.path.join(prot, cand)
            if os.path.exists(p):
                path = p
                break
        if path is None:
            missing += 1
            continue
        for name, seq in read_fasta(path):
            if name in wanted[gid]:
                fh.write(f">{gid}|{name}\n{seq}\n")
                n += 1
print(f"{n} sequences written to {out} ({missing} proteomes missing)")

# Fail rather than warn. Pointed at a directory holding one sequence, an earlier
# version wrote a one-sequence FASTA, said so on stderr and still exited 0, so a
# partial export could overwrite a complete one with nothing stopping the chain.
# These checks are blocking.
errs = []
if missing:
    errs.append(f"{missing} proteomes not found under {prot}")
if n != len(exp):
    errs.append(f"{n} sequences written for {len(exp)} expected")
got = set()
dup = set()
for _l in open(out):
    if _l.startswith(">"):
        _k = _l[1:].strip().split()[0]
        if _k in got:
            dup.add(_k)
        got.add(_k)
if dup:
    errs.append(f"{len(dup)} identifiant(s) en double : {sorted(dup)[:3]}")
if got != exp:
    errs.append(f"{len(exp - got)} manquant(s), {len(got - exp)} en trop "
                f"against the census")
if errs:
    print("EXPORT FAILED:", file=sys.stderr)
    for e in errs:
        print(f"  - {e}", file=sys.stderr)
    print(f"  The partial file is left at {out} for inspection; do not "
          f"deposit it.", file=sys.stderr)
    sys.exit(1)
print(f"[verified] {n} sequences, identifiers identical to the census, "
      f"no duplicates")
PYEOF

hmmalign --outformat A2M --amino results/hmm/Fluc_CrcB.hmm "$TMP_FAA" > "$TMP_ALN"
if [ ! -s "$TMP_ALN" ]; then
    echo "FAILED: hmmalign produced no alignment." >&2
    exit 1
fi
# Everything checks out, so publish now and not before.
mv "$TMP_FAA" "$OUT_FAA"
mv "$TMP_ALN" "$OUT_ALN"
trap - EXIT
echo "alignment written to $OUT_ALN"
echo
echo "Add these two files to the deposit, then regenerate the manifest:"
echo "  bash scripts/regenerate_checksums.sh"
echo
echo "Do not run the find/xargs/sha256sum one-liner directly: macOS has no"
echo "sha256sum, and the '>' redirection empties the manifest before the"
echo "command fails. The script picks the available tool and replaces the"
echo "manifest only on success."
