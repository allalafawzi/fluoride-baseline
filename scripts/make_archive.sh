#!/usr/bin/env bash
# Build the deposit archive.
#
# Packing by hand let derived bytecode into the tarball: files absent from the
# manifest and therefore not covered by any checksum. This script packs from the
# same exclusion list as the manifest (scripts/deposit_exclusions.py), and
# refuses to produce an archive containing a file that SHA256SUMS.txt does not
# cover -- the folder-to-manifest direction, which checking the manifest against
# the folder does not test.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NAME="$(basename "$ROOT")"
OUT="${1:-$ROOT/../${NAME}.tar.gz}"

cd "$ROOT"

# 1. the manifest must be current BEFORE packing
bash scripts/regenerate_checksums.sh

# 2. folder -> manifest coverage, from the same list the manifest uses
python3 - "$ROOT" <<'PY'
import sys, os, pathlib
sys.path.insert(0, str(pathlib.Path(sys.argv[1]) / "scripts"))
from deposit_exclusions import deposited_files   # the manifest's own list
root = pathlib.Path(sys.argv[1])
listed = set()
for line in (root / "SHA256SUMS.txt").read_text().splitlines():
    if line.strip():
        listed.add(line.split(None, 1)[1].strip().removeprefix("./"))
present = deposited_files(root)
uncovered = sorted(present - listed)
orphan = sorted(listed - present)
if uncovered or orphan:
    for f in uncovered:
        print(f"NOT COVERED by the manifest: {f}", file=sys.stderr)
    for f in orphan:
        print(f"ORPHAN (in the manifest, absent from the folder): {f}", file=sys.stderr)
    sys.exit(1)
print(f"coverage: {len(present)} files, {len(listed)} entries, 0 discrepancies")
PY

# 3. pack. The exclusions are read into an ARRAY, one argument per line, and
#    not through unquoted substitution: they contain '*', which the shell would
#    expand against the real files.
TAR_ARGS=()
while IFS= read -r _a; do TAR_ARGS+=("$_a"); done < <(
    python3 scripts/deposit_exclusions.py --tar-args0)
tar "${TAR_ARGS[@]}" \
    -czf "$OUT" -C "$(dirname "$ROOT")" "$NAME"

# 4. check afterwards that no excluded file survived into the archive
if tar -tzf "$OUT" | grep -E '__pycache__|\.pyc$|\.pyo$|\.DS_Store|\.ok$' ; then
    echo "FAILED: an excluded file survived packing" >&2
    rm -f "$OUT"
    exit 1
fi
echo "archive written: $OUT ($(tar -tzf "$OUT" | grep -vc '/$') files)"
