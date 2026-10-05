#!/usr/bin/env bash
# Regenerate SHA256SUMS.txt on Linux and on macOS.
#
# macOS ships `shasum -a 256` and not `sha256sum`, so the Linux one-liner
#     find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum
# fails on a Mac -- but only after the `>` redirection has already truncated
# SHA256SUMS.txt, leaving the deposit with an empty manifest. This script picks
# whichever tool is present, writes to a temporary file, and replaces the
# manifest only on success, so a failure leaves the existing one untouched.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if command -v sha256sum >/dev/null 2>&1; then
    HASH=(sha256sum)
elif command -v shasum >/dev/null 2>&1; then
    HASH=(shasum -a 256)
else
    echo "ERROR: neither sha256sum nor shasum is available." >&2
    exit 1
fi

TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT

# LC_ALL=C keeps the sort order the same on every machine; otherwise two correct
# manifests differ by line order alone.
#
# The manifest must cover deposited content and not files derived from a run, so
# the exclusions come from scripts/deposit_exclusions.py -- the single definition
# of what belongs to the deposit, shared with the coverage check and the tar
# step. Bytecode is the case that matters: a .pyc compiled from an older source
# would be checksummed as if it were deposited content while contradicting the
# .py beside it.
#
# The predicates are read into an array, one argument per line, and not through
# unquoted substitution: they contain '*', which the shell would expand against
# the real files.
FIND_ARGS=()
while IFS= read -r _a; do FIND_ARGS+=("$_a"); done < <(
    python3 "$(dirname "$0")/deposit_exclusions.py" --find-args0)
find . -type f "${FIND_ARGS[@]}" \
     -print0 \
  | LC_ALL=C sort -z \
  | xargs -0 "${HASH[@]}" > "$TMP"

if [ ! -s "$TMP" ]; then
    echo "ERROR: the computed manifest is empty; SHA256SUMS.txt is unchanged." >&2
    exit 1
fi

mv "$TMP" SHA256SUMS.txt
trap - EXIT
echo "SHA256SUMS.txt regenerated: $(wc -l < SHA256SUMS.txt | tr -d ' ') files"
echo
echo "To verify:"
echo "  ${HASH[*]} -c SHA256SUMS.txt | grep -v ': OK$'"
echo "  (should print nothing)"
