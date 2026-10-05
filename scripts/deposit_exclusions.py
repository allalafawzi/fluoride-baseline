#!/usr/bin/env python3
"""Single definition of which files belong to the deposit.

Three steps have to agree on what counts as part of the deposit:

  1. regenerate_checksums.sh  -> what enters the manifest
  2. the coverage check (make_archive.sh and check_deposit.py)
  3. the tar command in make_archive.sh

They all import this module, so the three answers cannot drift apart. When each
step carried its own list, the manifest excluded the two markers Snakemake
writes (.firth_selftest.ok and .check_deposit.ok) while the coverage check and
tar did not: once the workflow had actually been run and the markers existed,
make_archive.sh exited 1 and reported them as not covered by the manifest.

What is excluded, and why:

  - bytecode (.pyc, .pyo, __pycache__): derived from the sources, and a .pyc
    compiled from an older source gets timestamped and checksummed as deposited
    content while contradicting the .py beside it.
  - Snakemake markers (.ok): records that a step ran locally, with no content.
  - tooling files (.git, editor caches, .DS_Store): not part of the work.
  - SHA256SUMS.txt: it cannot contain its own digest.
"""

# directories we never descend into
EXCLUDED_DIRS = {
    "__pycache__", ".git", ".ipynb_checkpoints", ".pytest_cache", ".mypy_cache",
}

# suffixes of derived files
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}

# exact names (the basename is enough: each of these exists in one place only)
EXCLUDED_NAMES = {
    "SHA256SUMS.txt",
    ".DS_Store",
    ".firth_selftest.ok",   # written by the firth_selftest rule of the Snakefile
    ".check_deposit.ok",    # written by the check_deposit rule of the Snakefile
}

# Exact paths, relative to the root. These are kept apart from EXCLUDED_NAMES
# because a name like "figure.pdf" is too common to exclude everywhere: it is
# this file, in this place, that is derived.
#
# The deposit ships no figure file, but command 2 of the README writes these
# two. Without them here, a reader who follows the README and then reruns the
# self-test, the order the README recommends, sees one check fail on a deposit
# that has done nothing but follow its own documentation. Note that the failure
# is invisible to anyone who only checks the untouched archive.
EXCLUDED_PATHS = {
    "results/gtdb/figure.pdf",   # both written by command 2 of the README
    "results/gtdb/figure.png",
}

# Directories whose contents are produced locally and never deposited.
#
# results/power/ ships empty: it is the default destination of scripts/power.py,
# not a content directory. Running power.py writes six files into it, and
# without the exclusion the self-test then fails on a deposit that has only run
# one of its own scripts. power.py is called by neither the README nor the
# Snakefile, so the exclusion of the README figures above does not cover it.
#
# A directory belongs here only if it is an output destination with no deposited
# content. seeds/ is also a default destination (facd_seed.py) but it ships
# content: regenerating it overwrites files listed in the manifest, which the
# checksums detect, and that is the intended behaviour.
EXCLUDED_OUTPUT_DIRS = {
    "results/power",        # default destination of scripts/power.py
    "data/gtdb/proteomes",  # default destination of scripts/predict_genes.py
}


def is_excluded(relpath):
    """relpath: path relative to the deposit root, '/' separators."""
    # lstrip("./") would strip every leading '.' or '/' character, including the
    # initial dot of ".check_deposit.ok". Only the "./" prefix is removed here.
    rel = relpath.replace("\\", "/")
    if rel.startswith("./"):
        rel = rel[2:]
    if rel in EXCLUDED_PATHS:
        return True
    if any(rel.startswith(d + "/") for d in EXCLUDED_OUTPUT_DIRS):
        return True
    parts = rel.split("/")
    if any(p in EXCLUDED_DIRS for p in parts[:-1]):
        return True
    name = parts[-1]
    if name in EXCLUDED_NAMES:
        return True
    return any(name.endswith(s) for s in EXCLUDED_SUFFIXES)


def deposited_files(root):
    """Set of the relative paths that make up the deposit."""
    import os
    out = set()
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if d not in EXCLUDED_DIRS]
        for fn in fns:
            rel = os.path.relpath(os.path.join(dp, fn), root).replace("\\", "/")
            if not is_excluded(rel):
                out.add(rel)
    return out


if __name__ == "__main__":
    # Output meant to be consumed from a shell script.
    #
    # The '0' variants emit one argument per line, and that is the only usable
    # form: the patterns contain '*', so an unquoted command substitution would
    # have the shell expand them against the real files instead of passing them
    # through to find or tar.
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "--list"
    if mode in ("--find-args", "--find-args0"):
        args = []
        for d in sorted(EXCLUDED_DIRS):
            args += ["!", "-path", f"*/{d}/*"]
        for suf in sorted(EXCLUDED_SUFFIXES):
            args += ["!", "-name", f"*{suf}"]
        for n in sorted(EXCLUDED_NAMES):
            args += ["!", "-name", n]
        for pth in sorted(EXCLUDED_PATHS):
            args += ["!", "-path", f"./{pth}"]
        for d in sorted(EXCLUDED_OUTPUT_DIRS):
            args += ["!", "-path", f"./{d}/*"]
        print("\n".join(args) if mode.endswith("0") else " ".join(args))
    elif mode in ("--tar-args", "--tar-args0"):
        pats = sorted(EXCLUDED_DIRS) + [f"*{x}" for x in sorted(EXCLUDED_SUFFIXES)]
        pats += sorted(n for n in EXCLUDED_NAMES if n != "SHA256SUMS.txt")
        pats += sorted(EXCLUDED_PATHS)
        pats += [f"{d}/*" for d in sorted(EXCLUDED_OUTPUT_DIRS)]
        args = [f"--exclude={x}" for x in pats]
        print("\n".join(args) if mode.endswith("0") else " ".join(args))
    else:
        for n in sorted(EXCLUDED_DIRS | EXCLUDED_SUFFIXES | EXCLUDED_NAMES
                        | EXCLUDED_PATHS | EXCLUDED_OUTPUT_DIRS):
            print(n)
