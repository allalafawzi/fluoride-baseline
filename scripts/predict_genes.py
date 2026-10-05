#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Predict proteins from nucleotide assemblies with Prodigal.

  python3 scripts/predict_genes.py --genomes data/gtdb/genomes \
      --out data/gtdb/proteomes --manifest data/gtdb/log/proteomes_predicted.tsv

Two reasons to predict rather than download.

Coverage: 90 % of the archaeal species representatives in the balanced plan
are GenBank assemblies (GCA_ accessions) and many carry no protein annotation
at NCBI. The phyla affected throughout are Nanobdellota, Micrarchaeota,
Aenigmatarchaeota, Asgardarchaeota, Iainarchaeota, Hydrothermarchaeota,
Altiarchaeota and Undinarchaeota, which are the uncultivated lineages of
interest. Using only the downloaded proteomes would re-select the cultured
lineages, the same bias as a completeness filter by another route.

Comparability: where an annotation does exist it does not come from one
source. RefSeq entries go through PGAP, GenBank entries carry whatever the
depositor supplied. If gene calling differs between clades then the presence
of a gene differs for a reason that is not biological. Calling every gene with
the same tool makes the census comparable across genomes. The origin of each
proteome, downloaded or predicted, is recorded in the manifest and is meant to
be used as a covariate.

Prodigal has two modes. `single` trains its model on the genome itself and
needs a sequence that is long and complete enough; `meta` uses pre-computed
models and works on fragmented assemblies. MAGs are often fragmented, so the
mode switches to `meta` below --min-single (default: 100 kb of total
sequence), and the mode used is recorded for each genome.
"""
import argparse, csv, gzip, os, sys, tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# Two ways of reaching Prodigal, tried in this order:
#   1. pyrodigal, the Python library: nothing to compile;
#   2. the `prodigal` binary on the PATH, which helps when the Python
#      installation is locked (an "externally managed" Debian or Ubuntu
#      environment).
# The engine is the same in both cases, so the result is the same.
import shutil as _shutil
import subprocess as _sp

try:
    import pyrodigal
    METHOD = "pyrodigal"
except ImportError:
    pyrodigal = None
    METHOD = "binary" if _shutil.which("prodigal") else None

MISSING_ENGINE = (
    "[FAIL] neither pyrodigal nor the prodigal binary is available.\n"
    "  Install one of the following:\n"
    "    conda install -c conda-forge -c bioconda pyrodigal\n"
    "    python3 -m pip install pyrodigal\n"
    "    conda install -c bioconda prodigal        (installs the binary)\n"
    "  Then check with:  python3 -c \"import pyrodigal\"  or  which prodigal")


class _SerialExecutor:
    """Stand-in for ProcessPoolExecutor that runs in this process.

    Constructing a pool queries SC_SEM_NSEMS_MAX, which is refused in
    containers and on HPC login nodes without POSIX semaphores, so gene
    prediction would fail there before the first assembly. At --jobs 1 no pool
    is built, and a pool that cannot be constructed falls back here. A serial
    loop changes no prediction: the same work runs one assembly after another
    rather than side by side.
    """

    def __init__(self, initializer=None, initargs=()):
        if initializer is not None:
            initializer(*initargs)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, fn, it, chunksize=None):
        return map(fn, it)


def make_executor(jobs, initializer=None, initargs=()):
    if jobs <= 1:
        return _SerialExecutor(initializer=initializer, initargs=initargs)
    try:
        return ProcessPoolExecutor(max_workers=jobs, initializer=initializer,
                                   initargs=initargs)
    except (OSError, NotImplementedError, ValueError) as e:
        print(f"[predict] process pool unavailable ({e}); falling back to serial.",
              file=sys.stderr)
        return _SerialExecutor(initializer=initializer, initargs=initargs)


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


_CFG = {}


def _init(outdir, min_single):
    _CFG["outdir"], _CFG["min_single"] = outdir, min_single


def _prodigal_binary(fna_path, dest, meta):
    """Call the prodigal binary. Output identical to the pyrodigal route."""
    fd, tmp_faa = tempfile.mkstemp(suffix=".faa"); os.close(fd)
    try:
        cmd = ["prodigal", "-i", str(fna_path), "-a", tmp_faa,
               "-p", "meta" if meta else "single", "-q", "-o", os.devnull]
        r = _sp.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 and not meta:
            # Sequence too short to train on: fall back to meta mode.
            cmd[cmd.index("single")] = "meta"
            meta = True
            r = _sp.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return 0, meta, (r.stderr or "prodigal failed")[:300]
        n = 0
        tmp_gz = str(dest) + ".part"
        with open(tmp_faa) as fi, gzip.open(tmp_gz, "wt") as fo:
            for line in fi:
                if line.startswith(">"):
                    n += 1
                    fo.write(line.split(" # ")[0].rstrip() + "\n")
                else:
                    fo.write(line.replace("*", ""))
        os.replace(tmp_gz, dest)
        return n, meta, None
    finally:
        if os.path.exists(tmp_faa):
            os.unlink(tmp_faa)


def _one_genome(fna_str):
    fna = Path(fna_str)
    gid = fna.name.split(".fna")[0]
    dest = Path(_CFG["outdir"]) / f"{gid}.faa.gz"
    try:
        contigs = [(n, s) for n, s in read_fasta(fna) if s]
        if not contigs:
            return gid, 0, "", "empty assembly"
        total = sum(len(s) for _, s in contigs)
        meta = total < _CFG["min_single"]

        if METHOD == "binary":
            n, meta, err = _prodigal_binary(fna, dest, meta)
            if err:
                return gid, 0, "", err
            return gid, n, ("meta" if meta else "single") + "/binary", None

        gf = pyrodigal.GeneFinder(meta=True)
        if not meta:
            try:
                gf = pyrodigal.GeneFinder()
                gf.train(*[s for _, s in contigs])
            except Exception:                      # too short or atypical to train
                gf, meta = pyrodigal.GeneFinder(meta=True), True
        n = 0
        tmp = str(dest) + ".part"
        with gzip.open(tmp, "wt") as fo:
            for name, seq in contigs:
                genes = gf.find_genes(seq)
                genes.write_translations(fo, sequence_id=name)
                n += len(genes)
        os.replace(tmp, dest)
        return gid, n, ("meta" if meta else "single"), None
    except Exception as e:                          # noqa: BLE001
        p = str(dest) + ".part"
        if os.path.exists(p):
            os.unlink(p)
        return gid, 0, "", f"{type(e).__name__}: {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", default="data/gtdb/genomes",
                    help="directory of .fna/.fna.gz")
    ap.add_argument("--out", default="data/gtdb/proteomes",
                    help="where to write the predicted .faa.gz")
    ap.add_argument("--manifest", default="data/gtdb/log/proteomes_predicted.tsv")
    ap.add_argument("--min-single", type=int, default=100_000,
                    help="below this total size, meta mode")
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--force", action="store_true",
                    help="re-predict even if the .faa.gz already exists")
    a = ap.parse_args()

    if METHOD is None:
        sys.exit(MISSING_ENGINE)

    src = Path(a.genomes)
    fnas = sorted(src.glob("*.fna")) + sorted(src.glob("*.fna.gz"))
    if not fnas:
        sys.exit(f"[FAIL] no .fna in {a.genomes}")
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.dirname(a.manifest) or ".", exist_ok=True)

    if not a.force:
        remaining = [str(p) for p in fnas
                     if not (Path(a.out) / f"{p.name.split('.fna')[0]}.faa.gz").exists()]
    else:
        remaining = [str(p) for p in fnas]
    jobs = a.jobs or max(1, (os.cpu_count() or 4) // 2)
    print(f"[predict] engine: Prodigal via {METHOD}", file=sys.stderr)
    print(f"[predict] {len(fnas)} assemblies, {len(remaining)} to process, "
          f"{jobs} in parallel", file=sys.stderr)
    if not remaining:
        print("[predict] nothing to do.", file=sys.stderr)
        return

    new = not os.path.exists(a.manifest)
    n_err = 0
    with open(a.manifest, "a", newline="") as mf:
        w = csv.writer(mf, delimiter="\t", lineterminator="\n")
        if new:
            w.writerow(["genome_id", "n_proteins", "mode", "source"])
        with make_executor(jobs, initializer=_init,
                           initargs=(a.out, a.min_single)) as ex:
            for i, (gid, n, mode, err) in enumerate(
                    ex.map(_one_genome, remaining, chunksize=2), 1):
                if err:
                    n_err += 1
                    print(f"  [error] {gid}: {err}", file=sys.stderr)
                    continue
                w.writerow([gid, n, mode, "prodigal"])
                if i % 50 == 0:
                    mf.flush()
                    print(f"  {i}/{len(remaining)}", file=sys.stderr)
    print(f"[OK] {len(remaining)-n_err} proteomes predicted -> {a.out}"
          + (f"  ({n_err} in error)" if n_err else ""), file=sys.stderr)
    print(f"[OK] provenance: {a.manifest}  (to be used as a covariate)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
