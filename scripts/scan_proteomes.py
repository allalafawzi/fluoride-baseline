#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Census of the protein markers present in each genome.

  python3 scripts/scan_proteomes.py --proteomes data/gtdb/proteomes \
      --hmm results/hmm/Fluc_CrcB.hmm --hmm results/hmm/CLC_F.hmm \
      --out results/gtdb/scan_archaea.tsv --jobs 6 --cpu 2

Output: a TSV of genome_id, marker, n_copies, best_score, hit_ids. A genome
without the marker is written out with n_copies=0; those rows carry the
information about absence and are not padding.

By default the scan uses hmmsearch --cut_ga, which reports nothing below the
gathering threshold. In that table a genome scoring just under the threshold
looks exactly like a genome with no signal at all, so no threshold-sensitivity
analysis can be done from it. --all-scores exists for that: it adds a second
hmmsearch pass per genome and records the best sequence and domain scores
whatever the threshold.

Three choices that matter on a desktop machine:
  1. the models are concatenated into a single HMM file, so hmmsearch reads
     each proteome once instead of once per marker;
  2. genomes are processed in parallel (--jobs), each hmmsearch using only a
     few threads (--cpu). On a 6-core, 12-thread processor, --jobs 6 --cpu 2
     is considerably faster than --jobs 1 --cpu 12, because hmmsearch scales
     poorly with internal threads and well across processes;
  3. the run resumes: if the output file exists, genomes already listed in it
     are skipped, so an interrupted scan can be restarted without loss.
"""
import argparse, csv, gzip, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


class _SerialExecutor:
    """Drop-in stand-in for ProcessPoolExecutor that runs in this process.

    A process pool is constructed even when a single worker is requested, and
    constructing one queries SC_SEM_NSEMS_MAX. That query is refused in
    containers and on HPC login nodes that have no POSIX semaphores, where the
    scan would then fail before the first hmmsearch. At --jobs 1 no pool is
    built at all, and a pool that cannot be constructed falls back here.

    A serial loop changes no score, no threshold and no output line: the same
    _one_genome runs on the same inputs, one after another rather than side by
    side.
    """

    def __init__(self, max_workers=1, initializer=None, initargs=()):
        if initializer is not None:
            initializer(*initargs)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, fn, it, chunksize=1):
        return map(fn, it)


def _make_executor(jobs, initializer, initargs):
    """Parallel when possible, serial when the platform refuses a pool."""
    if jobs <= 1:
        return _SerialExecutor(initializer=initializer, initargs=initargs)
    try:
        return ProcessPoolExecutor(max_workers=jobs, initializer=initializer,
                                   initargs=initargs)
    except (OSError, NotImplementedError, ValueError) as e:
        print(f"[scan] process pool unavailable ({e}); falling back to serial.",
              file=sys.stderr)
        return _SerialExecutor(initializer=initializer, initargs=initargs)
from pathlib import Path


def read_names(hmm_paths):
    """Model names, in order, without duplicates."""
    names, seen = [], set()
    for h in hmm_paths:
        with open(h) as fh:
            for line in fh:
                if line.startswith("NAME"):
                    n = line.split(None, 1)[1].strip()
                    if n not in seen:
                        seen.add(n); names.append(n)
    return names


def merge(hmm_paths):
    """Concatenate the models: one read of the proteome instead of one per model."""
    fd, path = tempfile.mkstemp(suffix=".hmm"); os.close(fd)
    with open(path, "w") as fo:
        for h in hmm_paths:
            with open(h) as fi:
                shutil.copyfileobj(fi, fo)
            fo.write("\n")
    return path


_CFG = {}


def _init(hmm, use_ga, cpu, raw_dir=None, all_scores=False):
    _CFG["hmm"], _CFG["ga"], _CFG["cpu"] = hmm, use_ga, cpu
    _CFG["raw_dir"] = raw_dir
    _CFG["all_scores"] = all_scores


def _one_genome(faa_str):
    """Scans a proteome. Returns (genome_id, {model: {target: score}}, error)."""
    faa = Path(faa_str)
    gid = faa.name.split(".faa")[0]
    tmp_faa = dtbl = None
    try:
        if faa.name.endswith(".gz"):
            fd, tmp_faa = tempfile.mkstemp(suffix=".faa"); os.close(fd)
            with gzip.open(faa, "rb") as fi, open(tmp_faa, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            target = tmp_faa
        else:
            target = str(faa)
        fd, dtbl = tempfile.mkstemp(suffix=".dtbl"); os.close(fd)
        threshold = "--cut_ga" if _CFG["ga"] else "-E 1e-5"
        # Native HMMER output, kept when --raw-dir is given. What matters
        # there is not the hits, which are already in the TSV, but the
        # provenance footer hmmsearch writes: its own version, the exact
        # options, the query and target files, and the date. Bit scores depend
        # on the HMMER version through its background model, so without that
        # footer a census cannot be traced back to the software that produced
        # it.
        raw_out = os.devnull
        if _CFG.get("raw_dir"):
            raw_out = os.path.join(_CFG["raw_dir"], f"{gid}.hmmsearch.txt")
        r = subprocess.run(
            ["hmmsearch", *threshold.split(), "--cpu", str(_CFG["cpu"]),
             "--domtblout", dtbl, "-o", raw_out, _CFG["hmm"], target],
            capture_output=True, text=True)
        if r.returncode != 0:
            return gid, None, (r.stderr or "")[:400]
        raw = defaultdict(dict)
        with open(dtbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                if len(f) < 14:
                    continue
                target_id, model, score = f[0], f[3], float(f[7])
                if score > raw[model].get(target_id, -1e9):
                    raw[model][target_id] = score

        # Second pass. With --cut_ga, hmmsearch reports nothing below the
        # threshold, so a genome scoring just under it is indistinguishable in
        # the output from one with no signal at all, and the census table on
        # its own cannot support a threshold-sensitivity analysis. This pass
        # records the best score per model regardless of the threshold, which
        # is what fills the uncensored columns.
        best_any, best_dom = {}, {}
        if _CFG.get("all_scores"):
            fd2, tbl2 = tempfile.mkstemp(suffix=".domtbl"); os.close(fd2)
            try:
                # --domtblout, not --tblout: --cut_ga applies the gathering
                # threshold to the sequence score and also requires at least
                # one domain to clear it. A genome whose best sequence score
                # clears GA while every domain falls just under is therefore
                # legitimately absent from the census. GCF_001488575.1 is that
                # case for CLC_F: sequence 167.5, domain 167.1, GA 167.15.
                # With the sequence score alone the census recomputes to 540
                # instead of 539, so both scores are stored and the count can
                # be re-derived from the table.
                r2 = subprocess.run(
                    ["hmmsearch", "--max", "-E", "1e9", "--domE", "1e9",
                     "--cpu", str(_CFG["cpu"]), "--domtblout", tbl2,
                     "-o", os.devnull, _CFG["hmm"], target],
                    capture_output=True, text=True)
                if r2.returncode == 0:
                    with open(tbl2) as fh:
                        for line in fh:
                            if line.startswith("#"):
                                continue
                            f = line.split()
                            if len(f) < 14:
                                continue
                            # domtblout columns (1-based): 4 = query name (the model),
                            # 8 = full-sequence score, 14 = this-domain score.
                            model = f[3]
                            sq, dm = float(f[7]), float(f[13])
                            if sq > best_any.get(model, -1e9):
                                best_any[model] = sq
                            if dm > best_dom.get(model, -1e9):
                                best_dom[model] = dm
            finally:
                Path(tbl2).unlink(missing_ok=True)
        return gid, (dict(raw), best_any, best_dom), None
    except Exception as e:                       # noqa: BLE001
        return gid, None, f"{type(e).__name__}: {e}"
    finally:
        for t in (tmp_faa, dtbl):
            if t and os.path.exists(t):
                os.unlink(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomes", required=True,
                    help="directory of .faa/.faa.gz, file name = genome_id")
    ap.add_argument("--hmm", required=True, action="append",
                    help="HMM path (repeatable); the marker takes the NAME of the model")
    ap.add_argument("--no-ga", action="store_true",
                    help="use E<1e-5 instead of the calibrated GA threshold")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpu", type=int, default=2, help="threads per hmmsearch")
    ap.add_argument("--all-scores", action="store_true",
                    help="also record the best score per model BELOW the threshold, "
                         "in a best_score_unthresholded column. Without it a genome "
                         "under threshold is indistinguishable from one with no hit, "
                         "and no threshold-sensitivity analysis is possible from the "
                         "table. Costs a second hmmsearch pass per genome.")
    ap.add_argument("--jobs", type=int, default=0,
                    help="genomes in parallel (0 = auto: physical cores)")
    ap.add_argument("--raw-dir", default=None,
                    help="directory in which to keep the native hmmsearch "
                         "output for every genome, one file each. Each file "
                         "carries HMMER's own provenance footer: version, "
                         "exact options, query and target files, date. Deposit "
                         "this directory alongside the TSV so a census can be "
                         "traced back to the software that produced it -- bit "
                         "scores depend on the HMMER version through its "
                         "background model. Costs disk, not time.")
    ap.add_argument("--restart", action="store_true",
                    help="ignore the existing output and redo everything")
    a = ap.parse_args()

    d = Path(a.proteomes)
    faas = sorted(d.glob("*.faa")) + sorted(d.glob("*.faa.gz"))
    if not faas:
        sys.exit(f"[FAIL] no .faa or .faa.gz in {a.proteomes}")
    models = read_names(a.hmm)
    if not models:
        sys.exit("[FAIL] no NAME line in the HMM files provided")
    if a.raw_dir:
        os.makedirs(a.raw_dir, exist_ok=True)
        print(f"[scan] native hmmsearch output kept in {a.raw_dir}",
              file=sys.stderr)


    jobs = a.jobs or max(1, (os.cpu_count() or 4) // 2)

    # Resume: skip the genomes already present in the output file.
    done = set()
    if os.path.exists(a.out) and not a.restart:
        with open(a.out) as fh:
            rd = csv.reader(fh, delimiter="\t")
            next(rd, None)
            for row in rd:
                if row:
                    done.add(row[0])
    remaining = [str(p) for p in faas if p.name.split(".faa")[0] not in done]

    print(f"[scan] {len(faas)} proteomes, {len(models)} markers "
          f"({', '.join(models)})", file=sys.stderr)
    if done:
        print(f"[scan] {len(done)} already scanned -> {len(remaining)} to do "
              f"(--restart to redo everything)", file=sys.stderr)
    print(f"[scan] {jobs} genomes in parallel x {a.cpu} threads per hmmsearch",
          file=sys.stderr)
    if not remaining:
        print("[scan] nothing to do.", file=sys.stderr)
        return

    hmm = merge(a.hmm)
    new = not done or a.restart
    mode = "w" if new else "a"
    n_err = 0
    try:
        # newline="\n", not "": csv.writer defaults to CRLF line endings,
        # which would leave mixed line endings across the deposit.
        with open(a.out, mode, newline="\n") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            if new:
                header = ["genome_id", "marker", "n_copies", "best_score", "hit_ids"]
                if a.all_scores:
                    header += ["best_score_unthresholded",
                               "best_domain_score_unthresholded"]
                w.writerow(header)
            with _make_executor(jobs, _init,
                                (hmm, not a.no_ga, a.cpu, a.raw_dir,
                                 a.all_scores)) as ex:
                for i, (gid, res, err) in enumerate(
                        ex.map(_one_genome, remaining, chunksize=4), 1):
                    if err is not None:
                        n_err += 1
                        print(f"  [error] {gid}: {err.splitlines()[0] if err else '?'}",
                              file=sys.stderr)
                        continue
                    raw_res, best_any, best_dom = res
                    for m in models:
                        hits = raw_res.get(m, {})
                        row = [gid, m, len(hits),
                               f"{max(hits.values()):.1f}" if hits else "",
                               ";".join(sorted(hits))]
                        if a.all_scores:
                            b, d = best_any.get(m), best_dom.get(m)
                            row.append(f"{b:.1f}" if b is not None else "")
                            row.append(f"{d:.1f}" if d is not None else "")
                        w.writerow(row)
                    if i % 100 == 0:
                        fh.flush()
                        print(f"  {i}/{len(remaining)}", file=sys.stderr)
    finally:
        os.unlink(hmm)

    print(f"[OK] {len(remaining) - n_err} genomes scanned x {len(models)} markers "
          f"-> {a.out}" + (f"  ({n_err} in error)" if n_err else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
