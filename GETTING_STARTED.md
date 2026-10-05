# Getting started

## 0. Verify integrity, then check the deposit

First, that the files are the ones we deposited:

```bash
# Linux
sha256sum -c SHA256SUMS.txt | grep -v ': OK$'
# macOS
shasum -a 256 -c SHA256SUMS.txt | grep -v ': OK$'
```

Neither should print anything. **macOS has `shasum -a 256`, not `sha256sum`** — this
bites often enough to be worth stating. To regenerate the manifest after adding a file,
always use the script rather than a hand-typed pipeline:

```bash
bash scripts/regenerate_checksums.sh
```

It picks whichever tool exists, sorts under `LC_ALL=C` so the order is identical on
every machine, and writes through a temporary file so a failure leaves the existing
manifest intact. A hand-typed `... > SHA256SUMS.txt` truncates the manifest *before*
the command runs, so on a machine without `sha256sum` it destroys the file it was
meant to produce.

Then, that the files agree with each other:

```bash
python3 scripts/check_deposit.py
```

Cross-checks between the deposited files; the script prints how many it ran. Exit code 0
means they agree. Run this before anything else: every previous review round found
defects that were invisible in any single file and only appeared on cross-check.

## 1. Environment

```bash
conda env create -f envs/environment.yml
conda activate fluoride-baseline
hmmsearch -h | head -1     # must report HMMER 3.4
```

Bit scores depend on the HMMER version through its background model. A different
version will shift every score and therefore every count.

There is no lock file. `envs/README.md` explains why: the one shipped earlier was an
export of the development machine's conda `base` environment, not the project's, and
following it left a third party unable to run anything. We removed it rather than
fabricate a replacement we could not verify.

### A trap worth naming before you replay anything

This is the most dangerous failure mode in a reproduction attempt: a replay that does
not run leaves the output file unchanged, so a value-by-value comparison reports
"identical", the most favourable verdict possible. Two ways it happens in practice are
an unresolved `{input.scan}` placeholder and a read-only extracted copy on which the
script dies on permissions. Check the exit code, and check that the output file's
modification time moved.

Before comparing any value, check that the computation actually happened: verify the
exit code, or delete the output first and confirm it was recreated. `check_deposit.py`
now runs the recipe's commands itself for this reason, rather than reading them.

## 2. Use the models as they are

```bash
hmmsearch --cut_ga results/hmm/Fluc_CrcB.hmm your_proteome.faa
```

`--cut_ga` uses the GA line in the header, which `check_deposit.py` test 1 has
verified against the calibration record.

## 3. Reproduce the census

The commands are in README §3. The `workflow/Snakefile` covers the part that replays
from the DEPOSITED tables (statistics and controls); it does NOT cover model building or
the proteome scan, which need HMMER, MAFFT and proteomes that are not in this archive.

If anything ever disagrees, **`check_deposit.py` is authoritative**: it re-derives the
published numbers from the deposited tables and fails loudly. An earlier version of this
file declared the Snakefile authoritative while that Snakefile described a pipeline that
had never been run — it omitted FAcD, used the wrong negative set for CLC_F, and called
two scripts absent from the archive. A recipe that is wrong and declared authoritative is
worse than no recipe, so the rule is now the check, not the recipe.

```bash
snakemake -s workflow/Snakefile --cores 4
```

## 4. Rebuild the models from the seeds

```bash
bash scripts/regenerate_models.sh
```

Expect identical thresholds: 40.95, 167.15, 227.15. If one differs, the seed FASTA or
the HMMER version differs from ours; check both before anything else.

## 5. Running in a container

Supported. `scan_proteomes.py` falls back to a serial executor when the platform
refuses a process pool, which is the normal case in a container and on HPC login
nodes. You will see one notice on stderr and no other difference.
