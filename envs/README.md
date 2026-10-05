# Environment

## What this directory contains

- `environment.yml` — the project dependencies, with their version constraints.

## What it does not contain, and why

There is no locked environment file (`environment.lock.yml`), and none is going to be
added.

The reason is that the environment the figures came out of no longer exists in a form
we can export faithfully. A lock file written by hand after the fact would have the
appearance of a guarantee without being one, and a false lock is worse than a missing
one because it gets believed: it prescribes an exact dependency tree that nobody ever
ran. A lock file exported from the wrong conda environment, for instance the machine's
`base` environment instead of the project's, looks the same at a glance as a correct
one, while in fact listing none of HMMER, MAFFT, numpy, scipy or snakemake.

So this directory declares constraints and records versions, and says plainly what it
cannot pin.

## What that means for reproducibility

`environment.yml` declares the version constraints that matter, and
`results/versions.txt` records the versions actually used. Together they are enough to
rebuild an equivalent environment; they do not guarantee the same dependency tree
package by package.

The HMMER version is the one that really matters for the figures, because bit scores
depend on its background model. It is constrained in `environment.yml`
(`hmmer>=3.4`), recorded in `results/versions.txt`, and verified at run time by
`scripts/uncensor_scan.sh`, which refuses to run on another version without an
explicit override.

## Installation

```bash
conda env create -f envs/environment.yml
conda activate fluoride-baseline
hmmsearch -h | head -1     # must report HMMER 3.4
python3 scripts/check_deposit.py
```
