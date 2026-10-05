# data/gtdb: what is deposited and what is not

## Deposited

  ar53_plan.tsv    the sampling plan, 6 657 rows, ordered so that any prefix of the
                   file is itself a valid balanced sample. It is the only file the
                   published figures in the body of the paper depend on, and it
                   already carries the columns the analyses need (phylum, family,
                   completeness, genome size, number of proteins).

## Not deposited: the GTDB metadata file

  ar53_metadata_r226.tsv.gz   17 245 data rows, the archaeal genomes of GTDB r226,
                              and 113 columns; about 5.7 MB compressed.

We do not redistribute it: it belongs to GTDB, which publishes and versions it.
Three scripts require it and stop with an explicit message when it is missing:

  scripts/sample_design.py        produces data/gtdb/ar53_plan.tsv
      --meta PATH   the full path to the file; required argument, no detection

  scripts/plan_sizing.py          produces results/gtdb/plan_sizing.json
      --meta-dir DIR   a directory, default data/gtdb; detection by filename pattern

  scripts/completeness_audit.py   completeness audit by phylum
      --meta PATH or --meta-dir DIR   default data/gtdb

To get the file:

  https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_metadata_r226.tsv.gz

then put it in this directory. The two scripts that detect it look for the pattern
ar53_metadata*.tsv.gz (bac120_metadata*.tsv.gz for --domain Bacteria).

## Detection takes the last name in alphabetical order

  c = sorted(glob.glob(os.path.join(meta_dir, f"{pre}_metadata*.tsv.gz")))
  meta = c[-1]

If you leave two releases side by side, say ar53_metadata_r226.tsv.gz and
ar53_metadata_r232.tsv.gz, the scripts will silently take r232 and their figures will
no longer match the ones published here. Keep a single release in this directory.
Both scripts print the file they selected on standard error ([reading] ...), and
results/gtdb/plan_sizing.json records its basename in the "file" key: that is where to
check afterwards which one was used.

## How to check that you have the right file

After the two filters of the plan, completeness >= 50 % and contamination <= 5 %, there
must be 6 657 species representatives left in 591 families. results/gtdb/plan_sizing.json
carries those counts and plan_sizing.py reproduces them:

  python3 scripts/plan_sizing.py --meta-dir data/gtdb --out results/gtdb/plan_sizing.json

If your numbers differ, you do not have release 226.
