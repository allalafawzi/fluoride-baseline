# Changes

Release history of the deposit. Each entry says what changed and, where a number
moved, how to check the new one. `python3 scripts/check_deposit.py` verifies the
current state in one command.

Two published numbers have changed since the first release, both toward less
precision than was previously claimed. No point estimate has moved.

---

## Corrections that changed a published number

### The design effect was a lower bound reported as a measurement

The effective cluster size was derived from the most even possible allocation of
1259 genomes over 591 families, giving m_eff = 2.18, DEFF = 1.86, N_eff = 675.
The plan file records the real allocation: 293 families of 3, 82 of 2, 216 of 1.
The even allocation is a lower bound on m_eff, so using it overstated the
information available.

Corrected to m_eff = 2.527, DEFF = 2.115, N_eff = 595.

| interval | was | now |
|---|---|---|
| channel, 13.03 % | [10.7 ; 15.8] | [10.6 ; 16.0] |
| dehalogenase, 0 % | [0 ; 0.57] | [0 ; 0.64] |

Check: `python3 scripts/design_effect.py --plan data/gtdb/ar53_plan.tsv --analysis results/gtdb/baseline_archaea.json --prefix 1259`

### The clade-inclusion threshold was undocumented

The mean-by-phylum estimator restricts to phyla with at least 10 genomes, and the
estimate moves with that cut. The whole sweep is now published rather than the
single retained point:

| min-n | 1 | 2 | 3 | 5 | 8 | 10 | 12 | 15 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|
| estimate (%) | 11.9 | 13.2 | 13.2 | 13.2 | 14.7 | 15.6 | 16.7 | 16.1 | 13.8 | 13.1 |

The range is 13.1 to 16.7 % over cuts 3 to 30, an amplitude of 3.6 points, with
the published value high in the range. Two facts bound the consequence: the raw
proportion does not depend on this cut at all, and the qualitative conclusion,
that the per-clade estimator exceeds the raw proportion, holds from min-n 2 to 30.
It does not hold at min-n 1, where eight one- or two-genome phyla with no
positives enter the mean at equal weight and pull it to 11.9 %.

Check: `python3 scripts/clade_cut_sweep.py --analysis results/gtdb/baseline_archaea.json`

---

## Traceability of the calibration records

Each `.calib.json` now carries the filter that produced its own
`negative_max`, because that value is a maximum over a subset and the subset has
to be stated. `CLC_F` shipped `negative_max: 160.4` beside
`constraining_negatives: "prefix sp|"`, and those two statements are
inconsistent: the maximum over the `sp|` negatives is 523.2. The real filter is
the two-part rule, manually reviewed and established function name, which keeps
85 of 3351 reported sequences.

Each record also carries `n_negatives_above_GA`, `n_negatives_above_GA_reviewed`,
`unfiltered_negative_max` and a `separation_scope` string, so that
`separation: "OK"` is read as what it is. Separation is against the confirmed
negatives only; it is literally unconditional only for `Fluc_CrcB`.

| model | GA | negatives above GA | of which reviewed | absolute max |
|---|---|---|---|---|
| Fluc_CrcB | 40.95 | 0 | 0 | 13.8 |
| CLC_F | 167.15 | 39 | 2 | 605.5 |
| FAcD | 227.15 | 18 | 0 | 443.7 |

37 of the 39 CLC-F negatives above GA are automatic entries, so the
reviewed-status criterion excludes a real population.

Check: `check_deposit.py` re-scores each deposited negative set with `hmmsearch`
and counts the confirmed subset from the sequences, rather than from the record's
own `fraction_kept` arithmetic. Without `hmmsearch` on PATH it falls back to an
upper bound counted from the FASTA headers.

---

## Reproducibility of the census and the figures

- `baseline_archaea.json` and the three calibration records were written before
  the repository was standardised on English keys, so the deposited artefacts
  carried field names the deposited scripts do not write, and
  `figure_archaea.py` stopped on `KeyError: 'by_phylum'`. The keys were
  translated without recomputing anything: the multiset of 291 scalar values is
  unchanged, and the figure script's access pattern now resolves.
- `scan_proteomes.py` constructed a process pool even at `--jobs 1`, which
  queries `SC_SEM_NSEMS_MAX` and is refused in containers and on HPC login
  nodes, so the scan died before the first `hmmsearch`. It now runs serially at
  one job and falls back to serial if a pool cannot be built.
  `predict_genes.py` has the same guard.
- `hmmsearch --cut_ga` reports nothing below the threshold, so in
  `scan_archaea.tsv` a genome scoring just under it is indistinguishable from
  one with no signal: 5846 of 7224 rows (80.9 %). `scan_proteomes.py
  --all-scores` records the best score per genome and marker whatever it is, and
  the uncensored table for the archaeal panel is deposited
  (`scan_archaea_all_scores.tsv`), which is what makes the threshold-sensitivity
  analysis redoable.
- The 239 detected Fluc sequences and their alignment are deposited
  (`results/gtdb/fluc_hits.faa`, `.a2m`), so Table 10 can be redone without the
  700 GB of proteomes.
- `regenerate_checksums.sh` replaced the documented `find | xargs sha256sum`
  one-liner, which on macOS fails after the `>` redirection has already emptied
  the manifest.
- The manifest, the coverage check and the tar step read one exclusion list,
  `scripts/deposit_exclusions.py`, instead of three divergent copies.

---

## The N_eff ceiling

The deposit carried a table of N_eff ceilings for 1 000, 4 000, 20 000 and 40 000
genomes drawn, produced by `power.py` through `deff_at()`. `power.py` builds its
taxonomy with `load_taxo()`, which looks for the bacterial GTDB metadata and
falls back silently to `modelled_taxo()`, a synthetic taxonomy of 4 000 families
and 113 073 species, when that file is absent. The table therefore described the
model and not the archaeal domain, which is why two of its rows drew more genomes
than the domain contains. Its proportional column was also a single unseeded
multinomial draw: at N = 1 000 the value ranges from 259 to 387 across draws.

`scripts/neff_ceiling.py` replaces it. It names the taxonomy it used, computes
the archaeal row on the plan's own universe (591 families, 6657 representatives,
capacity 1259 at three layers), and reports a median over 200 replicates with the
range, from deterministic seeds.

| genomes drawn | proportional, median (range) | balanced, 3 layers |
|---|---|---|
| 1 000 | 199 (179 to 220) | N_eff 260 |
| 4 000 | 234 (221 to 248) | capacity exhausted at 1259 |
| 20 000 | 245 (239 to 255) | capacity exhausted |
| 40 000 | 247 (243 to 252) | capacity exhausted |

A ceiling computed under a two-level model with an assumed ICC is not comparable
to the single-level design effect `design_effect.py` measures on the drawn panel.

The plan itself was dimensioned with an assumed intra-family correlation of 0.15,
the script default, which predicts N_eff = 1024 for 1259 genomes. The measured
correlation is 0.730, giving N_eff = 595. The geometry of the plan was predicted
correctly (m_eff 2.53 predicted, 2.5266 realised); its information yield was
overstated by 72 %. `plan_sizing.json` records both values and which is which,
and the script warns when they differ.

Check: `python3 scripts/neff_ceiling.py --meta-dir data/gtdb --reps 200`

---

## Measurements added in response to specific objections

- **How far the Fluc model reaches into phyla absent from its seed**
  (`detector_reach.py`). For each phylum, the distribution of the best
  unthresholded score among undetected genomes: 1095 undetected panel genomes,
  median 8.0 bits, maximum 28.1, and 45 reaching the noise threshold of 13.8,
  against a minimum of 52.1 among the detected. The 24.0-bit empty band means the
  model is not half-seeing a divergent subfamily. It does not exclude a homolog
  divergent enough to score at background level; the decisive test needs the
  proteomes.
- **What the 95 % completeness filter moves** (`completeness_decomposition.py`).
  The 9.10-point gap decomposes into 0.15 points of phylum composition, 8.95
  points of residual at fixed phylum, and 3.61 points within the 113 mixed
  families. Genome size is not held fixed at any level and is correlated with
  completeness, so the within-family gap still contains it.
- **A penalised fit under complete separation** (`firth_phylum.py`). Three of the
  16 retained phyla have no positives, so three coefficients are unidentified in
  the unpenalised fit and the reported 15 degrees of freedom should be 12. The
  Firth fit gives finite estimates; the test reported with it is the penalised
  likelihood-ratio test.
- **A re-check of the 266 positive seeds against UniProt**
  (`seeds/uniprot_recheck_2026-10-01.tsv`), one row per seed, with the
  deposit-side columns verified against the FASTA.

---

## Declared gaps

- The UniProt and Pfam consultation dates were not recorded. `results/versions.txt`
  declares the gap rather than filling it; release 2026_01 is the only available
  bound.
- GTDB is recorded as R10-RS226, frozen; R11-RS232 has been released since.
- The bacterial column of the completeness table is incomplete at 70 % and 85 %,
  which needs `bac120_metadata.tsv.gz` (236 MB), not redistributed here.
- The standing scientific limits are in README section 7.

---

## Controls added after an independent re-reading of the deposit

Three fields were verified for their presence and never re-derived, which is the
same gap that let a confirmed-negative count one too high survive. All three are
now re-derived, and each check was tested by reintroducing the defect.

- **`n_negatives_above_GA`** is re-scored from the deposited negative set against
  the deposited model. A wrong value and a null value both fail now; previously
  either passed. The strict CLC-F record had shipped it as null behind a note
  claiming the number was not derivable from the deposit, which was wrong: the
  full negative set is deposited, and a rescore gives 10, one of them manually
  reviewed. The four values (Fluc 0, CLC-F 39, CLC-F strict 10, FAcD 18) now
  re-derive exactly.
- **`n_negatives_above_GA_reviewed`** is re-derived the same way.
- **`provenance_control`** in `baseline_archaea.json` had been null since it was
  introduced, because the documented command omitted `--predicted`. The manifest
  of Prodigal-predicted proteomes is now deposited
  (`data/gtdb/log/proteomes_predicted.tsv`, 372 rows), the documented command
  passes it, and the control reports a measurement: 124/887 = 13.98 %
  [11.85 ; 16.42] on downloaded proteomes against 40/372 = 10.75 %
  [8.00 ; 14.31] on predicted ones, with 11 phylum pairs. The intervals overlap,
  so gene prediction does not shift the prevalence detectably. A null value now
  fails the self-test and names the flag to pass.

`hmmsearch` is memoised across the checks that need it, so the suite runs in
about 45 seconds rather than twice that.

One figure quoted in the manuscript is declared in README section 7 as not
re-derivable here: the 49.9 % prevalence of the unbalanced sample depends on the
download order of the first 1510 proteomes, which was not recorded.

The self-test is at **127 checks**, and passes on a clean extraction of the archive.

---

## A language check that tested a list instead of the language

The check named "no French keys in any deposited JSON" compared every JSON
against a hand-written list of fourteen key names. That is a list of the
offenders already found, not a test of the property, and it reported green on
`seeds/facd_seed_report.json`, whose nine top-level keys were in French and on
no list. Three other deposited files were in the same state.

The files involved, all artefacts written before the repository switched to
English, whose producing scripts already emit English:

| file | was | now |
|---|---|---|
| `seeds/facd_seed_report.json` | `seuil`, `discontinuites`, `controles`, `partage_valide`, `effectifs`, `genres`, `requete`, `sans_hit`, and the nested `ecart`/`haut`/`bas`/`n_au_dessus`, `graine`/`non_etiquetes`/`negatifs`, `cote` | the names `facd_seed.py` writes: `threshold`, `discontinuities`, `controls`, `split_valid`, `counts`, `genera`, `query`, `no_hit`, `gap`/`above`/`below`/`n_above`, `seed`/`unlabelled`/`negatives`, `side` |
| `seeds/uniprot_recheck_2026-10-01.tsv` | `famille`, `statut_depot`, `nom_depot`, `PE_depot`, `longueur_aa`, `organisme`, `classe_libelle`, `alerte` | `family`, `status_deposit`, `name_deposit`, `PE_deposit`, `length_aa`, `organism`, `class_label`, `alert` |
| `results/hmm/FAcD.sensitivity.tsv` | `marge`, `rappel`, `faux_positifs`, `non_etiquetes_retenus`, `taux_FP` | the names `sensitivity_audit.py` writes: `margin`, `recall`, `false_positives`, `unlabelled_retained`, `FP_rate` |
| `results/gtdb/fluc_copy_number.json` | the group label `non classe` | `unclassified`, regenerated from the script |

No value moved in any of them: the row counts, scores, thresholds and totals
are unchanged, and the readers in `check_deposit.py` were updated with the
headers.

### The two checks that replace it

- **French by morphology, not by list.** Every deposited JSON key and TSV header
  is split on separators and tested against a list of French stems, so a French
  name nobody has seen yet still fails. UniProt accessions, which are keys in
  the seed report, are exempted by pattern rather than by name. Verified against
  48 French forms, including all fourteen from the old list and all the ones
  found here: 48 caught, 0 missed, and 0 false positives over 55 legitimate
  English names and accessions. `classes` is deliberately left out of the stem
  list, being a plausible English field name.
- **The seed report's names must be the ones its script writes.** They are read
  out of `facd_seed.py` by parsing its single `json.dump(dict(...))` call, so
  the file and its producer can no longer drift apart in any language, not only
  in French.

### Aliases removed

The plan readers carried French-to-English column aliases, `couche` to `layer`,
`famille` to `family`, `completude` to `completeness`, `taille_genome` to
`genome_size`, `n_proteines` to `n_proteins` and three more, so that a plan file
written before the rename would still load. The deposited plan has English headers, so
the aliases were dead code: removed from `analyse_archaea.py`, `baseline.py`,
`design_effect.py`, `facd_balanced.sh`, `rerun_analysis.sh`,
`download_genomes.sh` and a comment in `sample_design.py`. Re-running the census
afterwards gives the same 311 numeric values, none changed. The header check
that requires English column names stays, and it names the French spellings
because a test for French has to.

### The documentation had to be covered too

`seeds/README.md` still described two columns of the re-check table by their old
names, `classe_libelle` to `class_label` and `alerte` to `alert`, so the file
both carried French and named columns that no longer existed. A detector that
needs two French words on a line does not see a lone identifier, which is why
this survived the sweep over the data files.

Two more checks close that:

- the identifiers inside backticks in every deposited Markdown file go through
  the same morphology test as the fields themselves. `docs/CHANGES.md` is not
  exempted outright, since it has to quote what was renamed: a French name there
  is accepted only on a line that also names its English replacement, and an
  isolated one still fails;
- every column `seeds/README.md` names must exist in the table it documents,
  which catches a half-applied rename in any language.

Both were verified by reintroducing the defect. `proteome` was also removed from
the French stem list, being an English word that would have failed a legitimate
field name, and the last French fallback in `check_deposit.py`
(`"family" if "family" in rows[0] else "famille"`) is gone: the English header is
required, not preferred.

The self-test is at **130 checks**.
