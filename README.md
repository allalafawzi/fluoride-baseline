# fluoride-baseline

Profile HMMs for three fluoride-handling protein families, and a balanced-sampling
census of them across archaeal genomes.

Start here: `python3 scripts/check_deposit.py`. It runs 130 cross-checks between the
deposited files and exits non-zero if anything disagrees. Each one tests an agreement
between two files that reading either file alone would not reveal.

---

## 1. What is in this deposit

| Model | Family | Seed | Window (bits) | GA | Recall on unseen genera |
|---|---|---:|---|---:|---|
| `Fluc_CrcB` | fluoride channel (PF02537) | 30 genera | [13.8 ; 68.1] | **40.95** | 98.9 % (430/435) |
| `CLC_F` | F⁻/H⁺ antiporter | 22 genera | [160.4 ; 173.9] | **167.15** | 98 % (49/50) |
| `FAcD` | fluoroacetate dehalogenase | 214 genera | [202.2 ; 252.1] | **227.15** | 100 % (355/355) |

A fourth model, `CLC_F_strict`, is shipped as a documented result rather than as a
detector: see §5.

**Headline result.** Applied to 1259 archaeal genomes drawn one per family across
the 591 GTDB families, the fluoride channel is present in 164 genomes, 13.03 %
(cluster-corrected 95 % interval [10.6 ; 16.0]). Fluoroacetate dehalogenase is
found in none of them; upper bound 0.64 %.

---

## 2. The threshold rule

The gathering threshold is the **midpoint, in bits, of the window between the best
confirmed negative and the lowest leave-one-out positive**.

This rule uses only calibration data. Nothing from the held-out genus set enters it,
so that set stays an assessment rather than a tuning criterion. `check_deposit.py`
test 2 verifies for each model that the shipped GA really is that midpoint, and that
the calibration record says so.

A confirmed negative must satisfy two independent conditions: the entry is
manually reviewed (`--confirmed-prefix 'sp|'`) and its name states an established
function (`--confirmed-pattern`). Neither suffices alone, and CLC-F proves both
directions:

- `sp|Q4VFY6|SYCA_RHITR` — reviewed, but of unknown ionic specificity. Scores 526 bits.
  Curation establishes that a protein exists and what family it belongs to, not which
  ion it moves.
- `tr|U5MTF6|U5MTF6_CLOSA` — a precise name, but propagated automatically. Scores 524 bits.

Either criterion alone admits one of them; the conjunction excludes both. Before this
rule, CLC_F carried `NC = 605.2` against `GA = 160.0` — a noise floor four times the
collection threshold, because unlabelled CLC homologues were being counted as negatives.

### What `negative_max` is a maximum *of*

It is a maximum over the confirmed negatives, which is a subset. Each record carries the
filter that produced it, so the number is re-derivable from the deposited negative set:
`check_deposit.py` re-scores that set with `hmmsearch` and counts the subset itself
rather than trusting the record's arithmetic.

| model | scored | confirmed | `negative_max` (= NC) | negatives above GA | absolute max |
|---|---:|---:|---:|---:|---:|
| `Fluc_CrcB` | 3097 | 3097 (prefix) | 13.8 | **0** | 13.8 |
| `CLC_F` | 3351 | **85** (prefix AND pattern) | 160.4 | 39 (2 reviewed) | **605.5** |
| `FAcD` | 6537 | 772 (prefix) | 202.2 | 18 (0 reviewed) | 443.7 |

So `separation: "OK"` is literally true only for `Fluc_CrcB`. For the other two it means
separation against the confirmed negatives, and each record says so in
`separation_scope`, beside `n_negatives_above_GA` and `unfiltered_negative_max`. A CLC-F
model collects bacterial chloride channels because they are its homologues — that is the
subject of §5, not a hidden defect. Of the 39 CLC-F negatives above GA, **37 are
automatic entries**: the reviewed-status criterion excludes a real population.

Rebuild all three with one command:

```bash
bash scripts/regenerate_models.sh
```

---

## 3. Reproducing the headline figure

**Read this first, it costs people an afternoon.** The `layer` column of `ar53_plan.tsv`
is a *stratum index*, not a sample/census flag. The balanced sample is the **first 1259
rows of the plan** (`--prefix 1259`, recorded in `baseline_archaea.json`). Selecting rows
by `layer`, or taking all 2408 scanned genomes, gives a different set and a wildly
different prevalence. The plan is ordered so that **any prefix is a valid balanced
sample**; that ordering is the whole point of the design (§4).

```bash
# 1. the census result, from the plan and the scan table.
#    --predicted is what populates the provenance control: it compares the
#    prevalence on downloaded proteomes against the prevalence on proteomes
#    predicted here with Prodigal. Without the flag that control is null.
python3 scripts/analyse_archaea.py \
    --plan data/gtdb/ar53_plan.tsv \
    --scan results/gtdb/scan_archaea.tsv \
    --marker Fluc_CrcB --prefix 1259 \
    --predicted data/gtdb/log/proteomes_predicted.tsv \
    --out results/gtdb/baseline_archaea.json

# 2. the figure
python3 scripts/figure_archaea.py results/gtdb/baseline_archaea.json results/gtdb/figure

# 3. the design effect, counted from the plan.
#    --out regenerates the deposited file; without it the script only prints.
python3 scripts/design_effect.py \
    --plan data/gtdb/ar53_plan.tsv \
    --analysis results/gtdb/baseline_archaea.json --prefix 1259 \
    --out results/gtdb/design_effect.json

# 4. the clade-cut sweep
python3 scripts/clade_cut_sweep.py --analysis results/gtdb/baseline_archaea.json

# 6. how far the Fluc model reaches into phyla absent from its seed
python3 scripts/detector_reach.py \
    --plan data/gtdb/ar53_plan.tsv \
    --scan results/gtdb/scan_archaea_all_scores.tsv \
    --out results/gtdb/detector_reach.json

# 7. what the 95 % completeness filter actually moves
python3 scripts/completeness_decomposition.py \
    --plan data/gtdb/ar53_plan.tsv \
    --scan results/gtdb/scan_archaea.tsv \
    --out results/gtdb/completeness_decomposition.json

# 8. the sizing of the plan: families visible per completeness threshold, and the
#    N_eff the plan was dimensioned for. This one needs a file we do not
#    redistribute; see data/gtdb/README.md for where to get it. Without it the
#    script stops with a one-line message and writes nothing.
# REQUIRES-NON-DEPOSITED-INPUT: data/gtdb/ar53_metadata_r226.tsv.gz
python3 scripts/plan_sizing.py \
    --meta-dir data/gtdb --domain Archaea \
    --out results/gtdb/plan_sizing.json

# 9. the N_eff ceiling: where proportional drawing stops buying information.
#    Runs with nothing but the deposit. With the GTDB archaeal metadata also
#    present it adds the archaeal row, which is the one the manuscript uses;
#    without it, only the synthetic model row, and the script says so loudly.
python3 scripts/neff_ceiling.py \
    --meta-dir data/gtdb --reps 200 \
    --out results/gtdb/neff_ceiling.json
```

**Why the ceiling needs its own script.** `power.py` can compute an N_eff ceiling
through `deff_at()`, but it builds its taxonomy with `load_taxo()`, which looks for
`bac120_metadata*.tsv.gz` (the *bacterial* GTDB metadata) and falls back silently to
`modelled_taxo()` when that file is absent: a synthetic taxonomy of 4 000 families and
113 073 species. A ceiling computed that way describes the model, not any real domain,
and nothing in the output said so. Its outputs also live in `results/power/`, which is
not deposited, so no documented command reached the table at all.

`neff_ceiling.py` addresses three things:

- it names the taxonomy in its output, and prints a notice rather than letting the
  synthetic one pass for a measurement when the archaeal row cannot be computed;
- it computes the archaeal row on the plan's own universe: 591 families, 6 657
  representatives, capacity 1 259 at three layers, which is the size of the plan
  actually drawn;
- it reports a median over 200 replicates with the range, from deterministic seeds,
  instead of one draw. The proportional column is a random multinomial: at N = 1 000 it
  ranges from 259 to 387 across draws, so a single draw quoted to three significant
  figures is not reproducible. Two runs of `neff_ceiling.py` give the same table.

The archaeal ceiling is about 247, and is reached by 4 000 genomes drawn; going to
40 000 moves it to 247 from 234. That is the argument for the balanced plan.

A ceiling computed under a two-level (family and order) model with an assumed ICC is
not comparable to the single-level design effect `design_effect.py` measures on the
drawn panel. The two answer different questions: what a plan of a given size could
yield, and what the plan drawn did yield.

**What command 8 showed, and why the deposit now carries it.** The plan was
dimensioned with an *assumed* intra-family correlation of 0.15, which predicts
N_eff = 1024 for the 1259 genomes drawn. The correlation *measured* on the drawn
panel is 0.730 (`design_effect.json`), which gives N_eff = 595. The geometry of the
plan was predicted correctly — m_eff 2.53 predicted, 2.5266 realised — but its
information yield was overstated by 72 %. `plan_sizing.json` now records both
values, says which is the assumption and which the measurement, and
`check_deposit.py` verifies that its corrected column lands on `design_effect.py`'s
number by an independent path. Published intervals use the measured value
throughout; nothing in the manuscript rests on 1024.

### The census table censors scores below the threshold

`hmmsearch --cut_ga` reports nothing below the gathering threshold, so in
`scan_archaea.tsv` a genome scoring just under it is indistinguishable from one with
no signal: both show an empty cell. That is 5846 of 7224 rows (80.9 %) — Fluc 65 %,
CLC_F 78 %, FAcD 100 %.

That would block the check this work most invites — a threshold-sensitivity analysis — so
the uncensored table is deposited alongside it, as
`results/gtdb/scan_archaea_all_scores.tsv`. It carries a `best_score_unthresholded`
column holding the best score per (genome, marker) whatever it is. Verified against the
census: same 7224 rows, **0 differences on the GA column, 0 on `n_copies`, 0 empty
cells**. (Reproducing the GA column exactly also proves the rescan ran under the same
HMMER version, since bit scores depend on it.)

**Two score columns, and why both are needed.** `hmmsearch --cut_ga` requires the
sequence score *and* at least one domain score to clear GA. A genome can therefore
clear GA on sequence score alone and still be legitimately absent from the census.
`GCF_001488575.1` is exactly that: CLC_F sequence score 167.5, best domain score 167.1,
GA 167.15. The table carries `best_score_unthresholded` and
`best_domain_score_unthresholded` so the census count re-derives from the scores alone —
`check_deposit.py` test 13 does exactly that, per marker, and fails if any count differs.

What it buys, on the 1259-genome balanced sample:

| threshold (bits) | positives | prevalence | |
|---:|---:|---:|---|
| 13.8 | 195 | 15.49 % | NC, lower bound |
| 27.2 | 165 | 13.11 % | old threshold |
| **40.95** | **164** | **13.03 %** | **published (midpoint)** |
| 61.3 | 158 | 12.55 % | old margin rule |
| 68.1 | 158 | 12.55 % | TC, upper bound |

From the midpoint to the upper bound of the window the proportion moves half a point,
13.03 % to 12.55 %; across the whole window it spans 12.55 % to 15.49 %. The comparison
with the confidence interval has to be made on bounds and not on widths: [12.55; 15.49]
is not contained in the naive interval [11.28; 15.00], since 15.49 falls outside it, and
is contained in the cluster-corrected [10.56; 15.97], which is the one reported
throughout.

Note also that the lower bound of the window is its least defensible point: at
NC = 13.8 bits the threshold coincides with the best verified negative, i.e. it admits
that negative.

The count applies the census rule: `--cut_ga` requires the sequence score and at least
one domain score to clear the threshold. It makes no difference at the published threshold
(164 either way) and a large one at the edges: counting on sequence score alone gives 209
at NC instead of 195. The sweep is produced by `scripts/threshold_sweep.py`, which writes
the rule it used into its own output, and `check_deposit.py` re-derives all eleven points
from the raw scores. It regenerates from the deposited tables alone:

```bash
# 5. the threshold sweep
python3 scripts/threshold_sweep.py \
    --scan results/gtdb/scan_archaea_all_scores.tsv \
    --plan data/gtdb/ar53_plan.tsv \
    --prefix 1259 --marker Fluc_CrcB \
    --out results/gtdb/threshold_sensitivity.json
```

Regenerate it yourself where the proteomes live (they are too large to deposit):

```bash
bash scripts/uncensor_scan.sh /path/to/proteomes
```

The script refuses to start if `hmmsearch` is missing, and warns if it is not HMMER 3.4.


### Clade effect beyond genome size: Firth

Three retained phyla (*B1Sed10-29*, *Korarchaeota*, *Nanohalarchaeota*) have zero
positives, so their coefficients are completely separated and the ordinary maximum
likelihood estimate does not exist. The manuscript therefore reports a **Firth penalised
logistic regression with a penalised** likelihood-ratio test:

| effect | penalised chi2 | df | p |
|---|---:|---:|---|
| phylum, size held constant | 139.11 | 16 | 1.1e-21 |
| size, phylum held constant | 43.77 | 1 | 3.7e-11 |

```bash
python3 scripts/firth_phylum.py        # writes results/gtdb/firth_phylum.json
python3 scripts/firth_selftest.py      # four validation checks, exits non-zero on failure
```

**The trap this implementation avoids.** The Firth penalty (1/2 log det I) depends on the
*design matrix*. Fitting the full and reduced models separately and subtracting their
penalised log-likelihoods compares two quantities carrying different penalties; the
penalty difference — which measures no signal at all — enters the statistic. A null
simulation of that version rejected at 100 % where it should reject at 5 %. The correct
test constrains the tested coefficients to zero without changing the design matrix, so
both models are evaluated under the same penalty. `firth_selftest.py` check 4 is the
simulation that catches this; it is the reason the file exists.

---

## 4. Confidence intervals: read this before quoting one

Genomes drawn one per family are not independent observations. The intraclass
correlation of the fluoride-channel trait, measured on the 2408 scanned genomes, is
0.730 at family rank and 0.449 at phylum rank. A Wilson interval on N = 1259
pretends otherwise.

The plan's composition is counted in `ar53_plan.tsv`, not assumed: 293 families
contribute 3 genomes, 82 contribute 2, and 216 contribute 1.

|  | value |
|---|---:|
| Kish effective cluster size m_eff | **2.527** |
| Design effect DEFF = 1 + (m_eff − 1)·ICC | **2.115** |
| Effective information N_eff | **595** (not 1259) |

| Estimate | naive | cluster-corrected |
|---|---|---|
| Channel, 13.03 % | [11.3 ; 15.0] | **[10.6 ; 16.0]** |
| Dehalogenase, 0 % | [0 ; 0.30] | **[0 ; 0.64]** |

An earlier version of this analysis derived m_eff from the most even possible
allocation of 1259 genomes over 591 families, giving 2.18 and N_eff = 675. That is a
lower bound, and reporting it as a measurement overstated the available
information. `design_effect.py` prints both, and `check_deposit.py` test 7 fails if
the recorded value is at or below the bound — i.e. if it was assumed rather than counted.

---

## 5. CLC_F does not separate its target family, and that is a result

CLC-F is a subfamily buried among far more numerous chloride antiporters. Its
admissible window is 13.5 bits wide, against 54.3 for Fluc/CrcB — the width is
itself a measure of how hard the problem is.

Splitting the seed on the three published fluoride signature motifs (`GNNLI/GMGLI`,
`GREGT/V`, `GEVTP`) shows how much the labelling decision is worth:

| seed | window | archaeal genomes called positive |
|---|---:|---:|
| relaxed (22 sequences) | 13.5 bits | **539** |
| strict (14 sequences carrying ≥ 2 of 3 signatures) | 249.2 bits | **9** |

**A factor of 60 on the count, from a decision about which sequences deserve to be in
the seed.** The eight signature-free members were pulling the model toward the
chloride side.

The decisive control: among the 598 archaeal hits of the relaxed model, the fraction
carrying all three signatures is 0 % — against 18.2 % in the seed and 0.2 % in the
chloride antiporter family used as a negative control.

Consequently no CLC-F prevalence is reported anywhere, archaeal or bacterial.
`CLC_F.hmm` is deposited so the result is checkable, not so the count is used.

---

## 6. What the FAcD seed is, and the circularity objection

`seed_audit_facd.json`, produced by `scripts/seed_audit.py`:

| | count | share |
|---|---:|---:|
| manually reviewed (`sp|`) | 3 | 1.4 % |
| automatic (`tr|`) | 211 | 98.6 % |
| evidence at protein level (PE 1) | 4 | 1.9 % |
| inferred from homology (PE 3) | 10 | 4.7 % |
| predicted (PE 4) | 200 | **93.5 %** |

So the seed is a characterised core of four sequences extended by homology to 214
genera. That is normal practice. Describing it as 214 characterised enzymes would not
be, and this deposit does not.

**The four PE 1 sequences are `Q6NAM1`, `Q1JU72`, `Q01398` and `A9BLX5`. Only the
first three are anchors.** `A9BLX5` (*Delftia acidovorans* SPH-1) does carry
protein-level evidence, but it is unreviewed and its UniProt name is the generic fold
label "Alpha/beta hydrolase fold" — no EC number, no named function. It meets the
evidence criterion and fails the established-function criterion. This is the same
conjunction used for confirmed negatives, applied on the positive side: *observing* a
protein does not establish *what it does*.

**The circularity objection.** This work argues that automatic annotations are not
reliable labels, while the FAcD seed is 98.6 % automatic. That would be contradictory
if the model had never been confronted with anything but its own source labels. Three
things prevent it, none of which is an annotation:

1. Construction — the seed/non-seed split is made on score against a preliminary
   model built from three bench-characterised enzymes, not on any label.
2. Structure — the `HDRG` nucleophile elbow, read from the sequence, correctly
   separates the twelve manually checked cases.
3. External validation — the five defluorinases characterised experimentally by
   Ji *et al.* (2026), whose public labels mention no defluorination activity, are
   recovered at 404–417 bits. A model that had only learned the annotator's habits
   would have no reason to recover precisely the enzymes that annotator missed.

---

## 7. Known limits

- CLC-F is not quantified. No archaeal prevalence is reported; the seed contains
  no archaeon.
- The detection floor is unmeasured. The seven characterised enzymes recovered
  were cryptic *by label*, not divergent *by sequence*. How far the model reaches on
  genuinely distant homologues is not known.
- A balanced plan samples a rare, clustered trait badly. All eight FAcD-positive
  genomes fall outside the balanced prefix. The plan protects a global prevalence
  against clade over-representation; it does not protect detection of a rare trait.
  For such traits, either stratify toward candidate clades or accept that the upper
  bound is the only accessible result.
- The ICC is measured on 2408 genomes and applied to the 1259 subsample, which
  assumes it is homogeneous between the two sets. Untested.
- Presence of a protein family is not physiological tolerance. Other fluoride
  tolerance routes exist (e.g. BenE-I, Ets *et al.* 2026).
- One figure quoted in the manuscript is not re-derivable from this deposit. The
  49.9 % prevalence of the unbalanced sample, the first 1510 proteomes downloaded,
  depends on the order in which those genomes were fetched, and that order was not
  recorded. It appears in no deposited file and cannot be recomputed here. It is
  given in the manuscript as the historical starting point of the work and is
  labelled there as not reproducible. Everything the conclusions rest on uses the
  balanced sample, the first 1259 rows of `ar53_plan.tsv`.

---

## 8. Repository layout

```
scripts/     analysis code, English only, one language throughout
seeds/       seed and truth-set FASTA, seed composition audit
results/hmm/ models, calibration records, validation notes
results/gtdb/census output, design effect, clade-cut sweep
data/gtdb/   the sampling plan
workflow/    Snakefile
envs/        locked conda environment
docs/        release notes
```

## 9. Software

HMMER 3.4 (Aug 2023), MAFFT, Prodigal, Python 3.12. Bit scores depend on the HMMER
version through its background model, so that version is required for any reproduction.
The seed alignment runs single-threaded: multi-threaded MAFFT is not deterministic and
drifted 1.2 bits across three runs on CLC-F, a tenth of that model's 13.5-bit window.
