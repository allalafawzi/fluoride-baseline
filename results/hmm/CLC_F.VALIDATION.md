# CLC-F detector: usable model, with an explicit reservation

27 August 2026. The contrast with the Fluc/CrcB model belongs in the paper.

## The difficulty, quantified

UniProt Swiss-Prot (the manually verified part) contains 192 proteins of the CLC family
and only one mentions fluoride: `Q87WD2 ERIC_PSESM`, and even that one is annotated
ambiguously as "Chloride/fluoride channel protein".

A seed cannot be built from curated data. We had to go down into TrEMBL (the unverified
part, automatic annotation): 98 sequences named "Chloride/fluoride channel protein",
96 usable after removing 2 fragments, 48 genera.

## First attempt: failure declared by the pipeline

Taking as negatives the whole CLC family:

```
lowest leave-one-out score (CLC-F positives)   173.9 bits
highest score (negatives)                      605.2 bits
separation : "SEPARATION IMPOSSIBLE"
```

The best-scoring negative scored 3.5 times higher than the weakest true positive.
No threshold separates the two groups.

## Diagnosis: the problem is in the labels, not in the biology

Decomposing the negatives by the precision of their annotation:

| Group of negatives | n | median | max | above 152.5 |
|---|---|---|---|---|
| "chloride channel protein", vague annotation | 157 | 110.6 | 605.2 | 43 |
| confirmed ClcA ("H(+)/Cl(-) exchange transporter") | 250 | 128.4 | 158.7 | 7 |
| confirmed ClcB | 28 | 142.8 | 152.4 | 0 |
| YfeO (function unresolved) | 38 | 69.6 | 74.6 | 0 |
| miscellaneous *E. coli* proteins | 2,303 | −2.5 | 230.4 | 2 |

The chloride transporters whose annotation is precise cap at **160.4 bits**, over the
86 confirmed negatives selected by the reproducible rule set out in the next section.
CLC-F has a median of 517.5. The biological separation exists.

All of the overlap comes from the 157 vaguely annotated entries, which are in all
likelihood, for the most part, true CLC-F that TrEMBL has named generically.
That is an annotation defect, not a signal defect.

## What counts as a confirmed negative: two conditions, not one

The 316 ClcA/ClcB/YfeO of the first partition were selected by hand. The reproducible
rule that replaces that selection requires both conditions at once: the entry must be
manually reviewed (`--confirmed-prefix 'sp|'`) and carry the name of an established
function (`--confirmed-pattern`, matching "H(+)/Cl(-) exchange transporter",
"Chloride channel protein ClcB", "putative ion channel protein YfeO" or
"Chloride channel protein EriC"). Two entries show that each condition alone is
insufficient, in opposite directions:

- `sp|Q4VFY6|SYCA_RHITR`, "Symbiosis-assisting ClC homolog", *Rhizobium tropici*,
  455 aa, Pfam Voltage_CLC. Manually reviewed, but of unknown ionic specificity; it
  scores 526.5 bits. Swiss-Prot curation establishes that a protein exists and which
  family it belongs to, not which ion it moves. The prefix alone would have made it a
  negative and destroyed the separation.
- `tr|U5MTF6|U5MTF6_CLOSA`, a precise name, but propagated automatically (UniRule),
  scoring 524.8 bits. The name pattern alone would have made it a negative.

The conjunction excludes both. It yields 86 confirmed negatives, capping at 160.4 bits,
against 316 entries capping at 158.7 for the hand-made partition: a set that is
smaller, slightly higher-scoring, and above all reproducible by anyone from the
deposited command.

## Threshold adopted

Counting as negatives only the precisely annotated chloride transporters, and treating
the vague "chloride channel" entries as unlabelled:

| Threshold | CLC-F recall (50 seq., 26 unseen genera) | FP on the 316 ClcA/ClcB/YfeO of the first partition |
|---|---|---|
| 152.5 | 100.0% | 7 |
| 160 | 98.0% | 0 |
| 200 | 98.0% | 0 |
| 250 | 94.0% | 0 |

This sweep predates the reproducible rule: its false-positive column is counted against
the 316 hand-selected entries, which capped at 158.7 bits. Against the 86 confirmed
negatives now in force, which cap at 160.4, a threshold of 160 would retain one of them.

GA = **167.15 bits**, the midpoint, in bits, of the window between the highest confirmed
negative (160.4) and the lowest leave-one-out positive (173.9), that is the window
[160.4 ; 173.9]. NC = 160.40, TC = 173.90. Recall 98.0%, and zero false positives among
the 86 confirmed negatives.

The model also retains 41 vaguely annotated sequences. That is a prediction: these are
misnamed CLC-F, and correcting the annotation is the objective.

## Reservation to be written in the paper

The Fluc/CrcB seed comes from Swiss-Prot, where function is established
experimentally. The CLC-F seed comes from TrEMBL, where function is inferred by
automatic annotation.

The two prevalence figures therefore do not have the same status:

- Fluc/CrcB prevalence: a measurement, with a detection floor of 1.1% at the
  adopted threshold
- CLC-F prevalence: a lower bound, calibrated on annotation, to be validated against
  the list of characterised CLC-F of McIlwain, Ruprecht & Stockbridge (2021),
  *Annu Rev Biochem* 90:559, and Dodge, O'Connor & Wackett (2026), *Environ Microbiol*

This asymmetry is why the ">85%" of Stockbridge & Wackett 2024 cannot be reproduced:
half of their definition rests on a family that the curated databases do not yet
distinguish from the chloride transporters.

Usage: `hmmsearch --cut_ga CLC_F.hmm proteome.faa`
