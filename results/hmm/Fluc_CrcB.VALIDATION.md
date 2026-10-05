# Fluc/CrcB detector: validated model

27 August 2026. Built and validated on real UniProt/Swiss-Prot data.

## Data

| Set | Content | N |
|---|---|---|
| Seed | 1 sequence per genus, 30 genera | 30 |
| Validation (*held-out*) | all the sequences of the 164 genera never seen by the model | 435 |
| Negatives | *E. coli* K-12 proteome (minus its own FluC) + 191 CLC transporters | 4,593 |

Source: UniProt Swiss-Prot, `(xref:pfam-PF02537) AND (reviewed:true)`, 522 sequences in
194 genera. All are annotated "Fluoride-specific ion channel FluC" or "Fluoride export
protein", so their function is established by manual curation with references.

## Calibration

```
lowest leave-one-out score (positives)        68.1 bits
highest score (confirmed negatives)           13.8 bits
gap                                           54.3 bits
threshold window [confirmed negative max ; LOO min]    [13.80 ; 68.10]
GA = midpoint of the window                   40.95 bits
NC 13.80    GA 40.95    TC 68.10
```

## Recall / false positive trade-off

| Margin | GA (bits) | Recall on 435 seq. from unseen genera | False positives / 4,593 |
|---|---|---|---|
| 1.20 | 81.7 | 90.3% | 0 |
| 1.00 | 68.1 | 95.2% | 0 |
| 0.90 | 61.3 | 97.2% | 0 |
| 0.80 | 54.5 | 98.2% | 0 |
| **0.60** | **40.9** | 98.9% | 0 |
| 0.40 | 27.2 | 100.0% | 0 |
| 0.20 | 13.6 | 100.0% | 1 |

Threshold adopted: GA = 40.95 bits, the midpoint, in bits, of the window between the
highest confirmed negative (13.8) and the lowest leave-one-out positive (68.1).

Detection floor: 1.1%. Five of the 435 held-out sequences score below the threshold
(recall 98.9% at the nearest swept point, GA 40.9), with zero false positives. A
threshold of 27.2 bits misses none of them and still admits no false positive, so the
adopted value is conservative: it is the one the documented midpoint rule gives, not the
one that maximises recall.

The single false positive at the extreme threshold (13.6) is `P0A9R2 ESSD_ECOLI`, a
prophage lysis protein, score 13.8, hence below even the naive threshold itself. It has
no effect on the adopted threshold.

## Controls

Positive control. The FluC of *E. coli* (`P37002`), removed from the negative set, is
recovered at 137.7 bits (E = 1.1 × 10⁻⁴⁴), more than 3 times the threshold.

Specificity on real environmental data. On 24,718 ORFs from a metagenome (HaloCycDB
example set):

| Detector | ORFs accepted | Rate |
|---|---|---|
| HaloCycDB `DRG` regex (the "fluoroacetate dehalogenase" filter) | 684 | 2.77% |
| This model, GA = 27.2 (threshold in force when the test was run) | **3** | **0.012%** |

A factor of 228 in specificity, on the same data. The test has not been rerun at
GA = 40.95; since raising the threshold can only remove accepted ORFs, that factor is a
lower bound.

## Files

- `Fluc_CrcB.hmm`: model, GA/TC/NC thresholds written into the header
- `Fluc_CrcB.calib.json`: full calibration log, per-sequence LOO scores
- `Fluc_CrcB.sensitivity.tsv`: the table above
- `Fluc_CrcB.seed.aln`: seed alignment

Usage: `hmmsearch --cut_ga Fluc_CrcB.hmm proteome.faa`
