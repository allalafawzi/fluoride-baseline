# FAcD detector: validation record

27 August 2026. Equivalog of fluoroacetate dehalogenase (EC 3.8.1.3), α/β-hydrolase
fold (Pfam PF00561). Third model of the project, after Fluc/CrcB and CLC-F.

## Why the EC number could not serve as a label

Two documented annotation errors, in two opposite directions, both of them
verifiable by anyone.

Same label, different biologies. `Q01399` (DehH2 of *Moraxella* sp. B) carries
EC 3.8.1.3 in the manually reviewed part of UniProt. Yet its own entry states:
"*dehH2 acts on chloro-, bromo- and iodoacetate, but not on fluoroacetate*". It
belongs moreover to the HAD superfamily (PF00702), not to the α/β-hydrolases. Against
the preliminary model built on the three characterised enzymes, it obtains no hit at
all, even at E = 1000.

Different labels, same biology. `A0A0C4Y433` and `A0A0C4Y4A6` come from the same
organism (*Cupriavidus basilensis*), same family PF00561, 299 and 296 amino acids,
consecutive accession numbers. The first is annotated EC 3.8.1.5 and therefore ended up
among the negatives; the second EC 3.8.1.3 and ended up among the positives. Neither
has been characterised. The split rested on the number that an automatic pipeline had
assigned.

## Building the seed

Three experimentally characterised anchors, and nothing else:

| Accession | Enzyme | Organism |
|---|---|---|
| `Q6NAM1` | RPA1163 | *Rhodopseudomonas palustris* |
| `Q1JU72` | FAc-DEX FA1 | *Burkholderia* sp. |
| `Q01398` | DehH1 | *Moraxella* sp. B |

A preliminary model built on those alone scored 738 candidates. The splitting threshold
was not chosen: it was read off the data, at the largest discontinuity lying below the
weakest anchor score. The discontinuity retained is 28.1 bits (193.1 → 165.0), against
15.9 / 12.0 / 11.8 for the next ones.

| Control | Score | Expected side | Side obtained |
|---|---:|---|---|
| Q01398 DehH1 | 586.3 | high | high |
| Q1JU72 FAc-DEX | 575.8 | high | high |
| Q6NAM1 RPA1163 | 540.0 | high | high |
| Q01399 DehH2 | no hit | low | low |

Independent validation of the split. Twelve sequences had been identified by hand
as lacking the `HDRG` nucleophile elbow shared by the three anchors
(epoxide hydrolases and other subfamilies). All twelve fall between 38 and 88 bits,
hence below the threshold. The eight identified as canonical are between 301 and 381.
The two criteria are independent: neither was used to build the other.

Result: a seed of 214 sequences over 214 genera, a validation set of 355
sequences from genera absent from the seed, 171 unlabelled sequences.

## Calibration

The first attempt declares SEPARATION IMPOSSIBLE (highest negative 443.7 > GA 227.15).
The decomposition of the negatives by annotation precision explains it:

| Source | n | max | median |
|---|---:|---:|---:|
| Swiss-Prot (manually reviewed) | 771 | 202.2 | 25.2 |
| TrEMBL (automatic) | 5,765 | 443.7 | 91.1 |

The negatives verified by a human cap at 202.2; the minimum of the positives in
leave-one-out cross-validation is at 252.1. Threshold window: **[202.2 ; 252.1]**,
49.9 bits wide. The ten best-scoring Swiss-Prot negatives are epoxide hydrolases
(human, mouse, pig, *Xenopus*, *C. elegans*, *Bacillus subtilis*) between 103 and
202 bits: genuinely different enzymes, correctly excluded.

The 5,765 TrEMBL entries are therefore treated as unlabelled, never as
negatives. It is the same correction that rescued CLC-F, applied here to an
independent family.

```
GA    227.15
TC    252.10
NC    202.20
separation : OK
```

## Performance

| Test | Result |
|---|---|
| "Leave-one-out" cross-validation, 214 sequences | 100% recall at the GA threshold |
| Genera unknown to the seed, 355 sequences | 355 / 355 = 100% |
| False positives on 771 verified negatives | 0 |

Threshold sweep:

| Margin | GA | Recall | False positives | Unlabelled retained |
|---:|---:|---:|---:|---:|
| 1.00 | 252.1 | 99.7% | 0 | 13 |
| **0.90** | **226.9** | 100% | 0 | 18 |
| 0.80 | 201.7 | 100% | 1 | 26 |
| 0.50 | 126.1 | 100% | 2 | 221 |
| 0.30 | 75.6 | 100% | 35 | 4,175 |

Threshold adopted: GA = 227.15 bits, the midpoint, in bits, of the window between the
highest confirmed negative (202.2) and the lowest leave-one-out positive (252.1).
NC = 202.20, TC = 252.10.

Above the window one loses a positive; below it one gains a false positive, and
that one is `O52866` (HYES_CORS2, epoxide hydrolase) at 202.2 bits, exactly the
announced bound. Any threshold in [202.2 ; 252.1] gives 100% recall and 0 false
positives: the choice of the margin inside the window has no effect.

## A by-product: 18 public annotations that are probably wrong

At the adopted threshold, the model retains 18 TrEMBL sequences annotated as
haloalkane dehalogenases or epoxide hydrolases. Given the *Cupriavidus* case, these are
in all likelihood mislabelled fluoroacetate dehalogenases. The list is to be attached as
supplementary material, as a proposed correction to the public annotation.

## Reservations to be written in the paper

The detection floor is not zero, it is unmeasured. The sweep gives 0%
loss, but on a validation set drawn from the same UniProt pool, that is sequences
already annotated as FAcD. The loss on genuinely divergent environmental
homologues cannot be measured, for lack of characterised divergent enzymes.
That is a limitation of the field, not of the model, and it has to be declared as such.

Taxonomic composition of the seed. The 214 genera are mostly
bacteria, with a share of fungi (*Fusarium*, *Trichoderma*, *Metarhizium*,
*Colletotrichum*…), of cyanobacteria and a few eukaryotes. The archaeal domain is
represented only by *Haloarcula* (plus one Euryarchaeota entry, `A0A830GNL7`, at
398.3 bits). That is more than CLC-F, where there was no archaeon at all, but it remains
thin: any archaeal figure produced by this model has to carry that reservation.

The splitting discontinuity is a dip, not a chasm: 28.1 bits against 15.9 for
the next one. What makes it acceptable is the difference in density on either side.
Above, the lowest scores are at 193–219; below, only 14 sequences out of
171 exceed 100 bits and the other 157 are between 23 and 100.
