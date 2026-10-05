# Seeds: to be assembled before any computation

The principle that sets this census apart from that of Stockbridge & Wackett 2024:
the thresholds come from characterised enzymes, not from a generic Pfam model.
Each seed contains only proteins whose function has been demonstrated
experimentally, together with the reference that demonstrates it.

## How a threshold is set, for every marker

`build_hmm.py` applies one rule, the same in all three cases. It scores the
positives in leave-one-out cross-validation, scores the negatives, and places the
gathering threshold (GA) at the midpoint, in bits, of the window between the highest
confirmed negative and the lowest leave-one-out positive. NC is the lower edge of
that window, TC the upper edge.

Separation is evaluated against the confirmed negatives only, never against the whole
family. A confirmed negative has to meet two conditions at once: the entry is manually
reviewed and it carries the name of an established function. Entries whose annotation
is vague are treated as unlabelled, neither positive nor negative. If the window
closes, that is, if a confirmed negative outscores the weakest leave-one-out positive,
the pipeline says so (`"separation": "SEPARATION IMPOSSIBLE"`) rather than picking a
threshold anyway.

## `fluc_pos.faa`: characterised Fluc / CrcB channels

Target: >= 8 sequences, maximum taxonomic diversity.
Sources to comb through for the accessions:

| Protein | Organism | Reference |
|---|---|---|
| CrcB / Fluc-Ec | *Escherichia coli* | Baker, Sudarsan, Weinberg, Roth, Stockbridge, Breaker (2012) *Science* 335:233 |
| Fluc-Bpe | *Bordetella pertussis* | Stockbridge et al. (2013) *Nature*; McIlwain, Ruprecht, Stockbridge (2021) *Annu Rev Biochem* 90:559 |
| CrcB | *Pseudomonas putida* KT2440 | Calero, Gurdo, Nikel (2022) *Environ Microbiol* 24:5082 |
| Fluc / EriC | *Streptococcus mutans* | Men et al. (2016) *PLoS ONE* — PMID 27824896 |
| *frm* operon | *Enterobacter cloacae* FRM | Liu et al. (2017) *Sci Rep* — PMID 28754999 |
| CrcB | *Caballeronia* sp. S22 (CABS22_g0192) | Badel et al. (2024) *MRA* — assembly GCA_964261745 |
| FEX (eukaryote, outside the prokaryotic seed) | *S. cerevisiae*, *A. thaliana* | Li et al. (2013) *PNAS* — PMID 24173035; Tausta et al. (2021) *Plant Physiol* |

## `clcf_pos.faa`: characterised CLC-F F⁻/H⁺ antiporters

This is the difficult seed. CLC-F is a subfamily of PF00654 (Voltage_CLC),
drowned among the Cl⁻/H⁺ antiporters, which are far more numerous. A PF00654 hit says
nothing on its own.

| Protein | Organism | Reference |
|---|---|---|
| CLC^F | *Pseudomonas putida* ATCC 12633 | Dodge, O'Connor, Wackett (2026) *Environ Microbiol* — PMID 42469183 |
| CLC^F (subfamily lacking the chloride-binding motif) | various | McIlwain, Ruprecht, Stockbridge (2021) *Annu Rev Biochem* 90:559 |
| CLC^F | *Enterococcus casseliflavus* | Stockbridge lab, CLC^F literature |

## `clc_antiporter_neg.faa`: the negatives CLC-F cannot do without

Without them the threshold means nothing. They must contain:
- **ClC-ec1** (*E. coli*), the canonical Cl⁻/H⁺ antiporter, GSGIP motif intact
- eukaryotic CLC-0 / CLC-1
- a broad sample of uncharacterised PF00654

If separation fails against this set, the next step is a residue filter on the
alignment columns: the Ser of the GSGIP motif, which is substituted in CLC-F.

## `negatives.faa`: generic negatives

A random sample of prokaryotic proteomes (>= 2000 sequences) plus the related
paralogues. Used to measure the false positive rate.

## `fluc_truth.faa` / `clcf_truth.faa`: validation set, not included in the seed

Homologues whose function is demonstrated but deliberately held out of the seed.
This is what allows `sensitivity_audit.py` to measure the detection floor, that is the
percentage of divergent homologues the threshold misses. That figure belongs in the
paper: it is the quantity the original ">85%" leaves undeclared.

## `uniprot_recheck_2026-10-01.tsv`

An exhaustive re-check of the 266 positive seeds (214 FAcD, 22 CLC-F, 30 Fluc) against
UniProt on 1 October 2026. No sampling.

Result: 266/266 still exist, 0 changes of reviewed / unreviewed status, 0 changes of
evidence level, 1 change of label, `A0A7D8YNN0`, from which UniProt has removed the
"(Fragment)" mention that the deposited header still carries.

The table has to be read in two halves.

  `*_deposit` columns        verifiable here, permanently. check_deposit.py compares them
                           against the deposited FASTA on every run: status, label,
                           evidence level, length. If a seed changes without the
                           re-check being redone, the check fails.

  `*_uniprot_*` columns    not verifiable here. They come from a network query and will
                           age. That is why the date is in the filename.

The `class_label` column summarises what the FAcD label names: fluoroacetate (88),
haloacetate with no substrate specified (114), no dehalogenase at all (12). The
`alert` column flags lengths outside the distribution and the modified label.
