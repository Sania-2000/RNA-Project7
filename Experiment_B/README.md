# Experiment B — Chemical-Probing Association

## Project 7: RNA Secondary-Structure Probabilities

This folder contains the Student 2 contribution to Project 7.

The objective is to evaluate how predicted nucleotide-level
unpaired probabilities relate to experimental RNA chemical
reactivity measured by DMS-MaP and 2A3-MaP.

## Dataset

- Source: Ribonanza RNA chemical-probing dataset
- RNA sequences analyzed: 20,739
- Experimental profiles: 41,478
- Source groups: 17
- CD-HIT-EST sequence clusters: 16,520
- Sequence clustering threshold: 90% identity

The analyzed sequences passed provisional completeness-based
quality-control filters. This is not equivalent to official
signal-to-noise quality certification.

## Computational predictors

Three RNA secondary-structure probability predictors were evaluated:

1. ViennaRNA
2. LinearPartition-V
3. EternaFold

Each predictor generated nucleotide-level unpaired probabilities.

## Spearman correlation results

Median per-RNA Spearman correlations:

| Predictor | DMS-MaP | 2A3-MaP |
|---|---:|---:|
| ViennaRNA | 0.6264 | 0.4244 |
| LinearPartition-V | 0.6268 | 0.4239 |
| EternaFold | 0.6409 | 0.4616 |

Cluster-bootstrap confidence intervals and paired comparisons
are available in `results/tables/`.

## Linear regression

Ordinary least-squares regression modeled chemical reactivity
as a function of predicted unpaired probability.

An 80/20 CD-HIT-cluster-separated training/test split was used.

Held-out test performance:

| Predictor | DMS R2 | DMS RMSE | 2A3 R2 | 2A3 RMSE |
|---|---:|---:|---:|---:|
| ViennaRNA | 0.2191 | 0.7357 | 0.1255 | 0.5625 |
| LinearPartition-V | 0.2191 | 0.7357 | 0.1256 | 0.5625 |
| EternaFold | 0.2360 | 0.7277 | 0.1630 | 0.5503 |

## Source-to-source transfer

Regression models were trained on one RNA experimental source
and evaluated on another, excluding target RNAs belonging to
training-source CD-HIT clusters.

Two source pairs were evaluated in both directions:

- OpenKnot EternaPlayers and OpenKnot Round 2
- 15k and RNAmake designs

EternaFold achieved the highest transfer R2 in the eight
source-direction/probe comparisons tested.

However, transfer performance varied between sources.
In particular, RNAmake-to-15k DMS transfer produced
negative R2 values for all three predictors.

## Results and scripts

- `scripts/`: Python analysis and prediction scripts
- `results/tables/`: final statistical results
- `results/figures/`: comparison figures in PNG/PDF
- `environment/`: Conda environment and software versions
- `docs/`: Student 2 report

## Scientific limitations

1. QC used reactivity completeness thresholds rather than
   independently verified signal-to-noise certification.
2. Experimental sources have unequal sample sizes.
3. CD-HIT clusters reduce close-sequence overlap but do not
   establish independence between biological RNA families.
4. Minor numerical probability-bound violations were
   corrected for analysis; raw values were preserved in
   the original computational project.
5. The source-transfer tests cover selected source pairs
   and should not be interpreted as exhaustive cross-source
   validation.
6. Association with chemical reactivity is distinct from
   calibration against experimentally known RNA structures.

## Data availability

The original Ribonanza dataset and large intermediate prediction
files are not included in this repository because of their size.

The full computational outputs are retained separately.
See the analysis scripts and environment files for methodology.

## Relationship to Experiment A

Experiment A evaluates predicted probabilities against
reference RNA structures.

Experiment B evaluates associations between predicted
unpaired probabilities and chemical-probing measurements.

These are complementary analyses with different evaluation
targets.

## Large-scale dataset preparation

The large-scale analysis used RNA sequences and DMS-MaP and
2A3-MaP reactivity profiles from the Ribonanza QUICK_START dataset.

The dataset was prepared in two selection stages:

1. `select_large_ribonanza.py` selected 20,000 candidate RNAs.
2. `qc_large_ribonanza.py` retained 13,826 RNAs passing
   the project completeness checks for both probes.
3. `select_additional_ribonanza.py` selected another 10,000
   candidates from the remaining sequences (random seed 43).
4. `qc_additional_ribonanza.py` retained 6,913 additional RNAs
   passing the checks for both probes.
5. `combine_large_ribonanza.py` combined the two groups into
   `combined_20739_qc_passed.csv`.

The final dataset contained 20,739 unique RNA sequences
(13,826 + 6,913). The combination script checks for duplicate
sequence identifiers and duplicate RNA sequences.

DMS analyses used eligible A/C positions, while 2A3 analyses
used all measured nucleotide types. The project completeness
checks required at least 10 valid measurements and at least
50% coverage of eligible positions. These checks are separate
from any filtering applied to the original Ribonanza dataset.

### Sequence clustering

Sequence similarity was assessed using CD-HIT-EST at 90%
sequence identity.

`parse_large_clusters.py` processes the CD-HIT cluster output
and produces cluster membership and summary tables.

The final dataset contained 16,520 sequence clusters.

### Prediction batches

`create_large_batches.py` reads
`large_scale/prediction_inventory_20739.csv`, sorts RNAs by
sequence identifier, and divides the dataset into batches
of 500 sequences.

This produced 42 batches: 41 containing 500 RNAs and a final
batch containing 239 RNAs.

The script also saves `large_scale/batch_manifest.csv`.

The combined QC dataset and prediction inventory were checked
to contain the same 20,739 sequence identifiers, with no
missing or extra identifiers.

### Reproducibility notes

The original scripts were run from the RNA_Project7 working
directory, with input and intermediate files under `large_scale/`.

The original Ribonanza dataset and large intermediate prediction
files are not stored in this GitHub repository.

The exact procedure that added cluster metadata to produce
`prediction_inventory_20739.csv` still needs to be documented.
The final inventory itself has been verified against the
combined QC-passing dataset.

