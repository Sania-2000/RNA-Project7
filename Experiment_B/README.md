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
