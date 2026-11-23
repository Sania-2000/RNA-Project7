# Experiment B — Protocol and Change Log

## Objective

Experiment B examined the association between predicted RNA unpaired probabilities and chemical-probing measurements from the Ribonanza dataset. Three predictors were compared: ViennaRNA, LinearPartition-V and EternaFold.

## Experimental protocol

### 1. Dataset preparation

The analysis used the Ribonanza QUICK_START training dataset.

The workflow began with a small pilot and a development dataset before moving to the final analysis.

For the large-scale dataset, 20,000 initial RNA candidates were selected, of which 13,826 passed the project quality-control checks for both probes. Another 10,000 candidates were sampled from the remaining sequences using random seed 43, and 6,913 passed the same checks.

The final dataset contained 20,739 unique RNA sequences, each with DMS-MaP and 2A3-MaP measurements.

### 2. Sequence clustering

CD-HIT-EST was used at 90% sequence identity. The final dataset contained 16,520 sequence clusters.

Cluster assignments were used for bootstrap resampling and regression train/test separation.

### 3. Structure prediction

Base-pair probabilities were calculated using ViennaRNA, LinearPartition-V and EternaFold.

For nucleotide i, the unpaired probability was calculated as one minus the sum of its pairing probabilities.

Small numerical violations outside the interval [0, 1] were handled using a documented tolerance and clipping procedure. Raw values were preserved separately.

### 4. Large-scale processing

The 20,739 RNA sequences were processed in 42 batches: 41 batches of 500 sequences and one final batch of 239.

The three predictors were evaluated using matching RNA identifiers and nucleotide positions.

### 5. Statistical analysis

Spearman correlations were calculated separately for DMS-MaP and 2A3-MaP. DMS comparisons used adenine and cytosine positions; 2A3 comparisons used all measured nucleotide types.

The analysis included paired predictor comparisons, 500 sequence-cluster bootstrap replicates, source-specific comparisons and length-stratified results.

Linear regression was evaluated using cluster-separated training and testing datasets. Source-to-source transfer was tested in both directions for two selected RNA-source pairs.

## Change log

| Stage | Change and reason |
|---|---|
| Initial pilot | Tested data extraction, nucleotide alignment and prediction calculations on 10 RNA sequences. |
| Development analysis | Expanded testing to a length-diverse selection; 139 RNAs passed the development checks. |
| Probability validation | Checked predictor output formats and numerical probability-bound deviations before large-scale processing. |
| Full analysis | Expanded the study to 20,739 RNAs using 42 prediction batches. |
| Statistical evaluation | Added cluster-bootstrap uncertainty estimates, regression and source-to-source transfer testing. |
| Length analysis | Added final length-stratified Spearman results and paired comparisons. |
| Reproducibility update | Added the selection, QC, clustering, batching and inventory-reconstruction scripts to GitHub. |
| Software documentation | Recorded the source revisions for LinearPartition-V and EternaFold. |

## Reproducibility

The scripts are located in `Experiment_B/scripts/`. Statistical output tables and figures are located in `Experiment_B/results/`.

The software versions and Conda environment exports are stored in `Experiment_B/environment/`.

The original Ribonanza dataset and large intermediate prediction files are not included in GitHub because of their size.

The workflow and input-generation procedure are documented in `Experiment_B/README.md`.
