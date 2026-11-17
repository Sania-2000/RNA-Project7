# Experiment A — Protocol and Change Log

## Original Analysis

Experiment A evaluates the calibration of RNA base-pairing probabilities from ViennaRNA, LinearPartition-V, and EternaFold.

The initial reference collection contained 89 RNA structures. After excluding structures with crossing base pairs, 58 RNAs remained. The analysis was subsequently reduced to exploratory 50-RNA and 45-RNA datasets.

Calibration was evaluated using expected calibration error (ECE), Brier score, log loss, and reliability diagrams.

## Additional Validation

To check the robustness of the original findings, I carried out the following analyses:

- **27-group bootstrap:** Repeated the uncertainty analysis using 5,000 resamples and seed 42.
- **Sequence-similarity screening:** Used BLASTN to identify highly similar RNA sequences.
- **RNA3DB comparison:** Examined whether the supplementary-section groups could be treated as independent.
- **26-group sensitivity analysis:** Combined sections 2.4 and 4.5 as an exploratory grouping test.
- **Post-hoc calibration:** Evaluated temperature scaling and isotonic regression using leave-one-section-out fitting.

## Main Findings

EternaFold had the lowest pooled ECE in the 45-RNA analysis. However, its advantage was less certain when supplementary sections were weighted equally.

The 26-group analysis gave a similar conclusion. Post-hoc calibration improved some pooled metrics for ViennaRNA and LinearPartition-V, but did not consistently improve all metrics or models.

## Limitations

Supplementary sections are not confirmed independent evolutionary families. The original prediction environment and some reference-data processing steps could not be fully reconstructed.

For this reason, the results are presented as calibration analyses with exploratory grouping-sensitivity checks, rather than definitive evidence of performance on independent RNA families.
