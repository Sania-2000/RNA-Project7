# Experiment A — Structural Calibration and 45-RNA Sensitivity Analysis

**RNA Bioinformatics Project 7 | Université Paris-Saclay**

## 1. Objective

This analysis evaluates the calibration of RNA nucleotide pairing probabilities predicted by three RNA secondary-structure methods:

- ViennaRNA
- LinearPartition-V
- EternaFold

Predicted probabilities are compared with experimentally derived paired/unpaired annotations from RNA structural references.

## 2. Dataset Preparation

The structural-reference analysis began with 89 PDB-derived RNA structures.

After excluding structures containing crossing annotated base pairs, 58 RNAs were retained for the primary calibration analysis.

A 50-RNA subset was subsequently examined to assess potential training-data overlap.

A more conservative 45-RNA sensitivity analysis was performed to assess sequence duplication, potential training-related similarity, and conflicting structural annotations.

The 45-RNA analysis contains:

- 45 RNA records
- 2,708 nucleotide positions
- 27 supplementary-section groups

Two alternative datasets were evaluated, retaining either `7r6q_6` or `7u0h_6`, which have identical sequences but differing structural annotations.

## 3. Calibration Metrics

Three main evaluation metrics were used:

- **Expected Calibration Error (ECE):** Measures the agreement between predicted pairing probabilities and observed pairing frequencies.
- **Brier score:** Measures mean squared probability prediction error.
- **Log loss:** Evaluates the probability assigned to observed structural labels.

ECE was calculated using 10 equal-width probability bins.

## 4. Statistical Uncertainty

Grouped bootstrap resampling was performed using the 27 supplementary-section labels.

The analysis used:

- 5,000 bootstrap resamples
- Random seed 42
- Paired resampling across prediction methods
- 95% percentile confidence intervals

Both nucleotide-pooled and equal-section-weighted results were examined.

The supplementary-section labels were used as operational grouping units; they have not been confirmed to represent independent evolutionary RNA families.

## 5. Main Findings

EternaFold achieved the lowest nucleotide-pooled ECE, Brier score, and log loss in both 45-RNA sensitivity scenarios.

For the dataset retaining `7r6q_6`:

| Predictor | ECE | Brier score | Log loss |
|---|---:|---:|---:|
| ViennaRNA | 0.1213 | 0.1384 | 0.5203 |
| LinearPartition-V | 0.1213 | 0.1384 | 0.5176 |
| EternaFold | 0.0217 | 0.1035 | 0.3374 |

The nucleotide-pooled ECE advantage of EternaFold over ViennaRNA was supported by a positive paired bootstrap interval.

However, equal-section-weighted ECE comparisons were less conclusive, with confidence intervals including zero.

EternaFold showed lower within-section ECE than ViennaRNA in 13 of the 27 supplementary sections.

These results indicate that conclusions about calibration depend partly on the unit of analysis and how RNA groups are weighted.

## 6. Sensitivity Analyses

The analysis examined:

- Potential training-related sequence similarity
- Duplicate RNA sequences
- Conflicting structural annotations
- Alternative representatives of identical sequences
- Nucleotide-pooled versus group-weighted evaluation

The exclusions were used to test the robustness of the results and should not be interpreted as confirmed evidence of training-data contamination.

## 7. Reproducibility

The following files support this analysis:

- `45RNA_matched_keep_7r6q_6.csv`
- `45RNA_matched_keep_7u0h_6.csv`
- `45RNA_group_bootstrap_per_predictor.csv`
- `45RNA_group_bootstrap_paired_differences.csv`
- `45RNA_each_supp_section_metrics.csv`
- `reproduce_45RNA_grouped_bootstrap.py`

The bootstrap analysis can be reproduced using:

```bash
python reproduce_45RNA_grouped_bootstrap.py
```

The script requires Python, NumPy, and pandas.

The script reproduces the grouped-bootstrap result tables from the saved nucleotide-level prediction and reference data. Reproducing the original prediction probabilities requires the separate RNA folding software and structural-reference inputs.

## 8. Limitations

- The structural dataset is relatively small.
- Supplementary-section groups are not verified independent RNA families.
- Complete separation from predictor training data has not been established.
- The 45-RNA subset was developed through exploratory sensitivity analyses.
- Results depend partly on nucleotide versus group weighting.
- Structural-reference calibration and chemical-probing association are different evaluation tasks.

## 9. Conclusion

The structural-reference evaluation found that EternaFold had the lowest nucleotide-pooled calibration errors among the three methods on the selected RNA datasets.

The grouped sensitivity analysis also demonstrated important limitations: the observed advantage was less consistent when supplementary sections were given equal weight.

These findings motivate careful reporting of both pooled and group-level calibration performance, together with transparent documentation of sequence exclusions and potential training overlap.
