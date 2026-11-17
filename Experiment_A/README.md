# Experiment A — RNA Pairing Probability Calibration

**RNA Bioinformatics Project 7 | M2 GENIOMHE-AI | Université Paris-Saclay**

## Overview

Experiment A evaluates how well predicted RNA nucleotide pairing probabilities agree with experimentally derived secondary-structure annotations.

We compare three RNA secondary-structure prediction methods:

- **ViennaRNA**
- **LinearPartition-V**
- **EternaFold**

The aim is to assess probability calibration, compare model performance, and examine how dataset selection and RNA grouping influence the results.

## Dataset

The analysis uses structural references from the published Archive-NoFam benchmark.

The dataset was progressively examined through several stages:

| Dataset | Description |
|---|---|
| 89 RNAs | Initial PDB-derived structural records |
| 58 RNAs | Primary analysis after excluding crossing annotated base pairs |
| 50 RNAs | Subset examined for possible training-data overlap |
| 45 RNAs | Sensitivity analysis addressing duplication, similarity, and conflicting annotations |

The final 45-RNA sensitivity analysis contains **2,708 nucleotide positions across 27 supplementary-section groups**.

Two alternative versions were evaluated, retaining either `7r6q_6` or `7u0h_6`. These records have identical sequences but different structural annotations.

## Methods

Predicted nucleotide pairing probabilities were compared with binary paired/unpaired reference labels.

Calibration was assessed using:

- **Expected Calibration Error (ECE)** with 10 equal-width probability bins
- **Brier score**
- **Log loss**

Lower values indicate better performance for all three metrics.

Statistical uncertainty was assessed using 5,000 paired bootstrap resamples, a fixed random seed of 42, and 95% percentile confidence intervals.

We compared nucleotide-pooled results with analyses that give supplementary sections equal weight.

Additional analyses examined sequence duplication, potential training-data overlap, Rfam matches, alternative structural annotations, and post-hoc calibration using temperature scaling and isotonic regression.

## Main Results

For the 45-RNA dataset retaining `7r6q_6`:

| Predictor | ECE | Brier score | Log loss |
|---|---:|---:|---:|
| ViennaRNA | 0.1213 | 0.1384 | 0.5203 |
| LinearPartition-V | 0.1213 | 0.1384 | 0.5176 |
| EternaFold | **0.0217** | **0.1035** | **0.3374** |

EternaFold achieved the lowest nucleotide-pooled calibration errors.

The paired bootstrap supported its nucleotide-pooled ECE advantage over ViennaRNA. However, the equal-section-weighted ECE comparisons were less conclusive, with confidence intervals including zero.

EternaFold also had lower within-section ECE than ViennaRNA in 13 of the 27 supplementary sections.

These findings show that the interpretation of calibration performance depends partly on how RNA observations are weighted and grouped.

## Repository Structure

- `data/reference_58/` — Reference sequences, labels, and dataset manifest
- `analysis/45_RNA_bootstrap/` — 45-RNA datasets, grouped bootstrap analysis, and reproduction script
- `analysis/Rfam/` — Rfam screening scripts and results
- `analysis/posthoc_calibration/` — Temperature scaling, isotonic regression, and held-out-section calibration analyses
- `results/` — Prediction outputs, calibration metrics, and sensitivity analyses
- `figures/` — Reliability diagrams and calibration comparisons
- `scripts/REPRODUCIBILITY_STATUS.md` — Reproducibility details and known missing components

## Reproducibility

The 45-RNA grouped bootstrap analysis can be reproduced from the saved prediction and reference data:

```bash
cd Experiment_A/analysis/45_RNA_bootstrap
python reproduce_45RNA_grouped_bootstrap.py
```

This requires Python, NumPy, and pandas.

The post-hoc calibration analysis can be reproduced from the `Experiment_A` directory:

```bash
python analysis/posthoc_calibration/reproduce_tasks5_6.py --source results/50_RNAs/ExperimentA_50RNA_matched.csv --output analysis/posthoc_calibration/reproduced
```

This additionally requires SciPy and scikit-learn.

These scripts reproduce selected downstream analyses. Complete reproduction of the original prediction generation requires additional software configurations, original structural inputs, and computational records.

## Limitations

The 45-RNA dataset is relatively small and was developed through exploratory sensitivity analyses.

Although the Rfam screening found no matches passing the model-specific gathering thresholds, this does not establish evolutionary independence.

The 27 supplementary-section labels used for grouped bootstrap and held-out calibration have not been verified as independent evolutionary families. Therefore, these analyses should not be described as fully family-disjoint.

Complete independence from the predictors' training datasets has also not been established.

Finally, experimentally derived structural annotations are imperfect reference labels, and results from this structural-calibration experiment do not directly establish agreement with chemical-probing measurements.

## Conclusion

EternaFold showed the lowest nucleotide-pooled calibration errors among the three evaluated methods on the selected 45-RNA datasets.

However, its advantage was less consistent under equal-section weighting, highlighting the importance of RNA grouping, uncertainty estimation, and transparent dataset selection.

Further validation using independently established RNA group assignments would strengthen the interpretation of these findings.


