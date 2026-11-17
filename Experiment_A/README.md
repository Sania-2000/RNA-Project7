# Experiment A — RNA Pairing Probability Calibration

**RNA Bioinformatics Project 7 | M2 GENIOMHE-AI | Université Paris-Saclay**

## Overview

Experiment A investigates the calibration of RNA nucleotide pairing probabilities predicted by three secondary-structure methods:

- ViennaRNA
- LinearPartition-V
- EternaFold

Predicted probabilities are compared with experimentally derived paired/unpaired annotations to evaluate how reliably each method represents structural uncertainty.

## Dataset and Methods

The analysis uses RNA structural references from the Archive-NoFam benchmark.

The initial dataset contained 89 RNA structures, of which 58 were retained for the primary analysis after excluding crossing annotated base pairs. Further sensitivity analyses examined subsets of 50 and 45 RNAs.

The final 45-RNA sensitivity dataset contains **2,708 nucleotide positions across 27 supplementary-section groups**.

Calibration was assessed using:

- Expected Calibration Error (ECE; 10 bins)
- Brier score
- Log loss
- Grouped bootstrap confidence intervals (5,000 resamples)

Additional analyses include sequence-similarity screening, RNA3DB-informed grouping sensitivity, and post-hoc calibration using temperature scaling and isotonic regression.

## Key Findings

EternaFold achieved the lowest nucleotide-pooled calibration errors among the three methods on the 45-RNA sensitivity datasets.

However, its advantage was less conclusive when supplementary sections were weighted equally.

An exploratory RNA3DB-informed sensitivity analysis using 26 operational groups preserved this overall conclusion.

The supplementary-section groups have not been established as independent evolutionary families.

## Repository Contents

| Folder | Description |
|---|---|
| `data/` | RNA structural-reference data |
| `analysis/` | Calibration and statistical analysis scripts |
| `results/` | Prediction outputs and result tables |
| `results/validation_2026/` | Bootstrap reproduction, sequence screening, RNA3DB grouping, and post-hoc calibration verification |
| `figures/` | Calibration plots and comparison figures |
| `scripts/` | Supporting scripts and reproducibility documentation |

## Reproducibility

The 45-RNA bootstrap analysis can be reproduced from the saved nucleotide-level datasets:

```bash
cd analysis/45_RNA_bootstrap
python reproduce_45RNA_grouped_bootstrap.py
```

Python, NumPy, and pandas are required.

Additional reproduction scripts and validation evidence are provided in the corresponding analysis and results folders.

The original probability-generation workflow is not yet fully reproducible from the repository alone. See `scripts/REPRODUCIBILITY_STATUS.md` for details.

## Limitations

The analysis uses a relatively small structural-reference dataset. Although sequence-similarity and RNA3DB checks were performed, evolutionary-family independence and complete separation from predictor training data remain unverified.

The 45-RNA and 26-group analyses are exploratory sensitivity evaluations.

## Conclusion

Experiment A shows that EternaFold has lower nucleotide-pooled calibration error on the evaluated structural-reference datasets, while group-weighted comparisons provide a more cautious interpretation.

The analysis highlights the importance of dataset selection, grouping assumptions, and statistical uncertainty when evaluating RNA pairing probabilities.

---

