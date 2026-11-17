# Project 7 — Calibration of Base-Pair Probabilities and Their Relation to Chemical Probing

**M2 GENIOMHE-AI | Université Paris-Saclay**

**RNA Bioinformatics — Group Project**

**Supervisor:** Guillaume Postic

## Project Overview

This project investigates the reliability of RNA secondary-structure pairing probabilities and their relationship with experimental chemical-probing measurements.

Three RNA secondary-structure prediction methods are evaluated:

- ViennaRNA
- LinearPartition-V
- EternaFold

The study consists of two complementary experiments.

## Experiment A — Structural-Reference Calibration

Predicted nucleotide pairing probabilities are compared against experimentally determined RNA secondary-structure annotations derived from PDB structures.

The analysis includes:

- Expected Calibration Error (ECE)
- Brier score and log loss
- Reliability diagrams
- Bootstrap confidence intervals
- Sequence and structural sensitivity analyses
- Rfam family screening
- Temperature scaling and isotonic regression

## Experiment B — Chemical-Probing Association

Predicted nucleotide unpaired probabilities are compared with Ribonanza chemical-probing measurements.

The analysis includes:

- DMS-MaP and 2A3-MaP reactivity
- Per-RNA Spearman correlations
- Sequence-cluster bootstrap analysis
- Source-specific comparisons
- Held-out regression
- Cross-source regression transfer

## Repository Organization

The repository will contain:

- `Experiment_A/` — Structural calibration analysis
- `Experiment_B/` — Chemical-probing analysis
- `reports/` — Final research report
- `presentation/` — Final project presentation
- `docs/` — Protocol and reproducibility documentation

## Reproducibility

Analysis scripts, computational settings, dataset provenance, results, and figures will be documented to support reproducibility.

## Important Limitations

Structural calibration and chemical-probing association measure different properties and are interpreted separately.

Sequence-clustering and Rfam screening do not establish complete evolutionary-family independence or prove the absence of overlap with predictor training data.
