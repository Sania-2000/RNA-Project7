# Software Versions — Experiment A

This file summarizes the software and computational settings used in Experiment A of RNA Bioinformatics Project 7.

## RNA Structure Prediction

The pairing probabilities were generated using ViennaRNA, LinearPartition-V, and EternaFold.

| Software | Version / Revision |
|---|---|
| ViennaRNA | 2.8.0 |
| LinearPartition-V | `b450fb3e63189073b68d385589035f992080aa3a` |
| EternaFold | `87b9aac55cee14fd562049d08f7b92d3131f10ce` |
| EternaFold parameters | `EternaFoldParams.v1` |

These versions and revisions were recorded in the original prediction notebook.

LinearPartition-V was configured with a beam size of 100, while EternaFold used a posterior probability reporting threshold of `0.00001`.

## Sequence and Homology Analysis

The following versions are suitable for a fixed, reproducible RNA homology-screening environment.

| Software / Database | Version |
|---|---|
| Infernal | 1.1.5 |
| Rfam | 15.1 |
| NCBI BLASTN | 2.12.0+ |
| RNA3DB metadata | `2026-01-05-full-release` |

BLASTN 2.12.0+ was used for the additional sequence-similarity analysis. RNA3DB metadata was used to investigate possible relationships between the supplementary-section groups.

The original Rfam screening used Infernal with model-specific gathering thresholds (`--cut_ga`). Infernal 1.1.5 and Rfam 15.1 are specified here as versions for a reproducible screening setup; the versions used in the earlier screening have not been confirmed.

## Python Environment

The following versions were recorded in Google Colab while reproducing the bootstrap and post-hoc calibration analyses.

| Package | Version |
|---|---|
| Python | 3.13.18 |
| NumPy | 2.4.6 |
| pandas | 2.2.3 |
| SciPy | 1.16.3 |
| scikit-learn | 1.6.1 |

## Reproducibility

The original 27-group bootstrap and the additional 26-group sensitivity analysis were successfully reproduced using the saved prediction datasets.

The post-hoc calibration results were also checked against the original tables. All 18 metric records agreed within an absolute numerical tolerance of `1e-8`.

These checks confirm the reproducibility of the downstream statistical calculations. They do not establish complete reproduction of the original prediction-generation workflow or full independence between evolutionary RNA families.

The analysis scripts, result tables, and additional validation files are available in the Experiment A repository.
