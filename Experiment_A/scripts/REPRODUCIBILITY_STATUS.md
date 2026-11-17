# Reproducibility status — Experiment A

## Computations supported by scripts
1. **45-RNA paired/grouped bootstrap:** From inside `Experiment_A/analysis/45_RNA_bootstrap/`, run `python reproduce_45RNA_grouped_bootstrap.py`. Requires `numpy`, `pandas`; writes `reproduced_tables/`. Includes two 45-RNA matched datasets and saved summaries.
2. **Held-out supplementary-section recalibration:** From `Experiment_A/`, run `python analysis/posthoc_calibration/reproduce_tasks5_6.py --source source_inputs/ExperimentA_50RNA_results.zip --output analysis/posthoc_calibration/reproduced`. Requires `numpy`, `pandas`, `scipy`, `scikit-learn`; outputs temperature-scaled and isotonic predictions, and metric summaries. The grouping is **not verified evolutionary-family holdout**.
3. **58/50 numerical research audit:** From `Experiment_A/`, run `python scripts/reproduce_research_audit.py`. Requires `numpy`, `pandas`, the archives in `source_inputs/` and included manifest. Produces tables in `results/58_50_audit/reproduced/`. Note this script was made path-relative from its original archive version and has not been newly end-to-end reverified here.
4. **Rfam cmscan:** `analysis/Rfam/ExperimentA_Rfam_screening.py` and `.ipynb` implement the Rfam scan via a download of the then-current Rfam.cm from EMBL-EBI. Repeating the download later may use a different Rfam release. This is a reconstructed script, not an exact export of the original Colab notebook. `Rfam_45_results.tbl` is the originally supplied output with no reporting hits and an Infernal `[ok]` marker.

## Missing for complete end-to-end reproducibility
- Original production scripts and command logs for generating **all 58 RNAs'** ViennaRNA, LinearPartition-V, and EternaFold marginal probability outputs, from original sequence inputs; a limited ViennaRNA pilot is not equivalent.
- Original benchmark CT structures (and clear licensing/redistribution terms); derived labels and sequences are provided.
- Frozen dependency lockfile/environment, machine settings, predictor build commands and logs; commit hashes and major parameters appear in README.
- Exact checksum/date of Rfam covariance models and full per-query cmscan output.
- Frozen complete predictor-specific training manifests; family-disjointness has not been proven.

## Source responsibilities and publication review
The source structural benchmark and method code belong to their respective publishers. Derived data, third-party sequence content, code and figures should be reviewed for redistribution rights before publication in a **public** GitHub repository. Do not upload third-party archives or unpublished material unless authorized. No private user data or secret credentials are intentionally included, but file contents should be checked before publication.

## Interpretation
Observed improvement in pooled ECE does not establish universal group-level superiority, evolutionary-family generalization, or a causal effect. Experiment B (chemical probing) belongs in a separate sibling directory in the group repository and should not be presented as direct structural calibration.
