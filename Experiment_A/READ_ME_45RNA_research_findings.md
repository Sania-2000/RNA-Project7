# Experiment A: combined 45-RNA sensitivity and grouped uncertainty (2026-11-16)

Source: user-provided `ExperimentA_50RNA_results.zip`, nucleotide-level matched probability and experimental-label table. No predictor probabilities were rerun. 50-RNA original has 2,969 nucleotides and 29 supplementary-section labels.

## Two alternative 45-RNA scenarios

Both remove `6t2c_V` (plausibly related to EternaFold training influenza promoter RNA) and `3egz_B` (20-nt exact segment shared with a released training sequence), plus duplicate records `3p22_E` and `8hzj_A`. Then one of the disagreeing identical-sequence pair `7r6q_6` / `7u0h_6` is removed, retaining the other. Each scenario has 45 records, 2,708 nucleotide positions and 27 supplementary-section labels. These exclusions are *sensitivity analyses*, not confirmed contamination removals. No ground-truth labels have been edited.

## Measurements

* Pooled 10-bin ECE weights every nucleotide equally. Bins are [0,.1), ..., [.9,1] and p=1 falls in the last bin. Brier score and Bernoulli log loss are averaged per nucleotide, using numerical clipping 1e-12 only for log loss.
* Supplementary-section paired bootstrap: 5,000 draws, random seed 42, sample the 27 *published supplementary-section labels* with replacement, using the same drawn sections for all three predictors. The intervals are percentile intervals, not formal hypothesis tests. These are not independently verified evolutionary families nor identical to the original paper's 47 structural-group labels.
* `ECE_nt_pooled`: reaggregate all sampled nucleotide positions for each resample, so long sections have more influence.
* `ECE_equal_section_mixture`: each published section has equal weight when mixing the section-normalized probability-bin counts, then absolute bin deviations are computed. Equal-group pooled ECE does NOT equal the average per-section ECE.
* `ECE_mean_section`: calculate ECE separately within every section, then average the 27 scores equally. This is a separate estimand, with different bin behavior for short sections.
* Brier/log-loss `mean_section` variants weight each section equally rather than each nucleotide.

## Primary interpretation

For both scenarios, EternaFold has the lowest nucleotide-pooled ECE, Brier score, and log loss. The paired bootstrap 95% percentile interval for ViennaRNA-minus-EternaFold pooled ECE is positive: [0.0050, 0.1367] (retain 7r6q_6), [0.0059, 0.1366] (retain 7u0h_6). HOWEVER both equal-section-mixture ECE intervals include zero: [-0.0219, 0.0895] and [-0.0212, 0.0900]. Mean per-section ECE actually has a slightly higher point estimate for EternaFold than ViennaRNA; the difference interval includes zero. EternaFold has lower within-section ECE than ViennaRNA in only 13/27 supplementary sections in either scenario. The model ranking is sensitive to how calibration and sampling units are defined.

## Remaining unresolved limitations

- RNA-family-independent generalization cannot be claimed: the supplementary-section groups are not proven independent evolutionary families, and the exact frozen parameter-file training manifest was not established.
- There may be synthetic and structurally related training RNA data not represented in the released sources.
- The candidate 50/45 sets were chosen through successive sensitivity screens, so the sensitivity subsets are not a prospectively defined independent test.
- Ribonanza chemical-probing comparisons belong to Experiment B, and should not be treated as experimental base-pair labels.

## Contents

`45RNA_matched_keep_7r6q_6.csv`, `45RNA_matched_keep_7u0h_6.csv`: fully matched per-nucleotide data.
`45RNA_group_bootstrap_per_predictor.csv`: 95% intervals for each method/scenario/metric.
`45RNA_group_bootstrap_paired_differences.csv`: paired contrasts with 95% intervals.
`45RNA_each_supp_section_metrics.csv`: per-section ECE, Brier, and log loss.
`reproduce_45RNA_grouped_bootstrap.py`: portable script that reproduces the three summary tables from the two matched 45-RNA CSVs in this package. Run `python reproduce_45RNA_grouped_bootstrap.py` from a Python environment containing `numpy` and `pandas`. Its outputs appear in `reproduced_tables/`. It does not recompute the original predictor probabilities, which must be reproduced separately from predictor software and structural inputs.

