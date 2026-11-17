# Dataset Selection Audit — Experiment A

This directory documents the progression of the structural-reference evaluation from 89 initial records to the 58-RNA primary dataset and the 50- and 45-RNA sensitivity subsets.

The 58-RNA dataset was obtained after excluding structures containing crossing annotated base pairs.

Eight RNA records were subsequently excluded for the 50-RNA sensitivity analysis: `1f1t_A`, `4k27_U`, `4pcj_A`, `4r4v_A`, `5m0h_A`, `8d28_A`, `8d29_C`, and `8dk7_C`.

These exclusions were based on published overlap evidence, reported sequence similarity, and conservative subsection-level filtering. They should not all be interpreted as individually confirmed training contamination.

The 45-RNA analysis excluded four additional records in both scenarios: `3egz_B`, `3p22_E`, `6t2c_V`, and `8hzj_A`. The two scenarios additionally retained either `7r6q_6` or `7u0h_6`, excluding the alternative structural annotation.

The final sensitivity datasets each contain 45 RNAs and 2,708 nucleotide positions.

The dataset-selection audit is descriptive. Complete independence from predictor training data and evolutionary-family independence remain unverified.

The complete original 89-to-58 exclusion record and individual justifications for all 50-to-45 exclusions require further documentation.
