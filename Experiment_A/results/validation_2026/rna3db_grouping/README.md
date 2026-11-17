# RNA3DB Grouping Sensitivity Analysis

This directory documents an exploratory grouping assessment of the 45-RNA Experiment A sensitivity dataset.

RNA chains were compared against the RNA3DB metadata release dated 2026-01-05.

All 45 RNA chains were represented in the RNA3DB metadata. The initial mapping required additional checks to identify nested sequence-cluster members and resolve an ambiguous assignment for `6w2t_A`.

No RNA3DB sequence cluster was shared across supplementary sections in the evaluated mapping. However, `6aay_B` (section 2.4) and `7dlz_Y` (section 4.5) were assigned to the same nonzero RNA3DB component.

As an exploratory sensitivity analysis, supplementary sections 2.4 and 4.5 were merged into one operational resampling group, reducing the number of groups from 27 to 26.

The 26-group bootstrap retained the same overall interpretation as the original 27-group analysis: EternaFold showed lower nucleotide-pooled calibration error, while equal-group-weighted ECE comparisons remained less conclusive.

RNA3DB component membership should not be interpreted automatically as direct evolutionary-family membership. This analysis does not establish fully family-independent validation.
