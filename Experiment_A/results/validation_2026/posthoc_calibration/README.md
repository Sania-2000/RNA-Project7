# Post-hoc Calibration Reproducibility

This directory contains reproduction and verification evidence for the Experiment A post-hoc calibration analysis.

The analysis compares raw pairing probabilities with temperature scaling and isotonic regression for ViennaRNA, LinearPartition-V, and EternaFold.

Calibration parameters were fitted using other supplementary sections, with the evaluated section held out during fitting.

All 18 calibration metric records were independently regenerated from the saved prediction and reference data. The maximum absolute difference between the original and regenerated metric tables was approximately 1.20 × 10⁻⁹, and all numerical metrics agreed within an absolute tolerance of 10⁻⁸.

Temperature scaling and isotonic regression reduced pooled ECE for ViennaRNA and LinearPartition-V. For EternaFold, recalibration produced only small ECE changes and did not consistently improve Brier score or log loss.

These results establish reproducibility of the reported downstream metric table, but not complete end-to-end reproduction of the original predictor outputs or evolutionary-family-independent validation.
