# Structural Reference Dataset — Experiment A

## Data source

We used PDB-derived RNA secondary-structure annotations from the supplementary benchmark of Szikszai et al. (2026).

The starting collection contained 89 RNA structures.

## Dataset selection

We examined the annotated base pairs and excluded 31 structures containing crossing base pairs (pseudoknots).

The resulting dataset contained:

- 58 RNA structures
- 3,426 nucleotide positions
- 35 supplementary-section groups

This dataset was used for the initial comparison of ViennaRNA, LinearPartition-V and EternaFold.

## Reference files

The structural reference package contains:

- `Project7_ExperimentA_58_sequences.fasta` — RNA sequences
- `Project7_ExperimentA_58_reference_labels.csv` — nucleotide-level paired/unpaired annotations
- `Project7_ExperimentA_58_manifest.csv` — RNA identifiers and grouping information

The original CT structure files are retained separately as part of the source-data archive.

## Limitations

The supplementary-section groups are not confirmed independent evolutionary RNA families. The screening therefore does not establish complete family independence or exclude all possible overlap with predictor training data.

The 50-RNA and 45-RNA sensitivity analyses are documented separately in the repository.

## Reproducibility

The three-predictor workflow is available in the `Experiment_A/analysis/` directory. The reference input files and saved output tables support verification of the structural calibration analysis.

A complete fresh-run verification of the original prediction workflow remains to be performed.
