
import pandas as pd
from pathlib import Path

# Original experimental dataset
input_file = Path(
    "data/raw/train_data_QUICK_START.csv"
)

# Output files
output_file = Path(
    "data/processed/ribonanza_sample.csv"
)

summary_file = Path(
    "data/processed/ribonanza_sample_summary.csv"
)

output_file.parent.mkdir(
    parents=True, exist_ok=True
)

# Read the dataset
print("Loading Ribonanza dataset...")

df = pd.read_csv(input_file)

print("Total experimental profiles:", len(df))
print("Unique RNA sequences:", df["sequence_id"].nunique())

# Keep the two required experiment types
df = df[
    df["experiment_type"].isin(
        ["DMS_MaP", "2A3_MaP"]
    )
].copy()

# Keep only sequences with standard RNA nucleotides
df["sequence"] = df["sequence"].str.upper()

df = df[
    df["sequence"].str.fullmatch("[ACGU]+")
].copy()

# Only keep sequence IDs having both probes
probe_counts = (
    df.groupby("sequence_id")["experiment_type"]
    .nunique()
)

matched_ids = probe_counts[
    probe_counts == 2
].index

matched = df[
    df["sequence_id"].isin(matched_ids)
].copy()

# Select 10 matched RNAs reproducibly
selected_ids = (
    matched["sequence_id"]
    .drop_duplicates()
    .sample(
        n=min(10, len(matched_ids)),
        random_state=42
    )
)

sample = matched[
    matched["sequence_id"].isin(selected_ids)
].copy()

# Save the experimental profiles
sample.to_csv(output_file, index=False)

# Save a compact dataset summary
summary = sample[
    ["sequence_id", "sequence", "experiment_type"]
].copy()

summary["length"] = summary["sequence"].str.len()

summary.to_csv(summary_file, index=False)

print("\nSelected RNA sequences:", len(selected_ids))
print("Experimental profiles:", len(sample))

print("\nProbe distribution:")
print(sample["experiment_type"].value_counts())

print("\nSaved files:")
print(output_file)
print(summary_file)
