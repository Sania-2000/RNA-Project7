from pathlib import Path
import pandas as pd

SOURCE = Path("data/raw/train_data_QUICK_START.csv")
OUTPUT = Path("results/tables")
OUTPUT.mkdir(parents=True, exist_ok=True)

print("Reading Ribonanza metadata...")

columns = [
    "sequence_id",
    "sequence",
    "experiment_type",
    "dataset_name"
]

df = pd.read_csv(
    SOURCE,
    usecols=columns,
    dtype=str
)

df["sequence"] = df["sequence"].str.upper()
df["length"] = df["sequence"].str.len()

print("\n=== DATASET OVERVIEW ===")
print("Total profiles:", len(df))
print("Unique RNA IDs:", df.sequence_id.nunique())
print("Unique sequences:", df.sequence.nunique())

print("\n=== PROBE DISTRIBUTION ===")
print(df.experiment_type.value_counts())

print("\n=== RNA LENGTHS ===")
print(df["length"].describe())

print("\n=== DATASET SOURCES ===")
print(df.dataset_name.value_counts().head(20))

# Check consistency of sequence IDs
id_counts = df.groupby("sequence_id")["sequence"].nunique()
print("\nIDs mapping to multiple sequences:",
      (id_counts > 1).sum())

# Check matched DMS and 2A3 profiles
required = {"DMS_MaP", "2A3_MaP"}

probe_sets = df.groupby("sequence_id")[
    "experiment_type"
].agg(lambda x: set(x))

matched_ids = probe_sets[
    probe_sets.apply(lambda x: required.issubset(x))
].index

matched = df[df.sequence_id.isin(matched_ids)].copy()

print("\n=== MATCHED DATA ===")
print("RNA IDs with both probes:", len(matched_ids))
print("Matched profiles:", len(matched))
print("Unique matched sequences:",
      matched.sequence.nunique())

print("\n=== LENGTH DISTRIBUTION ===")
lengths = (
    matched.drop_duplicates("sequence_id")
    .groupby("length")
    .size()
    .sort_values(ascending=False)
)

print(lengths.head(15))

# Save reusable inspection tables
lengths.rename("n_RNAs").to_csv(
    OUTPUT / "ribonanza_length_distribution.csv"
)

matched[
    ["sequence_id", "sequence",
     "length", "dataset_name", "experiment_type"]
].to_csv(
    OUTPUT / "ribonanza_matched_inventory.csv",
    index=False
)

print("\nSUCCESS: Dataset inspection completed!")
print("Files saved in results/tables/")
