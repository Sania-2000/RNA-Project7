from pathlib import Path
import pandas as pd

inventory_file = Path(
    "results/tables/ribonanza_matched_inventory.csv"
)

output_dir = Path("data/processed")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(inventory_file)

# Select one row per RNA sequence
unique = df.drop_duplicates("sequence_id").copy()

# Proposed balanced pilot across RNA lengths
targets = {
    177: 50,
    170: 50,
    115: 35,
    206: 35,
    155: 30
}

selected = []

for length, target in targets.items():

    group = unique[
        unique["length"] == length
    ]

    n = min(target, len(group))

    sample = group.sample(
        n=n,
        random_state=42
    )

    selected.append(sample)

selected = pd.concat(
    selected, ignore_index=True
)

# Preserve selected RNA identifiers
selected_ids = set(selected["sequence_id"])

expanded = df[
    df["sequence_id"].isin(selected_ids)
].copy()

selected.to_csv(
    output_dir / "expanded_rna_inventory.csv",
    index=False
)

expanded.to_csv(
    output_dir / "expanded_ribonanza_metadata.csv",
    index=False
)

print("\n=== EXPANDED DATASET ===")
print("Unique RNAs:", selected.sequence_id.nunique())
print("Metadata profiles:", len(expanded))

print("\nLength distribution:")
print(selected["length"].value_counts().sort_index())

print("\nSource distribution:")
print(selected["dataset_name"].value_counts())

print("\nSaved dataset inventory and metadata.")
