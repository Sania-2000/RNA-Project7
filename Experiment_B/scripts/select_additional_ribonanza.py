from pathlib import Path
import pandas as pd

OUT = Path("large_scale")

remaining = pd.read_csv(
    OUT / "remaining_rna_inventory.csv",
    dtype={"sequence_id": str}
)

TARGET = 10000
SEED = 43

# Draw candidates proportionally across source groups.
sample = remaining.sample(
    n=TARGET,
    random_state=SEED
).copy()

assert len(sample) == TARGET
assert sample.sequence_id.nunique() == TARGET
assert sample.sequence.nunique() == TARGET

sample.to_csv(
    OUT / "additional_10000_inventory.csv",
    index=False
)

print("Additional candidates:", len(sample))
print("Unique sequences:", sample.sequence.nunique())

print("\nLength distribution:")
print(sample.length.value_counts().sort_index())

print("\nTop 10 source groups:")
print(sample.source_group.value_counts().head(10))

print("\nSaved additional candidate inventory.")
