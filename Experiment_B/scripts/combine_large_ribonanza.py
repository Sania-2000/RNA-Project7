from pathlib import Path
import pandas as pd

OUT = Path("large_scale")

first = pd.read_csv(
    OUT / "candidate_20000_qc_passed.csv",
    dtype={"sequence_id": str}
)

second = pd.read_csv(
    OUT / "additional_10000_qc_passed.csv",
    dtype={"sequence_id": str}
)

combined = pd.concat(
    [first, second],
    ignore_index=True
)

print("First batch:", len(first))
print("Second batch:", len(second))

assert combined.sequence_id.is_unique, \
    "Duplicate sequence IDs found"

assert combined.sequence.is_unique, \
    "Duplicate RNA sequences found"

assert len(combined) == 20739, \
    "Unexpected combined count"

combined = combined.sort_values(
    "sequence_id"
).reset_index(drop=True)

combined.to_csv(
    OUT / "combined_20739_qc_passed.csv",
    index=False
)

print("\nCombined QC-passing RNAs:", len(combined))
print("Unique sequences:", combined.sequence.nunique())

print("\nRNA length distribution:")
print(combined.length.value_counts().sort_index())

print("\nTop 10 source groups:")
print(combined.source_group.value_counts().head(10))

print("\nSaved: large_scale/combined_20739_qc_passed.csv")
