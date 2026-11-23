from pathlib import Path
import sys
import pandas as pd

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("large_scale")

combined = pd.read_csv(
    root / "combined_20739_qc_passed.csv",
    dtype=str
)

clusters = pd.read_csv(
    root / "rna_cluster_membership.csv",
    dtype=str
)

inventory = combined.merge(
    clusters[
        ["sequence_id", "cluster_id",
         "is_representative", "cluster_size"]
    ],
    on="sequence_id",
    how="left",
    validate="one_to_one",
    indicator=True
)

assert inventory["_merge"].eq("both").all(), "Missing cluster assignments"

inventory = (
    inventory.drop(columns="_merge")
    .sort_values("sequence_id")
    .reset_index(drop=True)
)

assert len(inventory) == 20739
assert inventory["sequence_id"].is_unique

output = root / "prediction_inventory_20739.csv"

if output.exists():
    original = pd.read_csv(output, dtype=str)
    original = original.sort_values("sequence_id").reset_index(drop=True)

    assert list(inventory.columns) == list(original.columns)
    assert inventory.equals(original), "Inventory values differ"
    print("PASS: reconstructed inventory matches existing file exactly")
else:
    inventory.to_csv(output, index=False)
    print("Created:", output)

print("RNA sequences:", len(inventory))
