from pathlib import Path
import pandas as pd

ROOT = Path("large_scale")
OUT = ROOT / "results/tables"
OUT.mkdir(parents=True, exist_ok=True)

files = [
    OUT / f"probing_batch_{i:03d}.csv"
    for i in range(1, 43)
]

for file in files:
    assert file.exists(), f"Missing file: {file}"

results = pd.concat(
    [
        pd.read_csv(file, dtype={"sequence_id": str})
        for file in files
    ],
    ignore_index=True
)

inventory = pd.read_csv(
    ROOT / "prediction_inventory_20739.csv",
    dtype={"sequence_id": str}
)

assert len(inventory) == 20739
assert inventory.sequence_id.is_unique

assert len(results) == 124434
assert results.sequence_id.nunique() == 20739

assert not results.duplicated(
    ["sequence_id", "predictor", "probe"]
).any()

assert results.spearman_rho.notna().all()
assert results.spearman_rho.between(-1, 1).all()

assert set(results.sequence_id) == set(inventory.sequence_id)

metadata = inventory[
    ["sequence_id", "length", "source_group", "cluster_id"]
]

combined = results.merge(
    metadata,
    on="sequence_id",
    how="left",
    validate="many_to_one"
)

assert combined.cluster_id.notna().all()
assert combined.source_group.notna().all()
assert combined.cluster_id.nunique() == 16520

combined_file = OUT / "large_probing_correlations_20739.csv"
combined.to_csv(combined_file, index=False)

summary = (
    combined.groupby(["probe", "predictor"])
    .agg(
        n_rnas=("sequence_id", "nunique"),
        median_rho=("spearman_rho", "median"),
        mean_rho=("spearman_rho", "mean")
    )
    .reset_index()
)

summary_file = OUT / "large_probing_summary_20739.csv"
summary.to_csv(summary_file, index=False)

print("RNA sequences:", combined.sequence_id.nunique())
print("Correlation records:", len(combined))
print("Sequence clusters:", combined.cluster_id.nunique())

print("\nMedian Spearman correlations:")
print(
    summary.pivot(
        index="predictor",
        columns="probe",
        values="median_rho"
    ).round(4).to_string()
)

print("\nSaved:", combined_file)
print("Saved:", summary_file)
print("\nSUCCESS: All 42 batches combined!")
