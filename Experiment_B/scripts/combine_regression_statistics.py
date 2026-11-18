from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("large_scale/results/tables")

files = [
    ROOT / f"regression_statistics_batch_{i:03d}.csv"
    for i in range(1, 43)
]

assert all(p.is_file() for p in files), "Missing batch files"

df = pd.concat(
    [
        pd.read_csv(p, dtype={"sequence_id": str})
        for p in files
    ],
    ignore_index=True
)

assert len(df) == 124434
assert df.sequence_id.nunique() == 20739

assert not df.duplicated(
    ["sequence_id", "probe", "predictor"]
).any()

assert df.cluster_id.notna().all()
assert df.source_group.notna().all()

assert df.cluster_id.nunique() == 16520

assert set(df.predictor) == {
    "vienna", "linearpartition", "eternafold"
}

assert set(df.probe) == {"DMS_MaP", "2A3_MaP"}

numeric = [
    "n", "sum_x", "sum_y",
    "sum_xx", "sum_xy", "sum_yy"
]

assert np.isfinite(df[numeric].to_numpy()).all()
assert df.n.ge(10).all()

counts = df.pivot(
    index=["sequence_id", "probe"],
    columns="predictor",
    values="n"
)

assert len(counts) == 41478
assert counts.notna().all().all()
assert (counts.nunique(axis=1) == 1).all()

output = ROOT / "full_regression_statistics_20739.csv"
df.to_csv(output, index=False)

print("RNA sequences:", df.sequence_id.nunique())
print("Regression records:", len(df))
print("Sequence clusters:", df.cluster_id.nunique())
print("Experimental sources:", df.source_group.nunique())

print("\nPositions per probe and predictor:")
print(
    df.groupby(["probe", "predictor"])["n"]
    .sum()
    .to_string()
)

print("\nSaved:", output)
print("SUCCESS: Full regression statistics combined!")
