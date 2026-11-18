from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path("large_scale/results/tables")
INPUT = ROOT / "full_regression_statistics_20739.csv"
OUTPUT = ROOT / "full_regression_results_20739.csv"
SPLIT_OUTPUT = ROOT / "regression_cluster_split_20739.csv"

df = pd.read_csv(
    INPUT,
    dtype={"sequence_id": str, "cluster_id": str}
)

assert len(df) == 124434
assert df.sequence_id.nunique() == 20739

# Every RNA must belong to one cluster.
rnas = df[
    ["sequence_id", "cluster_id"]
].drop_duplicates()

assert rnas.sequence_id.is_unique
assert len(rnas) == 20739

# Split complete RNA clusters, not nucleotides.
splitter = GroupShuffleSplit(
    n_splits=1,
    test_size=0.20,
    random_state=108
)

train_index, test_index = next(
    splitter.split(
        rnas,
        groups=rnas.cluster_id
    )
)

rnas["split"] = ""

rnas.loc[rnas.index[train_index], "split"] = "train"
rnas.loc[rnas.index[test_index], "split"] = "test"

assert set(rnas.split) == {"train", "test"}

train_clusters = set(
    rnas.loc[rnas.split == "train", "cluster_id"]
)

test_clusters = set(
    rnas.loc[rnas.split == "test", "cluster_id"]
)

assert train_clusters.isdisjoint(test_clusters)

rnas.to_csv(SPLIT_OUTPUT, index=False)

df = df.merge(
    rnas[["sequence_id", "split"]],
    on="sequence_id",
    validate="many_to_one"
)

def totals(frame):
    cols = [
        "n", "sum_x", "sum_y",
        "sum_xx", "sum_xy", "sum_yy"
    ]
    return frame[cols].sum()

def fit_model(t):
    n = float(t["n"])
    sx = float(t["sum_x"])
    sy = float(t["sum_y"])
    sxx = float(t["sum_xx"])
    sxy = float(t["sum_xy"])

    mean_x = sx / n
    mean_y = sy / n

    var_x = sxx - sx * sx / n
    cov_xy = sxy - sx * sy / n

    assert var_x > 0, "Predictor has no variation"

    slope = cov_xy / var_x
    intercept = mean_y - slope * mean_x

    return intercept, slope

def evaluate(t, intercept, slope):
    n = float(t["n"])
    sx = float(t["sum_x"])
    sy = float(t["sum_y"])
    sxx = float(t["sum_xx"])
    sxy = float(t["sum_xy"])
    syy = float(t["sum_yy"])

    # Sum of squared prediction errors on the test set.
    sse = (
        syy
        - 2 * intercept * sy
        - 2 * slope * sxy
        + n * intercept**2
        + 2 * intercept * slope * sx
        + slope**2 * sxx
    )

    # Total sum of squares relative to test-set mean.
    sst = syy - sy**2 / n

    assert sse >= -1e-6
    assert sst > 0

    sse = max(0.0, sse)

    rmse = np.sqrt(sse / n)
    r2 = 1 - sse / sst

    return r2, rmse

records = []

for probe in ("DMS_MaP", "2A3_MaP"):
    for predictor in (
        "vienna", "linearpartition", "eternafold"
    ):
        part = df[
            (df.probe == probe) &
            (df.predictor == predictor)
        ]

        train = part[part.split == "train"]
        test = part[part.split == "test"]

        assert len(train) > 0
        assert len(test) > 0

        train_totals = totals(train)
        test_totals = totals(test)

        intercept, slope = fit_model(train_totals)

        test_r2, test_rmse = evaluate(
            test_totals,
            intercept,
            slope
        )

        records.append({
            "probe": probe,
            "predictor": predictor,
            "train_rnas": train.sequence_id.nunique(),
            "test_rnas": test.sequence_id.nunique(),
            "train_clusters": train.cluster_id.nunique(),
            "test_clusters": test.cluster_id.nunique(),
            "train_positions": int(train_totals["n"]),
            "test_positions": int(test_totals["n"]),
            "intercept": intercept,
            "slope": slope,
            "test_r2": test_r2,
            "test_rmse": test_rmse
        })

results = pd.DataFrame(records)

assert len(results) == 6
assert np.isfinite(
    results[["intercept", "slope", "test_r2", "test_rmse"]]
).all().all()

results.to_csv(OUTPUT, index=False)

print("Total RNA sequences:", len(rnas))
print("Total clusters:", rnas.cluster_id.nunique())
print("Training RNAs:", (rnas.split == "train").sum())
print("Testing RNAs:", (rnas.split == "test").sum())
print("Shared train/test clusters:", len(
    train_clusters & test_clusters
))

print("\nFULL-DATASET REGRESSION RESULTS:")
print(results.round(4).to_string(index=False))

print("\nSaved:", OUTPUT)
print("Saved:", SPLIT_OUTPUT)

print("\nSUCCESS: Full regression completed!")
