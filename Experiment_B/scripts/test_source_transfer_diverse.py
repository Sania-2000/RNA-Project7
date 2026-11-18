from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("large_scale/results/tables")

df = pd.read_csv(
    ROOT / "full_regression_statistics_20739.csv",
    dtype={"sequence_id": str, "cluster_id": str}
)

sources = [
    "15k",
    "DasLabBigLib_OneMil_RNAmake_designs"
]

assert all(
    (df.source_group == source).any()
    for source in sources
)

print("Testing source transfer between:")
for s in sources:
    print(" -", s)

def totals(frame):
    return frame[
        ["n", "sum_x", "sum_y",
         "sum_xx", "sum_xy", "sum_yy"]
    ].sum()

def fit(t):
    n = t["n"]
    slope = (
        t["sum_xy"] - t["sum_x"] * t["sum_y"] / n
    ) / (
        t["sum_xx"] - t["sum_x"] ** 2 / n
    )
    intercept = (t["sum_y"] - slope * t["sum_x"]) / n
    return intercept, slope

def evaluate(t, a, b):
    n = t["n"]
    sse = (
        t["sum_yy"]
        - 2 * a * t["sum_y"]
        - 2 * b * t["sum_xy"]
        + n * a**2
        + 2 * a * b * t["sum_x"]
        + b**2 * t["sum_xx"]
    )
    sst = t["sum_yy"] - t["sum_y"]**2 / n
    assert sst > 0
    return 1 - sse / sst, np.sqrt(max(0, sse) / n)

records = []

for train_source, test_source in [
    (sources[0], sources[1]),
    (sources[1], sources[0])
]:
    train_all = df[df.source_group == train_source]
    test_all = df[df.source_group == test_source]

    train_clusters = set(train_all.cluster_id)

    # Exclude target RNAs from any training cluster.
    test_all = test_all[
        ~test_all.cluster_id.isin(train_clusters)
    ]

    assert not (
        set(train_all.cluster_id) & set(test_all.cluster_id)
    )

    for probe in ("DMS_MaP", "2A3_MaP"):
        for predictor in (
            "vienna", "linearpartition", "eternafold"
        ):
            train = train_all[
                (train_all.probe == probe) &
                (train_all.predictor == predictor)
            ]
            test = test_all[
                (test_all.probe == probe) &
                (test_all.predictor == predictor)
            ]

            assert len(train) > 0 and len(test) > 0

            a, b = fit(totals(train))
            r2, rmse = evaluate(totals(test), a, b)

            records.append({
                "train_source": train_source,
                "test_source": test_source,
                "probe": probe,
                "predictor": predictor,
                "train_rnas": train.sequence_id.nunique(),
                "test_rnas": test.sequence_id.nunique(),
                "train_clusters": train.cluster_id.nunique(),
                "test_clusters": test.cluster_id.nunique(),
                "slope": b,
                "intercept": a,
                "transfer_r2": r2,
                "transfer_rmse": rmse
            })

results = pd.DataFrame(records)

assert len(results) == 12
assert np.isfinite(
    results[["transfer_r2", "transfer_rmse"]]
).all().all()

output = ROOT / "source_transfer_diverse.csv"
results.to_csv(output, index=False)

print("\nSOURCE-TO-SOURCE TRANSFER RESULTS:")
print(
    results[
        ["train_source", "test_source", "probe",
         "predictor", "test_rnas",
         "transfer_r2", "transfer_rmse"]
    ].round(4).to_string(index=False)
)

print("\nSaved:", output)
print("SUCCESS: Source transfer pilot completed!")
