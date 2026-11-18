from pathlib import Path
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--batch", type=int, required=True)
args = parser.parse_args()

assert 1 <= args.batch <= 42

ROOT = Path("large_scale")
OUT = ROOT / "results/tables"
OUT.mkdir(parents=True, exist_ok=True)

batch = pd.read_csv(
    ROOT / f"batches/batch_{args.batch:03d}.csv",
    dtype={"sequence_id": str}
)

assert batch.sequence_id.is_unique
ids = set(batch.sequence_id)
sequences = batch.set_index("sequence_id").sequence.to_dict()

react = pd.read_csv(
    ROOT / "reactivities_20739.csv.gz",
    dtype={"sequence_id": str}
)
react = react[react.sequence_id.isin(ids)]

assert len(react) == 2 * len(batch)

columns = sorted(
    [
        c for c in react.columns
        if c.startswith("reactivity_")
        and c.split("_")[-1].isdigit()
    ],
    key=lambda c: int(c.split("_")[-1])
)

profiles = {}

for row in react.itertuples(index=False):
    seq_id = row.sequence_id
    sequence = sequences[seq_id]
    probe = row.experiment_type

    assert row.sequence == sequence

    values = react.loc[
        (react.sequence_id == seq_id) &
        (react.experiment_type == probe),
        columns[:len(sequence)]
    ].iloc[0].to_numpy(dtype=float)

    profiles[(seq_id, probe)] = values

records = []

for predictor in ("vienna", "linearpartition", "eternafold"):
    prediction = pd.read_csv(
        ROOT / (
            f"full_predictions/"
            f"{predictor}_batch_{args.batch:03d}.csv"
        ),
        dtype={"sequence_id": str}
    )

    assert prediction.sequence_id.nunique() == len(batch)

    for seq_id, group in prediction.groupby("sequence_id"):
        group = group.sort_values("position")
        sequence = sequences[seq_id]

        assert group.position.tolist() == list(
            range(1, len(sequence) + 1)
        )
        assert "".join(group.nucleotide) == sequence

        x_all = group.p_unpaired.to_numpy(dtype=float)

        for probe in ("DMS_MaP", "2A3_MaP"):
            y_all = profiles[(seq_id, probe)]

            eligible = (
                np.array([base in "AC" for base in sequence])
                if probe == "DMS_MaP"
                else np.ones(len(sequence), dtype=bool)
            )

            mask = (
                eligible
                & np.isfinite(x_all)
                & np.isfinite(y_all)
            )

            x = x_all[mask]
            y = y_all[mask]

            assert len(x) >= 10

            records.append({
                "sequence_id": seq_id,
                "probe": probe,
                "predictor": predictor,
                "n": len(x),
                "sum_x": x.sum(),
                "sum_y": y.sum(),
                "sum_xx": np.dot(x, x),
                "sum_xy": np.dot(x, y),
                "sum_yy": np.dot(y, y)
            })

stats = pd.DataFrame(records)

assert len(stats) == 6 * len(batch)

metadata = pd.read_csv(
    ROOT / "prediction_inventory_20739.csv",
    dtype={"sequence_id": str}
)

stats = stats.merge(
    metadata[
        ["sequence_id", "cluster_id", "source_group", "length"]
    ],
    on="sequence_id",
    validate="many_to_one"
)

assert stats.cluster_id.notna().all()

output = OUT / f"regression_statistics_batch_{args.batch:03d}.csv"
stats.to_csv(output, index=False)

print("RNAs:", stats.sequence_id.nunique())
print("Regression records:", len(stats))
print("Total nucleotide observations:", stats.n.sum())
print("Saved:", output)
print("SUCCESS: Regression batch prepared!")
