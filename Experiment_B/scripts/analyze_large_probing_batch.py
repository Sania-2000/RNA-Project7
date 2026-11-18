from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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

ids = set(batch.sequence_id)
sequences = batch.set_index("sequence_id").sequence.to_dict()

react = pd.read_csv(
    ROOT / "reactivities_20739.csv.gz",
    dtype={"sequence_id": str}
)
react = react[react.sequence_id.isin(ids)].copy()

cols = sorted(
    [c for c in react if c.startswith("reactivity_")
     and c[len("reactivity_"):].isdigit()],
    key=lambda c: int(c.split("_")[-1])
)

assert len(react) == 2 * len(batch)
assert len(cols) == 206
assert not react.duplicated(
    ["sequence_id", "experiment_type"]
).any()

profiles = {}

for row in react.itertuples(index=False):
    seq_id = str(row.sequence_id)
    probe = row.experiment_type

    assert row.sequence == sequences[seq_id]
    assert probe in ("DMS_MaP", "2A3_MaP")

for _, row in react.iterrows():
    seq_id = row["sequence_id"]
    probe = row["experiment_type"]
    sequence = sequences[seq_id]

    values = pd.to_numeric(
        row[cols[:len(sequence)]],
        errors="coerce"
    ).to_numpy(dtype=float)

    profiles[(seq_id, probe)] = values

records = []

for predictor in ("vienna", "linearpartition", "eternafold"):
    pred = pd.read_csv(
        ROOT / f"full_predictions/{predictor}_batch_{args.batch:03d}.csv",
        dtype={"sequence_id": str}
    )

    assert pred.sequence_id.nunique() == len(batch)

    for seq_id, group in pred.groupby("sequence_id"):
        group = group.sort_values("position")
        sequence = sequences[seq_id]

        assert group.position.tolist() == list(
            range(1, len(sequence) + 1)
        )
        assert "".join(group.nucleotide) == sequence

        p = group.p_unpaired.to_numpy(dtype=float)

        for probe in ("DMS_MaP", "2A3_MaP"):
            values = profiles[(seq_id, probe)]

            eligible = np.array(
                [base in "AC" for base in sequence]
            ) if probe == "DMS_MaP" else np.ones(
                len(sequence), dtype=bool
            )

            mask = eligible & np.isfinite(values) & np.isfinite(p)
            n_valid = int(mask.sum())

            rho = np.nan
            if n_valid >= 10:
                if (np.ptp(values[mask]) > 0 and
                        np.ptp(p[mask]) > 0):
                    rho = float(spearmanr(
                        p[mask], values[mask]
                    ).statistic)

            records.append({
                "sequence_id": seq_id,
                "predictor": predictor,
                "probe": probe,
                "n_valid": n_valid,
                "spearman_rho": rho
            })

results = pd.DataFrame(records)

assert len(results) == 6 * len(batch)

path = OUT / f"probing_batch_{args.batch:03d}.csv"
results.to_csv(path, index=False)

print("RNAs:", len(batch))
print("Correlation records:", len(results))
print("Undefined correlations:", results.spearman_rho.isna().sum())

print("\nMedian Spearman correlations:")
print(
    results.groupby(["probe", "predictor"])
    .spearman_rho.median().unstack()
    .round(4).to_string()
)

print("\nSaved:", path)
print("SUCCESS: First 500 RNA probing correlations calculated!")
