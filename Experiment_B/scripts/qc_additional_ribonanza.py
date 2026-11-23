from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("large_scale")
RAW = Path("data/raw/train_data_QUICK_START.csv")
CANDIDATES = OUT / "additional_10000_inventory.csv"

sample = pd.read_csv(
    CANDIDATES,
    dtype={"sequence_id": str}
)

ids = set(sample.sequence_id)

MIN_VALID_POSITIONS = 10
MIN_VALID_FRACTION = 0.50

reactivity_cols = pd.read_csv(
    RAW, nrows=0
).filter(regex=r"^reactivity_\d+$").columns.tolist()

usecols = [
    "sequence_id",
    "experiment_type"
] + reactivity_cols

results = []

for chunk in pd.read_csv(
    RAW,
    usecols=usecols,
    dtype={"sequence_id": str},
    chunksize=2000
):
    chunk = chunk[chunk.sequence_id.isin(ids)]

    for row in chunk.itertuples(index=False, name=None):
        seq_id = row[0]
        probe = row[1]
        values = np.asarray(row[2:], dtype=float)

        seq = sample.loc[
            sample.sequence_id == seq_id,
            "sequence"
        ].iloc[0]

        if probe == "DMS_MaP":
            valid_bases = np.array(
                [base in "AC" for base in seq]
            )
        else:
            valid_bases = np.ones(len(seq), dtype=bool)

        measured = values[:len(seq)]
        usable = valid_bases & np.isfinite(measured)

        n_eligible = int(valid_bases.sum())
        n_valid = int(usable.sum())

        fraction = (
            n_valid / n_eligible
            if n_eligible else 0
        )

        passes = (
            n_valid >= MIN_VALID_POSITIONS
            and fraction >= MIN_VALID_FRACTION
        )

        results.append({
            "sequence_id": seq_id,
            "probe": probe,
            "n_eligible": n_eligible,
            "n_valid": n_valid,
            "valid_fraction": fraction,
            "passes_qc": passes
        })

qc = pd.DataFrame(results)

qc.to_csv(
    OUT / "additional_10000_reactivity_qc.csv",
    index=False
)

passed = (
    qc[qc.passes_qc]
    .groupby("sequence_id")["probe"]
    .nunique()
)

passed_ids = set(passed[passed == 2].index)

final = sample[
    sample.sequence_id.isin(passed_ids)
].copy()

final.to_csv(
    OUT / "additional_10000_qc_passed.csv",
    index=False
)

print("Candidate RNAs:", len(sample))
print("QC records:", len(qc))
print("RNAs passing both probes:", len(final))

print("\nPassing profiles by probe:")
print(qc.groupby("probe").passes_qc.sum())

print("\nSaved QC results in large_scale/")
