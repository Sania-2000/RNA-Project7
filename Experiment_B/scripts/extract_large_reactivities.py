from pathlib import Path
import pandas as pd

ROOT = Path("large_scale")
RAW = Path("data/raw/train_data_QUICK_START.csv")
INVENTORY = ROOT / "prediction_inventory_20739.csv"
OUTPUT = ROOT / "reactivities_20739.csv.gz"
TEMP = ROOT / "reactivities_20739.partial.csv.gz"

inventory = pd.read_csv(
    INVENTORY,
    dtype={"sequence_id": str}
)

assert len(inventory) == 20739
assert inventory.sequence_id.is_unique

sequences = inventory.set_index("sequence_id")["sequence"].to_dict()
wanted = set(sequences)

header = pd.read_csv(RAW, nrows=0).columns.tolist()
reactivity_cols = [
    col for col in header
    if col.startswith("reactivity_")
    and col[len("reactivity_"):].isdigit()
]

usecols = [
    "sequence_id",
    "sequence",
    "experiment_type",
    "dataset_name",
] + reactivity_cols

seen = set()
counts = {"DMS_MaP": 0, "2A3_MaP": 0}

with TEMP.open("wb") as binary:
    pass

import gzip

with gzip.open(TEMP, "wt", newline="") as output:
    first = True

    for chunk in pd.read_csv(
        RAW,
        usecols=usecols,
        dtype={"sequence_id": str},
        chunksize=2000
    ):
        part = chunk[chunk.sequence_id.isin(wanted)].copy()

        if part.empty:
            continue

        for row in part[
            ["sequence_id", "sequence", "experiment_type"]
        ].itertuples(index=False, name=None):
            seq_id, sequence, probe = row

            assert sequence == sequences[seq_id]
            assert probe in counts

            key = (seq_id, probe)
            assert key not in seen, f"Duplicate profile: {key}"

            seen.add(key)
            counts[probe] += 1

        part.to_csv(
            output,
            index=False,
            header=first
        )
        first = False

assert counts["DMS_MaP"] == 20739
assert counts["2A3_MaP"] == 20739
assert len(seen) == 41478

TEMP.replace(OUTPUT)

print("RNAs:", len(wanted))
print("DMS profiles:", counts["DMS_MaP"])
print("2A3 profiles:", counts["2A3_MaP"])
print("Total profiles:", len(seen))
print("Reactivity columns:", len(reactivity_cols))
print("Saved:", OUTPUT)
print("\nSUCCESS: Large-scale reactivities extracted!")
