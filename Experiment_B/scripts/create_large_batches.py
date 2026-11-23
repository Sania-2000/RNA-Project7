from pathlib import Path
import pandas as pd

ROOT = Path("large_scale")
BATCHES = ROOT / "batches"
BATCHES.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    ROOT / "prediction_inventory_20739.csv",
    dtype={"sequence_id": str}
)

assert len(df) == 20739
assert df.sequence_id.is_unique

df = df.sort_values("sequence_id").reset_index(drop=True)

BATCH_SIZE = 500
summary = []

for start in range(0, len(df), BATCH_SIZE):
    batch_number = start // BATCH_SIZE + 1
    batch = df.iloc[start:start + BATCH_SIZE]

    filename = BATCHES / f"batch_{batch_number:03d}.csv"

    batch.to_csv(filename, index=False)

    summary.append({
        "batch": batch_number,
        "n_sequences": len(batch),
        "filename": str(filename)
    })

pd.DataFrame(summary).to_csv(
    ROOT / "batch_manifest.csv",
    index=False
)

print("Total RNAs:", len(df))
print("Batches created:", len(summary))
print("Batch size:", BATCH_SIZE)
print("Final batch size:", summary[-1]["n_sequences"])
print("\nSaved: large_scale/batch_manifest.csv")
