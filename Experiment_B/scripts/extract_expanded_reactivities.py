from pathlib import Path
import pandas as pd
import numpy as np

RAW = Path("data/raw/train_data_QUICK_START.csv")
SELECTED = Path("data/processed/expanded_rna_nonredundant.csv")
OUT = Path("data/processed")

OUT.mkdir(parents=True, exist_ok=True)

# Load selected RNA IDs and sequences
selected = pd.read_csv(SELECTED, dtype={"sequence_id": str})

selected_ids = set(selected["sequence_id"])
expected = dict(zip(selected["sequence_id"], selected["sequence"]))

print("Selected RNA sequences:", len(selected_ids))

profiles = []

# Read the large original CSV in manageable chunks
for chunk in pd.read_csv(
    RAW,
    chunksize=2000,
    dtype={"sequence_id": str},
    low_memory=False
):
    subset = chunk[
        chunk["sequence_id"].isin(selected_ids)
    ].copy()

    if not subset.empty:
        profiles.append(subset)

if not profiles:
    raise RuntimeError("No matching profiles found.")

df = pd.concat(profiles, ignore_index=True)

# Check that IDs correspond to correct sequences
for row in df.itertuples(index=False):
    if row.sequence.upper() != expected[row.sequence_id].upper():
        raise ValueError(
            f"Sequence mismatch: {row.sequence_id}"
        )

# Keep required probes
df = df[
    df["experiment_type"].isin(["DMS_MaP", "2A3_MaP"])
].copy()

# Check for duplicate profiles
if df.duplicated(
    ["sequence_id", "experiment_type"]
).any():
    raise ValueError("Duplicate RNA/probe profiles found.")

# Count usable nucleotide measurements
qc_rows = []

for _, row in df.iterrows():
    seq = row["sequence"].upper()
    probe = row["experiment_type"]

    valid = 0
    eligible = 0

    for pos, nt in enumerate(seq, start=1):
        if probe == "DMS_MaP" and nt not in "AC":
            continue

        eligible += 1

        col = f"reactivity_{pos:04d}"

        if col in df.columns:
            value = pd.to_numeric(row[col], errors="coerce")
            if np.isfinite(value):
                valid += 1

    qc_rows.append({
        "sequence_id": row["sequence_id"],
        "probe": probe,
        "length": len(seq),
        "eligible_positions": eligible,
        "valid_measurements": valid,
        "valid_fraction": valid / eligible if eligible else 0
    })

qc = pd.DataFrame(qc_rows)

# Initial minimum coverage rule for the expanded analysis
qc["passes_qc"] = (
    (qc["valid_measurements"] >= 10) &
    (qc["valid_fraction"] >= 0.5)
)

# Retain RNAs that pass for BOTH probes
passing_ids = (
    qc[qc["passes_qc"]]
    .groupby("sequence_id")["probe"]
    .nunique()
)

passing_ids = set(
    passing_ids[passing_ids == 2].index
)

final = df[
    df["sequence_id"].isin(passing_ids)
].copy()

qc.to_csv(OUT / "expanded_reactivity_qc.csv", index=False)
final.to_csv(OUT / "expanded_reactivities.csv", index=False)

print("\n=== EXTRACTION RESULTS ===")
print("Profiles extracted:", len(df))
print("RNA IDs extracted:", df.sequence_id.nunique())

print("\n=== QUALITY FILTERING ===")
print("RNAs passing both probes:", len(passing_ids))
print("Final profiles:", len(final))

print("\nProbe distribution:")
print(final.experiment_type.value_counts())

print("\nSaved:")
print(OUT / "expanded_reactivity_qc.csv")
print(OUT / "expanded_reactivities.csv")
