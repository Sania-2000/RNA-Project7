from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import time

import pandas as pd

PROJECT = Path.home() / "RNA_Project7"

EXE = (
    PROJECT / "tools/LinearPartition"
    / "linearpartition"
)

INPUT = (
    PROJECT / "expanded/expanded_139_sequences.csv"
)

OUTPUT = PROJECT / "expanded/predictions"
OUTPUT.mkdir(exist_ok=True)

if not EXE.exists():
    raise FileNotFoundError(EXE)

df = pd.read_csv(INPUT)

unique_rnas = df[
    ["sequence_id", "sequence"]
].drop_duplicates()

pair_rows = []
nucleotide_rows = []
summary_rows = []

for record in unique_rnas.itertuples(index=False):

    rna_id = record.sequence_id
    sequence = record.sequence.upper()
    n = len(sequence)

    start = time.perf_counter()

    # Run one RNA inside a temporary folder
    with TemporaryDirectory() as temp:
        command = [
            str(EXE),
            "-V",
            "--prefix", "pairs"
        ]

        run = subprocess.run(
            command,
            input=sequence + "\n",
            text=True,
            capture_output=True,
            cwd=temp,
            check=True
        )

        # One input sequence produces pairs_1
        pair_file = Path(temp) / "pairs_1"

        if not pair_file.exists():
            raise FileNotFoundError(
                f"No probability file for {rna_id}"
            )

        paired = [0.0] * (n + 1)

        # LinearPartition output is sparse.
        # Unlisted pairs have probability zero.
        observed_pairs = {}

        for line in pair_file.read_text().splitlines():
            parts = line.split()

            if len(parts) != 3:
                continue

            i = int(parts[0])
            j = int(parts[1])
            p = float(parts[2])

            if not (1 <= i < j <= n):
                raise ValueError("Invalid pair index")

            if not (0 <= p <= 1):
                raise ValueError("Invalid probability")

            if (i, j) in observed_pairs:
                raise ValueError("Duplicate pair")

            observed_pairs[(i, j)] = p

        # Save the complete candidate-pair matrix
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):

                p = observed_pairs.get((i, j), 0.0)

                pair_rows.append({
                    "sequence_id": rna_id,
                    "predictor": "LinearPartition-V",
                    "i": i,
                    "j": j,
                    "probability": p
                })

                paired[i] += p
                paired[j] += p

    runtime = time.perf_counter() - start

    for i in range(1, n + 1):
        p_paired = paired[i]

        if not -1e-3 <= p_paired <= 1.001:
            raise ValueError(
                f"Invalid paired probability "
                f"for {rna_id}, position {i}: p_paired={p_paired:.12f}"
            )

        nucleotide_rows.append({
            "sequence_id": rna_id,
            "predictor": "LinearPartition-V",
            "position": i,
            "nucleotide": sequence[i - 1],
            "p_paired": p_paired,
            "p_unpaired": 1 - p_paired
        })

    summary_rows.append({
        "sequence_id": rna_id,
        "length": n,
        "runtime_seconds": runtime,
        "reported_pairs": len(observed_pairs)
    })

    print(f"Processed {rna_id}: {n} nt")

pd.DataFrame(pair_rows).to_csv(
    OUTPUT / "linearpartition_pairs.csv",
    index=False
)

pd.DataFrame(nucleotide_rows).to_csv(
    OUTPUT / "linearpartition_unpaired.csv",
    index=False
)

pd.DataFrame(summary_rows).to_csv(
    OUTPUT / "linearpartition_summary.csv",
    index=False
)

print("\nSUCCESS: 10 RNA predictions completed!")
print("All CSV files saved.")
