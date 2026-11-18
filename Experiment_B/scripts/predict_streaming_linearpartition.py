from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import argparse
import csv

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

PROJECT = Path.cwd()
EXE = PROJECT / "tools/LinearPartition/linearpartition"

if not EXE.exists():
    raise FileNotFoundError(EXE)

df = pd.read_csv(args.input, dtype={"sequence_id": str})
assert df.sequence_id.is_unique

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "sequence_id", "predictor", "position",
        "nucleotide", "p_paired",
        "p_unpaired_raw", "p_unpaired",
        "was_clipped"
    ])

    for count, row in enumerate(df.itertuples(index=False), 1):
        sequence = row.sequence.upper()
        n = len(sequence)

        with TemporaryDirectory() as temp:
            subprocess.run(
                [str(EXE), "-V", "--prefix", "pairs"],
                input=sequence + "\n",
                text=True,
                capture_output=True,
                cwd=temp,
                check=True
            )

            pair_file = Path(temp) / "pairs_1"
            if not pair_file.exists():
                raise FileNotFoundError(pair_file)

            paired = [0.0] * (n + 1)
            seen = set()

            for line in pair_file.read_text().splitlines():
                parts = line.split()
                if len(parts) != 3:
                    continue

                i, j = int(parts[0]), int(parts[1])
                p = float(parts[2])

                assert 1 <= i < j <= n
                assert 0 <= p <= 1
                assert (i, j) not in seen
                seen.add((i, j))

                paired[i] += p
                paired[j] += p

        # Compute and validate the entire RNA before writing.
        records = []
        for i in range(1, n + 1):
            raw = 1.0 - paired[i]

            if not (-0.001 <= raw <= 1.001):
                raise ValueError(
                    f"Excessive numerical violation: "
                    f"{row.sequence_id}, position {i}: {raw}"
                )

            corrected = min(1.0, max(0.0, raw))

            records.append([
                row.sequence_id,
                "LinearPartition-V",
                i,
                sequence[i - 1],
                paired[i],
                raw,
                corrected,
                raw != corrected
            ])

        writer.writerows(records)

        if count % 25 == 0 or count == len(df):
            print(
                f"Completed {count}/{len(df)}",
                flush=True
            )

print("\nSaved:", output)
print("RNAs processed:", len(df))
