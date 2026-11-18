from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import csv
import subprocess

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

ROOT = Path.cwd()
EXE = ROOT / "tools/EternaFold/src/contrafold"
PARAMS = ROOT / "tools/EternaFold/parameters/EternaFoldParams.v1"

for path in (EXE, PARAMS):
    if not path.exists():
        raise FileNotFoundError(path)

df = pd.read_csv(
    args.input,
    dtype={"sequence_id": str}
)

assert df.sequence_id.is_unique

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "sequence_id",
        "predictor",
        "position",
        "nucleotide",
        "p_paired",
        "p_unpaired_raw",
        "p_unpaired",
        "was_clipped"
    ])

    for count, row in enumerate(df.itertuples(index=False), 1):
        sequence = row.sequence.upper()
        n = len(sequence)

        with TemporaryDirectory() as temp:
            temp = Path(temp)
            fasta = temp / "input.fa"
            fasta.write_text(
                f">{row.sequence_id}\n{sequence}\n"
            )

            command = [
                str(EXE), "predict", str(fasta),
                "--params", str(PARAMS),
                "--posteriors", "0.00001", "bps.txt"
            ]

            subprocess.run(
                command,
                cwd=temp,
                capture_output=True,
                text=True,
                check=True
            )

            posterior = temp / "bps.txt"

            if not posterior.exists():
                raise FileNotFoundError(posterior)

            paired = [0.0] * (n + 1)
            seen = set()

            for line in posterior.read_text().splitlines():
                parts = line.split()

                if len(parts) < 2:
                    continue

                i = int(parts[0])

                for item in parts[2:]:
                    j_text, p_text = item.split(":")
                    j = int(j_text)
                    p = float(p_text)

                    if not 1 <= i < j <= n:
                        raise ValueError("Invalid pair index")

                    if not 0 <= p <= 1:
                        raise ValueError("Invalid pair probability")

                    if (i, j) in seen:
                        raise ValueError("Duplicate pair")

                    seen.add((i, j))
                    paired[i] += p
                    paired[j] += p

        # Validate all nucleotides before saving this RNA.
        records = []

        for i in range(1, n + 1):
            p_unpaired_raw = 1.0 - paired[i]

            if not -0.001 <= p_unpaired_raw <= 1.001:
                raise ValueError(
                    f"Excessive numerical violation: "
                    f"{row.sequence_id}, position {i}, "
                    f"value={p_unpaired_raw}"
                )

            p_unpaired = min(
                1.0, max(0.0, p_unpaired_raw)
            )

            records.append([
                row.sequence_id,
                "EternaFold",
                i,
                sequence[i - 1],
                paired[i],
                p_unpaired_raw,
                p_unpaired,
                p_unpaired_raw != p_unpaired
            ])

        writer.writerows(records)

        if count % 25 == 0 or count == len(df):
            print(
                f"Completed {count}/{len(df)}",
                flush=True
            )

print("\nSaved:", output)
print("RNAs processed:", len(df))
