from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import time
import pandas as pd

PROJECT = Path.home() / "RNA_Project7"
ROOT = PROJECT / "tools/EternaFold"

EXE = ROOT / "src/contrafold"
PARAMS = ROOT / "parameters/EternaFoldParams.v1"

INPUT = PROJECT / "expanded/expanded_139_sequences.csv"
OUTPUT = PROJECT / "expanded/predictions"
OUTPUT.mkdir(exist_ok=True)

for path in [EXE, PARAMS, INPUT]:
    if not path.exists():
        raise FileNotFoundError(path)

df = pd.read_csv(INPUT)

unique_rnas = df[
    ["sequence_id", "sequence"]
].drop_duplicates()

pair_rows = []
nucleotide_rows = []
summary_rows = []

for record in unique_rnas.itertuples(index=False):

    rna_id = str(record.sequence_id)
    sequence = str(record.sequence).upper()
    n = len(sequence)

    start = time.perf_counter()

    with TemporaryDirectory() as temp:
        temp = Path(temp)

        fasta = temp / "input.fa"
        fasta.write_text(
            f">{rna_id}\n{sequence}\n"
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

        observed_pairs = {}

        for line in posterior.read_text().splitlines():
            parts = line.split()

            if len(parts) < 2:
                continue

            i = int(parts[0])

            for item in parts[2:]:
                j_text, p_text = item.split(":")
                j = int(j_text)
                p = float(p_text)

                if not (1 <= i < j <= n):
                    raise ValueError(
                        f"Invalid pair {i}, {j}"
                    )

                if not (0 <= p <= 1):
                    raise ValueError(
                        f"Invalid probability {p}"
                    )

                if (i, j) in observed_pairs:
                    raise ValueError("Duplicate pair")

                observed_pairs[(i, j)] = p

        paired = [0.0] * (n + 1)

        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                p = observed_pairs.get((i, j), 0.0)

                pair_rows.append({
                    "sequence_id": rna_id,
                    "predictor": "EternaFold",
                    "i": i,
                    "j": j,
                    "probability": p
                })

                paired[i] += p
                paired[j] += p

    runtime = time.perf_counter() - start

    for i in range(1, n + 1):
        if not -0.001 <= paired[i] <= 1.001:
            raise ValueError(
                f"Invalid paired probability "
                f"at {rna_id}, position {i}"
            )

        nucleotide_rows.append({
            "sequence_id": rna_id,
            "predictor": "EternaFold",
            "position": i,
            "nucleotide": sequence[i - 1],
            "p_paired": paired[i],
            "p_unpaired": 1 - paired[i]
        })

    summary_rows.append({
        "sequence_id": rna_id,
        "length": n,
        "runtime_seconds": runtime,
        "reported_pairs": len(observed_pairs)
    })

    print(f"Processed {rna_id}: {n} nt")

pd.DataFrame(pair_rows).to_csv(
    OUTPUT / "eternafold_pairs.csv",
    index=False
)

pd.DataFrame(nucleotide_rows).to_csv(
    OUTPUT / "eternafold_unpaired.csv",
    index=False
)

pd.DataFrame(summary_rows).to_csv(
    OUTPUT / "eternafold_summary.csv",
    index=False
)

print("\nSUCCESS: EternaFold predictions completed!")
print("All three CSV files saved.")
