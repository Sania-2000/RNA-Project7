from pathlib import Path
from time import perf_counter
import argparse
import csv
import pandas as pd
import RNA

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--limit", type=int, default=0)
args = parser.parse_args()

df = pd.read_csv(
    args.input,
    dtype={"sequence_id": str}
)

if args.limit:
    df = df.head(args.limit)

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
        "p_unpaired"
    ])

    for count, row in enumerate(df.itertuples(index=False), 1):
        sequence = row.sequence.upper()
        n = len(sequence)

        start = perf_counter()

        fc = RNA.fold_compound(sequence)
        _, mfe = fc.mfe()

        fc.exp_params_rescale(mfe)
        fc.pf()

        bpp = fc.bpp()
        paired = [0.0] * (n + 1)

        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                p = float(bpp[i][j])
                paired[i] += p
                paired[j] += p

        for i in range(1, n + 1):
            p_paired = paired[i]
            p_unpaired = 1.0 - p_paired

            if not (-1e-6 <= p_unpaired <= 1 + 1e-6):
                raise ValueError(
                    f"Invalid probability: {row.sequence_id}, {i}"
                )

            writer.writerow([
                row.sequence_id,
                "ViennaRNA",
                i,
                sequence[i - 1],
                p_paired,
                p_unpaired
            ])

        print(
            f"Completed {count}/{len(df)} | "
            f"{n} nt | "
            f"{perf_counter() - start:.3f}s",
            flush=True
        )

print("\nSaved:", output)
print("RNAs processed:", len(df))
