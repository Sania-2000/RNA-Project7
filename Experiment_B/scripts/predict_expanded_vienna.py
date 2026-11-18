
import time
from pathlib import Path

import pandas as pd
import RNA

input_file = "expanded/expanded_139_sequences.csv"

output_dir = Path("expanded/predictions")
output_dir.mkdir(exist_ok=True)

df = pd.read_csv(input_file)

# Both probes use the same sequence,
# so each unique sequence is predicted once.
unique_rnas = df[
    ["sequence_id", "sequence"]
].drop_duplicates()

pair_rows = []
nucleotide_rows = []
summary_rows = []

for record in unique_rnas.itertuples(index=False):

    rna_id = record.sequence_id
    sequence = record.sequence.upper()
    length = len(sequence)

    start = time.perf_counter()

    fc = RNA.fold_compound(sequence)
    structure, mfe = fc.mfe()

    fc.exp_params_rescale(mfe)
    fc.pf()

    bpp = fc.bpp()

    paired = [0.0] * (length + 1)

    for i in range(1, length + 1):
        for j in range(i + 1, length + 1):

            probability = float(bpp[i][j])

            # Save every possible pair
            pair_rows.append({
                "sequence_id": rna_id,
                "predictor": "ViennaRNA",
                "i": i,
                "j": j,
                "probability": probability
            })

            paired[i] += probability
            paired[j] += probability

    runtime = time.perf_counter() - start

    for i in range(1, length + 1):

        p_paired = paired[i]

        if not -1e-6 <= p_paired <= 1 + 1e-6:
            raise ValueError(
                f"Invalid probability at "
                f"{rna_id}, position {i}"
            )

        nucleotide_rows.append({
            "sequence_id": rna_id,
            "predictor": "ViennaRNA",
            "position": i,
            "nucleotide": sequence[i - 1],
            "p_paired": p_paired,
            "p_unpaired": 1 - p_paired
        })

    summary_rows.append({
        "sequence_id": rna_id,
        "length": length,
        "mfe_structure": structure,
        "mfe_kcal_mol": mfe,
        "runtime_seconds": runtime
    })

    print(f"Processed {rna_id}: {length} nt")

pd.DataFrame(pair_rows).to_csv(
    output_dir / "vienna_pairs.csv",
    index=False
)

pd.DataFrame(nucleotide_rows).to_csv(
    output_dir / "vienna_unpaired.csv",
    index=False
)

pd.DataFrame(summary_rows).to_csv(
    output_dir / "vienna_summary.csv",
    index=False
)

print("\nAll predictions completed!")
print("Saved three files in prediction_outputs/")
