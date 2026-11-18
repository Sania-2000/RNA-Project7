from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path("large_scale")
BATCHES = ROOT / "batches"
OUTPUT = ROOT / "full_predictions"
LOGS = ROOT / "logs"

OUTPUT.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(ROOT / "batch_manifest.csv")


def validate(batch_file, output_file):
    batch = pd.read_csv(
        batch_file,
        dtype={"sequence_id": str}
    )

    pred = pd.read_csv(
        output_file,
        dtype={"sequence_id": str}
    )

    assert batch.sequence_id.is_unique
    assert set(pred.sequence_id) == set(batch.sequence_id)
    assert len(pred) == int(batch.length.sum())

    assert not pred.duplicated(
        ["sequence_id", "position"]
    ).any()

    assert np.isfinite(
        pred[["p_paired", "p_unpaired"]].to_numpy()
    ).all()

    assert pred.p_unpaired.between(0, 1).all()
    # Raw paired sums can exceed 1 slightly due to
    # numerical approximation. Do not alter these values.
    assert pred.p_paired.between(-0.001, 1.001).all()

    if "p_unpaired_raw" in pred.columns:
        # New eight-column format.
        assert "was_clipped" in pred.columns
        assert np.isfinite(pred.p_unpaired_raw).all()
        assert pred.p_unpaired_raw.between(-0.001, 1.001).all()

        assert np.allclose(
            pred.p_unpaired_raw,
            1 - pred.p_paired,
            atol=1e-10,
            rtol=0
        )

        assert np.allclose(
            pred.p_unpaired,
            pred.p_unpaired_raw.clip(0, 1),
            atol=1e-12,
            rtol=0
        )

        expected_flags = (
            (pred.p_unpaired_raw < 0) |
            (pred.p_unpaired_raw > 1)
        )
        saved_flags = (
            pred.was_clipped.astype(str).str.lower().eq("true")
        )

        assert (expected_flags == saved_flags).all()

    else:
        # Only the already-validated batch 001 uses
        # the original six-column output.
        assert output_file.name == "eternafold_batch_001.csv"

        assert np.allclose(
            pred.p_unpaired,
            1 - pred.p_paired,
            atol=1e-10,
            rtol=0
        )

    reference = batch.set_index("sequence_id")
    counts = pred.groupby("sequence_id").size()

    assert counts.sort_index().equals(
        reference["length"].sort_index()
    )

    for seq_id, group in pred.groupby("sequence_id"):
        sequence = reference.loc[seq_id, "sequence"]
        group = group.sort_values("position")

        assert group.position.tolist() == list(
            range(1, len(sequence) + 1)
        )

        assert "".join(group.nucleotide) == sequence

    return len(batch), len(pred)


for row in manifest.itertuples(index=False):
    number = int(row.batch)

    batch_file = BATCHES / f"batch_{number:03d}.csv"
    output_file = OUTPUT / f"eternafold_batch_{number:03d}.csv"
    log_file = LOGS / f"eternafold_batch_{number:03d}.log"

    if output_file.exists():
        n, positions = validate(batch_file, output_file)

        print(
            f"Batch {number:03d}: already valid, skipping",
            flush=True
        )
        continue

    print(f"Running batch {number:03d}...", flush=True)

    command = [
        sys.executable,
        "scripts/predict_streaming_eternafold.py",
        "--input", str(batch_file),
        "--output", str(output_file)
    ]

    with log_file.open("w") as log:
        subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True
        )

    n, positions = validate(batch_file, output_file)

    print(
        f"PASS batch {number:03d}: "
        f"{n} RNAs, {positions} positions",
        flush=True
    )

print("\nSUCCESS: All EternaFold batches validated!")
