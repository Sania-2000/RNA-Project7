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
        batch_file, dtype={"sequence_id": str}
    )
    pred = pd.read_csv(
        output_file, dtype={"sequence_id": str}
    )

    assert batch.sequence_id.is_unique
    assert pred.sequence_id.nunique() == len(batch)
    assert set(pred.sequence_id) == set(batch.sequence_id)
    assert len(pred) == int(batch.length.sum())

    assert not pred.duplicated(
        ["sequence_id", "position"]
    ).any()

    assert np.isfinite(
        pred[["p_paired", "p_unpaired_raw", "p_unpaired"]]
    ).all().all()

    assert pred.p_unpaired.between(0, 1).all()

    assert np.allclose(
        pred.p_unpaired_raw,
        1 - pred.p_paired,
        atol=1e-10
    )

    correction = (
        pred.p_unpaired - pred.p_unpaired_raw
    ).abs()

    assert correction.max() <= 0.001000001

    expected_clipped = (
        (pred.p_unpaired_raw < 0) |
        (pred.p_unpaired_raw > 1)
    )

    assert np.allclose(
        pred.p_unpaired,
        pred.p_unpaired_raw.clip(0, 1),
        atol=1e-12,
        rtol=0
    )
    assert (
        pred.was_clipped.astype(str).str.lower().eq("true")
        == expected_clipped
    ).all()

    counts = pred.groupby("sequence_id").size()
    lengths = batch.set_index("sequence_id").length

    assert counts.sort_index().equals(
        lengths.sort_index()
    )

    expected_positions = batch.set_index("sequence_id")[
        "sequence"
    ].to_dict()

    for seq_id, group in pred.groupby("sequence_id"):
        seq = expected_positions[seq_id]
        group = group.sort_values("position")
        assert group.position.tolist() == list(
            range(1, len(seq) + 1)
        )
        assert "".join(group.nucleotide) == seq

    return (
        len(batch),
        len(pred),
        int(expected_clipped.sum()),
        float(correction.max())
    )


for row in manifest.itertuples(index=False):
    number = int(row.batch)

    batch_file = BATCHES / f"batch_{number:03d}.csv"
    output_file = OUTPUT / f"linearpartition_batch_{number:03d}.csv"
    log_file = LOGS / f"linearpartition_batch_{number:03d}.log"

    if output_file.exists():
        n, positions, corrected, maximum = validate(
            batch_file, output_file
        )
        print(
            f"Batch {number:03d}: already valid, skipping",
            flush=True
        )
        continue

    print(f"Running batch {number:03d}...", flush=True)

    command = [
        sys.executable,
        "scripts/predict_streaming_linearpartition.py",
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

    n, positions, corrected, maximum = validate(
        batch_file, output_file
    )

    print(
        f"PASS batch {number:03d}: "
        f"{n} RNAs, {positions} positions, "
        f"{corrected} corrections, "
        f"maximum={maximum:.6f}",
        flush=True
    )

print("\nSUCCESS: All LinearPartition batches validated!")
