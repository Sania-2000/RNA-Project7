from pathlib import Path
import subprocess
import sys
import pandas as pd

ROOT = Path("large_scale")
TABLES = ROOT / "results/tables"
LOGS = ROOT / "logs"

TABLES.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

required = {
    "sequence_id", "probe", "predictor", "n",
    "sum_x", "sum_y", "sum_xx", "sum_xy",
    "sum_yy", "cluster_id", "source_group"
}

def validate(number, output):
    batch = pd.read_csv(
        ROOT / f"batches/batch_{number:03d}.csv",
        dtype={"sequence_id": str}
    )

    df = pd.read_csv(
        output,
        dtype={"sequence_id": str}
    )

    assert required.issubset(df.columns)
    assert len(df) == 6 * len(batch)
    assert df.sequence_id.nunique() == len(batch)

    assert set(df.sequence_id) == set(batch.sequence_id)

    assert not df.duplicated(
        ["sequence_id", "probe", "predictor"]
    ).any()

    assert set(df.probe) == {"DMS_MaP", "2A3_MaP"}
    assert set(df.predictor) == {
        "vienna", "linearpartition", "eternafold"
    }

    assert df.n.ge(10).all()
    assert df.cluster_id.notna().all()
    assert df.source_group.notna().all()

    # All three predictors must have identical numbers
    # of measured positions for each RNA and probe.
    counts = df.pivot(
        index=["sequence_id", "probe"],
        columns="predictor",
        values="n"
    )

    assert counts.notna().all().all()
    assert (counts.nunique(axis=1) == 1).all()

    return len(batch), len(df)


for number in range(1, 43):
    output = (
        TABLES /
        f"regression_statistics_batch_{number:03d}.csv"
    )

    if output.exists():
        n, records = validate(number, output)
        print(
            f"Batch {number:03d}: already valid, skipping",
            flush=True
        )
        continue

    print(f"Running regression batch {number:03d}...",
          flush=True)

    log = LOGS / f"regression_batch_{number:03d}.log"

    with log.open("w") as handle:
        subprocess.run(
            [
                sys.executable,
                "scripts/prepare_regression_batch.py",
                "--batch", str(number)
            ],
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True
        )

    n, records = validate(number, output)

    print(
        f"PASS batch {number:03d}: "
        f"{n} RNAs, {records} regression records",
        flush=True
    )

print("\nSUCCESS: All regression batches validated!")
