from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("expanded")
TABLES = ROOT / "results/tables"
TABLES.mkdir(parents=True, exist_ok=True)

experimental = pd.read_csv(
    "data/processed/expanded_reactivities.csv",
    dtype={"sequence_id": str}
)

predictors = {
    "ViennaRNA": ROOT / "predictions/vienna_unpaired.csv",
    "LinearPartition-V": ROOT / "predictions/linearpartition_unpaired_corrected.csv",
    "EternaFold": ROOT / "predictions/eternafold_unpaired.csv"
}

# Convert experimental profiles to nucleotide-level data
rows = []

for record in experimental.itertuples(index=False):
    seq = record.sequence.upper()
    probe = record.experiment_type

    for pos, nt in enumerate(seq, start=1):
        if probe == "DMS_MaP" and nt not in "AC":
            continue

        col = f"reactivity_{pos:04d}"
        value = getattr(record, col, np.nan)

        if not pd.notna(value):
            continue

        rows.append({
            "sequence_id": str(record.sequence_id),
            "probe": probe,
            "position": pos,
            "nucleotide": nt,
            "reactivity": value
        })

observed = pd.DataFrame(rows)
observed["reactivity"] = pd.to_numeric(
    observed["reactivity"], errors="coerce"
)
observed = observed.replace([np.inf, -np.inf], np.nan)
observed = observed.dropna(subset=["reactivity"])

observed.to_csv(
    TABLES / "expanded_probing_observations.csv",
    index=False
)

all_results = []

for name, path in predictors.items():
    predicted = pd.read_csv(
        path, dtype={"sequence_id": str}
    )

    merged = observed.merge(
        predicted[["sequence_id", "position", "p_unpaired"]],
        on=["sequence_id", "position"],
        how="inner",
        validate="many_to_one"
    )

    assert len(merged) == len(observed), (
        f"{name}: missing experimental positions"
    )

    assert merged["p_unpaired"].between(0, 1).all()

    for (rna_id, probe), group in merged.groupby(
        ["sequence_id", "probe"]
    ):
        if len(group) < 10:
            continue

        if (
            group["p_unpaired"].nunique() < 2
            or group["reactivity"].nunique() < 2
        ):
            continue

        rho, _ = spearmanr(
            group["p_unpaired"], group["reactivity"]
        )

        all_results.append({
            "sequence_id": rna_id,
            "probe": probe,
            "predictor": name,
            "n_nucleotides": len(group),
            "spearman_rho": rho
        })

results = pd.DataFrame(all_results)
results.to_csv(
    TABLES / "expanded_three_predictor_correlations.csv",
    index=False
)

summary = (
    results.groupby(["probe", "predictor"])
    .agg(
        n_RNAs=("sequence_id", "nunique"),
        median_rho=("spearman_rho", "median"),
        mean_rho=("spearman_rho", "mean")
    )
    .reset_index()
)

summary.to_csv(
    TABLES / "expanded_three_predictor_summary.csv",
    index=False
)

print("\n=== EXPANDED THREE-PREDICTOR RESULTS ===")
print(summary.round(4).to_string(index=False))
print("\nSaved expanded correlation tables.")
