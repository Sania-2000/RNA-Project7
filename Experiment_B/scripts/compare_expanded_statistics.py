from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

TABLES = Path("expanded/results/tables")
TABLES.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    TABLES / "expanded_three_predictor_correlations.csv",
    dtype={"sequence_id": str}
)

comparisons = [
    ("EternaFold", "ViennaRNA"),
    ("EternaFold", "LinearPartition-V"),
    ("LinearPartition-V", "ViennaRNA")
]

rng = np.random.default_rng(42)
rows = []

for probe, group in df.groupby("probe"):

    wide = group.pivot(
        index="sequence_id",
        columns="predictor",
        values="spearman_rho"
    )

    for model_a, model_b in comparisons:

        paired = wide[[model_a, model_b]].dropna()

        differences = (
            paired[model_a] - paired[model_b]
        ).to_numpy()

        n = len(differences)

        if n == 0:
            raise ValueError("No matched RNA sequences")

        bootstrap = np.array([
            np.median(
                rng.choice(
                    differences,
                    size=n,
                    replace=True
                )
            )
            for _ in range(5000)
        ])

        lower, upper = np.percentile(
            bootstrap, [2.5, 97.5]
        )

        if np.allclose(differences, 0):
            p = 1.0
        else:
            p = wilcoxon(differences).pvalue

        rows.append({
            "probe": probe,
            "model_a": model_a,
            "model_b": model_b,
            "n_RNAs": n,
            "median_paired_difference": np.median(differences),
            "CI_lower": lower,
            "CI_upper": upper,
            "wilcoxon_p": p
        })

results = pd.DataFrame(rows)

results.to_csv(
    TABLES / "expanded_statistical_comparison.csv",
    index=False
)

print("\n=== EXPANDED STATISTICAL COMPARISON ===")
print(results.round(5).to_string(index=False))

print("\nSaved expanded_statistical_comparison.csv")
