from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("large_scale/results/tables")
DATA = ROOT / "large_probing_correlations_20739.csv"

BOOTSTRAPS = 500
SEED = 91

df = pd.read_csv(DATA, dtype={"sequence_id": str})
rng = np.random.default_rng(SEED)

predictors = ["vienna", "linearpartition", "eternafold"]

comparisons = [
    ("eternafold", "vienna"),
    ("eternafold", "linearpartition"),
    ("linearpartition", "vienna")
]

summary_rows = []
comparison_rows = []

for probe in ["DMS_MaP", "2A3_MaP"]:
    sub = df[df.probe == probe]

    wide = sub.pivot(
        index="sequence_id",
        columns="predictor",
        values="spearman_rho"
    )

    metadata = (
        sub[["sequence_id", "cluster_id"]]
        .drop_duplicates()
        .set_index("sequence_id")
    )

    wide = wide.join(metadata)

    assert len(wide) == 20739
    assert wide[predictors].notna().all().all()
    assert wide.cluster_id.notna().all()

    wide = wide.sort_values("cluster_id")

    clusters = wide.cluster_id.to_numpy()
    unique_clusters, starts, counts = np.unique(
        clusters,
        return_index=True,
        return_counts=True
    )

    assert len(unique_clusters) == 16520

    values = wide[predictors].to_numpy(dtype=float)
    n_clusters = len(unique_clusters)

    medians = np.median(values, axis=0)

    diff_pairs = [
        (
            predictors.index(a),
            predictors.index(b)
        )
        for a, b in comparisons
    ]

    observed_diffs = np.array([
        np.median(values[:, a] - values[:, b])
        for a, b in diff_pairs
    ])

    boot_medians = np.empty((BOOTSTRAPS, 3))
    boot_diffs = np.empty((BOOTSTRAPS, 3))

    for b in range(BOOTSTRAPS):
        sampled = rng.integers(
            0, n_clusters, size=n_clusters
        )

        # Sample whole clusters with replacement.
        # Repeated clusters retain their multiplicity.
        indices = np.concatenate([
            np.arange(starts[k], starts[k] + counts[k])
            for k in sampled
        ])

        sample_values = values[indices]

        boot_medians[b] = np.median(
            sample_values, axis=0
        )

        boot_diffs[b] = [
            np.median(
                sample_values[:, a] - sample_values[:, c]
            )
            for a, c in diff_pairs
        ]

    for j, name in enumerate(predictors):
        lo, hi = np.quantile(
            boot_medians[:, j], [0.025, 0.975]
        )

        summary_rows.append({
            "probe": probe,
            "predictor": name,
            "n_rnas": len(wide),
            "n_clusters": n_clusters,
            "median_rho": medians[j],
            "ci_lower": lo,
            "ci_upper": hi,
            "bootstrap_replicates": BOOTSTRAPS
        })

    for j, (a, c) in enumerate(comparisons):
        lo, hi = np.quantile(
            boot_diffs[:, j], [0.025, 0.975]
        )

        comparison_rows.append({
            "probe": probe,
            "predictor_a": a,
            "predictor_b": c,
            "median_paired_difference": observed_diffs[j],
            "ci_lower": lo,
            "ci_upper": hi,
            "bootstrap_replicates": BOOTSTRAPS
        })

    print(
        f"PASS: {probe}, {len(wide)} RNAs, "
        f"{n_clusters} clusters",
        flush=True
    )

summary = pd.DataFrame(summary_rows)
comparisons_df = pd.DataFrame(comparison_rows)

summary.to_csv(
    ROOT / "cluster_bootstrap_medians.csv",
    index=False
)

comparisons_df.to_csv(
    ROOT / "cluster_bootstrap_comparisons.csv",
    index=False
)

print("\nMEDIAN CORRELATIONS AND 95% CIs:")
print(summary.round(4).to_string(index=False))

print("\nPAIRED DIFFERENCES AND 95% CIs:")
print(comparisons_df.round(4).to_string(index=False))

print("\nSUCCESS: Cluster bootstrap completed!")
