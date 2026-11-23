"""Experiment B: length-stratified per-RNA Spearman correlations.

Run from anywhere:
    python scripts/analyze_length_stratification.py

Input: Experiment_B/results/tables/large_probing_correlations_20739.csv
Output: two CSVs and one PNG under Experiment_B/results/.

Resample complete 90%-CD-HIT sequence clusters *within each length*;
retain all their RNA observations and preserve sampled multiplicity.
Report medians and 95% percentile CIs over 500 cluster-bootstrap samples.
No RNA structure predictions or reactivity profiles are reprocessed.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

DATA = TABLES / "large_probing_correlations_20739.csv"
B = 500
SEED = 91
RNG = np.random.default_rng(SEED)
PREDICTORS = ("vienna", "linearpartition", "eternafold")
PROBES = ("DMS_MaP", "2A3_MaP")
PAIRS = (("eternafold", "vienna"),
         ("eternafold", "linearpartition"),
         ("linearpartition", "vienna"))

frame = pd.read_csv(DATA, dtype={"sequence_id": str, "cluster_id": str})
assert len(frame) == 20_739 * 6
assert frame.sequence_id.nunique() == 20_739
assert not frame.duplicated(["sequence_id", "predictor", "probe"]).any()
assert frame.spearman_rho.notna().all()
assert frame.length.notna().all() and frame.cluster_id.notna().all()

meta = frame[["sequence_id", "length", "cluster_id", "source_group"]].drop_duplicates()
assert len(meta) == 20_739
source_counts = meta.groupby("length").source_group.nunique().to_dict()

summaries, differences = [], []
for length in sorted(frame.length.unique()):
    for probe in PROBES:
        sub = frame[(frame.length == length) & (frame.probe == probe)]
        wide = sub.pivot(index="sequence_id", columns="predictor", values="spearman_rho")
        attrs = sub[["sequence_id", "cluster_id"]].drop_duplicates().set_index("sequence_id")
        wide = wide.join(attrs)
        assert set(PREDICTORS) == set(wide.columns) - {"cluster_id"}
        assert wide[list(PREDICTORS)].notna().all().all()
        assert len(wide) == meta.loc[meta.length == length, "sequence_id"].nunique()

        v = wide[list(PREDICTORS)].to_numpy(dtype=float)
        _, cluster_idx = np.unique(wide.cluster_id.to_numpy(), return_inverse=True)
        clusters = int(cluster_idx.max()) + 1
        n = len(v)
        rna_indices = np.arange(n)
        n_sources = source_counts[length]

        observed = np.median(v, axis=0)
        diff_matrix = np.column_stack([v[:, PREDICTORS.index(a)] - v[:, PREDICTORS.index(b)] for a, b in PAIRS])
        observed_differences = np.median(diff_matrix, axis=0)
        bootstrap_medians = np.empty((B, 3))
        bootstrap_diffs = np.empty((B, 3))

        for replicate in range(B):
            # Sampling cluster identifiers with replacement implements cluster bootstrap.
            draws = RNG.integers(0, clusters, size=clusters)
            multiplicities = np.bincount(draws, minlength=clusters)
            # RNA rows inherit multiplicity of their parent cluster.
            # Repetition exactly matches sampling all RNA members of each drawn cluster.
            sampled_rnas = np.repeat(rna_indices, multiplicities[cluster_idx])
            bootstrap_medians[replicate, :] = np.median(v[sampled_rnas], axis=0)
            bootstrap_diffs[replicate, :] = np.median(diff_matrix[sampled_rnas], axis=0)

        for col, predictor in enumerate(PREDICTORS):
            lower, upper = np.quantile(bootstrap_medians[:, col], (0.025, 0.975))
            summaries.append({
                "length_nt": int(length), "probe": probe, "predictor": predictor,
                "n_rnas": n, "n_clusters_within_length": clusters,
                "n_sources": n_sources, "median_rho": observed[col],
                "ci_lower": lower, "ci_upper": upper,
                "bootstrap_replicates": B, "bootstrap_seed": SEED,
            })

        for col, (a, b) in enumerate(PAIRS):
            lower, upper = np.quantile(bootstrap_diffs[:, col], (0.025, 0.975))
            differences.append({
                "length_nt": int(length), "probe": probe,
                "predictor_a": a, "predictor_b": b,
                "n_rnas": n, "n_clusters_within_length": clusters,
                "n_sources": n_sources,
                "median_paired_difference": observed_differences[col],
                "ci_lower": lower, "ci_upper": upper,
                "bootstrap_replicates": B, "bootstrap_seed": SEED,
            })
        print(f"Length {length:3}, {probe:7}: {n:5} RNAs, {clusters:5} clusters, {n_sources} sources", flush=True)

summary = pd.DataFrame(summaries)
diffs = pd.DataFrame(differences)
summary_file = TABLES / "length_stratified_spearman.csv"
diffs_file = TABLES / "length_stratified_paired_differences.csv"
summary.to_csv(summary_file, index=False)
diffs.to_csv(diffs_file, index=False)

# Two-panel visualization, displaying the original per-RNA correlation estimates.
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), sharey=True, constrained_layout=True)
labels = {"vienna": "ViennaRNA", "linearpartition": "LinearPartition-V", "eternafold": "EternaFold"}
lengths = sorted(summary.length_nt.unique())
for ax, probe in zip(axes, PROBES):
    for predictor in PREDICTORS:
        s = summary[(summary.probe == probe) & (summary.predictor == predictor)].set_index("length_nt").loc[lengths]
        ax.errorbar(lengths, s.median_rho, yerr=[s.median_rho - s.ci_lower, s.ci_upper - s.median_rho],
                    marker="o", capsize=3, label=labels[predictor], linewidth=1.6)
    ax.set_title("DMS-MaP (A/C positions)" if probe == "DMS_MaP" else "2A3-MaP (all measured bases)")
    ax.set_xlabel("RNA length (nt)")
    ax.set_xticks(lengths)
    ax.grid(axis="y", alpha=0.25)
axes[0].set_ylabel("Median per-RNA Spearman correlation")
axes[1].legend(loc="best", fontsize=8)
fig.suptitle("Experiment B: correlation stratified by RNA length (95% cluster-bootstrap CIs)")
fig.savefig(FIGURES / "length_stratified_spearman.png", dpi=220)
plt.close(fig)

print("Saved:", summary_file, diffs_file, FIGURES / "length_stratified_spearman.png", sep="\n - ")
print("CAUTION: Length is heavily confounded with source identity; these are descriptive strata.")
