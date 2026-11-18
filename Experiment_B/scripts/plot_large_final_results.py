from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path("large_scale")
TABLES = ROOT / "results/tables"
FIGURES = ROOT / "results/figures"
FIGURES.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    TABLES / "cluster_bootstrap_medians.csv"
)

predictors = ["vienna", "linearpartition", "eternafold"]

labels = {
    "vienna": "ViennaRNA",
    "linearpartition": "LinearPartition-V",
    "eternafold": "EternaFold"
}

probes = ["DMS_MaP", "2A3_MaP"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

for ax, probe in zip(axes, probes):
    part = df[df.probe == probe].set_index("predictor")
    part = part.loc[predictors]

    x = np.arange(len(predictors))
    medians = part.median_rho.to_numpy()

    lower = medians - part.ci_lower.to_numpy()
    upper = part.ci_upper.to_numpy() - medians

    ax.bar(
        x,
        medians,
        color=["#607d8b", "#78909c", "#4c956c"],
        width=0.65
    )

    ax.errorbar(
        x,
        medians,
        yerr=np.vstack([lower, upper]),
        fmt="none",
        ecolor="black",
        capsize=5,
        linewidth=1.4
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [labels[p] for p in predictors],
        rotation=15,
        ha="right"
    )

    ax.set_ylim(0, 0.75)
    ax.set_ylabel("Median Spearman correlation (ρ)")
    ax.set_title(probe.replace("_MaP", ""))
    ax.grid(axis="y", alpha=0.2)

fig.suptitle(
    "RNA unpaired-probability association with chemical probing",
    fontsize=13
)

fig.text(
    0.5, 0.01,
    "20,739 RNAs; 95% cluster-bootstrap confidence intervals",
    ha="center",
    fontsize=9
)

fig.tight_layout(rect=[0, 0.06, 1, 0.94])

png = FIGURES / "large_predictor_comparison.png"
pdf = FIGURES / "large_predictor_comparison.pdf"

fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

print("Saved:", png)
print("Saved:", pdf)
print("SUCCESS: Final comparison figure created!")
