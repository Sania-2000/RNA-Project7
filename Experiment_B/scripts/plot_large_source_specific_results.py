from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from textwrap import shorten

ROOT = Path("large_scale")
TABLE = ROOT / "results/tables/large_source_specific_differences.csv"
FIGURES = ROOT / "results/figures"
FIGURES.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(TABLE)

required = {
    "source_group",
    "probe",
    "n_rnas",
    "median_eterna_minus_vienna",
    "median_eterna_minus_linear",
}

assert required.issubset(df.columns)

probes = ["DMS_MaP", "2A3_MaP"]
comparisons = [
    ("median_eterna_minus_vienna", "EternaFold - ViennaRNA"),
    ("median_eterna_minus_linear", "EternaFold - LinearPartition-V"),
]

# Select the same 10 largest source groups for all panels.
top_sources = (
    df.groupby("source_group")["n_rnas"]
    .max()
    .nlargest(10)
    .index
    .tolist()
)

labels = {
    name: shorten(name, width=48, placeholder="...")
    for name in top_sources
}

fig, axes = plt.subplots(
    2, 2,
    figsize=(17, 12),
    sharex=True
)

for row, probe in enumerate(probes):
    part = (
        df[df["probe"] == probe]
        .set_index("source_group")
        .reindex(top_sources)
    )

    for col, (column, description) in enumerate(comparisons):
        ax = axes[row, col]
        values = part[column].to_numpy(dtype=float)
        y = np.arange(len(top_sources))

        colors = [
            "#4c956c" if v >= 0 else "#b65c5c"
            for v in values
        ]

        ax.barh(y, values, color=colors)
        ax.axvline(0, color="black", linewidth=0.8)

        ax.set_yticks(y)
        ax.set_yticklabels([labels[s] for s in top_sources], fontsize=8)
        ax.invert_yaxis()

        ax.set_title(f"{probe}: {description}")
        ax.set_xlabel("Median paired Spearman difference (Δρ)")
        ax.grid(axis="x", alpha=0.2)

        for i, source in enumerate(top_sources):
            n = part.loc[source, "n_rnas"]
            if pd.notna(n):
                ax.text(
                    0.98, i,
                    f"n={int(n)}",
                    transform=ax.get_yaxis_transform(),
                    ha="right",
                    va="center",
                    fontsize=7
                )

fig.suptitle(
    "Source-specific RNA probing performance",
    fontsize=16
)

fig.text(
    0.5, 0.015,
    "20,739-RNA dataset; positive Δρ favors EternaFold. "
    "Bars are descriptive source-specific medians.",
    ha="center",
    fontsize=10
)

fig.tight_layout(rect=[0, 0.045, 1, 0.95])

png = FIGURES / "large_source_specific_comparison.png"
pdf = FIGURES / "large_source_specific_comparison.pdf"

fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

print("Saved:", png)
print("Saved:", pdf)
print("SUCCESS: Source-specific comparison figure created!")
