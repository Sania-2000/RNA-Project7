from pathlib import Path
import re
import pandas as pd

ROOT = Path("large_scale")
cluster_file = ROOT / "clustering/ribonanza_90.clstr"

records = []
cluster_id = None

with open(cluster_file) as f:
    for line in f:
        line = line.strip()

        if line.startswith(">Cluster"):
            cluster_id = int(line.split()[1])
            continue

        match = re.search(r">(.+?)\.\.\.", line)

        if not match:
            raise ValueError(f"Cannot parse: {line}")

        sequence_id = match.group(1)
        representative = line.endswith("*")

        records.append({
            "sequence_id": sequence_id,
            "cluster_id": cluster_id,
            "is_representative": representative
        })

clusters = pd.DataFrame(records)

inventory = pd.read_csv(
    ROOT / "combined_20739_qc_passed.csv",
    dtype={"sequence_id": str}
)

clusters = clusters.merge(
    inventory[
        ["sequence_id", "length", "source_group"]
    ],
    on="sequence_id",
    validate="one_to_one"
)

assert len(clusters) == 20739
assert clusters.sequence_id.is_unique
assert clusters.cluster_id.nunique() == 16520

sizes = clusters.groupby("cluster_id").size()
clusters["cluster_size"] = clusters.cluster_id.map(sizes)

assert (
    clusters.groupby("cluster_id")["is_representative"]
    .sum() == 1
).all()

clusters.to_csv(
    ROOT / "rna_cluster_membership.csv",
    index=False
)

summary = (
    clusters.groupby("cluster_id")
    .agg(
        n_sequences=("sequence_id", "size"),
        n_sources=("source_group", "nunique"),
        representative_length=("length", "first")
    )
    .reset_index()
    .sort_values("n_sequences", ascending=False)
)

summary.to_csv(
    ROOT / "rna_cluster_summary.csv",
    index=False
)

print("RNAs assigned:", len(clusters))
print("Clusters identified:", clusters.cluster_id.nunique())
print("Largest cluster size:", sizes.max())
print("Singleton clusters:", (sizes == 1).sum())

print("\nLargest 10 clusters:")
print(summary.head(10).to_string(index=False))

print("\nSaved:")
print(ROOT / "rna_cluster_membership.csv")
print(ROOT / "rna_cluster_summary.csv")
