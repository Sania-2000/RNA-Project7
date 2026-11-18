from pathlib import Path
import pandas as pd

DATA = Path(
    "data/processed/expanded_rna_inventory.csv"
)

SIMILARITY = Path(
    "results/tables/expanded_sequence_similarity.csv"
)

OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)

THRESHOLD = 0.90

df = pd.read_csv(DATA, dtype={"sequence_id": str})
sim = pd.read_csv(
    SIMILARITY,
    dtype={
        "sequence_id_a": str,
        "sequence_id_b": str
    }
)

ids = df["sequence_id"].tolist()
parent = {x: x for x in ids}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

# Connect RNAs with high positional identity
similar = sim[sim["identity"] >= THRESHOLD]

for row in similar.itertuples(index=False):
    union(row.sequence_id_a, row.sequence_id_b)

# Assign clusters
df["cluster_id"] = df["sequence_id"].map(find)

# Choose one representative per cluster
representatives = (
    df.sort_values("sequence_id")
    .drop_duplicates("cluster_id")
    .copy()
)

df.to_csv(
    OUT / "expanded_rna_clusters.csv",
    index=False
)

representatives.to_csv(
    OUT / "expanded_rna_nonredundant.csv",
    index=False
)

print("\n=== RNA CLUSTERING ===")
print("Original RNA sequences:", len(df))
print("Clusters:", df.cluster_id.nunique())
print("Representative sequences:", len(representatives))

print("\nRepresentatives by RNA length:")
print(
    representatives["length"]
    .value_counts()
    .sort_index()
)

print("\nLargest clusters:")
print(
    df.cluster_id.value_counts().head(10)
)

print("\nSaved cluster assignments and representatives.")
