from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("large_scale")
OUT.mkdir(exist_ok=True)

df = pd.read_csv(
    OUT / "available_rna_inventory.csv",
    dtype={"sequence_id": str}
)

TARGET = 20000
SEED = 42

# Preserve original source labels.
df["source_group_original"] = df["source_group"]

# Reproducible source-stratified sampling.
groups = {
    name: part
    for name, part in df.groupby("source_group", sort=True)
}

# Allocate samples proportionally to source size.
sizes = pd.Series({
    name: len(part) for name, part in groups.items()
})

quotas_float = sizes / sizes.sum() * TARGET
quotas = np.floor(quotas_float).astype(int)

remaining = TARGET - quotas.sum()

fractional = (
    quotas_float - quotas
).sort_values(ascending=False, kind="stable")

for name in fractional.index[:remaining]:
    quotas.loc[name] += 1

samples = []

for name in sorted(groups):
    n = int(quotas.loc[name])
    if n:
        samples.append(
            groups[name].sample(
                n=n,
                random_state=SEED
            )
        )

sample = pd.concat(samples, ignore_index=True)
sample = sample.sort_values("sequence_id").reset_index(drop=True)

assert len(sample) == TARGET
assert sample.sequence_id.nunique() == TARGET
assert sample.sequence.nunique() == TARGET

sample.to_csv(
    OUT / "candidate_20000_inventory.csv",
    index=False
)

summary = (
    sample.groupby("source_group")
    .size()
    .reset_index(name="n_selected")
    .sort_values("n_selected", ascending=False)
)

summary.to_csv(
    OUT / "candidate_20000_by_source.csv",
    index=False
)

print("Candidate sequences:", len(sample))
print("Unique RNA sequences:", sample.sequence.nunique())
print("\nLength distribution:")
print(sample.length.value_counts().sort_index())

print("\nTop source groups:")
print(summary.head(15).to_string(index=False))

print("\nSaved candidate inventory.")
