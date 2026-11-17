
"""
Experiment A — RNA3DB-informed 26-group sensitivity analysis.

Reproduces the 26-group bootstrap by merging supplementary
sections 2.4 and 4.5.

This is an exploratory operational grouping, not verified
evolutionary-family-disjoint validation.

Inputs:
    45RNA_matched_keep_7r6q_6.csv
    45RNA_matched_keep_7u0h_6.csv

Run:
    python reproduce_26group_sensitivity.py --input-dir PATH --output-dir PATH
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


METHODS = {
    "ViennaRNA": "p_vienna",
    "LinearPartition-V": "p_linear",
    "EternaFold": "p_eterna",
}

SCENARIOS = {
    "Keep 7r6q_6": "45RNA_matched_keep_7r6q_6.csv",
    "Keep 7u0h_6": "45RNA_matched_keep_7u0h_6.csv",
}

METRICS = [
    "ECE_nt_pooled",
    "Brier_nt_pooled",
    "LogLoss_nt_pooled",
    "ECE_equal_section_mixture",
    "ECE_mean_section",
    "Brier_mean_section",
    "LogLoss_mean_section",
]


def calculate_metrics(counts, n_group, sum_p, sum_y,
                      brier_sums, loss_sums):
    m = np.asarray(counts, dtype=float)
    n = m @ n_group
    delta = sum_p - sum_y

    pooled_ece = np.abs(m @ delta).sum(axis=-1) / n
    pooled_brier = (m @ brier_sums) / n
    pooled_loss = (m @ loss_sums) / n

    equal_bins = m @ (delta / n_group[:, None])

    equal_ece = (
        np.abs(equal_bins).sum(axis=-1)
        / m.sum(axis=-1)
    )

    mean_ece = (
        m @ (np.abs(delta).sum(axis=1) / n_group)
    ) / m.sum(axis=-1)

    mean_brier = (
        m @ (brier_sums / n_group)
    ) / m.sum(axis=-1)

    mean_loss = (
        m @ (loss_sums / n_group)
    ) / m.sum(axis=-1)

    return np.stack([
        pooled_ece, pooled_brier, pooled_loss,
        equal_ece, mean_ece, mean_brier, mean_loss
    ], axis=-1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path,
                        default=Path("."))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("26_group_bootstrap_results"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    n_boot = 5000
    seed = 42
    rng = np.random.default_rng(seed)

    summary = []
    paired = []
    sectionrows = []
    group_mapping = []

    for scenario, filename in SCENARIOS.items():

        df = pd.read_csv(
            args.input_dir / filename,
            dtype={"supp_section": "string"}
        )

        required = [
            "structure_id", "supp_section",
            "reference_paired_label",
            *METHODS.values()
        ]
        assert not df[required].isna().any().any()
        assert set(df.reference_paired_label.unique()) <= {0, 1}

        # Create the conservative 26-group alternative.
        df["original_supp_section"] = df["supp_section"]
        df["sensitivity_group"] = df["supp_section"].astype(str)

        df.loc[
            df["supp_section"].isin(["2.4", "4.5"]),
            "sensitivity_group"
        ] = "merged_2.4_4.5"

        groups = sorted(df.sensitivity_group.unique())
        G = len(groups)

        assert df.structure_id.nunique() == 45
        assert len(df) == 2708
        assert G == 26

        for row in df[
            ["structure_id", "original_supp_section",
             "sensitivity_group"]
        ].drop_duplicates().itertuples(index=False):
            group_mapping.append({
                "scenario": scenario,
                "structure_id": row.structure_id,
                "original_supp_section": row.original_supp_section,
                "sensitivity_group": row.sensitivity_group
            })

        print(
            "\nScenario:", scenario,
            "| RNAs:", df.structure_id.nunique(),
            "| nt:", len(df),
            "| groups:", G
        )

        draws = rng.integers(0, G, size=(n_boot, G))
        counts = np.zeros((n_boot, G), dtype=np.int16)

        for k in range(n_boot):
            counts[k] = np.bincount(
                draws[k], minlength=G
            )

        n_group = np.zeros(G)
        predictions = {}

        for method, column in METHODS.items():

            sum_p = np.zeros((G, 10))
            sum_y = np.zeros((G, 10))
            brier_sums = np.zeros(G)
            loss_sums = np.zeros(G)

            for gi, group in enumerate(groups):

                sub = df[df.sensitivity_group == group]
                y = sub.reference_paired_label.to_numpy(float)
                p = sub[column].to_numpy(float)

                assert np.all(np.isfinite(p))
                assert np.all((p >= 0) & (p <= 1))

                bins = np.minimum(
                    np.floor(p * 10).astype(int), 9
                )

                n_group[gi] = len(sub)

                sum_p[gi] = np.bincount(
                    bins, weights=p, minlength=10
                )
                sum_y[gi] = np.bincount(
                    bins, weights=y, minlength=10
                )

                brier_sums[gi] = np.sum((p - y) ** 2)

                p_safe = np.clip(p, 1e-12, 1 - 1e-12)

                loss_sums[gi] = np.sum(
                    -y * np.log(p_safe)
                    - (1-y) * np.log(1-p_safe)
                )

                sectionrows.append({
                    "scenario": scenario,
                    "supp_section": group,
                    "predictor": method,
                    "nt": len(sub),
                    "RNAs": sub.structure_id.nunique(),
                    "ECE_section": np.sum(
                        np.abs(sum_p[gi] - sum_y[gi])
                    ) / n_group[gi],
                    "Brier_section":
                        brier_sums[gi] / n_group[gi],
                    "LogLoss_section":
                        loss_sums[gi] / n_group[gi]
                })

            observed = calculate_metrics(
                np.ones(G), n_group, sum_p, sum_y,
                brier_sums, loss_sums
            )

            boots = calculate_metrics(
                counts, n_group, sum_p, sum_y,
                brier_sums, loss_sums
            )

            predictions[method] = (observed, boots)

            for j, metric in enumerate(METRICS):
                lo, hi = np.quantile(
                    boots[:, j], [0.025, 0.975]
                )

                summary.append({
                    "scenario": scenario,
                    "predictor": method,
                    "metric": metric,
                    "estimate": observed[j],
                    "CI_lower": lo,
                    "CI_upper": hi,
                    "n_RNAs": 45,
                    "n_nucleotides": 2708,
                    "n_groups": G,
                    "bootstrap_iterations": n_boot,
                    "seed": seed
                })

        for method in ["ViennaRNA", "LinearPartition-V"]:

            improvement = (
                predictions[method][0]
                - predictions["EternaFold"][0]
            )

            boot_difference = (
                predictions[method][1]
                - predictions["EternaFold"][1]
            )

            for j, metric in enumerate(METRICS):
                lo, hi = np.quantile(
                    boot_difference[:, j],
                    [0.025, 0.975]
                )

                paired.append({
                    "scenario": scenario,
                    "contrast": method + " minus EternaFold",
                    "metric": metric,
                    "improvement_estimate": improvement[j],
                    "CI_lower": lo,
                    "CI_upper": hi,
                    "n_RNAs": 45,
                    "n_groups": G,
                    "bootstrap_iterations": n_boot,
                    "seed": seed
                })

            print(
                method,
                "| pooled ECE improvement:",
                round(float(improvement[0]), 6),
                "| 95% CI:",
                np.quantile(
                    boot_difference[:, 0],
                    [0.025, 0.975]
                )
            )

    pd.DataFrame(summary).to_csv(
        args.output_dir / "26group_bootstrap_per_predictor.csv",
        index=False
    )

    pd.DataFrame(paired).to_csv(
        args.output_dir / "26group_bootstrap_paired_differences.csv",
        index=False
    )

    pd.DataFrame(sectionrows).to_csv(
        args.output_dir / "26group_each_group_metrics.csv",
        index=False
    )

    pd.DataFrame(group_mapping).to_csv(
        args.output_dir / "26group_input_mapping.csv",
        index=False
    )

    print("\n26-GROUP BOOTSTRAP COMPLETED")
    print("Outputs:", args.output_dir.resolve())


if __name__ == "__main__":
    main()
