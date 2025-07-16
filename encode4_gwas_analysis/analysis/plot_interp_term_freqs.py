"""Visualize the frequency of interpretation terms across samples."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.data_loader import DataLoader
from encode4_gwas_analysis.utils.meta_helper import MetaHelper
from encode4_gwas_analysis.utils.path_helper import PathHelper


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = parser.parse_args()
    assert os.path.isdir(arg_dict.out_dir), f"Directory {arg_dict.out_dir} does not exist."

    return vars(arg_dict)


if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    label_caas_data = []

    for sample_id in meta_helper.iter_samples(show_progress=True):
        sample_label_caas = data_loader.load_sample_label_caas(sample_id=sample_id)
        sample_label_coverage = data_loader.load_sample_label_coverage(
            sample_id=sample_id, normalized=True
        )
        for label, caas in sample_label_caas.items():
            label_caas_data.append(
                {
                    "sample": sample_id,
                    "label": label,
                    "interp_term": label.split("_")[-1],
                    "caas": caas,
                    "coverage": sample_label_coverage[label],
                }
            )

    label_caas_df = pd.DataFrame(label_caas_data)
    label_caas_df["rank"] = label_caas_df.groupby(["sample", "interp_term"])["caas"].rank()
    label_caas_df["interp_term"] = pd.Categorical(label_caas_df["interp_term"])

    fig, ax = plt.subplots(figsize=(10, 10), dpi=400)
    interp_term_overall_counts = (
        label_caas_df.groupby(["sample"])["interp_term"].value_counts().value_counts()
    )
    ax.bar(
        interp_term_overall_counts.index,
        interp_term_overall_counts.values,
        linewidth=1.5,
        edgecolor="black",
    )
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_xlabel(
        "# occurrences of interpretation terms\n(for every interpretation term, sample pair)",
        fontsize=16,
        weight="bold",
    )
    ax.set_ylabel("Frequency", fontsize=16, weight="bold")
    ax.set_title("Interpretation term occurrences per sample", fontsize=20, weight="bold")
    fig.savefig(
        os.path.join(arg_dict["out_dir"], "snp_interp_term_counts_per_sample_overall.pdf"),
        format="pdf",
    )
    fig.savefig(os.path.join(arg_dict["out_dir"], "snp_interp_term_counts_per_sample_overall.png"))

    fig, axes = plt.subplots(3, 4, figsize=(40, 30), dpi=400)
    label_counts = label_caas_df.groupby(["sample"])["interp_term"].value_counts().reset_index()
    interp_coverages = (
        label_caas_df.groupby(["sample", "interp_term"])["coverage"].sum().reset_index()
    )
    i = 0
    for interp_term, color in const.LABEL_COLOR_MAP.items():
        count_sub = label_counts[label_counts["interp_term"] == interp_term]
        coverage_map = interp_coverages[interp_coverages["interp_term"] == interp_term]
        assert coverage_map["sample"].nunique() == coverage_map.shape[0] == 234
        coverage_map = coverage_map.set_index("sample")["coverage"].to_dict()
        count_sub["coverage"] = count_sub["sample"].map(coverage_map)
        assert count_sub["coverage"].notnull().all()
        count_val_counts = count_sub["count"].value_counts()
        count_coverages = count_sub.groupby("count")["coverage"].mean()
        ax = axes[np.unravel_index(i, axes.shape)]
        ax.bar(
            count_val_counts.index.values,
            count_val_counts.values,
            facecolor=color,
            label=interp_term,
            linewidth=1.5,
            edgecolor="black",
        )
        ax.set_xlim([-0.69, 5.69])
        ax.set_xticks(range(6))

        # twiny to show coverage as function of count
        ax_twin = ax.twinx()
        ax_twin.plot(
            count_coverages.index.values, count_coverages.values, "X", color="black", markersize=16
        )
        for h in count_coverages.values:
            ax_twin.axhline(h, linestyle="--", alpha=0.3, color="gray")
        ax_twin.tick_params(axis="both", which="major", labelsize=18)
        ax_twin.set_ylabel("Mean coverage", fontsize=16, weight="bold")

        ax.tick_params(axis="both", which="major", labelsize=18)
        ax.set_xlabel("# Associated Segway integer labels", fontsize=16, weight="bold")
        ax.set_ylabel("# Biosamples", fontsize=16, weight="bold")
        ax.set_title(interp_term, fontsize=20, weight="bold")
        # ax.legend()
        i += 1

    fig.tight_layout()
    fig.savefig(
        os.path.join(arg_dict["out_dir"], "snp_interp_term_counts_per_sample_per_term.pdf"),
        format="pdf",
    )
    fig.savefig(os.path.join(arg_dict["out_dir"], "snp_interp_term_counts_per_sample_per_term.png"))
