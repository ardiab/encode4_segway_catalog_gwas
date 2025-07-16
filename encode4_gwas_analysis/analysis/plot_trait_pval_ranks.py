import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from encode4_gwas_analysis.utils.data_loader import DataLoader
from encode4_gwas_analysis.utils.meta_helper import MetaHelper
from encode4_gwas_analysis.utils.path_helper import PathHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = parser.parse_args()
    assert os.path.isdir(arg_dict.out_dir), f"Directory {arg_dict.out_dir} does not exist."

    return vars(arg_dict)


def plot_pval_ranks(
    sample_pval_series, trait_name, sig_threshold, ax, meta_helper, cap_pvals_at_1=True
):
    rank_counter = 1
    trait_pvals = sample_pval_series.sort_values(ascending=True)
    if cap_pvals_at_1:
        trait_pvals = trait_pvals.where(trait_pvals < 1, 1)
    sig_traits = trait_pvals[trait_pvals < sig_threshold]
    insig_traits = trait_pvals[trait_pvals >= sig_threshold]
    if len(sig_traits) > 0:
        ax.scatter(
            range(rank_counter, len(sig_traits) + rank_counter), sig_traits.values, color="red"
        )
        rank_counter += len(sig_traits)
    ax.scatter(
        range(rank_counter, len(insig_traits) + rank_counter), insig_traits.values, color="blue"
    )
    ax.set_xlabel("Biosample rank (ascending P-value)", fontsize=16, weight="bold")
    ax.set_ylabel("Corrected P-value", fontsize=16, weight="bold")
    ax.set_title(
        " ".join([w.capitalize() for w in trait_name.split(" ")]), fontsize=18, weight="bold"
    )
    ax.tick_params(axis="both", which="major", labelsize=18)
    red_circle = plt.Line2D(
        [0],
        [0],
        marker="o",
        color="white",
        markerfacecolor="red",
        markersize=10,
        label="Significant association",
    )
    blue_circle = plt.Line2D(
        [0],
        [0],
        marker="o",
        color="white",
        markerfacecolor="blue",
        markersize=10,
        label="Insignificant association",
    )

    top_n_ranks = 10
    annotation = "Top 10 associations (* indicates significance):\n"
    max_line_length = 0
    for i in range(top_n_ranks):
        sample_name = meta_helper.get_sample_description(trait_pvals.index.values[i]).split(
            "Homo sapiens "
        )[-1]
        sig = trait_pvals.iloc[i] < sig_threshold
        annotation += f"{i + 1}. {sample_name} {' (*)' if sig else ''}\n"
        max_line_length = max(max_line_length, len(sample_name))
    print(trait_name, max_line_length)
    ax.annotate(
        annotation,
        color="black",
        xy=(150, 0.05),
        xytext=(233, 0.055),
        fontsize=12,
        horizontalalignment="right",
    )
    ax.set_xlim([1, 235])
    if not cap_pvals_at_1:
        ax.set_ylim([0, 234])
    else:
        ax.set_ylim([0, 1])
    ax.axhline(sig_threshold, linestyle="--", color="black")
    ax.legend(handles=[red_circle, blue_circle])


if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    ALL_TRAITS = False
    WINDOW_SIZE = 10_000
    SIG_THRESH = 0.05
    MIN_N_SNPS = 30
    pval_df = (
        data_loader.load_pval_df(
            metric_name="mean_caas",
            test_type="signedrank_per_sample_null",
            window_size=WINDOW_SIZE,
            fake=False,
        ).fillna(1)
        * 234
    )

    if ALL_TRAITS:
        trait_snp_count_df = data_loader.load_trait_snp_count_df(window_size=WINDOW_SIZE)
        sorted_trait_list = (
            trait_snp_count_df[trait_snp_count_df["n_snps"] >= MIN_N_SNPS]
            .sort_values(by=["n_snps"], ascending=False)["trait"]
            .unique()
            .tolist()
        )
        n_cols = 10
        n_rows = 10
        split_point = len(sorted_trait_list) % (n_cols * n_rows)
        trait_groups = list(
            np.array_split(
                sorted_trait_list[:-split_point],
                (len(sorted_trait_list) - split_point) / (n_cols * n_rows),
            )
        ) + [sorted_trait_list[-split_point:]]
        for fig_idx, trait_group in tqdm(enumerate(trait_groups), total=len(trait_groups)):
            trait_min_n_snps = trait_snp_count_df[trait_snp_count_df["trait"].isin(trait_group)][
                "n_snps"
            ].min()
            trait_max_n_snps = trait_snp_count_df[trait_snp_count_df["trait"].isin(trait_group)][
                "n_snps"
            ].max()
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 10, n_rows * 10), dpi=300)
            for idx, trait in enumerate(trait_group):
                ax = axes[np.unravel_index(idx, axes.shape)]
                rank_counter = 1
                trait_pvals = pval_df.loc[trait].sort_values(ascending=True)
                sig_traits = trait_pvals[trait_pvals < SIG_THRESH]
                insig_traits = trait_pvals[trait_pvals >= SIG_THRESH]
                if len(sig_traits) > 0:
                    ax.scatter(
                        range(rank_counter, len(sig_traits) + rank_counter),
                        sig_traits.values,
                        color="red",
                    )
                    rank_counter += len(sig_traits)
                ax.scatter(
                    range(rank_counter, len(insig_traits) + rank_counter),
                    insig_traits.values,
                    color="blue",
                )
                ax.set_xlabel("Biosample rank (ascending P-value)", fontsize=16, weight="bold")
                ax.set_ylabel("Corrected P-value", fontsize=16, weight="bold")
                ax.set_title(trait)
                ax.tick_params(axis="both", which="major", labelsize=18)
                red_circle = plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="white",
                    markerfacecolor="red",
                    markersize=10,
                    label="Significant association",
                )
                blue_circle = plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="white",
                    markerfacecolor="blue",
                    markersize=10,
                    label="Insignificant association",
                )

                top_n_ranks = 10
                rank_start = 1
                rank_end = rank_start + top_n_ranks
                max_pval = trait_pvals.iloc[rank_end - 1]
                ax.plot([rank_start, rank_end], [max_pval + 5, max_pval + 5], color="black")
                ax.plot(
                    [rank_start + 0.6, rank_start + 0.6],
                    [max_pval + 5, max_pval + 3],
                    color="black",
                )
                ax.plot([rank_end, rank_end], [max_pval + 5, max_pval + 3], color="black")
                annotation = ""
                sig_annotation = ""
                insig_annotation = ""
                for i in range(top_n_ranks):
                    sample_name = meta_helper.get_sample_description(
                        trait_pvals.index.values[i]
                    ).split("Homo sapiens ")[-1]
                    if trait_pvals.iloc[i] < SIG_THRESH:
                        sig_annotation += f"{i + 1}. {sample_name}\n"
                    else:
                        insig_annotation += f"{i + 1}. {sample_name}\n"
                ax.annotate(
                    "",
                    xy=((rank_start + rank_end) / 2, max_pval + 5),
                    xytext=((rank_start + rank_end) / 2, max_pval + 15),
                    arrowprops=dict(arrowstyle="->", color="black"),
                    color="black",
                )
                if len(insig_annotation) > 0:
                    n_insig = len([a for a in insig_annotation.split("\n") if len(a) > 0])
                    ax.annotate(
                        insig_annotation,
                        color="blue",
                        xy=((rank_start + rank_end) / 2 + 2, max_pval),
                        xytext=(rank_start + 3, max_pval + 10),
                    )
                else:
                    n_insig = 0
                ax.annotate(
                    sig_annotation,
                    color="red",
                    xy=((rank_start + rank_end) / 2 + 2, max_pval),
                    xytext=(rank_start + 3, max_pval + 10 + (7 * n_insig)),
                )
                ax.set_xlim([1, 235])
                ax.set_ylim([0, 234])
                ax.legend(handles=[red_circle, blue_circle])

            fig.savefig(
                os.path.join(
                    arg_dict["out_dir"], f"{fig_idx}_{trait_min_n_snps}_{trait_max_n_snps}.pdf"
                ),
                format="pdf",
            )
            fig.savefig(
                os.path.join(
                    arg_dict["out_dir"], f"{fig_idx}_{trait_min_n_snps}_{trait_max_n_snps}.png"
                )
            )
    else:
        trait_subset = ["colorectal cancer", "coronary artery disease", "bipolar disorder"]
        trait_pval_dfs = []
        for trait in trait_subset:
            fig, ax = plt.subplots(figsize=(10, 10), dpi=400)
            trait_pvals = pval_df.loc[trait]
            trait_pval_df = (
                trait_pvals.sort_values(ascending=True)
                .reset_index()
                .rename(columns={trait: "pval"})
            )
            trait_pval_df["sample"] = trait_pval_df["sample"].apply(
                lambda s: meta_helper.get_sample_description(s).split("Homo sapiens ")[-1]
            )
            # trait_pval_df['trait'] = trait
            trait_pval_df["sig"] = trait_pval_df["pval"] < 0.05
            trait_pval_dfs.append(trait_pval_df)
            trait_pval_df.to_csv(
                os.path.join(arg_dict["out_dir"], f"{trait}_sample_pvals.csv"), index=False
            )
            plot_pval_ranks(
                sample_pval_series=trait_pvals,
                trait_name=trait,
                sig_threshold=0.05,
                ax=ax,
                meta_helper=meta_helper,
                cap_pvals_at_1=True,
            )
            fig.savefig(os.path.join(arg_dict["out_dir"], f"{trait}_rank_plots.pdf"), format="pdf")
        trait_pval_df_full = pd.concat(trait_pval_dfs, ignore_index=True)
        trait_pval_df_full.to_csv(
            os.path.join(arg_dict["out_dir"], "selected_trait_sample_pvals.csv"), index=False
        )
