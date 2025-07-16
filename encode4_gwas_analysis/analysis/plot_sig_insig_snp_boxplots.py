import argparse
import math
import os

import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tqdm

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.data_loader import DataLoader
from encode4_gwas_analysis.utils.meta_helper import MetaHelper
from encode4_gwas_analysis.utils.path_helper import PathHelper

sns.set_style("whitegrid")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = parser.parse_args()
    assert os.path.isdir(arg_dict.out_dir), f"Directory {arg_dict.out_dir} does not exist."

    return vars(arg_dict)


def create_barplot(metric_df, label_col, metric_col, sig_col, ylabel):
    fig, ax = plt.subplots(figsize=(20, 9), dpi=400)
    # sns.barplot(data=metric_df, x='label', y='coverage', hue='sig')
    label_order = {k: idx for idx, k in enumerate(const.LABEL_COLOR_MAP.keys())}
    ordered_interp_terms = sorted(metric_df[label_col].unique(), key=lambda l: label_order[l])
    bin_width = 0.3
    hatches = [None, "/", "\\"]
    snp_types = ["significant_real", "insignificant_real", "all_fake"]
    total_width = len(snp_types) * bin_width
    start_idx = math.ceil(total_width / 2)
    idx_stride = math.ceil(total_width)
    end_idx = start_idx + (idx_stride * (len(ordered_interp_terms) - 1))
    interp_term_x_locs = np.arange(start_idx, end_idx + idx_stride, idx_stride)
    for interp_term_center_x, interp_term in zip(interp_term_x_locs, ordered_interp_terms):
        # Add (bin_width / 2) to the upper end so it is endpoint inclusive. Added (bin_width / 2) instead of (bin_width)
        # because adding full bin width can result in an additional bin depending on float representation (e.g. if the
        # lower end is 1.7, the upper end is 2.3000000000000003).
        bar_ticks = np.arange(
            interp_term_center_x - (total_width / 2) + (bin_width / 2),
            interp_term_center_x + (total_width / 2) - (bin_width / 2) + (bin_width / 2),
            bin_width,
        )
        assert len(bar_ticks) == len(snp_types) == len(hatches)
        label_color = const.LABEL_COLOR_MAP[interp_term]
        for snp_type_bar_x, snp_type, snp_hatch in zip(bar_ticks, snp_types, hatches):
            snp_freq_df = metric_df[
                (metric_df[label_col] == interp_term) & (metric_df[sig_col] == snp_type)
            ]
            snp_freq = snp_freq_df[metric_col].values[0]
            assert snp_freq_df.shape[0] == 1
            ax.bar(
                snp_type_bar_x,
                snp_freq,
                color=label_color,
                edgecolor="black",
                width=bin_width,
                hatch=snp_hatch,
            )

        legend_patches = []
        legend_color = (0.8313725490196079, 0.8313725490196079, 0.8313725490196079)
        legend_descs = {
            "significant_real": "Associated trait-sample pairs",
            "insignificant_real": "Non-associated trait-sample pairs",
            "all_fake": "All fake SNPs",
        }
        for snp_type, snp_hatch in zip(snp_types, hatches):
            legend_patches.append(
                matplotlib.patches.Patch(
                    facecolor=legend_color,
                    label=legend_descs[snp_type],
                    hatch=snp_hatch * 3 if snp_hatch is not None else None,
                    edgecolor="black",
                    linewidth=1.5,
                )
            )

    ax.legend(handles=legend_patches, prop={"size": 18})
    # assert len(interp_term_x_locs) == len(ax.get_xticks())
    ax.set_xticks(interp_term_x_locs)
    ax.set_xticklabels(ordered_interp_terms, rotation=90)
    ax.tick_params(axis="both", which="major", labelsize=18)
    # ax.set_xlabel('SAGA Label', fontsize=18, weight='bold')
    ax.set_ylabel(ylabel, fontsize=18, weight="bold")
    # ax.set_title(title)
    fig.tight_layout()

    return fig, ax


if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    WINDOW_SIZE = 10_000
    SNP_CUTOFF = 30
    SIG_THRESHOLD = 0.05

    trait_linkage_df = data_loader.load_trait_snp_linkage_df(window_size=WINDOW_SIZE)
    trait_snp_count_df = data_loader.load_trait_snp_count_df(window_size=WINDOW_SIZE)

    traits_to_keep = trait_snp_count_df[trait_snp_count_df["n_snps"] >= SNP_CUTOFF][
        "trait"
    ].unique()
    snps_to_keep = trait_linkage_df[trait_linkage_df["trait"].isin(traits_to_keep)]["snp"].unique()
    pval_df = (
        data_loader.load_pval_df(
            metric_name="mean_caas",
            test_type="signedrank_per_sample_null",
            window_size=WINDOW_SIZE,
            fake=False,
        ).fillna(1)
        * 234
    )
    assert pval_df.shape[0] == len(traits_to_keep)

    pval_df = pval_df.loc[traits_to_keep]

    all_sample_label_coverages_concat = {
        "significant": {interp_term: [] for interp_term in const.LABEL_COLOR_MAP.keys()},
        "insignificant": {interp_term: [] for interp_term in const.LABEL_COLOR_MAP.keys()},
    }
    all_sample_label_coverages_means = {
        "significant": {interp_term: [] for interp_term in const.LABEL_COLOR_MAP.keys()},
        "insignificant": {interp_term: [] for interp_term in const.LABEL_COLOR_MAP.keys()},
    }

    coverage_data_concat = []
    coverage_data_mean = []
    enrichment_data_concat = []
    enrichment_data_mean = []

    for sample_id in tqdm.tqdm(pval_df.columns):
        sample_label_distrs = data_loader.load_sample_snp_label_distribution(
            sample_id=sample_id, snp_region_side_length=WINDOW_SIZE, fake=False
        )
        sample_fake_label_distrs = data_loader.load_sample_snp_label_distribution(
            sample_id=sample_id, snp_region_side_length=WINDOW_SIZE, fake=True
        )
        sample_label_genome_coverage = data_loader.load_sample_label_coverage(
            sample_id=sample_id, normalized=True
        )

        sample_pvals = pval_df[sample_id]
        sample_sig_traits = sample_pvals[sample_pvals < SIG_THRESHOLD].index.unique()
        sample_insig_traits = sample_pvals[sample_pvals >= SIG_THRESHOLD].index.unique()
        assert (len(sample_sig_traits) + len(sample_insig_traits)) == 1274

        sample_sig_snps = (
            trait_linkage_df[trait_linkage_df["trait"].isin(sample_sig_traits)]["snp"]
            .unique()
            .tolist()
        )
        sample_insig_snps = trait_linkage_df[trait_linkage_df["trait"].isin(sample_insig_traits)][
            "snp"
        ].unique()
        # Some SNPs may be in both significant and insignificant traits. Filter them from insignificant (only inlcude
        # in significant).
        sample_insig_snps = list(set(sample_insig_snps).difference(set(sample_sig_snps)))

        # A couple of SNPs are missing from some samples, because they were in unannotated regions.
        sample_sig_snps_proc = list(
            set(sample_sig_snps).intersection(set(sample_label_distrs.index.values))
        )
        sample_insig_snps_proc = list(
            set(sample_insig_snps).intersection(set(sample_label_distrs.index.values))
        )
        assert (len(sample_sig_snps) - len(sample_sig_snps_proc)) < 10
        assert (len(sample_insig_snps) - len(sample_insig_snps_proc)) < 10
        sample_sig_snps_proc_fake = list(
            set([snp + "_fake" for snp in sample_sig_snps]).intersection(
                set(sample_fake_label_distrs.index.values)
            )
        )
        sample_insig_snps_proc_fake = list(
            set([snp + "_fake" for snp in sample_insig_snps]).intersection(
                set(sample_fake_label_distrs.index.values)
            )
        )

        sig_real_snp_coverages = sample_label_distrs.loc[sample_sig_snps_proc]
        insig_real_snp_coverages = sample_label_distrs.loc[sample_insig_snps_proc]
        sig_fake_snp_coverages = sample_fake_label_distrs.loc[sample_sig_snps_proc_fake]
        insig_fake_snp_coverages = sample_fake_label_distrs.loc[sample_insig_snps_proc_fake]

        sample_label_mapping = {
            interp_term: [] for interp_term in list(const.LABEL_COLOR_MAP.keys())
        }
        for c in sample_label_distrs.columns:
            label_interp_term = c.split("_")[-1]
            sample_label_mapping[label_interp_term].append(c)

        for interp_term in list(const.LABEL_COLOR_MAP.keys()):
            interp_term_labels = sample_label_mapping[interp_term]
            interp_term_genome_coverage = sum(
                [sample_label_genome_coverage[label] for label in interp_term_labels]
            )
            for coverage_df, sig_desc, real in zip(
                [
                    sig_real_snp_coverages,
                    insig_real_snp_coverages,
                    sig_fake_snp_coverages,
                    insig_fake_snp_coverages,
                    pd.concat([sig_fake_snp_coverages, insig_fake_snp_coverages]),
                ],
                ["significant", "insignificant", "significant", "insignificant", "all"],
                ["real", "real", "fake", "fake", "fake"],
            ):
                interp_term_snp_coverages = coverage_df[interp_term_labels].sum(
                    axis=1
                ) / coverage_df.sum(axis=1)
                interp_term_enrichments = interp_term_snp_coverages / interp_term_genome_coverage
                proc_sig_desc = sig_desc + "_" + real
                enrichment_data_mean.append(
                    {
                        "label": interp_term,
                        "enrichment": interp_term_enrichments.mean(),
                        "sig": proc_sig_desc,
                    }
                )

    df_enrichment_mean = (
        pd.DataFrame(enrichment_data_mean)
        .groupby(["label", "sig"])["enrichment"]
        .mean()
        .reset_index()
    )
    barplot_enrichment_df = (
        df_enrichment_mean.groupby(["label", "sig"])["enrichment"].mean().reset_index()
    )
    barplot_enrichment_df["enrichment"] = np.log2(barplot_enrichment_df["enrichment"])
    barplot_enrichment_fig, barplot_enrichment_ax = create_barplot(
        metric_df=barplot_enrichment_df,
        label_col="label",
        metric_col="enrichment",
        sig_col="sig",
        ylabel="Enrichment within SNP Region\n$log_2 (observed/expected)$",
    )
    barplot_enrichment_ax.axhline(0, color="black")

    barplot_enrichment_fig.savefig(
        os.path.join(arg_dict["out_dir"], "snp_label_enrichment_barplot_with_fake.pdf"),
        format="pdf",
    )
    barplot_enrichment_fig.savefig(
        os.path.join(arg_dict["out_dir"], "snp_label_enrichment_barplot_with_fake.png"),
        format="png",
    )
