import argparse
import itertools
import os

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from tqdm import tqdm

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.data_loader import DataLoader
from encode4_gwas_analysis.utils.path_helper import PathHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--snp_cutoff", default=30, type=int)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict["out_dir"]), f"Directory {arg_dict['out_dir']} does not exist."
    assert arg_dict["snp_cutoff"] > 0, "Got negative SNP cutoff."

    return arg_dict


def _get_linkage_idxs(snp_metric_df, trait_snp_linkage_df, fake):
    linkage_idxs = {}
    for trait, trait_linkage_df in trait_snp_linkage_df.groupby("trait"):
        trait_snp_names = [
            snp + "_fake" if fake else snp for snp in trait_linkage_df["snp"].unique()
        ]
        linkage_idxs[trait] = np.argwhere(snp_metric_df.index.isin(trait_snp_names)).flatten()

    return linkage_idxs


def _calculate_signedrank_pvals_uniform_null(snp_metric_df, trait_snp_linkage_df, fake):
    data = []
    linkage_idxs = _get_linkage_idxs(
        snp_metric_df=snp_metric_df, trait_snp_linkage_df=trait_snp_linkage_df, fake=fake
    )

    ranks_across_biosamples_df = snp_metric_df.rank(axis=1)
    rank_array = ranks_across_biosamples_df.values
    median_rank = (snp_metric_df.shape[1] + 1) / 2

    for trait_name, trait_snp_idxs in tqdm(linkage_idxs.items()):
        for sample_idx, sample_id in enumerate(snp_metric_df.columns):
            trait_sample_ranks = rank_array[trait_snp_idxs, sample_idx]
            exp_median = np.median(trait_sample_ranks)
            test_stat, pval = wilcoxon(trait_sample_ranks - median_rank, alternative="greater")
            data.append(
                {
                    "trait": trait_name,
                    "test_stat": test_stat,
                    "exp_median": exp_median,
                    "pval": pval,
                    "sample": sample_id,
                }
            )

    pval_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="pval")
    test_stat_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="test_stat")
    median_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="exp_median")

    return test_stat_pivot, pval_pivot, median_pivot


def _calculate_signedrank_pvals_sample_specific_null(snp_metric_df, trait_snp_linkage_df, fake):
    data = []
    linkage_idxs = _get_linkage_idxs(
        snp_metric_df=snp_metric_df, trait_snp_linkage_df=trait_snp_linkage_df, fake=fake
    )

    ranks_across_biosamples_df = snp_metric_df.rank(axis=1)
    rank_array = ranks_across_biosamples_df.values

    sample_medians = ranks_across_biosamples_df.median(axis=0)

    for trait_name, trait_snp_idxs in tqdm(linkage_idxs.items()):
        for sample_idx, sample_id in enumerate(snp_metric_df.columns):
            trait_sample_ranks = rank_array[trait_snp_idxs, sample_idx]
            exp_median = np.median(trait_sample_ranks)
            test_stat, pval = wilcoxon(
                trait_sample_ranks - sample_medians[sample_id], alternative="greater"
            )
            data.append(
                {
                    "trait": trait_name,
                    "test_stat": test_stat,
                    "exp_median": exp_median,
                    "pval": pval,
                    "sample": sample_id,
                }
            )

    pval_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="pval")
    test_stat_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="test_stat")
    median_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="exp_median")

    return test_stat_pivot, pval_pivot, median_pivot


def _calculate_signedrank_pvals_normalized_within_sample(snp_metric_df, trait_snp_linkage_df, fake):
    data = []
    linkage_idxs = _get_linkage_idxs(
        snp_metric_df=snp_metric_df, trait_snp_linkage_df=trait_snp_linkage_df, fake=fake
    )

    ranks_within_biosamples_df = snp_metric_df.rank(axis=0)
    ranks_across_biosamples_df = ranks_within_biosamples_df.rank(axis=1)
    rank_array = ranks_across_biosamples_df.values
    median_rank = (snp_metric_df.shape[1] + 1) / 2

    for trait_name, trait_snp_idxs in tqdm(linkage_idxs.items()):
        for sample_idx, sample_id in enumerate(snp_metric_df.columns):
            trait_sample_ranks = rank_array[trait_snp_idxs, sample_idx]
            exp_median = np.median(trait_sample_ranks)
            test_stat, pval = wilcoxon(trait_sample_ranks - median_rank, alternative="greater")
            data.append(
                {
                    "trait": trait_name,
                    "test_stat": test_stat,
                    "exp_median": exp_median,
                    "pval": pval,
                    "sample": sample_id,
                }
            )

    pval_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="pval")
    test_stat_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="test_stat")
    median_pivot = pd.DataFrame(data).pivot(index="trait", columns="sample", values="exp_median")

    return test_stat_pivot, pval_pivot, median_pivot


def calculate_signedrank_pvals(snp_metric_df, trait_snp_linkage_df, fake):
    assert fake in [True, False]
    test_stats_uniform_null, pvals_uniform_null, medians_uniform_null = (
        _calculate_signedrank_pvals_uniform_null(
            snp_metric_df=snp_metric_df, trait_snp_linkage_df=trait_snp_linkage_df, fake=fake
        )
    )

    return {
        "signedrank_uniform_null": {
            "pval": pvals_uniform_null,
            "test_stat": test_stats_uniform_null,
            "effect_size": medians_uniform_null,
        }
    }


def calculate_sample_pvals(
    metric_name, snp_region_side_length, fake_snps, data_loader, trait_snp_linkage_df
):
    metric_df = data_loader.load_snp_metric_df(
        metric_name=metric_name, snp_region_side_length=snp_region_side_length, fake=fake_snps
    )
    metric_df_snp_names = (
        metric_df.index.tolist()
        if not fake_snps
        else [i.replace("_fake", "") for i in metric_df.index.tolist()]
    )
    assert set(trait_snp_linkage_df["snp"].values).issubset(set(metric_df_snp_names))

    unique_snps_to_keep = trait_snp_linkage_df["snp"].unique().tolist()
    if fake_snps:
        unique_snps_to_keep = [snp + "_fake" for snp in unique_snps_to_keep]
    metric_df_sub = metric_df.loc[unique_snps_to_keep]
    test_result_dict = calculate_signedrank_pvals(
        snp_metric_df=metric_df_sub, trait_snp_linkage_df=trait_snp_linkage_df, fake=fake_snps
    )

    return test_result_dict


def get_filtered_trait_snps(snp_region_side_length, snp_cutoff, data_loader):
    trait_snp_linkage_df = data_loader.load_trait_snp_linkage_df(window_size=snp_region_side_length)
    trait_snp_count_df = data_loader.load_trait_snp_count_df(window_size=snp_region_side_length)
    traits_to_keep = (
        trait_snp_count_df[trait_snp_count_df["n_snps"] >= snp_cutoff]["trait"].unique().tolist()
    )
    trait_snp_linkage_df_sub = trait_snp_linkage_df[
        trait_snp_linkage_df["trait"].isin(traits_to_keep)
    ]

    sorted_traits_fwd = (
        trait_snp_count_df[trait_snp_count_df["trait"].isin(traits_to_keep)]
        .sort_values(by=["n_snps"], ascending=False)["trait"]
        .tolist()
    )
    sorted_traits_rev = (
        trait_snp_count_df[trait_snp_count_df["trait"].isin(traits_to_keep)]
        .sort_values(by=["n_snps"], ascending=True)["trait"]
        .tolist()
    )
    results = {}
    for direction_traits, direction_desc in [
        (sorted_traits_fwd, "fwd"),
        (sorted_traits_rev, "rev"),
    ]:
        added_traits = []
        added_snps = set()
        overlap_data = []
        for trait in tqdm(direction_traits):
            trait_df = trait_snp_linkage_df_sub[trait_snp_linkage_df_sub["trait"] == trait]
            trait_snps = set(trait_df["snp"].tolist())
            trait_overlap = len(trait_snps.intersection(added_snps)) / len(trait_snps)
            overlap_data.append({"trait": trait, "overlap": trait_overlap})
            added_snps = added_snps.union(trait_snps)
            added_traits.append(trait)
        results[direction_desc] = {
            "traits": added_traits,
            "snps": added_snps,
            "overlap_df": pd.DataFrame(overlap_data),
        }

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.hist(results["fwd"]["overlap_df"]["overlap"].values, label="Descending", alpha=0.6)
    ax.hist(results["rev"]["overlap_df"]["overlap"].values, label="Ascending", alpha=0.6)
    ax.set_xlabel("Proportion SNP overlap with previously added traits")
    ax.set_ylabel("Frequency")
    ax.set_title("Trait filtering by SNP overlap")
    ax.legend()
    fig.show()

    return trait_snp_linkage_df_sub


if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    for metric_name in ["mean_caas"]:
        for snp_region_side_length, fake_snps in itertools.product(
            const.SNP_REGION_SIDE_LENGTHS, [True, False]
        ):
            trait_snp_linkage_df = get_filtered_trait_snps(
                snp_region_side_length=snp_region_side_length,
                snp_cutoff=arg_dict["snp_cutoff"],
                data_loader=data_loader,
            )

            metric_test_result_dict = calculate_sample_pvals(
                metric_name=metric_name,
                snp_region_side_length=snp_region_side_length,
                fake_snps=fake_snps,
                data_loader=data_loader,
                trait_snp_linkage_df=trait_snp_linkage_df,
            )
            for test_type, test_df_dict in metric_test_result_dict.items():
                pval_df_fpath = path_helper.get_pval_df_fpath(
                    window_size=snp_region_side_length,
                    metric_name=metric_name,
                    test_type=test_type,
                    fake=fake_snps,
                )
                test_df_dict["pval"].to_pickle(pval_df_fpath)
                test_stat_df_fpath = path_helper.get_test_stat_df_fpath(
                    window_size=snp_region_side_length,
                    metric_name=metric_name,
                    test_type=test_type,
                    fake=fake_snps,
                )
                test_df_dict["test_stat"].to_pickle(test_stat_df_fpath)
                effect_size_df_fpath = path_helper.get_effect_size_df_fpath(
                    window_size=snp_region_side_length,
                    metric_name=metric_name,
                    test_type=test_type,
                    fake=fake_snps,
                )
                test_df_dict["effect_size"].to_pickle(effect_size_df_fpath)
