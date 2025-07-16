import argparse
import itertools
import os

import pandas as pd

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.data_loader import DataLoader
from encode4_gwas_analysis.utils.meta_helper import MetaHelper
from encode4_gwas_analysis.utils.path_helper import PathHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict["out_dir"]), f"Directory {arg_dict['out_dir']} does not exist."

    return arg_dict


def calculate_mean_caas(label_distr_df, label_caas_dict, sample_id):
    assert set(label_distr_df.columns) == set(label_caas_dict.keys())
    mean_caas_series = label_distr_df.mul(label_caas_dict).sum(axis=1) / label_distr_df.sum(axis=1)
    mean_caas_series.name = sample_id

    return mean_caas_series


def calculate_sample_metrics(sample_id, data_loader, snp_region_side_length, fake_snps):
    sample_metrics = {"mean_caas": None}
    sample_label_caas_dict = data_loader.load_sample_label_caas(sample_id=sample_id)
    sample_label_distr_df = data_loader.load_sample_snp_label_distribution(
        sample_id=sample_id, snp_region_side_length=snp_region_side_length, fake=fake_snps
    )
    mean_caas_series = calculate_mean_caas(
        label_distr_df=sample_label_distr_df,
        label_caas_dict=sample_label_caas_dict,
        sample_id=sample_id,
    )
    sample_metrics["mean_caas"] = mean_caas_series

    return sample_metrics


def create_metric_dfs(metric_series_dict):
    metric_df_dict = {}

    sample_ids = list(metric_series_dict.keys())
    metric_names = metric_series_dict[sample_ids[0]].keys()
    assert all(
        [set(metric_names) == set(metric_series_dict[sample_id]) for sample_id in sample_ids]
    ), "Not all samples have results for the same metrics."

    for metric_name in metric_names:
        metric_df = pd.concat(
            [metric_series_dict[sample_id][metric_name] for sample_id in sample_ids], axis=1
        )
        metric_df_dict[metric_name] = metric_df

    return metric_df_dict


if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    base_enrichment_dir = path_helper.enrichment_dir
    assert not os.path.isdir(base_enrichment_dir), (
        f"Directory {base_enrichment_dir} already exists."
    )
    os.mkdir(base_enrichment_dir)

    for snp_region_side_length, fake_snps in itertools.product(
        const.SNP_REGION_SIDE_LENGTHS, [True, False]
    ):
        sample_metrics = {}
        for sample_id in meta_helper.iter_samples(show_progress=True):
            sample_metrics[sample_id] = calculate_sample_metrics(
                sample_id=sample_id,
                data_loader=data_loader,
                snp_region_side_length=snp_region_side_length,
                fake_snps=fake_snps,
            )
        metric_df_dict = create_metric_dfs(metric_series_dict=sample_metrics)
        for metric_name, metric_df in metric_df_dict.items():
            metric_df_fpath = path_helper.get_snp_metric_df_fpath(
                metric_name=metric_name, window_size=snp_region_side_length, fake=fake_snps
            )
            metric_df.to_pickle(metric_df_fpath)
