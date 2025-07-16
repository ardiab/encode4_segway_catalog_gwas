"""Create a clustered heatmap of Segway interpretation term CAAS values."""

import argparse
import os

import encode_gwas.utils.constants as const
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
        for label, caas in sample_label_caas.items():
            label_caas_data.append(
                {
                    "sample": sample_id,
                    "label": label,
                    "interp_term": label.split("_")[-1],
                    "caas": caas,
                }
            )

    label_caas_df = pd.DataFrame(label_caas_data)
    label_caas_df["rank"] = label_caas_df.groupby(["sample", "interp_term"])["caas"].rank()
    label_caas_df["interp_term"] = pd.Categorical(label_caas_df["interp_term"])

    n_rows_per_interp_term = 60
    heatmap_data = []
    heatmap_array = None
    for interp_term_idx, interp_term in enumerate(const.LABEL_COLOR_MAP.keys()):
        interp_term_sub = label_caas_df[label_caas_df["interp_term"] == interp_term]
        max_n_labels = interp_term_sub.groupby("sample").size().max()
        for sample_id, sample_df in interp_term_sub.groupby("sample"):
            sample_n_interp_terms = sample_df.shape[0]
            if sample_n_interp_terms > 0:
                n_rows_per_integer_label = n_rows_per_interp_term / sample_n_interp_terms
                assert float(n_rows_per_integer_label).is_integer()
                n_rows_per_integer_label = int(n_rows_per_integer_label)
                n_times_to_repeat_label = max_n_labels / sample_df.shape[0]
                label_idx = 0
                row_idx = 0
                for idx, row in sample_df.sort_values(by=["caas"], ascending=False).iterrows():
                    for i in range(n_rows_per_integer_label):
                        heatmap_data.append(
                            {
                                "sample": sample_id,
                                "interp_term": interp_term,
                                "label": row["label"],
                                "caas": row["caas"],
                                "label_idx": label_idx,
                                "row": row_idx,
                            }
                        )
                        row_idx += 1
                    label_idx += 1
            else:
                n_rows_per_integer_label = n_rows_per_interp_term
                assert float(n_rows_per_integer_label).is_integer()
                n_rows_per_integer_label = int(n_rows_per_integer_label)
                for i in range(n_rows_per_integer_label):
                    heatmap_data.append(
                        {
                            "sample": sample_id,
                            "interp_term": interp_term,
                            "label": "N/A",
                            "caas": np.nan,
                            "label_idx": 0,
                            "row": i,
                        }
                    )

    heatmap_df = pd.DataFrame(heatmap_data)
    heatmap_pivot = heatmap_df.pivot(
        index=["interp_term", "row"], columns=["sample"], values=["caas"]
    )
    label_order = [k for k in const.LABEL_COLOR_MAP.keys()]
    heatmap_pivot = heatmap_pivot.fillna(-1).loc[label_order]
    min_caas = heatmap_df["caas"].min()
    max_caas = heatmap_df["caas"].max()

    clustermap = sns.clustermap(heatmap_pivot.values, dendrogram_ratio=(0.05, 0.05), cmap="magma")
    plt.close(clustermap.figure)
    row_order = clustermap.dendrogram_row.reordered_ind
    col_order = clustermap.dendrogram_col.reordered_ind

    fig, ax = plt.subplots(figsize=(40, 35), dpi=400)
    cmap = plt.get_cmap("magma").copy()
    cmap.set_under((0.5882352941176471, 0.5882352941176471, 0.5882352941176471))
    sns.heatmap(
        heatmap_pivot.iloc[:, col_order],
        ax=ax,
        vmin=min_caas,
        vmax=max_caas,
        cmap=cmap,
        rasterized=True,
        cbar_kws={"shrink": 0.5, "use_gridspec": False, "location": "top"},
    )
    ax.collections[0].colorbar.ax.tick_params(labelsize=17)
    ax.set_yticks(np.arange(30, 660, 60))
    ax.set_yticklabels(label_order, fontsize=25, weight="bold")
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=True)
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.savefig(os.path.join(arg_dict["out_dir"], "label_caas_heatmap_clustered.png"))
    fig.savefig(os.path.join(arg_dict["out_dir"], "label_caas_heatmap_clustered.pdf"), format="pdf")
