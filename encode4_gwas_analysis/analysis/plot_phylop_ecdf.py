import argparse
import json
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.data_loader import DataLoader


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = parser.parse_args()
    assert os.path.isdir(arg_dict.out_dir), f"Directory {arg_dict.out_dir} does not exist."

    return vars(arg_dict)


if __name__ == "__main__":
    arg_dict = parse_args()
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    ecdf_plot_json_data_fpath = os.path.join(
        arg_dict["out_dir"], "ENCSR067NPP_label_phylop_ecdf_plot_data.json"
    )
    with open(ecdf_plot_json_data_fpath, "r") as in_f:
        ecdf_data_dict = json.load(in_f)

    label_caas_vals = data_loader.load_sample_label_caas(sample_id="ENCSR067NPP")
    label_caas_vals = list(sorted([(k, v) for k, v in label_caas_vals.items()], key=lambda t: t[1]))

    barplot_fig, barplot_ax = plt.subplots(figsize=(7, 5), dpi=400)
    width = 0.5
    x_ticks = np.arange(0.5, 0.5 * (len(label_caas_vals) + 1), 0.5)
    for idx, (label, caas) in enumerate(label_caas_vals):
        interp_term = label.split("_")[-1]
        label_color = (
            const.LABEL_COLOR_MAP[interp_term] if interp_term != "Quiescent" else (0, 0, 0)
        )
        barplot_ax.bar(x_ticks[idx], caas, color=label_color, width=width, edgecolor="black")
    barplot_ax.set_xticks(x_ticks)
    barplot_ax.set_xticklabels([t[0] for t in label_caas_vals], rotation=90, weight="bold")
    barplot_ax.tick_params(axis="y", which="both", labelsize=14)
    barplot_ax.set_ylabel("Label CAAS", fontsize=16, weight="bold")
    barplot_ax.set_ylim([0.4, 1.1])
    barplot_fig.tight_layout()
    barplot_fig.show()

    barplot_fig.savefig(
        os.path.join(arg_dict["out_dir"], "ENCSR067NPP_caas_hist.pdf"), format="pdf"
    )
    barplot_fig.savefig(os.path.join(arg_dict["out_dir"], "ENCSR067NPP_caas_hist.png"))

    max_abs_phylop = max(max(label_dict["bins"]) for label_dict in ecdf_data_dict.values())
    fig, ax = plt.subplots(figsize=(7, 7), dpi=400)
    for label, label_dict in ecdf_data_dict.items():
        interp_term = label.split("_")[-1]
        label_color = (
            const.LABEL_COLOR_MAP[interp_term] if interp_term != "Quiescent" else (0.0, 0.0, 0.0)
        )
        x_vals = label_dict["bins"]
        y_vals = [0] + label_dict["counts"]
        y_vals[-1] = 1
        if max(x_vals) < max_abs_phylop:
            x_vals += [max_abs_phylop]
            y_vals += [1]
        ax.plot(x_vals, y_vals, label=label, color=label_color)
    ax.set_xlim([0, max_abs_phylop])
    ax.set_ylim([0, 1])
    ax.axhline(0.75, color="gray", linestyle="--")

    legend_handles = []
    for interp_term, interp_color in const.LABEL_COLOR_MAP.items():
        sample_labels = list(
            sorted(
                [k for k in ecdf_data_dict.keys() if k.split("_")[-1] == interp_term],
                key=lambda s: int(s.split("_")[0]),
            )
        )
        for l in sample_labels:
            legend_handles.append(
                matplotlib.patches.Patch(
                    color=interp_color if interp_term != "Quiescent" else (0, 0, 0), label=l
                )
            )
    ax.set_xlabel("Absolute PhyloP Score", fontsize=16, weight="bold")
    ax.set_ylabel("Empirical Cumulative Density", fontsize=16, weight="bold")
    ax.set_yticklabels(["", 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis="both", which="major", labelsize=18)
    fig.tight_layout()
    fig.savefig(os.path.join(arg_dict["out_dir"], "ENCSR067NPP_phylop_ecdf.pdf"), format="pdf")
    fig.savefig(os.path.join(arg_dict["out_dir"], "ENCSR067NPP_phylop_ecdf.png"), format="png")
