import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import tqdm

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.data_loader import DataLoader
from encode4_gwas_analysis.utils.path_helper import PathHelper

sns.set_style("whitegrid")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--sample_id_list_fpath", required=True)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict["out_dir"])
    assert os.path.isfile(arg_dict["sample_id_list_fpath"])

    return arg_dict


def compute_caas(sample_id, chrom_cons_array_dict, data_loader, path_helper):
    ann_bed_df = data_loader.load_biosample_annotation_bed(
        sample_id=sample_id, drop_alternate_loci=True
    )

    sample_out_dir = path_helper.get_sample_label_info_dir_path(sample_id)
    os.mkdir(sample_out_dir)

    label_caas = {}
    # Get the phyloP scores of each position annotated by each label.
    for label, label_df in ann_bed_df.groupby("label"):
        label_caas[label] = []
        for chrom, chrom_df in label_df.groupby("chrom"):
            assert chrom in chrom_cons_array_dict
            chrom_idxs = [
                (row["chrom_start"], row["chrom_end"]) for idx, row in chrom_df.iterrows()
            ]
            for start_idx, stop_idx in chrom_idxs:
                label_caas[label].append(chrom_cons_array_dict[chrom][start_idx:stop_idx])

    fig, ax = plt.subplots(figsize=(15, 15))
    # Get the CAAS for each label by taking the 75th percentile of the absolute values of all phyloP scores of positions
    # annotated by the label.
    label_percentiles = {}
    label_plot_info = {}
    caas_scores = {}
    for label, cons_arr in label_caas.items():
        abs_concat_array = np.abs((np.concatenate(cons_arr)))
        label_percentiles[label] = {}
        for p in [0, 50, 75, 99, 100]:
            label_percentiles[label][p] = np.nanpercentile(abs_concat_array, p)
        caas_scores[label] = label_percentiles[label][75]
        interp_term = label.split("_")[1]
        # Color quiescent black instead of white so it shows up on the plot (which has a white background).
        label_color = (
            const.LABEL_COLOR_MAP[interp_term]
            if interp_term.lower() != "quiescent"
            else (0.0, 0.0, 0.0)
        )
        counts, bins, bars = ax.hist(
            np.clip(abs_concat_array, a_min=None, a_max=label_percentiles[label][99]),
            density=True,
            cumulative=True,
            label=label,
            histtype="step",
            bins=30,
            color=label_color,
        )
        label_plot_info[label] = {"counts": list(counts), "bins": list(bins)}

    ax.axhline(0.75, color="black", linestyle="--", alpha=0.4)
    ax.set_yticks(np.arange(0, 1.05, 0.05))
    ax.set_xlabel(
        "Absolute PhyloP score\n(clipped to remove values >= 99th percentile)",
        fontsize=14,
        weight="bold",
    )
    ax.set_ylabel("eCDF", fontsize=14, weight="bold")
    ax.set_title("Label PhyloP eCDFs", fontsize=18, weight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(os.path.join(sample_out_dir, "label_phylop_distr_hist.png"))
    plt.close(fig)

    with open(path_helper.get_label_phylop_percentiles_json_fpath(sample_id), "w") as out_f:
        json.dump(label_percentiles, out_f)

    with open(path_helper.get_label_phylop_ecdf_plot_data_json_fpath(sample_id), "w") as out_f:
        json.dump(label_plot_info, out_f)

    # Save label CAAS to a JSON.
    with open(path_helper.get_label_caas_json_fpath(sample_id), "w") as out_f:
        json.dump(caas_scores, out_f)

    # Create a dictionary keyed by Segway label ID which indicates the proportion of annotations that the label forms.
    # Note that this uses the label distribution after dropping the "bad" chromosomes like 'chrUn_...' -- upon examining
    # the annotation beds, it seems that these "bad" chromosomes constitute less than 1% of total annotated bases.
    ann_bed_df["length"] = ann_bed_df["chrom_end"] - ann_bed_df["chrom_start"]
    coverage_unnormed = ann_bed_df.groupby("label")["length"].sum().to_dict()
    coverage_normed = (
        ann_bed_df.groupby("label")["length"].sum() / ann_bed_df["length"].sum()
    ).to_dict()
    with open(path_helper.get_label_coverage_json_fpath(sample_id, normalized=False), "w") as out_f:
        json.dump(coverage_unnormed, out_f)
    with open(path_helper.get_label_coverage_json_fpath(sample_id, normalized=True), "w") as out_f:
        json.dump(coverage_normed, out_f)


def main(out_dir, sample_id_list_fpath):
    data_loader = DataLoader(out_dir=out_dir)
    path_helper = PathHelper(out_dir=out_dir)

    with open(sample_id_list_fpath, "r") as in_f:
        sample_ids = [sid for sid in in_f.read().split("\n") if len(sid) > 0]

    # Load all the phyloP arrays once into a dict keyed by chromosome. This dict is reused for every annotation .bed
    # file in the specified list.
    chrom_cons_array_dict = {}
    for chrom in const.CHROM_LENGTHS.keys():
        assert chrom not in chrom_cons_array_dict
        chrom_phylop_array_fpath = path_helper.get_phylop_numpy_fpath(chrom)
        chrom_cons_array_dict[chrom] = data_loader.load_phylop_array(
            array_fpath=chrom_phylop_array_fpath, chrom=chrom
        )

    for sample_id in tqdm.tqdm(sample_ids, desc="Computing sample CAAS"):
        compute_caas(
            sample_id=sample_id,
            chrom_cons_array_dict=chrom_cons_array_dict,
            data_loader=data_loader,
            path_helper=path_helper,
        )


if __name__ == "__main__":
    arg_dict = parse_args()
    main(
        out_dir=arg_dict["out_dir"],
        sample_id_list_fpath=arg_dict["sample_id_list_fpath"],
    )
