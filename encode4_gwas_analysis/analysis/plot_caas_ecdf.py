"""Plot the ECDF of SNP CAAS values for significant and insignificant traits."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import tqdm

from encode4_gwas_analysis.utils.data_loader import DataLoader

sns.set_style("whitegrid")


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = parser.parse_args()
    assert os.path.isdir(arg_dict.out_dir), f"Directory {arg_dict.out_dir} does not exist."

    return vars(arg_dict)


if __name__ == "__main__":
    arg_dict = parse_args()
    data_loader = DataLoader(out_dir=arg_dict["out_dir"])

    WINDOW_SIZE = 10_000
    SNP_CUTOFF = 30
    SIG_THRESHOLD = 0.05

    trait_linkage_df = data_loader.load_trait_snp_linkage_df(window_size=WINDOW_SIZE)
    trait_snp_count_df = data_loader.load_trait_snp_count_df(window_size=WINDOW_SIZE)

    traits_to_keep = trait_snp_count_df[trait_snp_count_df["n_snps"] >= SNP_CUTOFF][
        "trait"
    ].unique()
    pval_df = (
        data_loader.load_pval_df(
            metric_name="mean_caas",
            test_type="signedrank_per_sample_null",
            window_size=WINDOW_SIZE,
            fake=False,
        ).fillna(1)
        * 234
    )
    real_caas_df = data_loader.load_snp_metric_df(
        metric_name="mean_caas", snp_region_side_length=WINDOW_SIZE, fake=False
    )
    fake_caas_df = data_loader.load_snp_metric_df(
        metric_name="mean_caas", snp_region_side_length=WINDOW_SIZE, fake=True
    )
    assert pval_df.shape[0] == len(traits_to_keep)

    pval_df = pval_df.loc[traits_to_keep]
    snp_caas_data = {"real_significant": [], "real_insignificant": [], "fake": []}
    percentiles_to_get = np.arange(0, 100.1, 1)
    for sample_id in tqdm.tqdm(pval_df.columns):
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
        # Some SNPs may be in both significant and insignificant traits. Remove them from traits
        # where they are insignificant, and include them in traits where they are significant.
        sample_insig_snps = list(set(sample_insig_snps).difference(set(sample_sig_snps)))

        sig_real_snp_caas = real_caas_df[sample_id].loc[sample_sig_snps].tolist()
        insig_real_snp_caas = real_caas_df[sample_id].loc[sample_insig_snps].tolist()
        all_sig_insig_fake_snps = [snp + "_fake" for snp in sample_sig_snps + sample_insig_snps]
        sig_insig_fake_snps = list(
            set(all_sig_insig_fake_snps).intersection(set(fake_caas_df.index.values))
        )
        fake_snp_caas = fake_caas_df[sample_id].loc[sig_insig_fake_snps].tolist()

        snp_caas_data["real_significant"] += sig_real_snp_caas
        snp_caas_data["real_insignificant"] += insig_real_snp_caas
        snp_caas_data["fake"] += fake_snp_caas

    fig, ax = plt.subplots(figsize=(10, 10), dpi=400)
    for snp_type, snp_caas in snp_caas_data.items():
        caas_subset = []
        for p in tqdm.tqdm(percentiles_to_get):
            caas_subset.append(np.nanpercentile(snp_caas, p))
        ax.plot(caas_subset, [100 - p for p in percentiles_to_get], label=snp_type)

    ax.legend()
    ax.set_xlabel("SNP Mean CAAS", fontsize=18, weight="bold")
    ax.set_ylabel("% of SNPs with mean CAAS >= x", fontsize=18, weight="bold")
    fig.show()
