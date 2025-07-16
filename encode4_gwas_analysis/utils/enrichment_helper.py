import pandas as pd
from scipy.stats import ranksums, wilcoxon
from tqdm import tqdm

tqdm.pandas()


class EnrichmentHelper:
    @classmethod
    def get_label_distr(cls, biosample_name):
        pass

    @classmethod
    def get_mean_caas_series(cls, label_distr_df, label_caas_series):
        """
        This method expects a label distribution DF, where:
        - The index specifies both the trait name and the SNP ID (i.e. every row is 1 SNP)
        - There is one column per label which contains the number of bases in the SNP region that received
          that label as an annotation.
        """
        assert set(label_distr_df.columns) == set(label_caas_series.index)
        # Check that all columns are numeric data types.
        all([pd.api.types.is_numeric_dtype(label_distr_df[c]) for c in label_distr_df.columns])
        # Each row contains the label distribution in the window surrounding one SNP.
        # Mean CAAS for this window can be calculated by taking a weighted average:
        # 1. Multiply the number of bases corresponding to each of the labels by the CAAS score for the label and sum
        #    up the results.
        # 2. Divide the sum by the total number of bases in the window.
        return label_distr_df.apply(lambda row: (row * label_caas_series) / row.sum(), axis=1)

    @classmethod
    def get_max_caas_series(cls, label_distr_df, label_caas_series):
        """
        This method expects a label distribution DF, where:
        - The index specifies both the trait name and the SNP ID (i.e. every row is 1 SNP)
        - There is one column per label which contains the number of bases in the SNP region that received
          that label as an annotation.
        """
        assert set(label_distr_df.columns) == set(label_caas_series.index)
        # Check that all columns are numeric data types.
        assert all(
            [pd.api.types.is_numeric_dtype(label_distr_df[c]) for c in label_distr_df.columns]
        )
        return (
            label_distr_df.astype(bool)
            .astype(int)
            .apply(lambda row: (row * label_caas_series).max(), axis=1)
        )

    @classmethod
    def rank_snps_across_biosamples(cls, metric_df):
        rank_df = metric_df.rank(axis=1)
        return rank_df

    @classmethod
    def rank_snps_within_biosamples(cls, metric_df):
        rank_df = metric_df.rank(axis=0)
        return rank_df

    @classmethod
    def _convert_snp_index_to_trait_snp_cols(cls, df):
        df_proc = df.reset_index()
        df_proc["trait"] = df_proc["index"].apply(lambda s: s.split("_snp_")[0])
        df_proc["snp"] = df_proc["index"].apply(lambda s: s.split("_snp_")[1])

        return df_proc.drop(columns=["index"])

    @classmethod
    def get_signedrank_pvals_per_trait_across_biosamples(cls, metric_df, print_progress=False):
        rank_df = cls.rank_snps_across_biosamples(metric_df).reset_index().drop(columns=["snp"])
        median_rank = (rank_df.shape[1] / 2) + 1  # Add 1 because ranks are 1-indexed
        if print_progress:
            return rank_df.groupby("trait").progress_apply(
                lambda trait_df: trait_df.drop(columns=["trait"]).apply(
                    lambda col: wilcoxon(col.values - median_rank, alternative="greater")[1]
                )
            )
        else:
            return rank_df.groupby("trait").apply(
                lambda trait_df: trait_df.drop(columns=["trait"]).apply(
                    lambda col: wilcoxon(col.values - median_rank, alternative="greater")[1]
                )
            )

    @classmethod
    def get_signedrank_pvals_per_biosample_across_traits(cls, metric_df, print_progress=False):
        rank_df = cls.rank_snps_within_biosamples(metric_df).reset_index().drop(columns=["snp"])
        median_rank = (rank_df.shape[0] / 2) + 1  # Add 1 because ranks are 1-indexed
        if print_progress:
            return rank_df.groupby("trait").progress_apply(
                lambda trait_df: trait_df.drop(columns=["trait"]).apply(
                    lambda col: wilcoxon(col.values - median_rank, alternative="greater")[1]
                )
            )
        else:
            return rank_df.groupby("trait").apply(
                lambda trait_df: trait_df.drop(columns=["trait"]).apply(
                    lambda col: wilcoxon(col.values - median_rank, alternative="greater")[1]
                )
            )

    @classmethod
    def get_signedrank_pvals_ranked_across_traits_biosamples(cls, metric_df, print_progress=False):
        rank_df = cls.rank_snps_within_biosamples(metric_df)
        rank_df = cls.rank_snps_across_biosamples(rank_df).reset_index().drop(columns=["snp"])
        median_rank = (rank_df.shape[1] / 2) + 1  # Add 1 because ranks are 1-indexed
        if print_progress:
            return rank_df.groupby("trait").progress_apply(
                lambda trait_df: trait_df.drop(columns=["trait"]).apply(
                    lambda col: wilcoxon(col.values - median_rank, alternative="greater")[1]
                )
            )
        else:
            return rank_df.groupby("trait").apply(
                lambda trait_df: trait_df.drop(columns=["trait"]).apply(
                    lambda col: wilcoxon(col.values - median_rank, alternative="greater")[1]
                )
            )

    @classmethod
    def get_ranksum_pvals_per_trait_across_biosamples(cls, metric_df, print_progress=False):
        data = []
        biosample_names = metric_df.columns
        grouper = (
            tqdm(metric_df.reset_index().groupby("trait"))
            if print_progress
            else metric_df.reset_index().groupby("trait")
        )
        for trait, trait_df in grouper:
            for biosample_name in biosample_names:
                other_biosamples = [b for b in biosample_names if b != biosample_name]
                _, pval = ranksums(
                    trait_df.loc[:, biosample_name],
                    trait_df.loc[:, other_biosamples],
                    alternative="greater",
                )
                data.append({"biosample": biosample_name, "trait": trait, "pval": pval})

        return pd.DataFrame(data).pivot(index="trait", columns="biosample", values="pval")

    @classmethod
    def get_ranksum_pvals_per_biosample_across_traits(cls, metric_df, print_progress=False):
        data = []
        rank_df = cls.rank_snps_across_biosamples(metric_df)
        trait_names = set(rank_df.reset_index()["trait"].values)

        cols = metric_df.columns if not print_progress else tqdm(metric_df.columns)
        for biosample_name in cols:
            biosample_series = metric_df[biosample_name]
            for trait in tqdm(trait_names):
                # Very rarely, we get null metric values/CAAS because a SNP is missing annotations in one of
                # the biosamples. When this happens, simply drop the null values and proceed as normal.
                trait_ranks = biosample_series.loc[trait].dropna().values
                other_ranks = (
                    biosample_series.loc[list(trait_names.difference(trait))].dropna().values
                )
                _, pval = ranksums(trait_ranks, other_ranks, alternative="greater")
                data.append({"biosample": biosample_name, "trait": trait, "pval": pval})

        return pd.DataFrame(data).pivot(index="trait", columns="biosample", values="pval")
