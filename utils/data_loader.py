import json

import encode_gwas.utils.constants as const
import numpy as np
import pandas as pd
from encode_gwas.utils.constants import CHROM_LENGTHS
from encode_gwas.utils.meta_helper import MetaHelper
from encode_gwas.utils.path_helper import PathHelper


class DataLoader:
    BED_COLS = ['chrom', 'chrom_start', 'chrom_end', 'label', 'score', 'strand', 'thickStart', 'thickEnd',
                'itemRgb', 'label_detailed']

    def __init__(self, out_dir):
        self.path_helper = PathHelper(out_dir=out_dir)
        self.meta_helper = MetaHelper.from_csv(self.path_helper.meta_df_fpath)

    @classmethod
    def validate_annotation_bed(cls, df, allow_alternate_loci=False):
        if not allow_alternate_loci:
            assert df[cls.BED_COLS[0]].isin([f'chr{i}' for i in range(23)] + ['chrX', 'chrY']).all()
        else:
            assert df[cls.BED_COLS[0]].apply(lambda c: c.split('_')[0]).isin(
                [f'chr{i}' for i in range(23)] + ['chrX', 'chrY']).all()
        assert df[cls.BED_COLS[1]].dtype == int
        assert df[cls.BED_COLS[2]].dtype == int

    @classmethod
    def load_annotation_bed(cls, fpath, drop_alternate_loci=True):
        """
        Load an annotation bed file as a Pandas DataFrame. Supports the standard bed9 format, and a bed9+ format where
        the first 9 columns are standard bed9 columns, and the tenth column is a detailed annotation term.
        In the context of Segway annotations, the tenth column would contain to the integer labels produced by the
        algorithm, whereas the fourth column would contain the interpretation terms for the interpretation labels.
        """
        df = pd.read_csv(fpath, sep='\t', header=None)
        assert df.shape[1] in [9, 10]
        if df.shape[1] == 10:
            assert df[9].dtype == int
            df[3] = df[9].astype(str) + '_' + df[3].astype(str)

        # Drop everything except the first four columns, because they are not needed for the analysis.
        df = df.drop(columns=list(range(4, df.shape[1])))
        df.columns = cls.BED_COLS[:4]

        if drop_alternate_loci:
            df = df[df[cls.BED_COLS[0]].isin(const.CHROM_LENGTHS.keys())]

        cls.validate_annotation_bed(df, allow_alternate_loci=not drop_alternate_loci)

        return df

    def load_biosample_annotation_bed(self, sample_id, drop_alternate_loci=True):
        bed_fpath = self.meta_helper.get_sample_bedpath(sample_id=sample_id)

        return self.load_annotation_bed(bed_fpath, drop_alternate_loci=drop_alternate_loci)

    @classmethod
    def load_phylop_array(cls, array_fpath, chrom):
        with open(array_fpath, 'rb') as in_f:
            phylop_array = np.load(in_f)
        assert phylop_array.shape == (CHROM_LENGTHS[chrom], )

        return phylop_array

    @classmethod
    def validate_intersection_bed(cls, df):
        assert (df['ann_chrom'] == df['snp_chrom']).all()
        assert (df[['ann_start', 'ann_end', 'snp_start', 'snp_end', 'intersection_length']].dtypes == np.int64).all()

    @classmethod
    def load_intersection_bed(cls, fpath):
        df = pd.read_csv(fpath, sep='\t', header=None)
        assert df.shape[1] == 15
        df[3] = df[9].astype(str) + '_' + df[3]
        df = df[[0, 1, 2, 3, 10, 11, 12, 13, 14]]
        df.columns = ['ann_chrom', 'ann_start', 'ann_end', 'ann_label', 'snp_chrom', 'snp_start', 'snp_end',
                      'snp_name', 'intersection_length']

        cls.validate_intersection_bed(df)

        return df

    def load_sample_intersection_bed(self, sample_id, window_size, fake):
        assert fake in [True, False]
        return pd.read_pickle(self.path_helper.get_sample_intersection_bed_fpath(sample_id=sample_id,
                                                                                 window_size=window_size,
                                                                                 fake=fake))

    def load_sample_snp_label_distribution(self, sample_id, snp_region_side_length, fake):
        distr_df_fpath = self.path_helper.get_label_distr_df_fpath(sample_id=sample_id,
                                                                   window_size=snp_region_side_length,
                                                                   fake=fake)
        return pd.read_pickle(distr_df_fpath)

    def load_sample_label_caas(self, sample_id):
        caas_json_fpath = self.path_helper.get_label_caas_json_fpath(sample_id)
        with open(caas_json_fpath, 'r') as in_f:
            return json.load(in_f)

    def load_sample_label_coverage(self, sample_id, normalized):
        label_coverage_fpath = self.path_helper.get_label_coverage_json_fpath(sample_id=sample_id,
                                                                              normalized=normalized)
        with open(label_coverage_fpath, 'r') as in_f:
            return json.load(in_f)

    def load_sample_label_phylop_percentiles(self, sample_id):
        percentile_fpath = self.path_helper.get_label_phylop_percentiles_json_fpath(sample_id=sample_id)

        with open(percentile_fpath, 'r') as in_f:
            return json.load(in_f)

    def load_snp_bed(self, snp_region_side_length, fake):
        snp_bed_fpath = self.path_helper.get_snp_bed_fpath(window_size=snp_region_side_length, fake=fake)
        snp_bed = pd.read_csv(snp_bed_fpath, sep='\t', header=None)
        assert snp_bed.shape[1] == 4
        snp_bed.columns = ['chrom', 'snp_start', 'snp_end', 'snp_name']

        return snp_bed

    def load_snp_metric_df(self, metric_name, snp_region_side_length, fake):
        metric_df_fpath = self.path_helper.get_snp_metric_df_fpath(metric_name=metric_name,
                                                                   window_size=snp_region_side_length,
                                                                   fake=fake)
        return pd.read_pickle(metric_df_fpath)

    def load_trait_snp_linkage_df(self, window_size):
        snp_linkage_df_fpath = self.path_helper.get_trait_snp_linkage_df_fpath(window_size=window_size)

        return pd.read_csv(snp_linkage_df_fpath).drop(columns=['Unnamed: 0'])

    def load_trait_snp_count_df(self, window_size):
        trait_snp_linkage_df = self.load_trait_snp_linkage_df(window_size=window_size)
        assert trait_snp_linkage_df.groupby('trait').apply(
            lambda trait_df: trait_df['snp'].nunique() == trait_df.shape[0]).all()

        return trait_snp_linkage_df.groupby('trait')['snp'].nunique().reset_index().rename(columns={'snp': 'n_snps'})

    def load_pval_df(self, metric_name, test_type, window_size, fake):
        pval_df_fpath = self.path_helper.get_pval_df_fpath(metric_name=metric_name, test_type=test_type,
                                                           window_size=window_size, fake=fake)

        return pd.read_pickle(pval_df_fpath)

    def load_test_stat_df(self, metric_name, test_type, window_size, fake):
        test_stat_df_fpath = self.path_helper.get_test_stat_df_fpath(metric_name=metric_name, test_type=test_type,
                                                                     window_size=window_size, fake=fake)

        return pd.read_pickle(test_stat_df_fpath)

    def load_effect_size_df(self, metric_name, test_type, window_size, fake):
        effect_size_df_fpath = self.path_helper.get_effect_size_df_fpath(metric_name=metric_name, test_type=test_type,
                                                                         window_size=window_size, fake=fake)

        return pd.read_pickle(effect_size_df_fpath)
