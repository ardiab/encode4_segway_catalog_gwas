import os

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.meta_helper import MetaHelper


class PathHelper:
    def __init__(self, out_dir, meta_helper_obj=None):
        assert os.path.isdir(out_dir), f"out_dir {out_dir} does not exist."
        self.base_out_dir = out_dir

        if meta_helper_obj is not None:
            assert isinstance(meta_helper_obj, MetaHelper)
            self.meta_helper = meta_helper_obj
        else:
            self.meta_helper = MetaHelper.from_csv(csv_fpath=self.meta_df_fpath)

    # region Property paths
    @property
    def base_output_dir(self):
        return self.base_out_dir

    @property
    def phylop_dir(self):
        return os.path.join(self.base_output_dir, "phylop_hg38")

    @property
    def processed_bed_dir(self):
        return os.path.join(self.base_output_dir, "processed_beds")

    @property
    def encode_submission_dir(self):
        return os.path.join(self.base_output_dir, "encode_submission")

    @property
    def processed_snp_dir(self):
        return os.path.join(self.base_output_dir, "processed_snps")

    @property
    def label_info_dir(self):
        return os.path.join(self.base_output_dir, "label_info")

    @property
    def snp_intersection_dir(self):
        return os.path.join(self.base_output_dir, "snp_ann_intersections")

    @property
    def enrichment_dir(self):
        return os.path.join(self.base_output_dir, "enrichment")

    @property
    def plot_dir(self):
        return os.path.join(self.base_output_dir, "plots")

    @property
    def meta_df_fpath(self):
        return os.path.join(self.base_output_dir, "metadata.csv")

    # endregion

    def get_submission_log_fpath(self, test, dry_run):
        assert test in [True, False]
        assert dry_run in [True, False]
        suffix = "test" if test else "production"
        if dry_run:
            suffix += "dry_run"
        return os.path.join(self.encode_submission_dir, f"submission_log_{suffix}.txt")

    def get_sample_bed_fpath(self, sample_id):
        if self.meta_helper is None:
            raise AttributeError("PathHelper was initialized without metadata information.")
        else:
            return self.meta_helper.get_sample_bedpath(sample_id)

    def get_phylop_wiggle_fpath(self, chrom):
        # Add a "." after the chromosome name so that e.g. "chr1" doesn't match files with "chr11".
        # The wiggle file names are formatted as {chrom name}.phyloP30way.wigFix.
        wig_fname = [fname for fname in const.PHYLOP_WIGS if f"{chrom}." in fname]
        assert len(wig_fname) == 1, f"Unable to find PhyloP Wiggle file for chromosome {chrom}"

        return os.path.join(self.phylop_dir, wig_fname[0])

    def get_phylop_numpy_fpath(self, chrom):
        phylop_wiggle_fpath = self.get_phylop_wiggle_fpath(chrom)
        assert phylop_wiggle_fpath.endswith(".wigFix")

        return phylop_wiggle_fpath.replace(".wigFix", ".npy")

    def get_sample_label_info_dir_path(self, sample_id):
        return os.path.join(self.label_info_dir, sample_id)

    def get_label_caas_json_fpath(self, sample_id):
        return os.path.join(
            self.get_sample_label_info_dir_path(sample_id), f"{sample_id}_label_caas.json"
        )

    def get_label_phylop_percentiles_json_fpath(self, sample_id):
        return os.path.join(
            self.get_sample_label_info_dir_path(sample_id),
            f"{sample_id}_label_phylop_percentiles.json",
        )

    def get_label_phylop_ecdf_plot_data_json_fpath(self, sample_id):
        return os.path.join(
            self.get_sample_label_info_dir_path(sample_id),
            f"{sample_id}_label_phylop_ecdf_plot_data.json",
        )

    def get_label_coverage_json_fpath(self, sample_id, normalized=True):
        if normalized:
            return os.path.join(
                self.get_sample_label_info_dir_path(sample_id), f"{sample_id}_label_coverage.json"
            )
        else:
            return os.path.join(
                self.get_sample_label_info_dir_path(sample_id), f"{sample_id}_label_counts.json"
            )

    def get_snp_bed_fpath(self, window_size, fake=False):
        assert fake in [True, False]
        return os.path.join(
            self.processed_snp_dir,
            f"{'fake' if fake else 'real'}_snp_bed_{window_size}_bidirectional_window.tsv",
        )

    def get_trait_snp_linkage_df_fpath(self, window_size):
        return os.path.join(
            self.processed_snp_dir, f"trait_snp_linkage_{window_size}_bidirectional_window.csv"
        )

    def get_sample_intersection_dir(self, sample_id):
        return os.path.join(self.snp_intersection_dir, sample_id)

    def get_sample_intersection_bed_fpath(self, sample_id, window_size, fake=False, zipped=False):
        assert fake in [True, False]
        sample_intersection_bed_fpath = os.path.join(
            self.get_sample_intersection_dir(sample_id),
            f"{'fake' if fake else 'real'}_snp_ann_intersections_{window_size}.bed{'.gz' if zipped else ''}",
        )

        return sample_intersection_bed_fpath

    def get_label_distr_df_fpath(self, sample_id, window_size, fake=False):
        assert fake in [True, False]
        return os.path.join(
            self.get_sample_intersection_dir(sample_id),
            f"{sample_id}_label_distr_{'fake' if fake else 'real'}_snps_{window_size}.pkl",
        )

    # ENCODE submission
    def get_sample_base_encode_submission_sheet_dir(self, sample_id):
        return os.path.join(self.encode_submission_dir, sample_id)

    def get_sample_test_encode_submission_sheet_dir(self, sample_id):
        return os.path.join(
            self.get_sample_base_encode_submission_sheet_dir(sample_id), "test_submission"
        )

    def get_sample_prod_encode_submission_sheet_dir(self, sample_id):
        return os.path.join(
            self.get_sample_base_encode_submission_sheet_dir(sample_id), "production_submission"
        )

    def get_sample_submission_pdf_figure_fpath(self, sample_id):
        return os.path.join(
            self.get_sample_base_encode_submission_sheet_dir(sample_id),
            f"Information for annotation {sample_id}.pdf",
        )

    def get_sample_file_sheet_fpath(self, sample_id, test):
        assert test in [True, False]
        fname = "file_sheet.tsv"
        if test:
            return os.path.join(self.get_sample_test_encode_submission_sheet_dir(sample_id), fname)
        else:
            return os.path.join(self.get_sample_prod_encode_submission_sheet_dir(sample_id), fname)

    def get_sample_analysis_sheet_fpath(self, sample_id, test):
        assert test in [True, False]
        fname = "analysis_sheet.tsv"
        if test:
            return os.path.join(self.get_sample_test_encode_submission_sheet_dir(sample_id), fname)
        else:
            return os.path.join(self.get_sample_prod_encode_submission_sheet_dir(sample_id), fname)

    def get_sample_document_sheet_fpath(self, sample_id, test):
        assert test in [True, False]
        fname = "document_sheet.tsv"
        if test:
            return os.path.join(self.get_sample_test_encode_submission_sheet_dir(sample_id), fname)
        else:
            return os.path.join(self.get_sample_prod_encode_submission_sheet_dir(sample_id), fname)

    def get_sample_annotation_patch_sheet_fpath(self, sample_id, test):
        assert test in [True, False]
        fname = "annotation_patch_sheet.tsv"
        if test:
            return os.path.join(self.get_sample_test_encode_submission_sheet_dir(sample_id), fname)
        else:
            return os.path.join(self.get_sample_prod_encode_submission_sheet_dir(sample_id), fname)

    def get_snp_metric_df_fpath(self, metric_name, window_size, fake):
        assert fake in [True, False]
        fname = f"snp_{metric_name}_{window_size}_{'real' if not fake else 'fake'}.pkl"
        return os.path.join(self.enrichment_dir, fname)

    def get_pval_df_fpath(self, metric_name, test_type, window_size, fake):
        assert fake in [True, False]
        fname = f"pvals_{metric_name}_{test_type}_{window_size}_{'fake' if fake else 'real'}.pkl"
        return os.path.join(self.enrichment_dir, fname)

    def get_test_stat_df_fpath(self, metric_name, test_type, window_size, fake):
        assert fake in [True, False]
        fname = (
            f"test_stats_{metric_name}_{test_type}_{window_size}_{'fake' if fake else 'real'}.pkl"
        )
        return os.path.join(self.enrichment_dir, fname)

    def get_effect_size_df_fpath(self, metric_name, test_type, window_size, fake):
        assert fake in [True, False]
        fname = (
            f"effect_size_{metric_name}_{test_type}_{window_size}_{'fake' if fake else 'real'}.pkl"
        )
        return os.path.join(self.enrichment_dir, fname)

    def get_pval_heatmap_fig_fpath(self, metric_name, test_type, window_size, fake, desc):
        assert fake in [True, False]
        fname = f"pval_heatmap_{metric_name}_{test_type}_{window_size}_{'fake' if fake else 'real'}_{desc}.pdf"
        return os.path.join(self.plot_dir, fname)

    def get_test_stat_heatmap_fig_fpath(self, metric_name, test_type, window_size, fake):
        assert fake in [True, False]
        fname = f"test_stat_heatmap_{metric_name}_{test_type}_{window_size}_{'fake' if fake else 'real'}.pdf"
        return os.path.join(self.plot_dir, fname)

    def get_effect_size_heatmap_fig_fpath(self, metric_name, test_type, window_size, fake):
        assert fake in [True, False]
        fname = f"effect_size_heatmap_{metric_name}_{test_type}_{window_size}_{'fake' if fake else 'real'}.pdf"
        return os.path.join(self.plot_dir, fname)
