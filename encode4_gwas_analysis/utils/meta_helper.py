import os
import warnings

import pandas as pd
from tqdm import tqdm

pd.options.mode.chained_assignment = None


class MetaHelper:
    SAMPLE_ID_COL = "sample_id"
    BEDPATH_COL = "bed_path"
    SAMPLE_DESCRIPTION_COL = "sample_description"
    BIGWIG_ACCESSIONS_COL = "bigwig_accessions"
    FIG_PATHS_COL = "figure_paths"
    MNEM_PATH_COL = "mnemonics_path"
    LIST_DELIMITER = "@@@\t@@@"

    def __init__(self, metadata_df) -> None:
        self.meta_df = metadata_df
        self.validate_meta_df(self.meta_df)

    @classmethod
    def validate_meta_df(cls, meta_df):
        assert meta_df.notnull().all().all(), "The metadata DF contains null values."
        assert not meta_df[cls.SAMPLE_ID_COL].duplicated().any()
        for fp in meta_df[cls.BEDPATH_COL].values:
            assert os.path.isfile(fp), f"Unable to find bed file {fp}."
            assert fp.endswith(".bed") or fp.endswith(".bed.gz"), (
                "All bed files must have an extension of '.bed' or '.bed.gz'."
            )
        if meta_df[cls.SAMPLE_DESCRIPTION_COL].duplicated().any():
            warnings.warn("Found duplicated sample descriptions.")

    @classmethod
    def from_csv(cls, csv_fpath):
        meta_df = pd.read_csv(csv_fpath)
        for idx, row in meta_df.iterrows():
            for col, val in row.items():
                if isinstance(val, str) and cls.LIST_DELIMITER in val:
                    meta_df.at[idx, col] = val.split(cls.LIST_DELIMITER)

        return MetaHelper(metadata_df=meta_df)

    @classmethod
    def save_meta_df(cls, meta_df, df_path):
        cls.validate_meta_df(meta_df)
        meta_df_cp = meta_df.copy(deep=True)
        for idx, row in meta_df_cp.iterrows():
            for col, val in row.items():
                if isinstance(val, list):
                    meta_df_cp.loc[idx, col] = cls.LIST_DELIMITER.join(val)
        meta_df_cp.to_csv(df_path, index=False)

    def get_all_sample_ids(self):
        return self.meta_df[self.SAMPLE_ID_COL].values

    def iter_samples(self, show_progress=False):
        iterable = (
            tqdm(self.meta_df[self.SAMPLE_ID_COL].values)
            if show_progress
            else self.meta_df[self.SAMPLE_ID_COL].values
        )
        for sample_id in iterable:
            yield sample_id

    def validate_sample_id(self, sample_id):
        assert (self.meta_df[self.SAMPLE_ID_COL] == sample_id).any(), (
            f"Unable to find metadata for sample ID {sample_id}"
        )

    # This function was put here instead of in PathHelper because these bed files are inputs to the analysis, and can
    # have varying locations (specified by the metadata CSV).
    # PathHelper implements a method by the same name, which calls this method, for convenience.
    def get_sample_bedpath(self, sample_id):
        self.validate_sample_id(sample_id)
        return self.meta_df.set_index(self.SAMPLE_ID_COL).loc[sample_id, self.BEDPATH_COL]

    def get_sample_description(self, sample_id):
        self.validate_sample_id(sample_id)
        return self.meta_df.set_index(self.SAMPLE_ID_COL).loc[
            sample_id, self.SAMPLE_DESCRIPTION_COL
        ]

    def get_sample_bigwigs(self, sample_id):
        self.validate_sample_id(sample_id)
        return self.meta_df.set_index(self.SAMPLE_ID_COL).loc[sample_id, self.BIGWIG_ACCESSIONS_COL]

    def get_sample_fig_paths(self, sample_id):
        self.validate_sample_id(sample_id)
        return self.meta_df.set_index(self.SAMPLE_ID_COL).loc[sample_id, self.FIG_PATHS_COL]

    def get_sample_mnemonics_path(self, sample_id):
        self.validate_sample_id(sample_id)
        return self.meta_df.set_index(self.SAMPLE_ID_COL).loc[sample_id, self.MNEM_PATH_COL]
