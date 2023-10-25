import argparse
import math
import os

import encode_gwas.utils.constants as const
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from encode_gwas.utils.data_loader import DataLoader
from encode_gwas.utils.meta_helper import MetaHelper
from encode_gwas.utils.path_helper import PathHelper
from tqdm import tqdm


# sns.set_style('whitegrid')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', required=True)

    arg_dict = vars(parser.parse_args())

    assert os.path.isdir(arg_dict['out_dir']), f'Directory {arg_dict["out_dir"]} does not exist.'

    return arg_dict


def plot_snp_appearance():
    pass


def get_snp_window_contents(snp_anns, window_size, label_idx_dict):
    total_length = 2*window_size + 1
    assert snp_anns['intersection_length'].sum() == total_length
    content_array = np.zeros((len(label_idx_dict), total_length), dtype=np.int32)
    pos = 0
    for idx, row in snp_anns.iterrows():
        length = row['intersection_length']
        interp_term = row['ann_label'].split('_')[-1]
        label_idx = label_idx_dict[interp_term]
        content_array[label_idx, pos:pos+length] += 1
        pos += length

    return content_array


def create_snp_window_heatmap(snp_contents, title, label_idx_dict, ax):
    rev_dict = {v: k for k, v in label_idx_dict.items()}
    window_labels = (snp_contents * np.array(range(len(label_idx_dict)))[:, None]).sum(axis=0)
    label_arr = np.array([const.LABEL_COLOR_MAP[rev_dict[val]] for val in window_labels])[None, :]
    ax.imshow(label_arr, aspect=5000)

    ax.legend(handles=[matplotlib.patches.Patch(color=const.LABEL_COLOR_MAP[label], label=label) for label in
                       const.LABEL_COLOR_MAP.keys()], loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title(title)
    ax.set_xticks(np.arange(0, 20_001, 2500))
    ax.set_xticklabels(np.arange(0, 20_001, 2500))


def main(out_dir):
    WINDOW_SIZE = 10_000

    path_helper = PathHelper(out_dir=out_dir)
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)
    data_loader = DataLoader(out_dir=out_dir)

    trait_snp_linkage_df = data_loader.load_trait_snp_linkage_df(window_size=WINDOW_SIZE)
    pval_df = data_loader.load_pval_df(metric_name='mean_caas', test_type='signedrank_per_sample_null',
                                       window_size=WINDOW_SIZE, fake=False) * 234
    label_idxs = {label: idx for idx, label in enumerate(const.LABEL_COLOR_MAP.keys())}
    mean_caas_df = data_loader.load_snp_metric_df(metric_name='mean_caas', snp_region_side_length=WINDOW_SIZE,
                                                  fake=False)
    n_interp_terms = len(label_idxs)

    # These beds are for sample ENCSR011EQM
    real_bed = data_loader.load_intersection_bed(os.path.join(out_dir, 'ENCSR111ABE_real_snp_ann_intersections_10000.bed.gz'))
    real_bed = real_bed.sort_values(by=['snp_name', 'ann_start'], ascending=True).set_index('snp_name')

    ctype = 'liver tissue male adult (31 years)'
    trait = 'urinary metabolite levels in chronic kidney disease'
    meta_df = meta_helper.meta_df
    accession = meta_df[meta_df['sample_description'].str.contains(ctype, regex=False)]
    assert accession.shape[0] == 1
    accession = accession['sample_id'].values[0]
    trait_pvals = pval_df[accession]
    trait_snps = trait_snp_linkage_df[trait_snp_linkage_df['trait'] == trait]['snp'].unique()

    n_cols = 8
    n_rows = math.ceil(len(trait_snps) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*10, n_rows*4))
    for idx, snp_name in tqdm(enumerate(trait_snps), total=len(trait_snps)):
        ax = axes[np.unravel_index(idx, axes.shape)]
        snp_window_contents = get_snp_window_contents(snp_anns=real_bed.loc[snp_name], window_size=WINDOW_SIZE,
                                                      label_idx_dict=label_idxs)
        assert all(snp_window_contents.sum(axis=0) > 0)
        title = f'{snp_name} (Mean CAAS = {mean_caas_df.loc[snp_name, accession]:.3f})'
        create_snp_window_heatmap(snp_window_contents, title=title, label_idx_dict=label_idxs, ax=ax)

    fig.tight_layout()
    fig.suptitle(f'Cell Type: {ctype}\nTrait: {trait} ({len(trait_snps):,} SNPs)')
    # fig.show()
    fig.savefig(os.path.join(out_dir, 'snp_region_labels.pdf'), format='pdf')


if __name__ == '__main__':
    arg_dict = parse_args()
    main(**arg_dict)
