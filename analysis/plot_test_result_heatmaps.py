import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from encode_gwas.utils.data_loader import DataLoader
from encode_gwas.utils.meta_helper import MetaHelper
from encode_gwas.utils.path_helper import PathHelper
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', required=True)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict['out_dir']), f'Directory {arg_dict["out_dir"]} does not exist.'

    return arg_dict


# region Overlap
def get_overlap_colormap():
    c1 = (80, 0, 0)
    c2 = (233, 71, 0)
    c3 = (255, 249, 192)
    normed_colors = [[c/255 for c in c1], [c/255 for c in c2], [c/255 for c in c3]]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('', normed_colors)
    cmap.set_under('black')

    return cmap


def get_trait_overlap(trait_snp_linkage_df, sorted_trait_list):
    trait_linkage_sub = trait_snp_linkage_df[trait_snp_linkage_df['trait'].isin(sorted_trait_list)]
    data = []
    for trait in tqdm(sorted_trait_list):
        trait_df = trait_linkage_sub[trait_linkage_sub['trait'] == trait]
        trait_snps = set(trait_df['snp'].values)
        other_traits_df = trait_linkage_sub

        for other_trait, other_trait_df in other_traits_df.groupby('trait'):
            other_trait_snps = set(other_trait_df['snp'].values)
            intersection = trait_snps.intersection(other_trait_snps)
            union = trait_snps.union(other_trait_snps)
            data.append({'trait': trait, 'other_trait': other_trait, 'overlap': len(intersection) / len(trait_snps),
                         'iou': len(intersection) / len(union)})

    overlap_df = pd.DataFrame(data)
    overlap_pivot = overlap_df.pivot(index='trait', columns='other_trait', values='overlap')
    overlap_pivot = overlap_pivot.loc[sorted_trait_list, sorted_trait_list]
    iou_pivot = overlap_df.pivot(index='trait', columns='other_trait', values='iou')
    iou_pivot = iou_pivot.loc[sorted_trait_list, sorted_trait_list]

    return overlap_pivot, iou_pivot
# endregion


# region P-value heatmaps
def get_pval_colormap():
    c1 = (252, 71, 84)
    c2 = (241, 133, 18)
    c3 = (255, 249, 192)
    normed_colors = [[c/255 for c in c1], [c/255 for c in c2], [c/255 for c in c3]]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('', normed_colors)
    cmap.set_under('black')

    return cmap


def get_neg_log_pval_heatmap(pval_df, fig_title, meta_helper, sig_threshold, trait_subset_to_display=None,
                             trait_snp_count_dict=None):
    epsilon = 1e-10
    pval_df_proc = -np.log10(pval_df + epsilon)
    proc_sig_threshold = -np.log10(sig_threshold)
    cbar_title = f'-log10 Bonferroni-corrected P-value (+{epsilon} to avoid log 0)'

    clustermap = sns.clustermap(pval_df_proc.values, dendrogram_ratio=[.05, .05], cmap='magma')
    plt.close(clustermap.figure)
    row_order = clustermap.dendrogram_row.reordered_ind
    col_order = clustermap.dendrogram_col.reordered_ind

    n_traits_to_display = len(trait_subset_to_display) if trait_subset_to_display is not None else pval_df.shape[0]
    fig, ax = plt.subplots(figsize=(50, int(n_traits_to_display/7)), dpi=200)
    cmap = get_pval_colormap()

    clustered_df = pval_df_proc.iloc[row_order, col_order]
    ordered_traits = clustered_df.index.tolist()

    if trait_subset_to_display is not None:
        clustered_df_sub = clustered_df.reset_index()
        clustered_df_sub = clustered_df_sub[clustered_df_sub['trait'].isin(trait_subset_to_display)].set_index('trait')
    else:
        clustered_df_sub = clustered_df

    sns.heatmap(clustered_df_sub.values, cmap=cmap, ax=ax, vmin=proc_sig_threshold,
                linewidths=0.0, rasterized=True, square=True, cbar=False)
    ax.set_xticks(np.arange(0.5, len(col_order), 1))
    xtick_labels = [meta_helper.get_sample_description(sample_id) for sample_id in
                    list(clustered_df_sub.columns)]
    xtick_labels = [l.split('Homo sapiens ')[-1] for l in xtick_labels]
    ax.set_xticklabels(xtick_labels)

    ax.set_yticks(np.arange(0.5, clustered_df_sub.shape[0], 1))
    ytick_labels = clustered_df_sub.index.values
    if trait_snp_count_dict is not None:
        assert set(ytick_labels).issubset(set(trait_snp_count_dict.keys()))
        ytick_labels = [l + f' ({trait_snp_count_dict[l]})' for l in ytick_labels]
    ax.set_yticklabels(ytick_labels)

    ax.tick_params(labelleft=True, labelright=True, labeltop=True, labelbottom=True)
    ax.tick_params(axis='x', rotation=90)
    ax.tick_params(axis='y', rotation=0)
    fig.tight_layout()

    return fig, ordered_traits
# endregion


def draw_mean_caas_boxplot(mean_caas_df, cell_name, trait_name, trait_snp_linkage_df):
    fig, ax = plt.subplots(figsize=(10, 10))

    trait_snps = trait_snp_linkage_df[trait_snp_linkage_df['trait'] == trait_name]['snp'].unique()
    assert len(trait_snps) > 0
    other_snps = trait_snp_linkage_df[~(trait_snp_linkage_df['snp'].isin(trait_snps))]['snp'].unique()

    snp_caas = mean_caas_df.loc[trait_snps][cell_name].tolist()
    other_caas = mean_caas_df.loc[other_snps][cell_name].tolist()
    sns.boxplot(x=['Trait SNPs'] * len(snp_caas) + ['Other SNPs'] * len(other_caas), y=snp_caas + other_caas, ax=ax)
    ax.set_xlabel('SNP Type', fontsize=16, weight='bold')
    ax.set_ylabel('Mean CAAS', fontsize=16, weight='bold')
    ax.set_title(f'{cell_name}, {trait_name}')

    return fig


def get_filtered_trait_snps(snp_region_side_length, snp_cutoff, data_loader):
    trait_snp_linkage_df = data_loader.load_trait_snp_linkage_df(window_size=snp_region_side_length)
    trait_snp_count_df = data_loader.load_trait_snp_count_df(window_size=snp_region_side_length)
    traits_to_keep = trait_snp_count_df[trait_snp_count_df['n_snps'] >= snp_cutoff]['trait'].unique().tolist()
    trait_snp_linkage_df_sub = trait_snp_linkage_df[trait_snp_linkage_df['trait'].isin(traits_to_keep)]

    sorted_traits_desc = trait_snp_count_df[trait_snp_count_df['trait'].isin(traits_to_keep)].sort_values(
        by=['n_snps'], ascending=False)['trait'].tolist()
    sorted_traits_asc = trait_snp_count_df[trait_snp_count_df['trait'].isin(traits_to_keep)].sort_values(
        by=['n_snps'], ascending=True)['trait'].tolist()
    results = {}
    for direction_traits, direction_desc in [(sorted_traits_desc, 'descending'), (sorted_traits_asc, 'ascending')]:
        added_traits = []
        added_snps_total = set()
        added_snps_pairwise = []
        overlap_data = []
        for trait in tqdm(direction_traits):
            trait_df = trait_snp_linkage_df_sub[trait_snp_linkage_df_sub['trait'] == trait]
            trait_snps = set(trait_df['snp'].tolist())
            overlap_union = len(trait_snps.intersection(added_snps_total)) / len(trait_snps)
            if len(added_snps_pairwise) > 0:
                overlap_pairwise = max([len(trait_snps.intersection(s)) / len(trait_snps)
                                        for s in added_snps_pairwise])
            else:
                overlap_pairwise = 0
            overlap_data.append({'trait': trait, 'overlap_union': overlap_union, 'overlap_pairwise': overlap_pairwise})
            added_snps_total = added_snps_total.union(trait_snps)
            added_snps_pairwise.append(trait_snps)
            added_traits.append(trait)

        results[direction_desc] = {
            'traits': added_traits,
            'snps': added_snps_total,
            'overlap_df': pd.DataFrame(overlap_data)
        }

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.hist(results['descending']['overlap_df']['overlap_union'].values, label='Descending -- Union', alpha=0.6)
    ax.hist(results['ascending']['overlap_df']['overlap_union'].values, label='Ascending -- Union', alpha=0.6)
    ax.hist(results['descending']['overlap_df']['overlap_pairwise'].values, label='Descending -- Pairwise', alpha=0.6)
    ax.hist(results['ascending']['overlap_df']['overlap_pairwise'].values, label='Ascending -- Pairwise', alpha=0.6)
    ax.set_xlabel('Proportion SNP overlap with previously added traits')
    ax.set_ylabel('Frequency')
    ax.set_title('Trait filtering by SNP overlap')
    ax.legend()
    fig.show()

    return results


def main(out_dir):
    path_helper = PathHelper(out_dir=out_dir)
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)
    data_loader = DataLoader(out_dir=out_dir)

    overlap_threshold = 0.2
    min_n_significant_req = 1
    window_size = 10_000
    sig_threshold = 0.05

    trait_snp_count_dict = data_loader.load_trait_snp_count_df(window_size=window_size).set_index('trait')['n_snps'].to_dict()
    trait_filter_dict = get_filtered_trait_snps(snp_region_side_length=10_000, snp_cutoff=30, data_loader=data_loader)

    for sort_direction, overlap_filtering_method in [('descending', 'overlap_pairwise')]:
        pval_df = data_loader.load_pval_df(metric_name='mean_caas', test_type='signedrank_per_sample_null',
                                           window_size=window_size, fake=False).fillna(1) * 234

        overlap_df = trait_filter_dict[sort_direction]['overlap_df']
        traits_to_keep_after_overlap_filtering = overlap_df[overlap_df[overlap_filtering_method] <
                                                            overlap_threshold]['trait'].unique()
        traits_to_keep_after_sig_count_filtering = pval_df[(pval_df < sig_threshold).sum(axis=1) >=
                                                           min_n_significant_req].index.unique()

        traits_to_cluster_on = traits_to_keep_after_sig_count_filtering
        trait_subset_to_display = list(set(traits_to_cluster_on).intersection(
            set(traits_to_keep_after_overlap_filtering)))

        pval_df = pval_df.loc[traits_to_cluster_on]
        full_fig_title = f'-log10 P-value (trait X biosample) heatmap\n' \
                    f'Test: mean_caas, signedrank_per_sample_null\n' \
                    f'Overlap filtering method: sorted {sort_direction}, metric {overlap_filtering_method}\n' \
                    f'# traits clustered on: {pval_df.shape[0]:,}, # traits displayed in this heatmap: {pval_df.shape[0]:,}' \
                    f'SNP window side length 10,000\n' \
                    f'Bonferroni-corrected P-value > {sig_threshold}:' \
                    f' black\nBonferroni-corected P-value <= {sig_threshold}: not black'
        fig_full, clustered_traits_full = get_neg_log_pval_heatmap(
            pval_df=pval_df, fig_title=full_fig_title, meta_helper=meta_helper, sig_threshold=sig_threshold,
            trait_subset_to_display=None, trait_snp_count_dict=trait_snp_count_dict)
        fig_full_desc = f'{sort_direction}_{overlap_filtering_method}_full'
        full_fig_out_path = path_helper.get_pval_heatmap_fig_fpath(
            metric_name='mean_caas', test_type='signedrank_per_sample_null', window_size=window_size, fake=False,
            desc=fig_full_desc)
        fig_full.savefig(os.path.join(out_dir, full_fig_out_path.split('/')[-1]), format='pdf')
        fig_full.savefig(os.path.join(out_dir, full_fig_out_path.split('/')[-1]) + '.png')
        plt.close(fig_full)
        partial_fig_title = f'-log10 P-value (trait X biosample) heatmap\n' \
                    f'Test: mean_caas, signedrank_per_sample_null\n' \
                    f'Overlap filtering method: sorted {sort_direction}, metric {overlap_filtering_method}\n' \
                    f'# traits clustered on: {pval_df.shape[0]:,}, # traits displayed in this heatmap: {len(trait_subset_to_display):,}' \
                    f'SNP window side length 10,000\n' \
                    f'Bonferroni-corrected P-value > {sig_threshold}:' \
                    f' black\nBonferroni-corected P-value <= {sig_threshold}: not black'
        fig_partial, clustered_traits_partial = get_neg_log_pval_heatmap(
            pval_df=pval_df, fig_title=partial_fig_title, meta_helper=meta_helper, sig_threshold=sig_threshold,
            trait_subset_to_display=trait_subset_to_display, trait_snp_count_dict=trait_snp_count_dict)
        fig_partial_desc = f'{sort_direction}_{overlap_filtering_method}_partial'
        partial_fig_out_path = path_helper.get_pval_heatmap_fig_fpath(
            metric_name='mean_caas', test_type='signedrank_per_sample_null', window_size=window_size, fake=False,
            desc=fig_partial_desc)
        fig_partial.savefig(os.path.join(out_dir,
                                         partial_fig_out_path.split('/')[-1]), format='pdf')
        fig_partial.savefig(os.path.join(out_dir,
                                         partial_fig_out_path.split('/')[-1]) + '.png')
        plt.close(fig_partial)


if __name__ == '__main__':
    arg_dict = parse_args()
    main(**arg_dict)
