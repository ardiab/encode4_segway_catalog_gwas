import argparse
import os

import numpy as np
from encode_gwas.caas_calculation import run_caas_computation_split
from encode_gwas.caas_calculation.phylop_processing import process_phylop
from encode_gwas.utils.meta_helper import MetaHelper
from encode_gwas.utils.path_helper import PathHelper
from slurm_utils.splits.split import run_slurm_splits, Split


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--enable_cluster_mode', action='store_true', default=False)
    parser.add_argument('--cluster_mode_split_count', type=int, default=40)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict['out_dir'])

    return arg_dict


def write_sample_id_list(sample_id_list, out_fpath):
    with open(out_fpath, 'w') as out_f:
        for sample_id in sample_id_list:
            out_f.write(sample_id + '\n')


def run_caas_computation_single_split(path_helper, meta_helper, out_dir):
    sample_ids = meta_helper.get_all_sample_ids()
    sample_id_list_fpath = os.path.join(path_helper.label_info_dir, 'single_split_sample_id_list.txt')
    write_sample_id_list(sample_ids, out_fpath=sample_id_list_fpath)
    run_caas_computation_split.main(out_dir=out_dir, sample_id_list_fpath=sample_id_list_fpath)


@run_slurm_splits(sbatch_args={'max_time': '06:00:00', 'memory': '64G'})
def run_caas_computation_multi_split(path_helper, meta_helper, out_dir, split_count):
    """
    Splits the Segway annotation .bed files into the specified number of subgroups, then runs the CAAS script on the
    .bed files in each subgroup as a separate batch job.
    """
    splits = []
    sample_id_splits = np.array_split(meta_helper.get_all_sample_ids(), split_count)

    for sample_id_split in sample_id_splits:
        split = Split(base_output_dir=path_helper.label_info_dir)
        split_dir = split.setup()
        split_sample_id_list_fpath = os.path.join(split_dir, 'split_sample_id_list.txt')
        write_sample_id_list(sample_id_split, out_fpath=split_sample_id_list_fpath)

        split_command = f'python {os.path.join(os.path.dirname(__file__), "run_caas_computation_split.py")} ' \
                        f'--sample_id_list_fpath {split_sample_id_list_fpath} --out_dir {out_dir}'
        split.set_command(split_command)
        splits.append(split)

    return path_helper.label_info_dir, splits


def main(out_dir, enable_cluster_mode, cluster_mode_split_count):
    path_helper = PathHelper(out_dir=out_dir)
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)

    process_phylop(base_out_dir=out_dir)

    label_info_dir = path_helper.label_info_dir
    assert not os.path.isdir(label_info_dir), f'Label info directory {label_info_dir} already exists'
    os.mkdir(label_info_dir)

    if enable_cluster_mode:
        run_caas_computation_multi_split(path_helper=path_helper, meta_helper=meta_helper, out_dir=out_dir,
                                         split_count=cluster_mode_split_count)
    else:
        run_caas_computation_single_split(path_helper=path_helper, meta_helper=meta_helper, out_dir=out_dir)


if __name__ == '__main__':
    arg_dict = parse_args()
    main(out_dir=arg_dict['out_dir'], enable_cluster_mode=arg_dict['enable_cluster_mode'],
         cluster_mode_split_count=arg_dict['cluster_mode_split_count'])
