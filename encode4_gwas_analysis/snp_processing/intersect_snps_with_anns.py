import argparse
import os

from encode4_gwas_analysis.snp_processing import run_snp_ann_intersection_split
from encode4_gwas_analysis.utils.meta_helper import MetaHelper
from encode4_gwas_analysis.utils.path_helper import PathHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", required=True)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict["out_dir"])

    return arg_dict


def write_sample_id_list(sample_id_list, out_fpath):
    with open(out_fpath, "w") as out_f:
        for sample_id in sample_id_list:
            out_f.write(sample_id + "\n")


def run_snp_ann_intersection(path_helper, meta_helper, out_dir):
    sample_ids = meta_helper.get_all_sample_ids()
    sample_id_list_fpath = os.path.join(
        path_helper.snp_intersection_dir, "single_split_sample_id_list.txt"
    )
    write_sample_id_list(sample_ids, out_fpath=sample_id_list_fpath)
    run_snp_ann_intersection_split.main(out_dir=out_dir, sample_id_list_fpath=sample_id_list_fpath)



if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    meta_helper = MetaHelper.from_csv(path_helper.meta_df_fpath)

    snp_intersection_dir = path_helper.snp_intersection_dir
    assert not os.path.isdir(snp_intersection_dir), (
        f"SNP/Annotation intersection directory {snp_intersection_dir} already exists"
    )
    os.mkdir(snp_intersection_dir)

    run_snp_ann_intersection(
        path_helper=path_helper, meta_helper=meta_helper, out_dir=arg_dict["out_dir"]
    )
