import argparse
import os

import encode_gwas.utils.constants as const
import tqdm
from encode_gwas.utils.data_loader import DataLoader
from encode_gwas.utils.path_helper import PathHelper


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--sample_id_list_fpath', required=True)

    arg_dict = vars(parser.parse_args())
    assert os.path.isdir(arg_dict['out_dir'])
    assert os.path.isfile(arg_dict['sample_id_list_fpath'])

    return arg_dict


def run_intersection(sample_ann_bed_fpath, real_snp_bed_fpath, fake_snp_bed_fpath, snp_region_side_length, sample_id,
                     path_helper):
    out_real_snp_intersection_bed_fpath = path_helper.get_sample_intersection_bed_fpath(
        sample_id=sample_id, window_size=snp_region_side_length, zipped=False, fake=False)
    # Intersection of SNPs with Segway annotations
    os.system(f'bedtools intersect -a {sample_ann_bed_fpath} -b {real_snp_bed_fpath}  '
              f'-wo > {out_real_snp_intersection_bed_fpath}')

    # Intersection of random locations in the genome with Segway annotations
    out_fake_snp_intersection_bed_fpath = path_helper.get_sample_intersection_bed_fpath(
        sample_id=sample_id, window_size=snp_region_side_length, zipped=False, fake=True)
    os.system(f'bedtools intersect -a {sample_ann_bed_fpath} -b {fake_snp_bed_fpath}  '
              f'-wo > {out_fake_snp_intersection_bed_fpath}')


def zip_intersection_beds(sample_id, snp_region_side_length, path_helper):
    real_snp_intersection_bed_fpath = path_helper.get_sample_intersection_bed_fpath(sample_id=sample_id,
                                                                                    window_size=snp_region_side_length,
                                                                                    zipped=False, fake=False)
    os.system(f'gzip {real_snp_intersection_bed_fpath}')

    fake_snp_intersection_bed_fpath = path_helper.get_sample_intersection_bed_fpath(sample_id=sample_id,
                                                                                    window_size=snp_region_side_length,
                                                                                    zipped=False, fake=True)
    os.system(f'gzip {fake_snp_intersection_bed_fpath}')


def get_snp_bed_label_distribution(bed_fpath):
    intersection_bed = DataLoader.load_intersection_bed(fpath=bed_fpath)
    snp_label_counts = intersection_bed.groupby(['snp_name', 'ann_label'])['intersection_length'].sum().reset_index()
    snp_label_counts = snp_label_counts.pivot(index='snp_name', columns='ann_label',
                                              values='intersection_length').fillna(0)

    return snp_label_counts


def get_snp_region_label_distributions(sample_id, snp_region_side_length, path_helper):
    out_real_snp_intersection_bed_fpath = path_helper.get_sample_intersection_bed_fpath(
        sample_id=sample_id, window_size=snp_region_side_length, zipped=False, fake=False)
    out_fake_snp_intersection_bed_fpath = path_helper.get_sample_intersection_bed_fpath(
        sample_id=sample_id, window_size=snp_region_side_length, zipped=False, fake=True)

    real_snps_label_distr = get_snp_bed_label_distribution(bed_fpath=out_real_snp_intersection_bed_fpath)
    fake_snps_label_distr = get_snp_bed_label_distribution(bed_fpath=out_fake_snp_intersection_bed_fpath)

    real_snps_label_distr_fpath = path_helper.get_label_distr_df_fpath(sample_id=sample_id,
                                                                       window_size=snp_region_side_length,
                                                                       fake=False)
    fake_snps_label_distr_fpath = path_helper.get_label_distr_df_fpath(sample_id=sample_id,
                                                                       window_size=snp_region_side_length,
                                                                       fake=True)
    real_snps_label_distr.to_pickle(real_snps_label_distr_fpath)
    fake_snps_label_distr.to_pickle(fake_snps_label_distr_fpath)


def main(out_dir, sample_id_list_fpath):
    path_helper = PathHelper(out_dir=out_dir)

    with open(sample_id_list_fpath, 'r') as in_f:
        sample_ids = [sid for sid in in_f.read().split('\n') if len(sid) > 0]

    for sample_id in tqdm.tqdm(sample_ids, desc='Intersecting SNPs'):
        sample_out_dir = path_helper.get_sample_intersection_dir(sample_id)
        assert not os.path.isdir(sample_out_dir)
        os.mkdir(sample_out_dir)

        sample_ann_bed_fpath = path_helper.get_sample_bed_fpath(sample_id=sample_id)

        for snp_region_side_length in const.SNP_REGION_SIDE_LENGTHS:
            fake_snp_bed_fpath = path_helper.get_snp_bed_fpath(window_size=snp_region_side_length, fake=True)
            real_snp_bed_fpath = path_helper.get_snp_bed_fpath(window_size=snp_region_side_length, fake=False)
            assert os.path.isfile(real_snp_bed_fpath)
            assert os.path.isfile(fake_snp_bed_fpath)

            run_intersection(sample_ann_bed_fpath=sample_ann_bed_fpath, real_snp_bed_fpath=real_snp_bed_fpath,
                             fake_snp_bed_fpath=fake_snp_bed_fpath, snp_region_side_length=snp_region_side_length,
                             sample_id=sample_id, path_helper=path_helper)
            get_snp_region_label_distributions(sample_id=sample_id, snp_region_side_length=snp_region_side_length,
                                               path_helper=path_helper)
            zip_intersection_beds(sample_id=sample_id, snp_region_side_length=snp_region_side_length,
                                  path_helper=path_helper)


if __name__ == '__main__':
    arg_dict = parse_args()
    main(**arg_dict)
