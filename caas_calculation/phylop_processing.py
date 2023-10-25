import os
import re
import subprocess
import warnings

import encode_gwas.utils.constants as const
import numpy as np
import tqdm
from encode_gwas.utils.path_helper import PathHelper


def download_and_unzip_phylop_wigfix_files(out_dir):
    # Download PhyloP wigFix files from UCSC website
    for wig in tqdm.tqdm(const.PHYLOP_WIGS, desc='Downloading PhyloP wigFix files'):
        subprocess.call(['wget', '-q', '-P', out_dir, f'{const.PHYLOP_BASE_URL}/{wig}.gz'], stderr=subprocess.STDOUT)
        subprocess.call(['gzip', '-d', f'{os.path.join(out_dir, wig)}.gz'], stderr=subprocess.STDOUT,
                        stdout=subprocess.PIPE)


def wiggle_to_numpy(wiggle_fpath):
    # Read the wigFix file as a text file.
    with open(wiggle_fpath, 'r') as in_f:
        wiggle_data = [l for l in in_f.read().split('\n') if len(l) > 0]

    # Isolate the "information" lines which have format "fixedStep chrom=chr1 start=10701 step=1", keeping track of the
    # index in the parsed list where the information line occurs (knowing the index simplifies the process of finding
    # all PhyloP values corresponding to the "block" specified by an information line -- for example, if I know that the
    # first information line occurs at index 0 and the next one occurs at index 10,000, then I know that the PhyloP
    # scores belonging to the block of the first information line start at index 1 and end at index 9,999.
    # Note that in general, wigFix files are structured as
    # fixedStep chrom=chr1 start=10701 step=1
    # val 1
    # val 2
    # ...
    # val N
    # fixedStep chrom=chr1 start=20841 step=1
    # val 1
    # val 2
    # ...
    info_lines = [(idx, l) for idx, l in enumerate(wiggle_data) if 'chrom' in l]
    # Check that all entries in the wigFix file belong to the same chromosome, and that all step sizes are 1. While this
    # might not be true for general wiggle files, it is true for the data that we're using, and the code relies on this
    # assumption.
    chrom = re.findall(r'chrom=(.+?)\b', info_lines[0][1])[0]
    assert chrom in const.CHROM_LENGTHS.keys()
    assert all(['step=1' in l[1] for l in info_lines])
    assert all([f'chrom={chrom}' in l[1] for l in info_lines])
    assert all(['start=' in l[1] for l in info_lines])
    # Create an array of null values which has the same length as the chromosome in question assuming the hg38 assembly.
    chrom_array = np.empty(const.CHROM_LENGTHS[chrom], dtype=np.float64)
    chrom_array[:] = np.nan

    # If we throw away the information lines, we should have at most as many PhyloP scores as the length of the
    # chromosome in question. In general, not all positions will have PhyloP scores (due to e.g. unmappability), and the
    # array will have some null values after adding all PhyloP scores.
    assert (len(wiggle_data) - len(info_lines)) <= const.CHROM_LENGTHS[chrom]

    # Go over pairs of info lines, grabbing all PhyloP scores in between and inserting them into the array at the
    # appropriate locations.
    for ix, wiggle_line_start_idx in enumerate([info_line[0] for info_line in info_lines]):
        if ix < len(info_lines) - 1:
            wiggle_line_end_idx = info_lines[ix + 1][0]
            phylop_vals = [float(val) for val in wiggle_data[wiggle_line_start_idx+1:wiggle_line_end_idx]]
        else:
            phylop_vals = [float(val) for val in wiggle_data[wiggle_line_start_idx+1:]]
        # The position in the chromosome where the block starts (obtained from the information line which has the format
        # "fixedStep chrom=chr1 start=20841 step=1"). This is one-indexed, so subtract 1 to get a 0-indexed coordinate.
        chrom_start_idx = int(re.findall(r'start=(.+?)\b', wiggle_data[wiggle_line_start_idx])[0]) - 1
        chrom_array[chrom_start_idx:chrom_start_idx + len(phylop_vals)] = phylop_vals

    print(f'{chrom}: {np.isnan(chrom_array).mean() * 100:.3f}% null values')

    return chrom, chrom_array


def run_phylop_processing(phylop_wiggle_fpath, phylop_numpy_fpath, chrom):
    """
    Convert a specified wigFix file into a numpy array, then save the numpy array to the specified output directory.
    """
    chrom_extracted, proc_phylop_array = wiggle_to_numpy(phylop_wiggle_fpath)
    assert chrom_extracted == chrom

    with open(phylop_numpy_fpath, 'wb') as out_f:
        np.save(out_f, proc_phylop_array)


def process_phylop(base_out_dir):
    path_helper = PathHelper(out_dir=base_out_dir)
    phylop_dir = path_helper.phylop_dir
    try:
        assert not os.path.isdir(phylop_dir)
        os.mkdir(phylop_dir)
        download_and_unzip_phylop_wigfix_files(phylop_dir)
        phylop_wiggle_files = {chrom: path_helper.get_phylop_wiggle_fpath(chrom)
                               for chrom in const.CHROM_LENGTHS.keys()}
        assert all([os.path.isfile(f) for f in phylop_wiggle_files.values()])
        for chrom, wiggle_path in tqdm.tqdm(phylop_wiggle_files.items(),
                                            desc='Converting wiggle files to numpy arrays'):
            run_phylop_processing(phylop_wiggle_fpath=wiggle_path,
                                  phylop_numpy_fpath=path_helper.get_phylop_numpy_fpath(chrom=chrom),
                                  chrom=chrom)
    except AssertionError:
        phylop_numpy_fpaths = {chrom: path_helper.get_phylop_numpy_fpath(chrom)
                               for chrom in const.CHROM_LENGTHS.keys()}
        if all([os.path.isfile(fp) for fp in phylop_numpy_fpaths.values()]):
            warnings.warn('All phyloP numpy files were found. Skipping phyloP processing step...')
        else:
            raise AssertionError(f'Unable to find all PhyloP beds in {phylop_dir}')
