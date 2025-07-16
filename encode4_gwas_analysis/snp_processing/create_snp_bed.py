import argparse
import itertools
import os
import re

import numpy as np
import pandas as pd
import tqdm

import encode4_gwas_analysis.utils.constants as const
from encode4_gwas_analysis.utils.path_helper import PathHelper


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas_snp_df_fpath", required=True)
    parser.add_argument("--out_dir", required=True)

    arg_dict = dict(vars(parser.parse_args()))
    assert os.path.isfile(arg_dict["gwas_snp_df_fpath"])
    assert os.path.isdir(arg_dict["out_dir"])

    return arg_dict


def format_gwas_df(gwas_df):
    """
    Process the raw GWAS SNP DF -- the main functionality here is splitting up rows that correspond to multiple SNPs and
    chromosomes so that a single row corresponds to a single SNP.
    """
    print(f"Raw GWAS DF shape: {gwas_df.shape}")
    assert (gwas_df["CHR_POS"].isnull() == gwas_df["CHR_ID"].isnull()).all()
    print(
        f'{gwas_df["CHR_POS"].isnull().mean() * 100:.3f}% of rows have null values for "CHR_POS" and "CHR_ID".'
    )
    gwas_df = gwas_df[gwas_df["CHR_POS"].notnull()]
    # Some rows contain multiple SNPs, which can show up as positions separated by a space, a semicolon, or an "x".
    # For example, 1x10 would specify one SNP that is at position 1 and another at position 10. A similar pattern can
    # be observed for "CHR_ID", which encodes the chromosome on which the SNP falls.
    gwas_df["CHR_POS"] = gwas_df["CHR_POS"].apply(
        lambda s: re.split(r"[ ;x]", str(s).replace(" ", ""))
    )
    gwas_df["CHR_ID"] = gwas_df["CHR_ID"].apply(
        lambda s: re.split(r"[ ;x]", str(s).replace(" ", ""))
    )
    for idx, row in gwas_df.iterrows():
        chr_pos, chr_id = row["CHR_POS"], row["CHR_ID"]
        assert type(chr_pos) == type(chr_id)
        if isinstance(chr_pos, list):
            assert len(chr_pos) == len(chr_id)
    gwas_df = gwas_df.explode(["CHR_POS", "CHR_ID"])

    assert gwas_df[["CHR_POS", "CHR_ID", "DISEASE/TRAIT", "P-VALUE"]].notnull().all().all()
    gwas_df["CHR_POS"] = gwas_df["CHR_POS"].astype(float).astype(int)
    gwas_df["DISEASE/TRAIT"] = gwas_df["DISEASE/TRAIT"].str.lower()

    # Check that all values for CHR_ID are either in the range [1, 22], or are X or Y chromosomes
    assert gwas_df["CHR_ID"].isin([str(c) for c in range(1, 23)] + ["X", "Y"]).all()
    assert gwas_df["CHR_POS"].apply(lambda p: str(p).isalnum()).all()

    print(f"Processed GWAS DF shape: {gwas_df.shape}")

    return gwas_df


def create_snp_region(snp_loc, snp_region_side_length, chrom):
    """
    Create a SNP region, which is a genomic region that is centered on the GWAS-identified SNP position.
    Note that the start position calculated is inclusive, and the endpoint is exclusive.

    For example, a SNP location of 100 means that the SNP occurs at the 101th base pair, i.e. there are 100
    base pairs before it. Using a side length of 100, this function would calculate:
    - Start position: (100 - 100) = 0 (first position using a 0-indexed coordinate system).
    - End position: (100 + 100 + 1) = 201 (202nd position using a 0-indexed coordinate system). However, since bed files
        are endpoint-exclusive, the specified end point will actually be the 201st position.
    """
    # Make sure the starting position of the region is non-negative.
    start_pos = max(0, snp_loc - snp_region_side_length)
    # Make sure the end position of the region does not exceed the chromosome length.
    end_pos = min(snp_loc + snp_region_side_length + 1, const.CHROM_LENGTHS[chrom])

    return start_pos, end_pos


def create_snp_bed_from_gwas_df(formatted_gwas_df, snp_region_side_length):
    """
    Convert the SNPs from the processed GWAS SNP DF into a DataFrame formatted like a bed file, where each SNP defines a
    window specified by the argument snp_window_side_length (i.e. the location of the SNP and a bidirectional buffer
    with the specified side length).
    """
    # MHC in hg38 -- from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
    exclude_ranges = {"6": range(28_510_120, 33_480_577)}

    dropped_mhc = []
    dropped_overlap = []
    snp_data = []
    trait_linkage_data = []
    # Given the specified SNP window side length, a given trait might have SNP windows that overlap. After discussion,
    # it was decided that one of the two possible types of overlap is okay.
    #   - If a given published SNP falls within the window of another published SNP (i.e. a published SNP + the LD
    #     window defined by snp_window_side_length), then only one of the two SNPs is kept (the one with the lower
    #     p-value).
    #     For example, letting X represent the published SNP and - represent the LD window:
    #      SNP 1: ----X----
    #      SNP 2:   ----X----
    #      Only the SNP with the lower p-value is kept.
    #   - If the windows of two published SNPs overlap, but the actual published SNPs do not fall within each other's
    #     windows, then both are kept.
    #     For example, letting X represent the published SNP and - represent the LD window:
    #      SNP 1: ----X----
    #      SNP 2:      ----X----
    #      Both SNPs are kept.
    #
    # This loop iterates over SNPs grouped by the trait they are associated with and the chromosome they fall on,
    # because we only want to do the filtering step described above for SNPs within a given trait. Note that this means
    # that a SNP will not necessarily get filtered out if it is associated with multiple traits, since we only do the
    # filtering for one trait at a time.
    assert formatted_gwas_df["P-VALUE"].dtype == np.float64
    trait_chrom_group = formatted_gwas_df.groupby(["DISEASE/TRAIT", "CHR_ID"])

    for (trait, chrom), group_df in tqdm.tqdm(
        trait_chrom_group,
        total=trait_chrom_group.ngroups,
        desc=f"Creating trait SNP regions -- {snp_region_side_length} BP window",
        miniters=trait_chrom_group.ngroups / 100,
    ):
        trait_chrom_added_snp_ranges = []
        # Since the following steps will iteratively add non-overlapping SNPs for the trait, SNPs are sorted in order
        # of ascending p-value so that lower p-value SNPs are always added first. This allows the loop to
        # (automatically) select the SNPs with lowest p-values whenever overlap occurs without having to keep track of
        # all overlapping SNPs.
        group_df_sorted = group_df.sort_values(by="P-VALUE", ascending=True)
        for idx, row in group_df_sorted.iterrows():
            # Important: The line below converts the SNP from a 1-indexed coordinate system to a 0-indexed coordinate
            # system.
            snp_loc = int(row["CHR_POS"]) - 1
            snp_name = f"snp_chr{chrom}_{snp_loc}"
            # Skip the SNP if it occurs within one of the specified exclusion ranges.
            if any(
                [((str(chrom) == str(k)) and (snp_loc in v)) for k, v in exclude_ranges.items()]
            ):
                dropped_mhc.append(snp_name)
                continue
            # Skip the SNP if it falls within the window of an already-added (lower p-value) SNP.
            if any([snp_loc in r for r in trait_chrom_added_snp_ranges]):
                dropped_overlap.append(snp_name)
                continue

            # If the SNP doesn't overlap with an already-added SNP, add it to the overall list of SNPs (and also keep
            # track of its window so that any subsequent overlapping higher-p-value SNP is skipped).
            start_pos, end_pos = create_snp_region(
                snp_loc=snp_loc, snp_region_side_length=snp_region_side_length, chrom=f"chr{chrom}"
            )
            trait_chrom_added_snp_ranges.append(range(start_pos, end_pos))

            snp_data.append(
                {
                    "name": snp_name,
                    "chrom": f"chr{chrom}",
                    "chromStart": start_pos,
                    "chromEnd": end_pos,
                }
            )
            trait_linkage_data.append({"trait": trait, "snp": snp_name})

    snp_bed_df = pd.DataFrame(snp_data)[
        ["chrom", "chromStart", "chromEnd", "name"]
    ].drop_duplicates()
    trait_linkage_df = pd.DataFrame(trait_linkage_data)
    assert set(snp_bed_df["name"]) == set(trait_linkage_df["snp"])
    print(f"\nSNP filtering stats for side length {snp_region_side_length}:")
    print(
        f"Dropped SNPs {len(dropped_mhc):,} times ({len(set(dropped_mhc)):,} unique SNPs) from a trait because they "
        f"occurred in the MHC."
    )
    print(
        f"Dropped SNPs {len(dropped_overlap):,} times ({len(set(dropped_overlap)):,} unique SNPs) from a trait "
        f"because over 50% of their region overlapped with a previously-added SNP for the trait which had a lower "
        f"P-value. Note that this does not necessarily mean that any of these SNPs were dropped from all traits."
    )
    n_snps_before_filtering = (
        "snp_chr" + formatted_gwas_df["CHR_ID"] + "_" + formatted_gwas_df["CHR_POS"].astype(str)
    ).nunique()
    n_snps_after_filtering = snp_bed_df["name"].nunique()
    print(
        f"The number of unique SNPs went from {n_snps_before_filtering:,} before filtering to "
        f"{n_snps_after_filtering:,} after filtering ({n_snps_after_filtering - n_snps_before_filtering:,})"
    )
    n_traits_before_filtering = formatted_gwas_df["DISEASE/TRAIT"].nunique()
    n_traits_after_filtering = trait_linkage_df["trait"].nunique()
    print(
        f"The number of unique traits went from {n_traits_before_filtering:,} before filtering to "
        f"{n_traits_after_filtering:,} after filtering ({n_traits_after_filtering - n_traits_before_filtering:,})"
    )

    return snp_bed_df, trait_linkage_df


def create_fake_snp_bed(
    snp_bed, snp_window_side_length, existing_fake_snp_bed=None, offset=100_000, seed=1234567
):
    # Note: This does not check for overlap in SNP regions within a trait. A trait with non-overlapping real SNP regions
    # can get overlapping fake SNP regions. If we were to filter overlapping fake SNP regions, this would potentially
    # result in a number of fake SNPs per trait that is smaller than the number of real SNPs. Since the fake SNPs are
    # generated in the proximity of real (filtered) SNPs, we will just use them as-is without filtering for overlap
    # again.
    np.random.seed(seed)
    fake_snp_data = []
    bed_cp = snp_bed.copy(deep=True)
    if existing_fake_snp_bed is not None:
        assert not existing_fake_snp_bed["name"].duplicated().any()
        existing_fake_snp_bed = existing_fake_snp_bed.set_index("name")
    for idx, row in tqdm.tqdm(bed_cp.iterrows(), desc="Creating fake SNPs", total=bed_cp.shape[0]):
        snp_name = row["name"]
        snp_loc = int(snp_name.split("_")[-1])
        snp_chrom = row["chrom"]

        if (existing_fake_snp_bed is not None) and (
            (snp_name + "_fake") in existing_fake_snp_bed.index
        ):
            existing_fake_snp_row = existing_fake_snp_bed.loc[snp_name + "_fake"]
            assert existing_fake_snp_row["chrom"] == snp_chrom
            fake_snp_loc = existing_fake_snp_row["fake_snp_loc"]
            assert abs(fake_snp_loc - snp_loc) == offset
        else:
            if (
                snp_loc <= offset
            ):  # Subtracting the offset would result in a negative position for the SNP
                options = [offset]
            elif (snp_loc + offset) >= const.CHROM_LENGTHS[snp_chrom]:
                options = [-offset]
            else:
                options = [offset, -offset]
            sampled_offset = np.random.choice(options)
            fake_snp_loc = snp_loc + sampled_offset

        fake_snp_start_pos, fake_snp_end_pos = create_snp_region(
            snp_loc=fake_snp_loc, snp_region_side_length=snp_window_side_length, chrom=snp_chrom
        )
        fake_snp_data.append(
            {
                "chrom": snp_chrom,
                "chromStart": fake_snp_start_pos,
                "chromEnd": fake_snp_end_pos,
                "name": snp_name + "_fake",
                "fake_snp_loc": fake_snp_loc,
            }
        )

    return pd.DataFrame(fake_snp_data)[["chrom", "chromStart", "chromEnd", "name", "fake_snp_loc"]]


def create_snp_beds(gwas_df, snp_region_side_lengths):
    gwas_df_formatted = format_gwas_df(gwas_df)

    real_beds = {}
    fake_beds = {}
    trait_linkage_dfs = {}

    sorted_snp_region_side_lengths = sorted(snp_region_side_lengths)
    for idx, side_length in enumerate(sorted_snp_region_side_lengths):
        snp_bed, trait_linkage_df = create_snp_bed_from_gwas_df(
            formatted_gwas_df=gwas_df_formatted.copy(deep=True), snp_region_side_length=side_length
        )
        if idx == 0:
            existing_fake_snp_bed = None
        else:
            # Get the locations of fake SNPs that have already been generated.
            existing_fake_snp_bed = pd.concat(
                [df.drop(columns=["chromStart", "chromEnd"]) for df in fake_beds.values()]
            ).drop_duplicates()
            assert not existing_fake_snp_bed["name"].duplicated().any()
        # The function create_fake_snp_bed reuses previously-generated fake SNPs when possible.
        fake_snp_bed = create_fake_snp_bed(
            snp_bed=snp_bed,
            snp_window_side_length=side_length,
            existing_fake_snp_bed=existing_fake_snp_bed,
        )

        real_beds[side_length] = snp_bed
        fake_beds[side_length] = fake_snp_bed
        trait_linkage_dfs[side_length] = trait_linkage_df

    # Check that common SNPs have the same fake location in every window size.
    for ix1, ix2 in itertools.product(range(len(fake_beds)), range(len(fake_beds))):
        bed_1 = fake_beds[sorted_snp_region_side_lengths[ix1]]
        bed_2 = fake_beds[sorted_snp_region_side_lengths[ix2]]
        common_snps = set(bed_1["name"]).intersection(set(bed_2["name"]))
        assert (
            bed_1.set_index("name").loc[list(common_snps)]["fake_snp_loc"]
            == bed_2.set_index("name").loc[list(common_snps)]["fake_snp_loc"]
        ).all()

    return real_beds, fake_beds, trait_linkage_dfs


if __name__ == "__main__":
    arg_dict = parse_args()
    path_helper = PathHelper(out_dir=arg_dict["out_dir"])
    gwas_snp_df = pd.read_csv(arg_dict["gwas_snp_df_fpath"], sep="\t")

    real_beds, fake_beds, trait_linkage_dfs = create_snp_beds(
        gwas_df=gwas_snp_df, snp_region_side_lengths=const.SNP_REGION_SIDE_LENGTHS
    )

    assert real_beds.keys() == fake_beds.keys() == trait_linkage_dfs.keys()
    processed_snp_dir = path_helper.processed_snp_dir
    os.mkdir(processed_snp_dir)
    for window_size in real_beds.keys():
        window_real_bed = real_beds[window_size][["chrom", "chromStart", "chromEnd", "name"]]
        window_real_bed["chr_proc"] = window_real_bed["chrom"].apply(
            lambda c: f"{int(c.lstrip('chr')):02}"
            if c.lstrip("chr").isnumeric()
            else c.lstrip("chr")
        )
        window_real_bed = window_real_bed.sort_values(
            by=["chr_proc", "chromStart"], ascending=True
        ).drop(columns=["chr_proc"])

        window_fake_bed = fake_beds[window_size][["chrom", "chromStart", "chromEnd", "name"]]
        window_fake_bed["chr_proc"] = window_fake_bed["chrom"].apply(
            lambda c: f"{int(c.lstrip('chr')):02}"
            if c.lstrip("chr").isnumeric()
            else c.lstrip("chr")
        )
        window_fake_bed = window_fake_bed.sort_values(
            by=["chr_proc", "chromStart"], ascending=True
        ).drop(columns=["chr_proc"])

        assert set(window_real_bed["name"].values) == set(
            window_fake_bed["name"].apply(lambda s: s.replace("_fake", "")).values
        )

        window_real_bed.to_csv(
            path_helper.get_snp_bed_fpath(window_size=window_size, fake=False),
            index=False,
            header=None,
            sep="\t",
        )
        window_fake_bed.to_csv(
            path_helper.get_snp_bed_fpath(window_size=window_size, fake=True),
            index=False,
            header=None,
            sep="\t",
        )
        trait_linkage_dfs[window_size].to_csv(
            path_helper.get_trait_snp_linkage_df_fpath(window_size=window_size)
        )
