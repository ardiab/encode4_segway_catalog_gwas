# ENCODE4 Segway Catalog GWAS Analysis

## Overview

This repository provides a comprehensive pipeline for analyzing the enrichment of GWAS (Genome-Wide Association Study) SNPs in functional genomic annotations derived from the ENCODE4 Segway Catalog. The analysis includes calculating the Conservation Annotation Absolute Score (CAAS), intersecting SNPs with genomic annotations, and performing statistical enrichment analyses, with a suite of scripts for downstream visualization and interpretation.

## Main Features
- **CAAS Calculation**: Quantifies the conservation of genomic regions annotated by Segway using phyloP scores.
- **SNP Processing**: Prepares GWAS SNP regions, generates control (fake) SNPs, and intersects SNPs with annotations.
- **Enrichment Analysis**: Calculates metrics and statistical significance for SNP enrichment in functional annotations.
- **Visualization**: Provides scripts for plotting results, including heatmaps, boxplots, and eCDFs.

## Directory Structure
- `encode4_gwas_analysis/caas_calculation/`: Scripts for CAAS computation and phyloP score processing.
- `encode4_gwas_analysis/snp_processing/`: Scripts for preparing SNP regions, intersecting with annotations, and calculating SNP metrics and p-values.
- `encode4_gwas_analysis/analysis/`: Scripts for visualizing results (e.g., heatmaps, boxplots, eCDFs).
- `encode4_gwas_analysis/utils/`: Utility modules for data loading, path management, constants, and enrichment helpers.

## Usage

### 1. CAAS Calculation
Compute CAAS for each biosample:
```bash
python caas_calculation/compute_caas.py --out_dir <output_directory>
```
- `--out_dir`: Output directory containing metadata and where results will be written.

### 2. SNP Processing
- **Create SNP BED files from GWAS data:**
  ```bash
  python snp_processing/create_snp_bed.py --gwas_snp_df_fpath <gwas_file> --out_dir <output_directory>
  ```
- **Intersect SNPs with annotations:**
  ```bash
  python snp_processing/intersect_snps_with_anns.py --out_dir <output_directory>
  ```
- **Calculate SNP metrics:**
  ```bash
  python snp_processing/calculate_snp_metrics.py --out_dir <output_directory>
  ```
- **Calculate SNP p-values:**
  ```bash
  python snp_processing/calculate_snp_pvals.py --out_dir <output_directory> [--snp_cutoff N]
  ```

### 3. Analysis & Visualization
Scripts in the `analysis/` directory generate plots and summary statistics:
- `plot_trait_pval_ranks.py`: Plots GWAS trait p-value ranks across biosamples.
- `plot_test_result_heatmaps.py`: Generates heatmaps of test results.
- `plot_snp_region_appearance.py`: Visualizes SNP region annotation appearance.
- `plot_sig_insig_snp_boxplots.py`: Boxplots for significant/insignificant SNPs.
- `plot_phylop_ecdf.py`, `plot_caas_ecdf.py`: eCDF plots for phyloP/CAAS scores.
- `plot_label_caas_heatmap.py`: Heatmap of CAAS by label.
- `plot_interp_term_freqs.py`: Frequency of interpretation terms.

## Workflow Summary
1. **Prepare metadata and input files.**
2. **Run CAAS calculation.**
3. **Process GWAS SNPs and create SNP regions.**
4. **Intersect SNPs with annotations.**
5. **Calculate enrichment metrics and p-values.**
6. **Visualize and interpret results using analysis scripts.**

## Notes
- MetaHelper (used by multiple scripts) expects a CSV which contains rows linking sample IDs to local paths that the scripts can use to fetch associated .bed files. This file must currently be generated manually, but we will be releasing a script that automatically downloads all the Segway annotations we used in our analysis and generates the corresponding metadata CSV.