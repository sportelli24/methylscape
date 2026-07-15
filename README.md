# Targeted region-level methylation analysis

This repository contains the R code used to summarise DNA methylation across selected promoter, gene-body, and PCR amplicon regions for four genes of interest: **ALX1**, **MACF1**, **IGSF3**, and **CLPTM1L**.

The workflow compares one tumour methylation dataset, generated from collapsed tumour replicates, against three normal breast tissue methylation datasets. For each predefined genomic region, the script calculates regional methylation beta values and reports tumour-versus-normal methylation differences.

All coordinates used in this analysis are based on **hg19 / GRCh37**.

## Repository structure

```text
.
├── dmr_beta_analysis_pcr.R
├── amplicons_pcr_hg19.bed
├── cancer_breast_tissue/
│   ├── regions_of_interest.bed
│   └── output/
│       └── maxi_1112_regions_merged_cpg.tsv
└── healthy_breast_tissue/
    └── outputs/
        ├── normal_1_GSM5652347_regions.tsv
        ├── normal_2_GSM5652348_regions.tsv
        └── normal_3_GSM5652349_regions.tsv
```

If your file or folder names differ, update the file paths at the top of `dmr_beta_analysis_pcr.R`.

## Requirements

The script requires R and the following packages:

```r
GenomicRanges
rtracklayer
dplyr
readr
tibble
tidyr
```

Install required packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GenomicRanges", "rtracklayer"))
install.packages(c("dplyr", "readr", "tibble", "tidyr"))
```

## Input files

| File | Description |
|---|---|
| `regions_of_interest.bed` | Promoter and gene-body coordinates for genes of interest |
| `amplicons_pcr_hg19.bed` | PCR amplicon coordinates in hg19 |
| `maxi_1112_regions_merged_cpg.tsv` | Tumour CpG methylation table |
| `normal_*_regions.tsv` | Normal breast tissue CpG methylation tables |

Each CpG methylation table should be tab-delimited and contain the following columns:

```text
chr    pos    beta    cov
```

where:

- `chr` is the chromosome name, e.g. `chr1`
- `pos` is the CpG coordinate
- `beta` is the methylation beta value
- `cov` is the CpG coverage

The BED files must use the same chromosome naming style and genome build as the CpG methylation tables.

## Running the analysis

From the repository root, run:

```r
source("dmr_beta_analysis_pcr.R")
```

or from the terminal:

```bash
Rscript dmr_beta_analysis_pcr.R
```

Before running, check the file paths defined near the top of the script:

```r
regions_bed_path
regions_bed_amplicon
tumour_tsv
normal1_tsv
normal2_tsv
normal3_tsv
```

## Outputs

The default output prefix is:

```r
dmr_beta_4genes
```

The script generates:

| Output file | Description |
|---|---|
| `dmr_beta_4genes_region_summary.csv` | Per-sample methylation summary for each promoter, gene-body, and amplicon region |
| `dmr_beta_4genes_summary_vs_normals.csv` | Tumour methylation compared with the mean of normal samples |
| `dmr_beta_4genes_summary_amplicon_results.csv` | Amplicon-specific methylation summary |

The main comparison values are reported as tumour minus mean-normal beta values. Positive values indicate higher methylation in tumour relative to normal; negative values indicate lower methylation in tumour relative to normal.

## Notes

This workflow is a targeted regional methylation summary for selected genes and amplicons. It is intended to support locus-specific analysis and validation, rather than replace genome-wide differential methylation analysis.

All input files should be checked for consistent genome build, chromosome naming, and beta-value scale before running the script.

## Citation
If you use this repository, please cite: Zhang, Z., Ahmed, E., Constantin, N., Lu, J., Korbie,D.A., Wuethrich, A., Sina, A., Trau, M. (2025). Tracking breast cancer progression using Methylscape. bioRxiv*. https://doi.org/10.1101/2025.07.09.664004

*The manuscript is currently undergoing advanced review. Once published, we will update the citation with the proper reference.
