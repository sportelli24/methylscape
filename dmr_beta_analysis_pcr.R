############################################################
## Region-level methylation summary for 4 genes of interest
## - Regions: promoters + gene bodies (from BED file)
## - Samples: 1 tumour (collapsed reps) + 3 normals
## - Outputs:
##   * region_summary.csv          (per sample and region)
##   * summary_vs_normals.csv     
##
## Two β summarisation mode:
##   - simple_beta  : plain mean across CpGs
############################################################

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)

# Regions BED (promoter + gene body for ALX1, MACF1, IGSF3, CLPTM1L)
regions_bed_path <- "/cancer_breast_tissue/regions_of_interest.bed"   
regions_bed_amplicon <- "/amplicons_pcr_hg19.bed"

# Tumour CpG table (collapsed replicates, already filtered to regions)
tumour_tsv  <- "/cancer_breast_tissue/output/maxi_1112_regions_merged_cpg.tsv"         

# Normal CpG tables (BigWig-derived, filtered to regions)
normal1_tsv <- "/healthy_breast_tissue/outputs/normal_1_GSM5652347_regions.tsv"          
normal2_tsv <- "/healthy_breast_tissue/outputs/normal_2_GSM5652348_regions.tsv"          
normal3_tsv <- "/healthy_breast_tissue/outputs/normal_3_GSM5652349_regions.tsv"          

# Optional: set prefix for output
out_prefix <- "dmr_beta_4genes"

# Sample metadata

sample_files <- tribble(
  ~sample_id,   ~group,    ~path,
  "tumour",     "tumour",  tumour_tsv,
  "normal_1",   "normal",  normal1_tsv,
  "normal_2",   "normal",  normal2_tsv,
  "normal_3",   "normal",  normal3_tsv
)


# Load regions from BED and parse metadata

message("Importing regions from BED: ", regions_bed_path)

regions_gr <- rtracklayer::import(regions_bed_path)
amplicon_gr <- rtracklayer::import(regions_bed_amplicon)

# promoter and gene bodies metadata
region_names <- mcols(regions_gr)$name

gene_symbol  <- sub("_(promoter|gene_body)$", "", region_names)
region_type  <- sub("^.*_(promoter|gene_body)$", "\\1", region_names)

mcols(regions_gr)$gene_symbol  <- gene_symbol
mcols(regions_gr)$gene_id      <- gene_symbol  # no separate ID, use symbol
mcols(regions_gr)$region_type  <- region_type  # "promoter" or "gene_body"
mcols(regions_gr)$region_label <- region_names # full name

# amplicon metadata
amp_names <- mcols(amplicon_gr)$name
amp_gene  <- sub("_.*$", "", amp_names)     

mcols(amplicon_gr)$gene_symbol  <- amp_gene
mcols(amplicon_gr)$gene_id      <- amp_gene
mcols(amplicon_gr)$region_type  <- "amplicon"
mcols(amplicon_gr)$region_label <- amp_names

# combine all regions
regions_gr_all <- c(regions_gr, amplicon_gr)

message("Regions loaded:")
print(data.frame(
  chr   = as.character(seqnames(regions_gr_all)),
  start = start(regions_gr_all),
  end   = end(regions_gr_all),
  name  = mcols(regions_gr_all)$region_label,
  type  = mcols(regions_gr_all)$region_type
))

# Read CpG tables and convert to GRanges

# Each TSV should have: chr, pos, beta, cov
read_cpg_table_as_gr <- function(path, sample_id) {
  message("Reading CpG table: ", path, " (sample: ", sample_id, ")")
  
  dat <- read_tsv(path, col_types = cols(
    chr  = col_character(),
    pos  = col_integer(),
    beta = col_double(),
    cov  = col_double()
  ))
  
  gr <- GRanges(
    seqnames = dat$chr,
    ranges   = IRanges(start = dat$pos, end = dat$pos),
    strand   = "*"
  )
  
  mcols(gr)$beta      <- dat$beta
  mcols(gr)$cov       <- dat$cov
  mcols(gr)$sample_id <- sample_id
  
  gr
}

# Build list of CpG GRanges by sample
cpg_gr_list <- lapply(seq_len(nrow(sample_files)), function(i) {
  read_cpg_table_as_gr(sample_files$path[i], sample_files$sample_id[i])
})
names(cpg_gr_list) <- sample_files$sample_id

# Summarise methylation per region for each sample

summarise_regions_for_sample <- function(cpg_gr, regions_gr, sample_id) {
  hits <- findOverlaps(regions_gr, cpg_gr)
  
  if (length(hits) == 0L) {
    message("No CpGs overlapping regions for sample ", sample_id)
    return(tibble())
  }
  
  region_idx <- queryHits(hits)
  cpg_idx    <- subjectHits(hits)
  
  region_meta <- as.data.frame(mcols(regions_gr))[region_idx, , drop = FALSE]
  cpg_meta    <- as.data.frame(mcols(cpg_gr))[cpg_idx, c("beta", "cov"), drop = FALSE]
  
  tmp <- tibble(
    sample_id    = sample_id,
    gene_symbol  = region_meta$gene_symbol,
    gene_id      = region_meta$gene_id,
    region_type  = region_meta$region_type,
    region_label = region_meta$region_label,
    beta         = cpg_meta$beta,
    cov          = cpg_meta$cov
  )
  
  tmp %>%
    group_by(sample_id, gene_symbol, gene_id, region_type, region_label) %>%
    summarise(
      n_cpg         = sum(!is.na(beta) & cov > 0),
      total_cov     = sum(cov, na.rm = TRUE),
      weighted_beta = ifelse(
        total_cov > 0,
        sum(beta * cov, na.rm = TRUE) / total_cov,
        NA_real_
      ),
      simple_beta   = ifelse(
        n_cpg > 0,
        mean(beta[cov > 0], na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )
}

message("Summarising regions for each sample...")

region_summary_list <- lapply(names(cpg_gr_list), function(sid) {
  summarise_regions_for_sample(cpg_gr_list[[sid]], regions_gr_all, sid)
})

region_summary <- bind_rows(region_summary_list) %>%
  left_join(sample_files %>% select(sample_id, group), by = "sample_id")

# Write out full per-sample region summary
out_region_summary <- paste0(out_prefix, "_region_summary.csv")
write_csv(region_summary, out_region_summary)
message("Wrote region summary to: ", out_region_summary)

# Compute simple beta from tumour vs average of healthy samples

# Split tumour vs normals
tumour_region <- region_summary %>%
  filter(group == "tumour") %>%
  select(gene_symbol, gene_id, region_type, region_label,
         tumour_weighted_beta = weighted_beta,
         tumour_simple_beta   = simple_beta)

normal_region <- region_summary %>%
  filter(group == "normal")

# Tumour vs mean of normals

summary_vs_normals <- normal_region %>%
  group_by(gene_symbol, gene_id, region_type, region_label) %>%
  summarise(
    normal_weighted_mean_beta = mean(weighted_beta, na.rm = TRUE),
    normal_weighted_sd_beta   = sd(weighted_beta,  na.rm = TRUE),
    normal_simple_mean_beta   = mean(simple_beta,   na.rm = TRUE),
    normal_simple_sd_beta     = sd(simple_beta,     na.rm = TRUE),
    n_normals                 = n(),
    .groups = "drop"
  ) %>%
  left_join(tumour_region,
            by = c("gene_symbol", "gene_id", "region_type", "region_label")) %>%
  mutate(
    delta_beta_mean_weighted = tumour_weighted_beta - normal_weighted_mean_beta,
    delta_beta_mean_simple   = tumour_simple_beta   - normal_simple_mean_beta
  ) %>%
  arrange(gene_symbol, region_type, region_label)

out_summary_normals <- paste0(out_prefix, "_summary_vs_normals.csv")
write_csv(summary_vs_normals, out_summary_normals)
message("Wrote tumour vs mean-normal summary to: ", out_summary_normals)

summary_amplicon_results <- summary_vs_normals %>%
  filter(region_type == "amplicon") %>%
  arrange(gene_symbol, region_label)

summary_amplicon_results

out_amplicon_results <- paste0(out_prefix,"_summary_amplicon_results.csv")
write_csv(summary_amplicon_results, out_amplicon_results)
message("Wrote tumour vs mean-normal summary to: ", out_amplicon_results)


message("All done.")

