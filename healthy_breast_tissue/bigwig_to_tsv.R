############################################################
## BigWig → CpG table converter (hg19 controls)
##
## For each BigWig:
##   - imports as GRanges
##   - extracts chr, pos, beta
##   - infers scale (0–1 vs 0–100) and rescales if needed
##   - sets cov = 1 (no coverage info available)
##   - write TSV: chr, pos, beta, cov
############################################################

library(rtracklayer)
library(dplyr)
library(readr)

bigwig_paths <- c(
  "/healthy_breast_tissue/GSM5652347_Breast-Luminal-Epithelial-Z000000V2.bigwig",  
  "/healthy_breast_tissue/GSM5652348_Breast-Luminal-Epithelial-Z000000VJ.bigwig",  
  "/healthy_breast_tissue/GSM5652349_Breast-Luminal-Epithelial-Z000000VN.bigwig"   
)

sample_ids <- c("normal_1", "normal_2", "normal_3") 

out_tsv_paths <- c(
  "/healthy_breast_tissue/outputs/normal_1_GSM5652347_from_bigwig.tsv",
  "/healthy_breast_tissue/outputs/normal_2_GSM5652348_from_bigwig.tsv",
  "/healthy_breast_tissue/outputs/normal_3_GSM5652349_from_bigwig.tsv"
)


# Convert one BigWig to CpG/position-level TSV
bigwig_to_cpg_tsv <- function(bw_path, out_tsv, sample_id = NA) {
  message("Importing BigWig: ", bw_path)
  
  # Import BigWig as GRanges
  bw <- rtracklayer::import(bw_path)
  
  df <- as.data.frame(bw)
  
  if (!"score" %in% names(df)) {
    stop("BigWig ", bw_path, " does not have a 'score' column. Names are: ",
         paste(names(df), collapse = ", "))
  }
  
  # Check score range to guess scale
  score_range <- range(df$score, na.rm = TRUE)
  message("  score range: [", paste(score_range, collapse = ", "), "]")
  
  if (score_range[2] <= 1.5) {
    # already scalced from 0–1 
    beta_vals <- df$score
  } else if (score_range[2] <= 100.5) {
    # scaled as a percentage
    beta_vals <- df$score / 100
    message("  Detected 0–100 scale; dividing by 100 to get 0–1 beta.")
  } else {
    warning("Scores go above 100; please double-check BigWig scale. ",
            "Proceeding by dividing by 100.")
    beta_vals <- df$score / 100
  }
  
  # Build CpG-like table
  cpg_tab <- tibble(
    chr  = as.character(df$seqnames),
    pos  = as.integer(df$start),   # representative coordinate
    beta = beta_vals,
    cov  = 1L                      # no coverage info available
  )
  
  # Optional filter for standard chromosomes
  cpg_tab <- cpg_tab %>%
    filter(grepl("^chr[0-9XYM]+$", chr))
  
  message("  Writing TSV: ", out_tsv)
  readr::write_tsv(cpg_tab, out_tsv)
  
  invisible(cpg_tab)
}


# Loop over all 3 controls
stopifnot(length(bigwig_paths) == length(out_tsv_paths),
          length(bigwig_paths) == length(sample_ids))

for (i in seq_along(bigwig_paths)) {
  bw  <- bigwig_paths[i]
  tsv <- out_tsv_paths[i]
  sid <- sample_ids[i]
  
  message("Processing control ", i, " (", sid, ")")
  bigwig_to_cpg_tsv(bw_path = bw, out_tsv = tsv, sample_id = sid)
}

message("All BigWigs converted to TSVs:")
print(out_tsv_paths)
