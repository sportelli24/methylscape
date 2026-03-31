########################################################################
## Builds BED file of promoters and gene bodies for 4 genes of interest
## Genome: hg19, annotation: gencode v19 GTF
##
## Output: regions_of_interest.bed
##   chr  start  end  name                    score  strand
##   chr1 12345  13000  ALX1_promoter        0      +
##   chr1 13001  20000  ALX1_gene_body       0      +
##   ...
########################################################################

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# Path to gencode v19 GTF (hg19)
gtf_path <- "/genomes/gencode.v19.annotation.gtf.gz"   

# Genes of interest
gene_symbols <- c("ALX1", "MACF1", "IGSF3", "CLPTM1L")

# Promoter window around Transcription Start Site
promoter_upstream   <- 1500  
promoter_downstream <- 500   

# Output BED file
out_bed <- "/cancer_breast_tissue/regions_of_interest.bed"

# Load GTF and subset genes

message("Importing GTF...")
gtf <- rtracklayer::import(gtf_path)

# Filter for gene features only
genes_gtf <- gtf[gtf$type == "gene"]

if (!"gene_name" %in% colnames(mcols(genes_gtf))) {
  stop("GTF does not contain a 'gene_name' column. Check mcols(genes_gtf).")
}

# Subset to genes of interest
genes_of_interest <- genes_gtf[mcols(genes_gtf)$gene_name %in% gene_symbols]

if (length(genes_of_interest) == 0L) {
  stop("None of the requested gene symbols were found in gene_name of the GTF.")
}

# Clean metadata
mcols(genes_of_interest)$gene_symbol <- mcols(genes_of_interest)$gene_name

message("Found genes:")
print(unique(mcols(genes_of_interest)$gene_symbol))

# Identify promoter regions

promoters_gr <- promoters(
  genes_of_interest,
  upstream   = promoter_upstream,
  downstream = promoter_downstream
)

mcols(promoters_gr)$region_type  <- "promoter"
mcols(promoters_gr)$region_name  <- paste0(mcols(genes_of_interest)$gene_symbol,
                                           "_promoter")

# Identify gene body regions, i.e. gene span - (promoter window)

# Split by gene_symbol
genes_by_symbol     <- split(genes_of_interest, mcols(genes_of_interest)$gene_symbol)
promoters_by_symbol <- split(promoters_gr,     mcols(genes_of_interest)$gene_symbol)

bodies_list <- mapply(
  FUN = function(gene_gr, prom_gr) {
    res <- GenomicRanges::setdiff(gene_gr, prom_gr)
    mcols(res)$gene_symbol <- mcols(gene_gr)$gene_symbol   
    res
  },
  gene_gr = genes_by_symbol,
  prom_gr = promoters_by_symbol,
  SIMPLIFY = FALSE
)

bodies_gr <- do.call(c, unname(bodies_list))

mcols(bodies_gr)$region_type <- "gene_body"
mcols(bodies_gr)$region_name <- paste0(
  mcols(bodies_gr)$gene_symbol, "_gene_body"
)

# Combine and export bed file

all_regions <- c(promoters_gr, bodies_gr)

bed_df <- tibble(
  chr    = as.character(seqnames(all_regions)),
  start  = pmax(start(all_regions) - 1L, 0L),
  end    = end(all_regions),
  name   = mcols(all_regions)$region_name,
  score  = 0,
  strand = as.character(strand(all_regions))
) %>%
  arrange(chr, start, end)

message("Writing BED to: ", out_bed)

write.table(
  bed_df,
  file      = out_bed,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Done. Regions written:")
print(bed_df)
