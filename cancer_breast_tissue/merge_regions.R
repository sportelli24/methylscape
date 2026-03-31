library(readr)
library(dplyr)
library(tidyr)

cov1_file <- "/cancer_breast_tissue/Maxi_1_1.regions.cov"
cov2_file <- "/cancer_breast_tissue/Maxi_1_2.regions.cov"
out_tsv   <- "/cancer_breast_tissue/output/maxi_1112_regions_merged_cpg.tsv"

colspec <- cols(
  chr       = col_character(),
  start     = col_integer(),
  end       = col_integer(),
  perc_meth = col_double(),
  meth      = col_integer(),
  unmeth    = col_integer()
)

rep1 <- read_tsv(cov1_file, col_names = c("chr","start","end","perc_meth","meth","unmeth"),
                 col_types = colspec) |>
  transmute(chr, pos = start, meth1 = meth, unmeth1 = unmeth)

rep2 <- read_tsv(cov2_file, col_names = c("chr","start","end","perc_meth","meth","unmeth"),
                 col_types = colspec) |>
  transmute(chr, pos = start, meth2 = meth, unmeth2 = unmeth)

tumour_merged <- full_join(rep1, rep2, by = c("chr", "pos")) |>
  mutate(
    meth1   = replace_na(meth1,   0L),
    unmeth1 = replace_na(unmeth1, 0L),
    meth2   = replace_na(meth2,   0L),
    unmeth2 = replace_na(unmeth2, 0L),
    meth_tot   = meth1 + meth2,
    unmeth_tot = unmeth1 + unmeth2,
    cov        = meth_tot + unmeth_tot,
    beta       = ifelse(cov > 0, meth_tot / cov, NA_real_)
  ) |>
  select(chr, pos, beta, cov)

write_tsv(tumour_merged, out_tsv)
