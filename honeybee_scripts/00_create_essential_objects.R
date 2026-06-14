###############################
## 00_create_essential_objects.R
## Create shared analysis objects used by downstream figure scripts.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("00_create_essential_objects.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(SNPRelate)
  library(GENESIS)
  library(dplyr)
  library(readr)
  library(tibble)
})

snpset_list <- snpgdsLDpruning(
  genofile,
  method = "corr",
  start.pos = "first",
  sample.id = analysis_ids_raw,
  maf = 0.10,
  ld.threshold = 0.10,
  slide.max.bp = 5000,
  autosome.only = FALSE,
  verbose = FALSE
)

snpset <- unlist(snpset_list, use.names = FALSE)
saveRDS(snpset, snpset_rds_path)
write.table(
  snpset,
  snpset_txt_path,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

snp_annot <- tibble(
  snp.id = read.gdsn(index.gdsn(genofile, "snp.id")),
  chrom  = read.gdsn(index.gdsn(genofile, "snp.chromosome")),
  pos    = read.gdsn(index.gdsn(genofile, "snp.position")),
  id     = read.gdsn(index.gdsn(genofile, "snp.rs.id"))
)

pruned_snps <- snp_annot %>%
  filter(snp.id %in% snpset) %>%
  arrange(chrom, pos)

write_tsv(pruned_snps, file.path(out_dir, "pruned_snps_annotation.tsv"))
write.table(
  pruned_snps[, c("chrom", "pos")],
  file = file.path(out_dir, "pruned_snps_chrom_pos.tsv"),
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

ibd_mom <- snpgdsIBDKING(
  genofile,
  snp.id = snpset,
  sample.id = analysis_ids_raw,
  autosome.only = FALSE,
  verbose = FALSE
)

ids <- trimws(as.character(ibd_mom$sample.id))
kin_mat <- as.matrix(ibd_mom$kinship)
ibs0_mat <- as.matrix(ibd_mom$IBS0)
rownames(kin_mat) <- colnames(kin_mat) <- ids
rownames(ibs0_mat) <- colnames(ibs0_mat) <- ids
ibd_mom$kinship <- kin_mat
ibd_mom$IBS0 <- ibs0_mat

KINGmat <- kingToMatrix(ibd_mom)

saveRDS(ibd_mom, ibd_mom_rds_path)
saveRDS(KINGmat, kingmat_rds_path)

write_csv(as.data.frame(kin_mat) %>% rownames_to_column("sample.id_raw"),
          file.path(out_dir, "king_kinship_matrix.csv"))
write_csv(as.data.frame(ibs0_mat) %>% rownames_to_column("sample.id_raw"),
          file.path(out_dir, "king_ibs0_matrix.csv"))

summary_tbl <- tibble(
  object = c("analysis_samples", "snpset", "king_matrix"),
  n = c(length(analysis_ids_raw), length(snpset), nrow(KINGmat)),
  file = c(
    file.path(out_dir, "sample_map.csv"),
    snpset_rds_path,
    kingmat_rds_path
  )
)
write_csv(summary_tbl, file.path(out_dir, "essential_objects_summary.csv"))

showfile.gds(closeall = TRUE)
