###############################
## 00_config.R
## Basal configuration + shared objects for CA honey bee genomics analyses
###############################

suppressPackageStartupMessages({
  library(gdsfmt)
  library(SeqArray)
  library(SNPRelate)
  library(dplyr)
  library(readxl)
  library(readr)
  library(tibble)
  library(stringr)
  library(ragg)
})

# -----------------------------
# Paths
# -----------------------------
parent_ofile <- sys.frame(1)$ofile
if (is.null(parent_ofile)) {
  parent_ofile <- file.path("honeybee_scripts", "00_config.R")
}
helper_file <- file.path(dirname(parent_ofile), "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)

config_file <- resolve_script_path("00_config.R")
default_base_path <- resolve_project_root(config_file)

base_path <- Sys.getenv("CA_HONEYBEE_ROOT", unset = default_base_path)
base_path <- normalizePath(base_path, mustWork = TRUE)

gds_path    <- file.path(base_path, "final_cleaned.gds")
excel_path  <- file.path(base_path, "HB_fam.xlsx")
out_dir     <- file.path(base_path, "analysis_outputs")
fig_dir     <- file.path(base_path, "figures")
script_dir  <- file.path(base_path, "honeybee_scripts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

snpset_rds_path   <- file.path(out_dir, "snpset_ld_pruned_ids.rds")
snpset_txt_path   <- file.path(out_dir, "snpset_ld_pruned_ids.txt")
ibd_mom_rds_path  <- file.path(out_dir, "king_ibd_mom.rds")
kingmat_rds_path  <- file.path(out_dir, "king_matrix.rds")

# -----------------------------
# Helper functions
# -----------------------------
clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\.g\\.vcf\\.gz$", "", x)
  x
}

save_plot <- function(filename, width = 12, height = 8, units = "in", res = 300) {
  out_path <- file.path(fig_dir, filename)
  ragg::agg_png(filename = out_path, width = width, height = height, units = units, res = res)
}

require_essential_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      sprintf(
        "%s was not found at %s. In RStudio, open and Source honeybee_scripts/00_create_essential_objects.R first.",
        label,
        path
      ),
      call. = FALSE
    )
  }
}

load_snpset <- function() {
  require_essential_file(snpset_rds_path, "LD-pruned SNP set")
  readRDS(snpset_rds_path)
}

load_ibd_mom <- function() {
  require_essential_file(ibd_mom_rds_path, "KING IBD object")
  readRDS(ibd_mom_rds_path)
}

load_king_matrix <- function() {
  require_essential_file(kingmat_rds_path, "KING kinship matrix")
  readRDS(kingmat_rds_path)
}

assert_required_cols <- function(df, cols, df_name = deparse(substitute(df))) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("%s is missing required columns: %s", df_name, paste(missing_cols, collapse = ", ")))
  }
}

# -----------------------------
# Population definitions
# -----------------------------
lvl <- c("NorCal","SoCal","Meso","A","C","M","O","L","U","Y")

pop_colors <- c(
  "NorCal" = "#7fbc41ff",
  "SoCal"  = "#ab2f5e",
  "Meso"   = "#A020F0",
  "A"      = "#FF1406",
  "C"      = "#EBCF00",
  "M"      = "grey45",
  "O"      = "#12C9D2",
  "L"      = "#B59BD9",
  "U"      = "#B6D7A8",
  "Y"      = "#CDAA74"
)

# -----------------------------
# Exclusions
# Store exclusions as CLEAN sample IDs only.
# -----------------------------
remove_ids_clean <- c(
  "CA0124","CA0121","CA0112","CA0126","CA0108","CA0119","CA0102","CA0128",
  "CA0201","CA0230","CA0205","CA0227","CA0217","CA0226","CA0225","CA0222",
  "CA0220","CA0330","CA0303","CA0306","CA0308","CA0329","CA0317","CA0325",
  "CA0321","CA0403","CA0401","CA0411","CA0412","CA0416","CA0413","CA0427",
  "CA0423","CA0502","CA0505","CA0527","CA0509","CA0522","CA0524","CA0519",
  "CA0512","CA0529","CA0825","CA0816","CA0813","CA0822","CA0805","CA0819",
  "CA0810","CA0801","CA0914","CA0912","CA0910","CA0918","CA0919","CA0906",
  "CA0904","CA0902","CA1010","CA1013","CA1015","CA1012","CA1007","CA1020",
  "CA1005","CA1004","CA1002","CA1115","CA1118","CA1111","CA1120","CA1109",
  "CA1122","CA1106","CA1103","CA1215","CA1218","CA1216","CA1213","CA1211",
  "CA1209","CA1207","CA1201","CA1302","CA1305","CA1307","CA1317","CA1315",
  "CA1314","CA1309","CA1313","CA1319","CA1410","CA1414","CA1408","CA1406",
  "CA1416","CA1402","CA1417","CA1419",
  "Ayem_SN1","Ayem_SN2","Ayem_SN4",
  "Ayem_SI2","Ayem_SI3","Ayem_SI4",
  "Alamark12_S26", "Alamark31_S29"
)

# -----------------------------
# Metadata import
# -----------------------------
metadata <- read_excel(excel_path, sheet = "raw_data") %>%
  mutate(
    ID       = as.character(ID),
    ID       = trimws(ID),
    ID_clean = clean_id(ID),
    ID_clean = as.character(ID_clean),  
    pop_raw  = trimws(as.character(pop)),
    pop      = pop_raw,
    pop      = factor(pop, levels = lvl)
  )

assert_required_cols(metadata, c("ID", "ID_clean", "pop"), "metadata")

# Check duplicated cleaned IDs in metadata
meta_dups <- metadata %>%
  dplyr::count(ID_clean, sort = TRUE) %>%
  dplyr::filter(n > 1)

if (nrow(meta_dups) > 0) {
  write_csv(meta_dups, file.path(out_dir, "metadata_duplicate_clean_ids.csv"))
  warning("Duplicated cleaned IDs found in metadata. See metadata_duplicate_clean_ids.csv")
}

# -----------------------------
# GDS sample map creation
# This is the canonical sample map used downstream.
# -----------------------------
showfile.gds(closeall = TRUE)
genofile <- snpgdsOpen(gds_path)

sample_ids_raw <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample_map <- tibble(
  sample.id_raw   = as.character(sample_ids_raw),
  sample.id_clean = clean_id(sample_ids_raw)
) %>%
  left_join(
    metadata %>% dplyr::select(ID_clean, pop) %>% distinct(),
    by = c("sample.id_clean" = "ID_clean")
  ) %>%
  mutate(
    excluded = sample.id_clean %in% remove_ids_clean,
    in_metadata = !is.na(pop),
    pop = factor(pop, levels = lvl)
  )

# Diagnostics
write_csv(sample_map, file.path(out_dir, "sample_map.csv"))

metadata_not_in_gds <- metadata %>%
  filter(!ID_clean %in% sample_map$sample.id_clean)
write_csv(metadata_not_in_gds, file.path(out_dir, "metadata_not_in_gds.csv"))

gds_not_in_metadata <- sample_map %>%
  filter(!in_metadata)
write_csv(gds_not_in_metadata, file.path(out_dir, "gds_not_in_metadata.csv"))

# Canonical analysis sample vectors
analysis_samples <- sample_map %>%
  filter(in_metadata, !excluded) %>%
  mutate(pop = droplevels(pop))

analysis_ids_raw   <- analysis_samples$sample.id_raw
analysis_ids_clean <- analysis_samples$sample.id_clean

# Final counts after filtering
final_counts <- analysis_samples %>%
  dplyr::count(pop) %>%
  arrange(factor(pop, levels = lvl))
write_csv(final_counts, file.path(out_dir, "final_population_counts.csv"))

# Leave genofile open for downstream scripts that source this file.
