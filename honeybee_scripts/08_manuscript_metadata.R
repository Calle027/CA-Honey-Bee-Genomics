###############################
## 08_manuscript_metadata.R
## Create manuscript sample metadata table.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("08_manuscript_metadata.R")
base_path <- resolve_project_root(script_file)
out_dir <- file.path(base_path, "analysis_outputs")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^_", "", x)
  x <- sub("[.]g[.]vcf[.]gz$", "", x)
  x
}

remove_ids <- c(
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

required_files <- c(
  file.path(base_path, "HB_fam.xlsx"),
  file.path(base_path, "genome_coverage.tsv"),
  file.path(out_dir, "inbreeding_coefficients.csv"),
  file.path(out_dir, "pcair_scores_with_mclust.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

raw_metadata <- read_excel(file.path(base_path, "HB_fam.xlsx"), sheet = "raw_data") %>%
  transmute(
    sample_id = clean_id(ID),
    species = if_else(subspecies == "A. cerana" | pop == "Outgroup", "A. cerana", "A. mellifera"),
    population = pop,
    subspecies = subspecies,
    city_region = `city/region`,
    country = country,
    lat = lat,
    long = long,
    western_hemisphere = Western_Hemisphere
  ) %>%
  filter(!sample_id %in% remove_ids)

coverage <- read_tsv(
  file.path(base_path, "genome_coverage.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    sample_id = clean_id(sample),
    coverage = as.numeric(avg_coverage)
  )

inbreeding <- read_csv(file.path(out_dir, "inbreeding_coefficients.csv"), show_col_types = FALSE) %>%
  transmute(
    sample_id = clean_id(sample.id_clean),
    inbreeding = inbreeding
  )

pcair_mclust <- read_csv(file.path(out_dir, "pcair_scores_with_mclust.csv"), show_col_types = FALSE) %>%
  transmute(
    sample_id = clean_id(sample.id_clean),
    relatedness = relatedness,
    PC1 = PC1,
    PC2 = PC2,
    cluster = cluster
  )

manuscript_metadata <- raw_metadata %>%
  left_join(coverage, by = "sample_id") %>%
  left_join(inbreeding, by = "sample_id") %>%
  left_join(pcair_mclust, by = "sample_id") %>%
  mutate(
    coverage = if_else(species == "A. cerana", NA_real_, coverage),
    city_region = if_else(city_region == "Nothern CA", "Northern CA", city_region)
  ) %>%
  select(
    sample_id,
    species,
    population,
    subspecies,
    city_region,
    country,
    lat,
    long,
    western_hemisphere,
    coverage,
    relatedness,
    inbreeding,
    PC1,
    PC2,
    cluster
  )

write_csv(
  manuscript_metadata,
  file.path(base_path, "manuscript_sample_metadata.csv"),
  na = ""
)
