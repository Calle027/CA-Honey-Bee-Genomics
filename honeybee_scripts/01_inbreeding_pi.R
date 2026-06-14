###############################
## 01_inbreeding_pi.R
## Legacy wrapper. Prefer running 01_inbreeding.R and 05_nucleotide_diversity.R.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("01_inbreeding_pi.R")
script_dir <- dirname(script_file)

source(file.path(script_dir, "01_inbreeding.R"))

source(file.path(script_dir, "05_nucleotide_diversity.R"))
