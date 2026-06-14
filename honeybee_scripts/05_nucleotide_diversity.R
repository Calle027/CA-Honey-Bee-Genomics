###############################
## 05_nucleotide_diversity.R
## Nucleotide diversity summary and figure.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("05_nucleotide_diversity.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

pi_files <- list.files(
  path = base_path,
  pattern = "maf05.*\\.windowed\\.pi$",
  full.names = TRUE
)

if (length(pi_files) == 0) {
  warning("No .windowed.pi files found. Skipping nucleotide diversity section.")
} else {
  pi_df <- pi_files %>%
    set_names(basename(.)) %>%
    purrr::map_dfr(~ readr::read_tsv(.x, show_col_types = FALSE), .id = "file") %>%
    mutate(
      pop_raw = stringr::str_match(file, "^(?:pi_)?([^_]+)_maf05\\.windowed\\.pi$")[, 2],
      pop = factor(pop_raw, levels = lvl)
    ) %>%
    relocate(pop, file, .before = 1)
  
  write_csv(pi_df, file.path(out_dir, "pi_all_windows_raw.csv"))
  
  shared_windows <- pi_df %>%
    distinct(pop, CHROM, BIN_START, BIN_END) %>%
    count(CHROM, BIN_START, BIN_END, name = "n_pops") %>%
    filter(n_pops == length(lvl)) %>%
    select(CHROM, BIN_START, BIN_END)
  
  pi_shared <- pi_df %>%
    semi_join(shared_windows, by = c("CHROM", "BIN_START", "BIN_END")) %>%
    mutate(pop = factor(pop, levels = lvl))
  
  write_csv(pi_shared, file.path(out_dir, "pi_shared_windows.csv"))
  
  pi_plot <- ggplot(pi_shared, aes(x = pop, y = PI, fill = pop)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, color = "black", linewidth = 0.4) +
    scale_y_continuous(limits = c(0, 0.006), expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = pop_colors, breaks = lvl, drop = FALSE) +
    labs(
      x = "Population",
      y = expression("Nucleotide diversity (" * pi * ")")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  save_plot("pi_plot.png", width = 10, height = 6)
  print(pi_plot)
  dev.off()
  
  pi_summary <- pi_shared %>%
    group_by(pop) %>%
    summarise(
      mean_pi = mean(PI, na.rm = TRUE),
      sd_pi   = sd(PI, na.rm = TRUE),
      n       = dplyr::n(),
      .groups = "drop"
    )
  write_csv(pi_summary, file.path(out_dir, "pi_summary.csv"))
}

showfile.gds(closeall = TRUE)
