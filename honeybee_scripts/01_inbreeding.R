###############################
## 01_inbreeding.R
## Inbreeding coefficient analysis and figure.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("01_inbreeding.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(SNPRelate)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(FSA)
  library(rcompanion)
  library(tibble)
})

snpset <- load_snpset()

samp.inb <- snpgdsIndInb(
  genofile,
  snp.id = snpset,
  sample.id = analysis_ids_raw,
  autosome.only = FALSE,
  method = "mle"
)

inbreed_df <- tibble(
  sample.id_raw = as.character(samp.inb$sample.id),
  inbreeding    = samp.inb$inbreeding
) %>%
  left_join(analysis_samples %>% dplyr::select(sample.id_raw, sample.id_clean, pop),
            by = "sample.id_raw") %>%
  mutate(pop = factor(pop, levels = lvl))

write_csv(inbreed_df, file.path(out_dir, "inbreeding_coefficients.csv"))

shapiro_df <- inbreed_df %>%
  group_by(pop) %>%
  summarise(
    n = dplyr::n(),
    p_value = ifelse(n >= 3, shapiro.test(inbreeding)$p.value, NA_real_),
    .groups = "drop"
  )
write_csv(shapiro_df, file.path(out_dir, "inbreeding_shapiro_by_population.csv"))

kw_inbreeding <- kruskal.test(inbreeding ~ pop, data = inbreed_df)
capture.output(kw_inbreeding,
               file = file.path(out_dir, "inbreeding_kruskal_test.txt"))

inbreed_df_test <- inbreed_df %>% filter(!is.na(pop)) %>% droplevels()
dunn_res <- dunnTest(
  x = inbreed_df_test$inbreeding,
  g = inbreed_df_test$pop,
  method = "bh"
)
write_csv(dunn_res$res, file.path(out_dir, "inbreeding_dunn_results.csv"))

dunn_df <- dunn_res$res %>%
  mutate(Comparison = gsub(" ", "", Comparison))

cld_letters <- cldList(P.adj ~ Comparison, data = dunn_df, threshold = 0.05)
letter_df <- cld_letters %>%
  rename(pop = Group, Letter = Letter) %>%
  mutate(pop = factor(pop, levels = levels(inbreed_df_test$pop)))
write_csv(letter_df, file.path(out_dir, "inbreeding_letters.csv"))

letter_y <- inbreed_df %>%
  summarise(y = max(inbreeding, na.rm = TRUE) + 0.04) %>%
  pull(y)

inbreeding_coef_plot <- ggplot(inbreed_df, aes(x = pop, y = inbreeding, fill = pop)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.2, alpha = 0.65, color = "black", size = 1) +
  geom_text(
    data = letter_df,
    aes(x = pop, y = letter_y, label = Letter),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_manual(values = pop_colors, drop = FALSE) +
  labs(x = "Population", y = "Inbreeding coefficient") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

save_plot("inbreeding_coef.png", width = 10, height = 6)
print(inbreeding_coef_plot)
dev.off()

inbreed_summary <- inbreed_df %>%
  group_by(pop) %>%
  summarise(
    mean_inbreed = mean(inbreeding, na.rm = TRUE),
    sd_inbreed   = sd(inbreeding, na.rm = TRUE),
    n            = dplyr::n(),
    .groups      = "drop"
  )
write_csv(inbreed_summary, file.path(out_dir, "inbreeding_summary.csv"))

showfile.gds(closeall = TRUE)
