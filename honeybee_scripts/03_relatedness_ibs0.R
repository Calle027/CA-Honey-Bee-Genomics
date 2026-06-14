###############################
## 03_relatedness_ibs0.R
## KING relatedness + IBS0 heatmaps
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("03_relatedness_ibs0.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(SNPRelate)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  library(readr)
})

# -----------------------------
# Shared KING objects
# -----------------------------
ibd_mom <- load_ibd_mom()

ids <- trimws(as.character(ibd_mom$sample.id))
ibs0_mat <- as.matrix(ibd_mom$IBS0)
kin_mat  <- as.matrix(ibd_mom$kinship)
rownames(ibs0_mat) <- colnames(ibs0_mat) <- ids
rownames(kin_mat)  <- colnames(kin_mat)  <- ids

write_csv(as.data.frame(kin_mat) %>% rownames_to_column("sample.id_raw"),
          file.path(out_dir, "king_kinship_matrix.csv"))
write_csv(as.data.frame(ibs0_mat) %>% rownames_to_column("sample.id_raw"),
          file.path(out_dir, "king_ibs0_matrix.csv"))

# -----------------------------
# Global IBS0 heatmap
# -----------------------------
ord_all <- analysis_samples %>%
  filter(sample.id_raw %in% ids) %>%
  arrange(pop) %>%
  pull(sample.id_raw)

ibs0_global <- ibs0_mat[ord_all, ord_all, drop = FALSE]

ann_df <- analysis_samples %>%
  filter(sample.id_raw %in% ord_all) %>%
  dplyr::select(sample.id_raw, pop) %>%
  distinct() %>%
  arrange(match(sample.id_raw, ord_all)) %>%
  tibble::column_to_rownames("sample.id_raw")

ann_pop <- factor(ann_df$pop, levels = lvl)

tmp_global <- ibs0_global
diag(tmp_global) <- NA_real_
rng_global <- range(tmp_global, na.rm = TRUE)
if (rng_global[1] == rng_global[2]) rng_global <- rng_global + c(-1e-6, 1e-6)

col_fun_global <- colorRamp2(
  seq(rng_global[1], rng_global[2], length.out = 200),
  colorRampPalette(brewer.pal(11, "RdYlBu"))(200)
)

ha_global <- HeatmapAnnotation(
  Population = ann_pop,
  col = list(Population = pop_colors),
  show_annotation_name = FALSE,
  show_legend = c(Population = FALSE)
)

ra_global <- rowAnnotation(
  Population = ann_pop,
  col = list(Population = pop_colors),
  show_annotation_name = FALSE,
  show_legend = c(Population = TRUE)
)

ht_global <- Heatmap(
  ibs0_global,
  name = "IBS0",
  col = col_fun_global,
  top_annotation = ha_global,
  left_annotation = ra_global,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)

save_plot("IBS0_global_heatmap.png", width = 8.5, height = 6)
draw(
  ht_global,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

# -----------------------------
# Population-specific relatedness panels
# -----------------------------
K <- as.matrix(ibd_mom$kinship)
rownames(K) <- colnames(K) <- ids
diag(K) <- NA_real_
R <- 2 * K
diag(R) <- NA_real_

grp <- analysis_samples$pop[match(ids, analysis_samples$sample.id_raw)] |> as.character()

bin_levels <- c("0.044-0.088", "0.088-0.177", "0.177-0.354", "0.354-0.75", ">=0.75")
bin_cols <- c(
  "0.044-0.088" = "#9ecae1",
  "0.088-0.177" = "#6baed6",
  "0.177-0.354" = "#4292c6",
  "0.354-0.75"  = "#2171b5",
  ">=0.75"      = "#08306B"
)

bin_rel <- function(r) {
  if (!is.finite(r)) return(NA_character_)
  if (r >= 0.75)  return(">=0.75")
  if (r >= 0.354) return("0.354-0.75")
  if (r >= 0.177) return("0.177-0.354")
  if (r >= 0.0884) return("0.088-0.177")
  if (r >= 0.0442) return("0.044-0.088")
  NA_character_
}

make_rel_panel <- function(pop) {
  ids_pop <- ids[grp == pop]
  n <- length(ids_pop)
  if (n < 2) {
    return(ggplot() + theme_void() + labs(title = paste0(pop, " (n<2)")))
  }
  
  M <- R[ids_pop, ids_pop, drop = FALSE]
  diag(M) <- NA_real_
  
  df <- expand.grid(i = seq_len(n), j = seq_len(n)) |>
    mutate(
      r   = as.numeric(M[cbind(i, j)]),
      bin = vapply(r, bin_rel, character(1)),
      bin = factor(bin, levels = bin_levels),
      x   = i,
      y   = n - j + 1
    )
  
  diag_df <- data.frame(x = seq_len(n), y = n:1)
  df_bins <- filter(df, !is.na(bin) & !(x == (n - y + 1)))
  
  ggplot(df, aes(x = x, y = y)) +
    geom_tile(fill = "white", color = "grey92", linewidth = 0.20) +
    geom_tile(data = df_bins, aes(fill = bin), color = NA) +
    geom_tile(data = diag_df, inherit.aes = FALSE, aes(x = x, y = y),
              fill = "#041E42", color = NA) +
    coord_fixed(expand = FALSE) +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_fill_manual(values = bin_cols, limits = bin_levels, drop = FALSE, na.value = "white") +
    labs(title = paste0(pop, " (n=", n, ")")) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6)
    )
}

make_bin_legend_panel <- function(bin_levels, bin_cols, title = "Relatedness (r)") {
  df <- data.frame(bin = rev(bin_levels), stringsAsFactors = FALSE)
  df$y <- seq_len(nrow(df)) * 0.75
  ggplot(df, aes(y = y)) +
    geom_tile(aes(x = 0, fill = bin), width = 0.3, height = 0.3,
              color = "white", linewidth = 0.6) +
    geom_text(aes(x = 0.45, label = bin), hjust = 0, size = 5.5) +
    scale_fill_manual(values = bin_cols, limits = bin_levels, drop = FALSE) +
    guides(fill = "none") +
    coord_fixed(ratio = 1,
                xlim = c(-0.40, 1.55),
                ylim = c(0.45, max(df$y) + 0.3),
                expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 18) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 20))
}

# Hybrids
hybrid_pops <- c("NorCal", "SoCal", "Meso")
hybrid_plots <- lapply(hybrid_pops, make_rel_panel)
Hybrids_grid <- wrap_plots(hybrid_plots, ncol = 3, byrow = TRUE)
legend_panel <- make_bin_legend_panel(bin_levels, bin_cols)

ggsave(file.path(fig_dir, "relatedness_Hybrids.png"), Hybrids_grid,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "relatedness_main_only_legend.png"), legend_panel,
       width = 4, height = 4.3, dpi = 300)

# Ancestors
anc_A <- "A"
anc_others <- c("C", "M", "O", "L", "U", "Y")
A_plot <- make_rel_panel(anc_A)
other_plots <- lapply(anc_others, make_rel_panel)
others_grid <- wrap_plots(other_plots, ncol = 2, byrow = TRUE)
Ancestors_split <- A_plot | others_grid + plot_layout(widths = c(1.15, 1))

ggsave(file.path(fig_dir, "relatedness_Ancestors_Aleft.png"), Ancestors_split,
       width = 10, height = 8, dpi = 300)

showfile.gds(closeall = TRUE)
