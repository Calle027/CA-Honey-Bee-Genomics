###############################
## 06_phylogeny_admixture_mclust.R
## ML tree with K=7 ADMIXTURE bars and Mclust cluster tiles at tips.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("06_phylogeny_admixture_mclust.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(ape)
  library(treeio)
  library(ggtree)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(readxl)
  library(readr)
  library(tidyr)
  library(tibble)
  library(colorspace)
  library(grid)
})

# -----------------------------
# Inputs
# -----------------------------
K <- 7
tree_path <- file.path(base_path, "gtr_ascR4_nonsyn_constrained_varsites.contree")
q_path <- file.path(base_path, sprintf("admixture_input.%d.Q", K))
pcair_mclust_path <- file.path(out_dir, "pcair_scores_with_mclust.csv")
min_backbone_descendant_tips <- 20

if (!file.exists(pcair_mclust_path)) {
  stop(
    sprintf(
      "Mclust assignments not found at %s. In RStudio, open and Source honeybee_scripts/02_pcair_mclust.R first.",
      pcair_mclust_path
    ),
    call. = FALSE
  )
}

# -----------------------------
# Helpers
# -----------------------------
strip_iid_suffix <- function(x) {
  sub("\\.g\\.vcf\\.gz$", "", as.character(x))
}

assign_clusters_to_clades <- function(df, K, clades, fid_col = "FID") {
  cluster_cols <- paste0("Q", seq_len(K))
  
  means <- df %>%
    filter(.data[[fid_col]] %in% clades) %>%
    group_by(Clade = .data[[fid_col]]) %>%
    summarise(across(all_of(cluster_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  prefs <- means %>%
    pivot_longer(all_of(cluster_cols), names_to = "AdmixCluster", values_to = "MeanQ") %>%
    group_by(Clade) %>%
    arrange(desc(MeanQ), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  top_strength <- prefs %>%
    filter(rank == 1) %>%
    arrange(desc(MeanQ))
  
  assigned <- tibble(Clade = character(), AdmixCluster = character(), MeanQ = double())
  used_clusters <- character()
  
  for (i in seq_len(nrow(top_strength))) {
    clade_i <- top_strength$Clade[i]
    pick <- prefs %>%
      filter(Clade == clade_i, !AdmixCluster %in% used_clusters) %>%
      slice_head(n = 1)
    
    if (nrow(pick) == 0) {
      pick <- prefs %>%
        filter(Clade == clade_i) %>%
        slice_head(n = 1)
    }
    
    assigned <- bind_rows(assigned, pick %>% select(Clade, AdmixCluster, MeanQ))
    used_clusters <- c(used_clusters, pick$AdmixCluster)
  }
  
  setNames(assigned$Clade, assigned$AdmixCluster)
}

# -----------------------------
# ADMIXTURE K=7 table
# Sheet1 defines the ADMIXTURE .Q row order.
# -----------------------------
fam_order <- read_excel(excel_path, sheet = "Sheet1", col_names = FALSE)
stopifnot(ncol(fam_order) >= 6)
colnames(fam_order)[1:6] <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

fam_order <- fam_order %>%
  mutate(
    FID = as.character(FID),
    IID = as.character(IID),
    label = strip_iid_suffix(IID)
  )

Q <- read.table(q_path)
stopifnot(nrow(Q) == nrow(fam_order), ncol(Q) == K)
colnames(Q) <- paste0("Q", seq_len(K))

admix_df <- bind_cols(fam_order, Q)

clade_cols <- c(
  A = "#E31A1C",
  C = "#FFD400",
  M = "#4D4D4D",
  O = "#00BFC4",
  L = unname(pop_colors["L"]),
  U = unname(pop_colors["U"]),
  Y = unname(pop_colors["Y"])
)
cluster_to_clade <- assign_clusters_to_clades(
  admix_df,
  K = K,
  clades = names(clade_cols),
  fid_col = "FID"
)

admix_long <- admix_df %>%
  select(label, starts_with("Q")) %>%
  pivot_longer(starts_with("Q"), names_to = "AdmixCluster", values_to = "Prop") %>%
  mutate(
    CladeComp = unname(cluster_to_clade[AdmixCluster]),
    CladeComp = factor(CladeComp, levels = names(clade_cols))
  ) %>%
  group_by(label) %>%
  arrange(AdmixCluster, .by_group = TRUE) %>%
  mutate(
    x0_prop = cumsum(lag(Prop, default = 0)),
    x1_prop = cumsum(Prop)
  ) %>%
  ungroup()

write_csv(
  admix_long %>% select(label, AdmixCluster, CladeComp, Prop),
  file.path(out_dir, "tree_k7_admixture_tip_annotations.csv")
)

# -----------------------------
# Mclust assignments
# -----------------------------
mclust_df <- read_csv(pcair_mclust_path, show_col_types = FALSE) %>%
  mutate(
    label = as.character(sample.id_clean),
    cluster = as.character(cluster)
  ) %>%
  select(label, cluster)

cluster_levels <- sort(unique(mclust_df$cluster))
cluster_palette_source <- c(
  "#17154FFF", "#2F357CFF", "#6C5D9EFF", "#9D9CD5FF",
  "#B0799AFF", "#F6B3B0FF", "#E48171FF", "#BF3729FF",
  "#E69B00FF", "#F5BB50FF", "#ADA43BFF", "#355828FF"
)

# Chosen from cluster_palette_source for contrast among clusters while avoiding
# close matches to the population/tip colors.
cluster_palette <- c(
  "1" = "#4E79A7",
  "2" = "#F6B3B0FF",
  "3" = "#355828FF",
  "4" = "#9D9CD5FF",
  "5" = "#BF3729FF",
  "6" = "#F5BB50FF",
  "7" = "#2F357CFF",
  "8" = "#E48171FF",
  "9" = "#ADA43BFF"
)
cluster_cols <- cluster_palette[cluster_levels]

# -----------------------------
# Tree and tip metadata
# -----------------------------
tr <- read.tree(tree_path)

id_map <- metadata %>%
  transmute(
    label = as.character(ID),
    group = as.character(pop_raw)
  ) %>%
  distinct(label, group)

tip_df_all <- tibble(label = tr$tip.label) %>%
  left_join(id_map, by = "label") %>%
  mutate(group = if_else(is.na(group), "Unknown", group))

out_tips <- tip_df_all$label[tip_df_all$group == "Outgroup"]
stopifnot(length(out_tips) > 0)

tr_root <- root(tr, outgroup = out_tips, resolve.root = TRUE)

rootnode <- Ntip(tr_root) + 1
kids <- tr_root$edge[tr_root$edge[, 1] == rootnode, 2]
stopifnot(length(kids) == 2)

out_child <- MRCA(tr_root, out_tips)
stopifnot(!is.na(out_child))

out_desc <- treeio::offspring(as.treedata(tr_root), out_child)
if (is.data.frame(out_desc) && "node" %in% names(out_desc)) {
  out_nodes <- c(out_child, out_desc$node)
} else {
  out_nodes <- out_child
}

tr_root2 <- tr_root
longest_edge <- which.max(tr_root2$edge.length)
tr_root2$edge.length[longest_edge] <- 0.15
td_root <- as.treedata(tr_root2)

tip_pop_cols <- c(pop_colors, Outgroup = "black", Unknown = "grey80")
legend_order <- c("NorCal", "SoCal", "Meso", "A", "C", "M", "O", "L", "U", "Y", "Outgroup", "Unknown")
tip_levels <- legend_order[legend_order %in% unique(tip_df_all$group)]
tip_df_all <- tip_df_all %>%
  mutate(group = factor(group, levels = tip_levels))
tip_pop_cols <- tip_pop_cols[tip_levels]

p_base <- suppressWarnings(
  ggtree(td_root, layout = "rectangular") +
    geom_tree(linewidth = 0.35, color = "grey65")
) %<+% tip_df_all

p_base$data$support <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*", "\\1", p_base$data$label)))

children_by_node <- split(tr_root2$edge[, 2], tr_root2$edge[, 1])
count_descendant_tips <- function(node) {
  if (node <= Ntip(tr_root2)) {
    return(1L)
  }
  
  children <- children_by_node[[as.character(node)]]
  if (length(children) == 0) {
    return(0L)
  }
  
  sum(vapply(children, count_descendant_tips, integer(1)))
}

support_symbol_df <- p_base$data %>%
  filter(!isTip, support >= 80, !(node %in% out_nodes)) %>%
  mutate(descendant_tips = vapply(node, count_descendant_tips, integer(1))) %>%
  mutate(
    support_class = case_when(
      support == 100 ~ "100",
      support >= 95  ~ "95-99",
      support >= 80  ~ "80-94",
      TRUE           ~ NA_character_
    ),
    support_class = factor(support_class, levels = c("100", "95-99", "80-94"))
  ) %>%
  filter(descendant_tips >= min_backbone_descendant_tips)

tip_y <- p_base$data %>%
  filter(isTip) %>%
  select(label, y)

xmax <- max(p_base$data$x, na.rm = TRUE)
bar_x0 <- xmax + 0.020 * xmax
bar_width <- 0.470 * xmax
tile_gap <- 0.020 * xmax
tile_width <- 0.105 * xmax

admix_plot_df <- admix_long %>%
  inner_join(tip_y, by = "label") %>%
  mutate(
    xmin = bar_x0 + x0_prop * bar_width,
    xmax = bar_x0 + x1_prop * bar_width,
    ymin = y - 0.42,
    ymax = y + 0.42
  )

cluster_plot_df <- tip_y %>%
  filter(!label %in% out_tips) %>%
  inner_join(mclust_df, by = "label") %>%
  mutate(
    xmin = bar_x0 + bar_width + tile_gap,
    xmax = xmin + tile_width,
    ymin = y - 0.48,
    ymax = y + 0.48
  )

write_csv(
  cluster_plot_df %>% select(label, cluster),
  file.path(out_dir, "tree_mclust_tip_annotations.csv")
)

scale_ticks <- tibble(
  prop = seq(0, 1, by = 0.25),
  x = bar_x0 + prop * bar_width,
  y = max(tip_y$y, na.rm = TRUE) + 9.0,
  label = scales::percent(prop, accuracy = 1)
)

scale_line <- tibble(
  x = bar_x0,
  xend = bar_x0 + bar_width,
  y = max(scale_ticks$y),
  yend = max(scale_ticks$y)
)

p <- p_base +
  geom_tippoint(aes(color = group), size = 1.2, alpha = 0.95) +
  scale_color_manual(values = tip_pop_cols, limits = names(tip_pop_cols), drop = FALSE, name = "Population") +
  geom_point(
    data = support_symbol_df,
    aes(x = x, y = y, shape = support_class, fill = support_class),
    inherit.aes = FALSE,
    size = 2.3,
    color = "black",
    stroke = 0.55
  ) +
  scale_shape_manual(
    values = c("100" = 21, "95-99" = 21, "80-94" = 21),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("100" = "black", "95-99" = "grey55", "80-94" = "white"),
    name = "Node support",
    guide = guide_legend(order = 2, override.aes = list(shape = 21, color = "black"))
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(
    data = admix_plot_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = CladeComp),
    inherit.aes = FALSE,
    color = NA
  ) +
  scale_fill_manual(
    values = clade_cols,
    drop = FALSE,
    name = "K=7 ancestry",
    guide = guide_legend(order = 3)
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(
    data = cluster_plot_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
    inherit.aes = FALSE,
    color = "grey10",
    linewidth = 0.12
  ) +
  scale_fill_manual(
    values = cluster_cols,
    drop = FALSE,
    name = "Mclust cluster",
    guide = guide_legend(order = 4)
  ) +
  geom_segment(
    data = scale_line,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "grey20"
  ) +
  geom_segment(
    data = scale_ticks,
    aes(x = x, xend = x, y = y - 0.35, yend = y + 0.35),
    inherit.aes = FALSE,
    linewidth = 0.30,
    color = "grey20"
  ) +
  geom_text(
    data = scale_ticks,
    aes(x = x, y = y - 2.25, label = label),
    inherit.aes = FALSE,
    size = 2.6,
    color = "grey20"
  ) +
  annotate(
    "text",
    x = bar_x0 + 0.5 * bar_width,
    y = max(scale_ticks$y) + 2.35,
    label = "Ancestry proportion",
    size = 3.0,
    fontface = "bold",
    color = "grey20"
  ) +
  theme_tree() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.90, 0.52),
    legend.justification = c(0, 0.5),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.30, "cm"),
    legend.box.spacing = unit(0, "pt"),
    plot.margin = margin(16, 8, 14, 30)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3))
  ) +
  coord_cartesian(
    ylim = c(min(tip_y$y, na.rm = TRUE) - 2, max(scale_ticks$y) + 7),
    clip = "off"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.24))) +
  geom_treescale(x = 0, y = 100, width = 0.1, fontsize = 3)

print(p)

ggsave(
  file.path(fig_dir, "phylogeny_K7_admixture_mclust.png"),
  p,
  width = 13.5,
  height = 16,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "phylogeny_K7_admixture_mclust.pdf"),
  p,
  width = 13.5,
  height = 16,
  bg = "white"
)

showfile.gds(closeall = TRUE)
