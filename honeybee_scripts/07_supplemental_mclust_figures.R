###############################
## 07_supplemental_mclust_figures.R
## Supplemental PC-AiR and phylogeny figures colored by Mclust clusters.
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("07_supplemental_mclust_figures.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(ape)
  library(treeio)
  library(ggtree)
  library(ggplot2)
  library(ggstar)
  library(dplyr)
  library(readr)
  library(tibble)
  library(grid)
})

# -----------------------------
# Inputs and shared colors
# -----------------------------
pcair_mclust_path <- file.path(out_dir, "pcair_scores_with_mclust.csv")
pcair_var_path <- file.path(out_dir, "pcair_variance_explained.csv")
kingmat_path <- file.path(out_dir, "king_matrix.rds")
tree_path <- file.path(base_path, "gtr_ascR4_nonsyn_constrained_varsites.contree")
min_backbone_descendant_tips <- 20

required_files <- c(pcair_mclust_path, pcair_var_path, kingmat_path, tree_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    sprintf("Missing required file(s): %s", paste(missing_files, collapse = ", ")),
    call. = FALSE
  )
}

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

pc_df_all <- read_csv(pcair_mclust_path, show_col_types = FALSE) %>%
  mutate(
    cluster = factor(as.character(cluster), levels = sort(unique(as.character(cluster)))),
    pop = factor(pop, levels = lvl),
    relatedness = factor(relatedness, levels = c("Unrelated", "Related"))
  )

pc_percent_tbl <- read_csv(pcair_var_path, show_col_types = FALSE)
pc_percent <- setNames(pc_percent_tbl$percent_variance, pc_percent_tbl$PC)

cluster_cols <- cluster_palette[levels(pc_df_all$cluster)]

# -----------------------------
# Supplemental PC-AiR, all samples colored by Mclust cluster
# -----------------------------
p_pcair_mclust <- ggplot(pc_df_all, aes(PC1, PC2)) +
  geom_star(
    data = subset(pc_df_all, relatedness == "Unrelated"),
    aes(starshape = relatedness, fill = cluster),
    colour = "black",
    starstroke = 0.3,
    size = 3
  ) +
  geom_star(
    data = subset(pc_df_all, relatedness == "Related"),
    aes(starshape = relatedness, fill = cluster),
    colour = "black",
    starstroke = 0.3,
    size = 4
  ) +
  scale_fill_manual(values = cluster_cols, drop = FALSE, name = "Mclust cluster") +
  scale_starshape_manual(
    values = c(Unrelated = 15, Related = 9),
    limits = c("Unrelated", "Related"),
    name = "Relatedness"
  ) +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc_percent["PC1"]),
    y = sprintf("PC2 (%.2f%%)", pc_percent["PC2"])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

save_plot("pcair_mclust_clusters.png", width = 12, height = 8)
print(p_pcair_mclust)
dev.off()

# -----------------------------
# Supplemental SoCal/NorCal PC-AiR with kinship edges, colored by Mclust cluster
# -----------------------------
kin_classify <- function(k) {
  dplyr::case_when(
    k >= 0.0884 ~ "1st/2nd",
    k >= 0.0221 ~ "3rd/4th",
    TRUE        ~ NA_character_
  )
}

pc_df_sub <- pc_df_all %>%
  filter(pop %in% c("SoCal", "NorCal")) %>%
  mutate(pop = factor(pop, levels = c("NorCal", "SoCal")))

KINGmat <- as.matrix(readRDS(kingmat_path))
ids_sub <- pc_df_sub$sample.id_raw
KING_sub <- as.matrix(KINGmat[ids_sub, ids_sub, drop = FALSE])
ut <- which(upper.tri(KING_sub), arr.ind = TRUE)

edges_facet_xy <- tibble(
  id1 = rownames(KING_sub)[ut[, "row"]],
  id2 = colnames(KING_sub)[ut[, "col"]],
  kinship = KING_sub[ut]
) %>%
  mutate(degree = kin_classify(kinship)) %>%
  filter(!is.na(degree)) %>%
  left_join(pc_df_sub %>% select(sample.id_raw, pop),
            by = c("id1" = "sample.id_raw")) %>%
  rename(pop1 = pop) %>%
  left_join(pc_df_sub %>% select(sample.id_raw, pop),
            by = c("id2" = "sample.id_raw")) %>%
  rename(pop2 = pop) %>%
  filter(pop1 == pop2) %>%
  mutate(pop = pop1) %>%
  select(-pop1, -pop2) %>%
  left_join(pc_df_sub %>% select(sample.id_raw, PC1, PC2),
            by = c("id1" = "sample.id_raw")) %>%
  rename(x = PC1, y = PC2) %>%
  left_join(pc_df_sub %>% select(sample.id_raw, PC1, PC2),
            by = c("id2" = "sample.id_raw")) %>%
  rename(xend = PC1, yend = PC2) %>%
  filter(!is.na(x), !is.na(xend), !is.na(y), !is.na(yend)) %>%
  mutate(degree = factor(degree, levels = c("1st/2nd", "3rd/4th")))

p_socal_norcal_mclust <- ggplot(pc_df_sub, aes(PC1, PC2)) +
  geom_segment(
    data = edges_facet_xy,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = degree),
    inherit.aes = FALSE,
    linewidth = 0.42,
    alpha = 0.70,
    color = "grey20"
  ) +
  geom_point(aes(fill = cluster), shape = 21, colour = "black", size = 3, stroke = 0.3) +
  facet_wrap(~ pop, ncol = 2, scales = "free") +
  scale_fill_manual(values = cluster_cols, drop = TRUE, name = "Mclust cluster") +
  scale_linetype_manual(
    values = c("1st/2nd" = "solid", "3rd/4th" = "dashed"),
    name = "Kinship degree"
  ) +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc_percent["PC1"]),
    y = sprintf("PC2 (%.2f%%)", pc_percent["PC2"])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

save_plot("pcair_SoCal_NorCal_mclust_clusters.png", width = 12, height = 6.5)
print(p_socal_norcal_mclust)
dev.off()

# -----------------------------
# Supplemental phylogenies with numeric support labels >75
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
  left_join(
    pc_df_all %>% transmute(label = sample.id_clean, cluster = as.character(cluster)),
    by = "label"
  ) %>%
  mutate(
    group = if_else(is.na(group), "Unknown", group),
    cluster_plot = if_else(group == "Outgroup", "Outgroup", cluster)
  )

out_tips <- tip_df_all$label[tip_df_all$group == "Outgroup"]
stopifnot(length(out_tips) > 0)

tr_root <- root(tr, outgroup = out_tips, resolve.root = TRUE)
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

tree_base <- suppressWarnings(
  ggtree(td_root, layout = "rectangular") +
    geom_tree(linewidth = 0.35, color = "grey65")
) %<+% tip_df_all

tree_base$data$support <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*", "\\1", tree_base$data$label)))

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

support_symbol_df <- tree_base$data %>%
  filter(!isTip, support >= 80) %>%
  mutate(
    descendant_tips = vapply(node, count_descendant_tips, integer(1)),
    support_class = case_when(
      support == 100 ~ "100",
      support >= 95  ~ "95-99",
      support >= 80  ~ "80-94",
      TRUE           ~ NA_character_
    ),
    support_class = factor(support_class, levels = c("100", "95-99", "80-94"))
  ) %>%
  filter(descendant_tips >= min_backbone_descendant_tips)

make_support_tree <- function(point_aes, color_scale, legend_title, output_stem) {
  p <- tree_base +
    geom_tippoint(point_aes, size = 1.25, alpha = 0.95) +
    color_scale +
    geom_point(
      data = support_symbol_df,
      aes(x = x, y = y, shape = support_class, fill = support_class),
      inherit.aes = FALSE,
      size = 2.3,
      color = "black",
      stroke = 0.45
    ) +
    scale_shape_manual(
      values = c("100" = 21, "95-99" = 21, "80-94" = 21),
      name = "Node support"
    ) +
    scale_fill_manual(
      values = c("100" = "black", "95-99" = "grey55", "80-94" = "white"),
      name = "Node support"
    ) +
    theme_tree() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.30, "cm"),
      plot.margin = margin(10, 16, 14, 30)
    ) +
    guides(
      color = guide_legend(title = legend_title, override.aes = list(size = 3)),
      shape = guide_legend(order = 2),
      fill = guide_legend(order = 2)
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.16))) +
    geom_treescale(x = 0, y = 100, width = 0.1, fontsize = 3)
  
  ggsave(file.path(fig_dir, paste0(output_stem, ".png")),
         p, width = 10.5, height = 16, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(output_stem, ".pdf")),
         p, width = 10.5, height = 16, bg = "white")
  p
}

p_tree_support_pop <- make_support_tree(
  point_aes = aes(color = group),
  color_scale = scale_color_manual(values = tip_pop_cols, limits = names(tip_pop_cols), drop = FALSE),
  legend_title = "Population",
  output_stem = "phylogeny_support_gt75_population"
)
print(p_tree_support_pop)

cluster_tree_levels <- c(levels(pc_df_all$cluster), "Outgroup")
cluster_tree_cols <- c(cluster_cols, Outgroup = "black")
tip_df_all <- tip_df_all %>%
  mutate(cluster_plot = factor(cluster_plot, levels = cluster_tree_levels))

tree_base <- suppressWarnings(
  ggtree(td_root, layout = "rectangular") +
    geom_tree(linewidth = 0.35, color = "grey65")
) %<+% tip_df_all
tree_base$data$support <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*", "\\1", tree_base$data$label)))
support_symbol_df <- tree_base$data %>%
  filter(!isTip, support >= 80) %>%
  mutate(
    descendant_tips = vapply(node, count_descendant_tips, integer(1)),
    support_class = case_when(
      support == 100 ~ "100",
      support >= 95  ~ "95-99",
      support >= 80  ~ "80-94",
      TRUE           ~ NA_character_
    ),
    support_class = factor(support_class, levels = c("100", "95-99", "80-94"))
  ) %>%
  filter(descendant_tips >= min_backbone_descendant_tips)

p_tree_support_mclust <- make_support_tree(
  point_aes = aes(color = cluster_plot),
  color_scale = scale_color_manual(values = cluster_tree_cols, limits = names(cluster_tree_cols), drop = FALSE),
  legend_title = "Mclust cluster",
  output_stem = "phylogeny_support_gt75_mclust"
)
print(p_tree_support_mclust)

showfile.gds(closeall = TRUE)
