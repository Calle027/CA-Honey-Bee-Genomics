###############################
## 02_pcair_mclust.R
## PC-AiR + Mclust clustering
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("02_pcair_mclust.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(SNPRelate)
  library(GENESIS)
  library(GWASTools)
  library(ggplot2)
  library(ggstar)
  library(dplyr)
  library(tibble)
  library(mclust)
  library(readr)
  library(tidyr)
})

# -----------------------------
# Shared objects
# -----------------------------
snpset <- load_snpset()
KINGmat <- load_king_matrix()

# Close SNPRelate connection before GWASTools/GENESIS reader
showfile.gds(closeall = TRUE)

geno_reader <- GdsGenotypeReader(gds_path)
geno_data   <- GenotypeData(geno_reader)

mypcair <- pcair(
  geno_data,
  kinobj = KINGmat,
  divobj = KINGmat,
  snp.include = snpset,
  autosome.only = FALSE,
  eigen.cnt = 3
)
pc.percent <- mypcair$varprop * 100

part <- pcairPartition(kinobj = KINGmat, divobj = KINGmat)
relatedness_vector <- ifelse(
  mypcair$sample.id %in% part$unrels,
  "Unrelated",
  "Related"
)

pc_df_all <- tibble(
  sample.id_raw = as.character(mypcair$sample.id),
  relatedness   = factor(relatedness_vector, levels = c("Unrelated", "Related")),
  PC1 = mypcair$vectors[, 1],
  PC2 = mypcair$vectors[, 2],
  PC3 = mypcair$vectors[, 3]
) %>%
  left_join(analysis_samples %>% select(sample.id_raw, sample.id_clean, pop), by = "sample.id_raw") %>%
  mutate(pop = factor(pop, levels = lvl))

write_csv(pc_df_all, file.path(out_dir, "pcair_scores.csv"))
write_csv(tibble(PC = paste0("PC", seq_along(pc.percent)), percent_variance = pc.percent),
          file.path(out_dir, "pcair_variance_explained.csv"))
write_csv(tibble(sample.id_raw = part$unrels), file.path(out_dir, "pcair_unrelated_samples.csv"))
write_csv(tibble(sample.id_raw = part$rels),   file.path(out_dir, "pcair_related_samples.csv"))

p_pcair <- ggplot(pc_df_all, aes(PC1, PC2)) +
  geom_star(
    data = subset(pc_df_all, relatedness == "Unrelated"),
    aes(starshape = relatedness, fill = pop),
    colour = "black", starstroke = 0.3, size = 3
  ) +
  geom_star(
    data = subset(pc_df_all, relatedness == "Related"),
    aes(starshape = relatedness, fill = pop),
    colour = "black", starstroke = 0.3, size = 4
  ) +
  scale_fill_manual(values = pop_colors, drop = FALSE) +
  scale_starshape_manual(
    values = c(Unrelated = 15, Related = 9),
    limits = c("Unrelated", "Related"),
    name = "Relatedness"
  ) +
  guides(
    starshape = guide_legend(
      title = "Relatedness",
      override.aes = list(fill = NA, size = 4, alpha = 1, colour = "black")
    ),
    fill = guide_legend(
      title = "Population",
      override.aes = list(starshape = 15, size = 4, colour = "black")
    )
  ) +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc.percent[1]),
    y = sprintf("PC2 (%.2f%%)", pc.percent[2])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

save_plot("pcair.png", width = 12, height = 8)
print(p_pcair)
dev.off()

# -----------------------------
# Optional SoCal vs NorCal panel
# -----------------------------
pc_df_sub <- pc_df_all %>%
  filter(pop %in% c("SoCal", "NorCal")) %>%
  mutate(pop = factor(pop, levels = c("NorCal", "SoCal")))

kin_classify <- function(k) {
  dplyr::case_when(
    k >= 0.0884 ~ "1st/2nd",
    k >= 0.0221 ~ "3rd/4th",
    TRUE        ~ NA_character_
  )
}

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

write_csv(edges_facet_xy, file.path(out_dir, "pcair_SoCal_NorCal_kinship_edges.csv"))

p_socal_norcal <- ggplot(pc_df_sub, aes(PC1, PC2)) +
  geom_segment(
    data = edges_facet_xy,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = degree),
    inherit.aes = FALSE,
    linewidth = 0.42,
    alpha = 0.70,
    color = "grey20"
  ) +
  geom_point(aes(fill = pop), shape = 21, colour = "black", size = 3, stroke = 0.3) +
  facet_wrap(~ pop, ncol = 2, scales = "free") +
  scale_fill_manual(values = pop_colors[c("NorCal", "SoCal")], name = "Population") +
  scale_linetype_manual(
    values = c("1st/2nd" = "solid", "3rd/4th" = "dashed"),
    name = "Kinship degree"
  ) +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc.percent[1]),
    y = sprintf("PC2 (%.2f%%)", pc.percent[2])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "right"
  )

save_plot("pcair_SoCal_NorCal.png", width = 12, height = 6.5)
print(p_socal_norcal)
dev.off()

# -----------------------------
# Mclust on PC1 and PC2
# -----------------------------
mclust_fit <- Mclust(pc_df_all[, c("PC1", "PC2")])

pc_df_all <- pc_df_all %>%
  mutate(
    cluster = factor(
      mclust_fit$classification,
      levels = sort(unique(mclust_fit$classification)),
      labels = sort(unique(mclust_fit$classification))
    )
  )

write_csv(pc_df_all, file.path(out_dir, "pcair_scores_with_mclust.csv"))
capture.output(summary(mclust_fit), file = file.path(out_dir, "mclust_summary.txt"))

counts <- pc_df_all %>%
  mutate(
    pop = factor(pop, levels = lvl),
    cluster = factor(cluster, levels = sort(unique(as.character(cluster))))
  ) %>%
  dplyr::count(pop, cluster, name = "n")

write_csv(counts, file.path(out_dir, "mclust_population_cluster_counts.csv"))

p_mclust <- ggplot(counts, aes(x = pop, y = cluster)) +
  geom_point(
    aes(size = n, fill = pop),
    shape = 22, color = "black", stroke = 0.35
  ) +
  scale_fill_manual(values = pop_colors, name = "Population", drop = FALSE) +
  scale_size_continuous(
    name = "Sample Count",
    range = c(9, 25),
    breaks = c(10, 20, 30),
    labels = c("10", "20", "+30")
  ) +
  labs(x = "Population", y = "Cluster") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  scale_y_discrete(limits = rev(sort(unique(as.character(counts$cluster))))) +
  guides(
    fill = guide_legend(override.aes = list(shape = 22, size = 6, color = "black")),
    size = guide_legend(order = 2)
  )

save_plot("pop_vs_cluster.png", width = 10, height = 7)
print(p_mclust)
dev.off()

close(geno_reader)
showfile.gds(closeall = TRUE)
