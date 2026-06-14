###############################
## 04_tree.R
## Rooted ML tree with shortened long branch, tip points, and tip-order color strip
###############################

helper_file <- file.path(getwd(), "honeybee_scripts", "00_rstudio_helpers.R")
if (!file.exists(helper_file)) {
  helper_file <- file.path(getwd(), "00_rstudio_helpers.R")
}
if (!file.exists(helper_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  helper_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "00_rstudio_helpers.R")
}
source(helper_file)
script_file <- resolve_script_path("04_tree.R")

source(file.path(dirname(script_file), "00_config.R"))

suppressPackageStartupMessages({
  library(ape)
  library(treeio)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(tibble)
})

# ---- inputs ----
tree_path <- file.path(base_path, "gtr_ascR4_nonsyn_constrained_varsites.contree")

# ---- read ----
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

# ---- root on outgroup ----
out_tips <- tip_df_all$label[tip_df_all$group == "Outgroup"]
stopifnot(length(out_tips) > 0)

tr_root <- root(tr, outgroup = out_tips, resolve.root = TRUE)

# ---- identify rooted outgroup subtree ----
rootnode <- Ntip(tr_root) + 1
kids     <- tr_root$edge[tr_root$edge[, 1] == rootnode, 2]
stopifnot(length(kids) == 2)

out_child <- MRCA(tr_root, out_tips)
stopifnot(!is.na(out_child))

# nodes to exclude from support plotting:
# outgroup MRCA + all descendants of that node
out_desc <- treeio::offspring(as.treedata(tr_root), out_child)

if (is.data.frame(out_desc) && "node" %in% names(out_desc)) {
  out_nodes <- c(out_child, out_desc$node)
} else {
  out_nodes <- out_child
}

ingroup_child <- kids[kids != out_child]
stopifnot(length(ingroup_child) == 1)

# ============================================================
# Shrink the actual longest edge in the entire tree
# ============================================================
tr_root2 <- tr_root

longest_edge <- which.max(tr_root2$edge.length)
target_len <- 0.15
tr_root2$edge.length[longest_edge] <- target_len

dash_parent <- tr_root2$edge[longest_edge, 1]
dash_child  <- tr_root2$edge[longest_edge, 2]

td_root <- as.treedata(tr_root2)

# ---- colors ----
my_cols <- c(
  pop_colors,
  Outgroup = "black",
  Unknown  = "grey80"
)

legend_order <- c("NorCal", "SoCal", "Meso", "A", "C", "M", "O", "L", "U", "Y", "Outgroup", "Unknown")
lvl <- legend_order[legend_order %in% unique(tip_df_all$group)]

tip_df_all <- tip_df_all %>%
  mutate(group = factor(group, levels = lvl))

my_cols <- my_cols[lvl]

# ---- base tree ----
p_base <- suppressWarnings(
  ggtree(td_root, layout = "rectangular") +
    geom_tree(linewidth = 0.35, color = "grey65")
) %<+% tip_df_all

# ---- main plot ----
p <- p_base +
  geom_tippoint(aes(color = group), size = 1.2, alpha = 0.95) +
  scale_color_manual(values = my_cols, limits = names(my_cols), drop = FALSE) +
  theme_tree() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.33, "cm"),
    plot.margin      = margin(10, 16, 14, 30)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3))) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.12))) +
  geom_treescale(x = 0, y = 100, width = 0.1, fontsize = 3)

# ---- 100% support as black dots ----
p$data$support <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*", "\\1", p$data$label)))

p <- p +
  geom_point2(
    aes(subset = !isTip & support == 100 & !(node %in% out_nodes)),
    size = 1.5,
    shape = 25,
    fill = "black",
    color = "white",
    stroke = 0.4
  )

# ============================================================
# Solid tip-order color strip
# ============================================================
tip_bar <- p_base$data %>%
  dplyr::filter(isTip) %>%
  dplyr::select(label, y) %>%
  dplyr::left_join(tip_df_all, by = "label")

xmax  <- max(p_base$data$x, na.rm = TRUE)
x_bar <- xmax + 0.02 * xmax

p <- p +
  geom_tile(
    data = tip_bar,
    aes(x = x_bar, y = y, fill = group),
    width  = 0.02 * xmax,
    height = 1,
    color  = NA
  ) +
  scale_fill_manual(values = my_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.18)))

print(p)

ggsave(
  file.path(fig_dir, "phylogeny.png"),
  p,
  width = 8.5,
  height = 11.5,
  dpi = 300,
  bg = "white"
)
