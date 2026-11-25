# ------------------------------------------------------------------------------
# 01_major_tick_dated_tree.R
#
# Dated tree of the major tick lineage with numeric 95% HPD labels

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(deeptime)
  library(treeio)
  library(ape)
  library(grid)
})

MODULE <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(MODULE, "scripts"))) {
  stop("Run this script from the module directory (the folder containing data/ and scripts/).")
}
dir.create(file.path(MODULE, "outputs"), showWarnings = FALSE, recursive = TRUE)

tree_file <- file.path(MODULE, "data/major_tick/FigTree.tre")
output_file <- file.path(MODULE, "outputs/Figure.pdf")
font_family <- "Arial"

# The input tree stores ages in units of 100 Myr.
time_multiplier <- 100

if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

# ------------------------------------------------------------------------------
# Read the annotated tree
#
# read.mcmctree() accepts UTREE blocks and keeps the true internal node that each
# HPD annotation belongs to, so the intervals are not re-associated by position.
# ------------------------------------------------------------------------------

annotated_tree <- treeio::read.mcmctree(tree_file)
tree <- annotated_tree@phylo
annotations <- as.data.frame(annotated_tree@data)

if (!"0.95HPD" %in% names(annotations)) {
  stop("No 95% HPD annotation found in the tree file.")
}

hpd_data <- data.frame(
  node = as.integer(annotations$node),
  hpd_low = vapply(
    annotations[["0.95HPD"]], function(x) as.numeric(x[1]), numeric(1)
  ),
  hpd_high = vapply(
    annotations[["0.95HPD"]], function(x) as.numeric(x[2]), numeric(1)
  )
)
hpd_data$hpd_low <- hpd_data$hpd_low * time_multiplier
hpd_data$hpd_high <- hpd_data$hpd_high * time_multiplier

stopifnot(
  nrow(hpd_data) == Nnode(tree),
  !anyDuplicated(hpd_data$node),
  all(is.finite(hpd_data$hpd_low)),
  all(hpd_data$hpd_low < hpd_data$hpd_high)
)

tree$edge.length <- tree$edge.length * time_multiplier

depth <- node.depth.edgelength(tree)
if (diff(range(depth[seq_len(Ntip(tree))])) > 0.01) {
  stop(
    "Tip ages differ by more than 0.01 Myr; ",
    "check the time unit and the ultrametricity of the tree."
  )
}

tip_data <- data.frame(
  label = tree$tip.label,
  Genus = sub("[ _].*", "", tree$tip.label),
  label_clean = gsub("_", " ", tree$tip.label, fixed = TRUE)
)

# ------------------------------------------------------------------------------
# Genus colours
#
# Colouring by genus is kept, but the red/green combinations of the original
# palette are avoided. Genera that were missing before are added, and the
# outgroups are included.
# ------------------------------------------------------------------------------

genus_colors <- c(
  Ixodes = "#0072B2", Haemaphysalis = "#CC79A7",
  Rhipicephalus = "#D55E00", Dermacentor = "#56B4E9",
  Amblyomma = "#8491B4", Hyalomma = "#E69F00",
  Ornithodoros = "#882255", Argas = "#332288",
  Alectorobius = "#B09C85", Alveonasus = "#3C5488",
  Bothriocroton = "#8C6D31", Nuttalliella = "#777777",
  Centruroides = "#555555", Parasteatoda = "#555555"
)
stopifnot(all(tip_data$Genus %in% names(genus_colors)))

# ------------------------------------------------------------------------------
# Tree geometry and HPD label positions
# ------------------------------------------------------------------------------

p <- ggtree(tree, linewidth = 0.28) |> revts()
p <- p %<+% tip_data
coordinates <- as.data.frame(p$data)

hpd_data <- merge(
  hpd_data, coordinates[, c("node", "y", "x")],
  by = "node", sort = FALSE
)
hpd_data$age_Ma <- abs(hpd_data$x)
# En dash between the interval bounds, as in the original figure label.
hpd_data$hpd_label <- sprintf("%.1f\u2013%.1f", hpd_data$hpd_low, hpd_data$hpd_high)

stopifnot(
  all(hpd_data$age_Ma >= hpd_data$hpd_low),
  all(hpd_data$age_Ma <= hpd_data$hpd_high)
)

# Major clade names are placed next to the parent node of each clade.
clade_specs <- list(
  Ixodidae = c("Ixodes_ricinus", "Rhipicephalus_microplus"),
  Metastriata = c("Bothriocroton_concolor", "Rhipicephalus_microplus"),
  Prostriata = c("Ixodes_ricinus", "Ixodes_scapularis"),
  Argasidae = c("Argas_persicus", "Alectorobius_microlophi"),
  Nuttalliellidae = "Nuttalliella_namaqua"
)

clade_data <- do.call(rbind, lapply(names(clade_specs), function(label) {
  tips <- clade_specs[[label]]
  node <- if (length(tips) == 1L) {
    match(tips, tree$tip.label)
  } else {
    getMRCA(tree, tips)
  }
  parent <- tree$edge[tree$edge[, 2] == node, 1]
  coord <- coordinates[coordinates$node == node, ]
  parent_x <- coordinates$x[coordinates$node == parent]
  below <- label %in% c("Ixodidae", "Metastriata")
  offset <- if (below) -0.55 else if (label == "Argasidae") 0.95 else 0.5
  data.frame(
    label = label,
    x = parent_x + 8,
    y = coord$y + offset,
    vjust = if (below) 1 else 0
  )
}))

legend_genera <- c(
  "Alectorobius", "Alveonasus", "Amblyomma", "Argas", "Dermacentor",
  "Haemaphysalis", "Hyalomma", "Ixodes", "Ornithodoros",
  "Rhipicephalus", "Nuttalliella"
)

# ggplot2 geometric text sizes are given in millimetres; convert to points.
pt <- 72.27 / 25.4

p <- p +
  aes(color = Genus) +
  geom_tiplab(
    aes(label = label_clean),
    size = 8 / pt, family = font_family, fontface = "italic",
    align = TRUE, linetype = "dotted", linesize = 0.25, offset = 3,
    show.legend = FALSE
  ) +
  geom_text(
    data = hpd_data, aes(x = x, y = y, label = hpd_label),
    inherit.aes = FALSE, size = 7 / pt, family = font_family,
    color = "black", vjust = -0.65, hjust = 0.5
  ) +
  geom_text(
    data = clade_data, aes(x = x, y = y, label = label, vjust = vjust),
    inherit.aes = FALSE, size = 7.5 / pt, family = font_family,
    color = "black", hjust = 0, fontface = "plain"
  ) +
  geom_nodepoint(color = "#4A82B0", size = 1.1, alpha = 0.9) +
  scale_color_manual(
    values = genus_colors, breaks = legend_genera,
    na.value = "grey50", na.translate = TRUE
  ) +
  theme_tree2()

# ------------------------------------------------------------------------------
# Geological time axes
#
# Bands outside 0-500 Myr are truncated. Very narrow epochs keep the band but
# omit the text that would not fit; short names are used where space is tight.
# ------------------------------------------------------------------------------

epochs <- subset(deeptime::epochs, min_age < 500 & max_age > 0)
periods <- subset(deeptime::periods, min_age < 500 & max_age > 0)
epochs$max_age <- pmin(epochs$max_age, 500)
periods$max_age <- pmin(periods$max_age, 500)
epochs$abbr[epochs$max_age - epochs$min_age < 15] <- ""

short_periods <- c(
  Cambrian = "Cm", Ordovician = "Ord.", Silurian = "Sil.",
  Neogene = "Neo.", Quaternary = "Q"
)
periods$name <- ifelse(
  periods$name %in% names(short_periods),
  short_periods[periods$name], periods$name
)
periods$name[periods$max_age - periods$min_age < 8] <- ""

p_final <- p +
  scale_x_continuous(
    breaks = seq(-500, 0, 50), labels = abs(seq(-500, 0, 50)),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_geo(
    dat = list(epochs, periods), abbrv = list(TRUE, FALSE),
    size = list(6.5 / pt, 6.5 / pt),
    height = list(unit(10, "bigpts"), unit(12, "bigpts")),
    family = font_family, color = "black", lwd = 0.25,
    bord = c("white", "white"), neg = TRUE, pos = "bottom",
    xlim = c(-500, 0), ylim = c(-1, Ntip(tree) + 0.8), clip = "off"
  ) +
  scale_y_continuous(expand = expansion(add = c(2.0, 0.8))) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(
    text = element_text(family = font_family, size = 8, face = "plain"),
    legend.position = "inside", legend.position.inside = c(0.015, 0.97),
    legend.justification = c(0, 1),
    legend.title = element_text(size = 8, face = "plain"),
    legend.text = element_text(size = 7.5, face = "italic"),
    legend.key.height = unit(10, "bigpts"),
    legend.key.width = unit(8, "bigpts"),
    legend.spacing.x = unit(2, "bigpts"),
    legend.background = element_blank(),
    axis.line.x = element_line(linewidth = 0.265),
    axis.ticks.x = element_line(linewidth = 0.265),
    axis.ticks.length.x = unit(2, "bigpts"),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.title.x = element_text(size = 8, margin = margin(t = 5)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 4, r = 150, b = 4, l = 8, unit = "pt")
  ) +
  labs(x = "Million years ago", color = "Tick Genus")

# ------------------------------------------------------------------------------
# Write the vector PDF
#
# The outer ruler follows the official template; bigpts = 1/72 inch, which
# matches the point unit of the PDF device.
# ------------------------------------------------------------------------------

grDevices::cairo_pdf(
  output_file,
  width = 8.5, height = 11,
  family = font_family, bg = "white", onefile = TRUE
)
grid.newpage()
pushViewport(viewport(
  x = unit(51.0234, "bigpts"), y = unit(53.4896534, "bigpts"),
  width = unit(510.2363858, "bigpts"), height = unit(687.4869466, "bigpts"),
  just = c("left", "bottom")
))
grid.draw(ggplotGrob(p_final))
popViewport()
invisible(dev.off())

# ------------------------------------------------------------------------------
# Optional audit record
#
# Only written when the caller sets MCMCTREE_FIGURE_AUDIT; no PNG is produced.
# ------------------------------------------------------------------------------

audit_file <- Sys.getenv("MCMCTREE_FIGURE_AUDIT", "")
if (nzchar(audit_file)) {
  write.csv(hpd_data, audit_file, row.names = FALSE)
}

soft_hard <- getMRCA(tree, c("Argas_persicus", "Ixodes_ricinus"))
hard_split <- getMRCA(tree, c("Ixodes_ricinus", "Rhipicephalus_microplus"))

message("Internal nodes / HPD rows: ", Nnode(tree), " / ", nrow(hpd_data))
print(hpd_data[
  hpd_data$node %in% c(soft_hard, hard_split),
  c("node", "age_Ma", "hpd_low", "hpd_high")
])
message("PDF written to: ", output_file)
