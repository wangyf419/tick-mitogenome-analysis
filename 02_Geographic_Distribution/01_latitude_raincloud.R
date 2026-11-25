# ------------------------------------------------------------------------------
# 01_latitude_raincloud.R
#
# Raincloud plot of species latitude centres for the six focal genera, with
# BH-adjusted Dunn brackets. The Kruskal-Wallis and Dunn results are read from
# stored tables; no test is rerun here.
#
# Run from the module directory:
#   Rscript scripts/01_latitude_raincloud.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

MODULE <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(MODULE, "scripts"))) {
  stop("Run this script from the module directory (the folder containing data/ and scripts/).")
}
dir.create(file.path(MODULE, "outputs"), showWarnings = FALSE, recursive = TRUE)

output_dir <- file.path(MODULE, "outputs")
input_dir <- file.path(MODULE, "data/latitude")

focal_genera <- c(
  "Amblyomma", "Dermacentor", "Haemaphysalis",
  "Hyalomma", "Ixodes", "Rhipicephalus"
)

species_centers <- read.csv(file.path(input_dir, "03_species_latitude_centers.csv"))
genus_summary <- read.csv(file.path(input_dir, "04_genus_summary_species_level.csv"))
kw_primary <- read.csv(
  file.path(input_dir, "05_kruskal_wallis_species_level.csv"),
  check.names = FALSE
)
dunn_primary <- read.csv(
  file.path(input_dir, "07_dunn_posthoc_BH_species_level.csv"),
  check.names = FALSE
)

stopifnot(nrow(species_centers) == 154, nrow(dunn_primary) == 15)

format_p <- function(p) {
  ifelse(p < .001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p))
}

write_csv_utf8 <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE)
}

colours <- c(
  Amblyomma = "#0072B2", Dermacentor = "#882255",
  Haemaphysalis = "#D55E00", Hyalomma = "#C39A00",
  Ixodes = "#009E73", Rhipicephalus = "#AA4499"
)

# ------------------------------------------------------------------------------
# Display geometry: density, box and point layers
# ------------------------------------------------------------------------------

densities <- bind_rows(lapply(seq_along(focal_genera), function(i) {
  v <- species_centers$latitude_center_median[
    as.character(species_centers$genus) == focal_genera[i]
  ]
  den <- density(v, from = -50, to = 70, n = 512, adjust = 0.9)
  data.frame(
    genus = focal_genera[i], x = den$x, base = i + .24,
    top = i + .24 + .28 * den$y / max(den$y)
  )
}))

boxes <- bind_rows(lapply(seq_along(focal_genera), function(i) {
  v <- species_centers$latitude_center_median[
    as.character(species_centers$genus) == focal_genera[i]
  ]
  b <- boxplot.stats(v, coef = 1.5)$stats
  data.frame(
    genus = focal_genera[i], y = i + .07,
    lo = b[1], q1 = b[2], med = b[3], q3 = b[4], hi = b[5]
  )
}))

set.seed(20260910)
points <- species_centers %>% mutate(
  y = match(as.character(genus), focal_genera) - .16 + runif(n(), -.08, .08)
)

comparisons <- dunn_primary %>%
  filter(p.adj < 0.05) %>%
  mutate(
    y1 = match(group1, focal_genera) + .07,
    y2 = match(group2, focal_genera) + .07
  ) %>%
  arrange(y1, y2)

comparisons$x <- 78 + (seq_len(nrow(comparisons)) - 1) * 5.3
comparisons$symbol <- case_when(
  comparisons$p.adj < .0001 ~ "****",
  comparisons$p.adj < .001 ~ "***",
  comparisons$p.adj < .01 ~ "**",
  comparisons$p.adj < .05 ~ "*",
  TRUE ~ "ns"
)

write_csv_utf8(comparisons, "19_plot_pairwise_annotation_positions.csv")

counts <- setNames(genus_summary$species_n, as.character(genus_summary$genus))
labels <- paste0(focal_genera, "\n(n = ", counts[focal_genera], ")")

# ------------------------------------------------------------------------------
# Plot
#
# exact_p = FALSE gives significance symbols, exact_p = TRUE gives the adjusted
# P values themselves. Both variants share the same geometry.
# ------------------------------------------------------------------------------

plot_raincloud <- function(exact_p = FALSE) {
  ann <- comparisons
  ann$label <- if (exact_p) {
    vapply(ann$p.adj, format_p, character(1))
  } else {
    ann$symbol
  }

  ggplot() +
    geom_ribbon(
      data = densities,
      aes(
        x = x, ymin = base, ymax = top, group = genus,
        fill = genus, colour = genus
      ),
      alpha = .36, linewidth = .55
    ) +
    geom_segment(
      data = boxes,
      aes(x = lo, xend = hi, y = y, yend = y, colour = genus),
      linewidth = .65
    ) +
    geom_segment(
      data = boxes,
      aes(x = lo, xend = lo, y = y - .045, yend = y + .045, colour = genus),
      linewidth = .55
    ) +
    geom_segment(
      data = boxes,
      aes(x = hi, xend = hi, y = y - .045, yend = y + .045, colour = genus),
      linewidth = .55
    ) +
    geom_rect(
      data = boxes,
      aes(
        xmin = q1, xmax = q3, ymin = y - .065, ymax = y + .065,
        fill = genus, colour = genus
      ),
      linewidth = .7, alpha = .45
    ) +
    geom_segment(
      data = boxes,
      aes(x = med, xend = med, y = y - .065, yend = y + .065, colour = genus),
      linewidth = .9
    ) +
    geom_point(
      data = points,
      aes(x = latitude_center_median, y = y, colour = genus),
      size = 1.55, alpha = .75
    ) +
    geom_segment(
      data = ann, aes(x = x, xend = x, y = y1, yend = y2), linewidth = .35
    ) +
    geom_segment(
      data = ann, aes(x = x - 1.7, xend = x, y = y1, yend = y1), linewidth = .35
    ) +
    geom_segment(
      data = ann, aes(x = x - 1.7, xend = x, y = y2, yend = y2), linewidth = .35
    ) +
    geom_text(
      data = ann, aes(x = x + 1.45, y = (y1 + y2) / 2, label = label),
      angle = 90, size = if (exact_p) 2.35 else 2.8, hjust = .5, vjust = .5
    ) +
    annotate(
      "text",
      x = 116, y = 6.66, label = "Pairwise Dunn tests (BH-adjusted)",
      size = 3.05
    ) +
    annotate(
      "text",
      x = -49, y = 6.7, label = "a",
      size = 5.4, fontface = "bold", hjust = 0
    ) +
    annotate(
      "text",
      x = 8, y = 6.69,
      label = paste0(
        "Kruskal-Wallis: H = ", sprintf("%.2f", kw_primary$statistic),
        ", df = ", kw_primary$df, ", P = ", format_p(kw_primary$p)
      ),
      size = 2.95, hjust = .5
    ) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_fill_manual(values = colours, guide = "none") +
    scale_x_continuous(
      limits = c(-50, 159), breaks = seq(-50, 70, 20),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(.40, 6.88), breaks = seq_along(focal_genera),
      labels = labels, expand = expansion(mult = c(0, 0))
    ) +
    # The latitude baseline stays under the latitude panel only; the brackets
    # carry no geographic units.
    annotate(
      "segment",
      x = -50, xend = 70, y = .4, yend = .4, linewidth = .45
    ) +
    labs(
      x = "Species latitude centre (degrees)", y = NULL,
      caption = if (exact_p) {
        "Significant comparisons only; labels are BH-adjusted P values."
      } else {
        "BH-adjusted P: * < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001."
      }
    ) +
    theme_classic(base_size = 11, base_family = "Arial") +
    theme(
      axis.text.y = element_text(face = "italic", colour = "black", size = 10),
      axis.text.x = element_text(colour = "black"),
      axis.ticks.y = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_text(hjust = .20),
      plot.caption = element_text(size = 8, hjust = 1),
      plot.margin = margin(10, 12, 8, 8)
    )
}

# ------------------------------------------------------------------------------
# Write both variants
# ------------------------------------------------------------------------------

for (exact in c(FALSE, TRUE)) {
  p <- plot_raincloud(exact)
  name <- if (exact) {
    "18b_six_genera_raincloud_adjusted_P"
  } else {
    "18_six_genera_raincloud_pairwise"
  }

  ggsave(
    file.path(output_dir, paste0(name, ".pdf")), p,
    width = 240, height = 185, units = "mm", device = cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(name, ".png")), p,
    width = 240, height = 185, units = "mm", dpi = 600, bg = "white"
  )

  if (!exact) {
    ggsave(
      file.path(output_dir, "preview.png"), p,
      width = 240, height = 185, units = "mm", dpi = 150, bg = "white"
    )
  }
}
