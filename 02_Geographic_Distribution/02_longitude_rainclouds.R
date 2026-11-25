# ------------------------------------------------------------------------------
# 02_longitude_rainclouds_A4.R
#
# Six-panel A4 page of within-genus longitude rainclouds. Pairwise PERMANOVA
# and PERMDISP results are read from stored tables; only the display is built.


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

MODULE <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(MODULE, "scripts"))) {
  stop("Run this script from the module directory (the folder containing data/ and scripts/).")
}
dir.create(file.path(MODULE, "outputs"), showWarnings = FALSE, recursive = TRUE)

out <- file.path(MODULE, "outputs")
input_dir <- file.path(MODULE, "data/longitude")

locations <- read.csv(
  file.path(input_dir, "03_species_localities.csv"),
  check.names = FALSE
)
species_summary <- read.csv(
  file.path(input_dir, "04_species_summary.csv"),
  check.names = FALSE
)
pairs <- read.csv(
  file.path(input_dir, "08_pairwise_PERMANOVA_BH.csv"),
  check.names = FALSE, na.strings = ""
)
overall <- read.csv(
  file.path(input_dir, "07_overall_PERMANOVA_BH.csv"),
  check.names = FALSE
)

# Displayed species per genus, revised list including Amblyomma variegatum.
selected <- list(
  Ixodes = c(
    "ovatus", "ricinus", "pacificus", "scapularis", "sinensis",
    "pavlovskyi", "persulcatus"
  ),
  Dermacentor = c(
    "silvarum", "nuttalli", "reticulatus", "variabilis",
    "albipictus", "steini"
  ),
  Amblyomma = c("americanum", "variegatum", "testudinarium"),
  Hyalomma = c("asiaticum", "rufipes", "marginatum"),
  Rhipicephalus = c("linnaei", "microplus", "turanicus"),
  Haemaphysalis = c("longicornis", "montgomeryi", "japonica")
)
selected <- lapply(names(selected), function(g) paste(g, selected[[g]])) |>
  setNames(names(selected))

# Shared panel builder; loaded, not a standalone entry point.
source(file.path(MODULE, "scripts/04_longitude_plot_functions.R"), encoding = "UTF-8")

panels <- lapply(seq_along(selected), function(i) {
  g <- names(selected)[i]
  sp <- rev(selected[[i]])

  z <- species_summary[match(sp, species_summary$species), ]
  labels <- paste0(
    substr(g, 1, 1), ". ", sub("^[^ ]+ ", "", sp), " (n=", z$n_localities, ")"
  )

  p <- build_plot(g)

  # Scale the marks for A4 while keeping only the brackets with BH-adjusted
  # P < 0.05.
  for (j in seq_along(p$layers)) {
    geom <- class(p$layers[[j]]$geom)[1]
    if (geom == "GeomPoint") p$layers[[j]]$aes_params$size <- .72
    if (geom == "GeomText") p$layers[[j]]$aes_params$size <- 2.2
    if (geom %in% c("GeomSegment", "GeomRibbon", "GeomRect")) {
      p$layers[[j]]$aes_params$linewidth <- .25
    }
  }

  ov <- overall[overall$genus == g, ]

  sy <- p$scales$get_scales("y")
  sy$labels <- labels

  sx <- p$scales$get_scales("x")
  sx$breaks <- c(-180, -90, 0, 90, 180)
  sx$labels <- c("180W", "90W", "0", "90E", "180E")

  p + labs(
    title = paste0(letters[i], "  ", g),
    subtitle = paste0(
      "R\u00b2 = ", sprintf("%.3f", ov$R2), "; adjusted P = ",
      if (ov$p_BH < .001) {
        formatC(ov$p_BH, format = "e", digits = 1)
      } else {
        sprintf("%.4f", ov$p_BH)
      }
    ),
    x = "Longitude (\u00b0)",
    caption = NULL
  ) +
    theme(
      axis.text.y = element_text(size = 6.5, face = "italic", colour = "black"),
      axis.text.x = element_text(size = 6, colour = "black"),
      axis.title.x = element_text(size = 7, hjust = .32),
      plot.title = element_text(size = 9, face = "bold.italic"),
      plot.subtitle = element_text(size = 6.4),
      plot.margin = margin(4, 3, 3, 3),
      axis.ticks.length = grid::unit(1, "mm")
    )
})

combined <- wrap_plots(panels, ncol = 2, byrow = TRUE, heights = c(1.5, 1, 1)) +
  plot_annotation(
    caption = paste0(
      "Points: distinct sampled localities. Shared-site permutations; ",
      "BH correction across ", sum(is.finite(pairs$p_BH)), " comparisons.\n",
      "Significant comparisons only: * < 0.05; ** < 0.01; *** < 0.001."
    ),
    theme = theme(
      plot.caption = element_text(size = 6.3, hjust = 0, colour = "black"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

ggsave(
  file.path(out, "14_six_rainclouds.pdf"), combined,
  width = 210, height = 297, units = "mm", device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(out, "14_six_rainclouds_preview.png"), combined,
  width = 210, height = 297, units = "mm", dpi = 180, bg = "white"
)

message("Written one A4 page: 210 x 297 mm.")
