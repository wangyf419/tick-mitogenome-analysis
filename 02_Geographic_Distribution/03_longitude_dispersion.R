# ------------------------------------------------------------------------------
# 03_longitude_dispersion.R
#
# Letter-page panel of longitude dispersion (distance to the species spatial
# median) for the six focal genera. Distances and PERMDISP statistics are read


suppressPackageStartupMessages({
  library(grid)
})

MODULE <- normalizePath(getwd(), winslash = "/")
stopifnot(dir.exists(file.path(MODULE, "scripts")))
dir.create(file.path(MODULE, "outputs"), showWarnings = FALSE)

lon_summary <- read.csv(
  file.path(MODULE, "data/longitude/07_overall_PERMANOVA_BH.csv")
)
distances <- read.csv(
  file.path(MODULE, "data/longitude/distances_to_spatial_median.csv")
)
disp <- split(distances, distances$genus)

lon_order <- c(
  "Ixodes", "Dermacentor", "Amblyomma",
  "Hyalomma", "Rhipicephalus", "Haemaphysalis"
)

species_order <- list(
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
species_order <- lapply(names(species_order), function(g) {
  paste(g, species_order[[g]])
}) |>
  setNames(names(species_order))

species_palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#AA4499", "#CCBB44", "#56B4E9", "#332288"
)

# ------------------------------------------------------------------------------
# Drawing primitives
#
# Panel viewports use native units, so the helpers below work in the coordinates
# of the current panel: x to the right, y upward from 0 to the panel height.
# ------------------------------------------------------------------------------

txt <- function(s, x, y, size = 8, font = 1, col = "#222222", just = "left",
                rot = 0) {
  grid.text(
    s,
    x = unit(x, "native"), y = unit(y, "native"),
    just = just, rot = rot,
    gp = gpar(fontfamily = "Arial", fontsize = size, fontface = font, col = col)
  )
}

ln <- function(x, y, col = "#333333", width = .65) {
  grid.lines(
    unit(x, "native"), unit(y, "native"),
    gp = gpar(col = col, lwd = width)
  )
}

rect <- function(x, y, w, h, fill = "white", border = NA, width = .65) {
  grid.rect(
    x = unit(x, "native"), y = unit(y, "native"),
    width = unit(w, "pt"), height = unit(h, "pt"),
    just = c("left", "top"),
    gp = gpar(fill = fill, col = border, lwd = width)
  )
}

pt <- function(x, y, col, size = 1.5) {
  grid.points(
    unit(x, "native"), unit(y, "native"),
    pch = 16,
    size = unit(size, "pt"), gp = gpar(col = col)
  )
}

with_panel <- function(x, y, w, h, fun) {
  pushViewport(viewport(
    x = unit(x, "native"), y = unit(y, "native"),
    width = unit(w, "pt"), height = unit(h, "pt"),
    just = c("left", "top"), yscale = c(h, 0), xscale = c(0, w)
  ))
  fun()
  popViewport()
}

# ------------------------------------------------------------------------------
# Box panel
#
# Whiskers extend to 1.5 x IQR, as in the original display.
# ------------------------------------------------------------------------------

draw_box <- function(vals, y, xmap, col, height = 8) {
  qs <- quantile(vals, c(.25, .5, .75), names = FALSE)
  span <- qs[3] - qs[1]
  bs <- c(
    min(vals[vals >= qs[1] - 1.5 * span]),
    qs,
    max(vals[vals <= qs[3] + 1.5 * span])
  )

  ln(xmap(bs[c(1, 5)]), c(y, y), col = col)
  ln(rep(xmap(bs[1]), 2), y + c(-2, 2), col = col)
  ln(rep(xmap(bs[5]), 2), y + c(-2, 2), col = col)
  rect(
    xmap(bs[2]), y - height / 2, xmap(bs[4]) - xmap(bs[2]), height,
    fill = adjustcolor(col, .4), border = col
  )
  ln(rep(xmap(bs[3]), 2), y + c(-height / 2, height / 2), col = col, width = .9)
}

# ------------------------------------------------------------------------------
# One dispersion panel
# ------------------------------------------------------------------------------

disp_panel <- function(g, tag, tall = FALSE) {
  panel_data <- disp[[g]]
  ref <- lon_summary[lon_summary$genus == g, ]
  species <- species_order[[g]]
  stopifnot(setequal(species, unique(panel_data$species)))

  spcols <- setNames(
    rev(species_palette[seq_along(species)]), species
  )

  left <- 90
  right <- 239
  top <- 43
  bottom <- if (tall) 178 else 123
  baseline <- if (tall) 193 else 138

  txt(paste0("(", tag, ")"), .5, 5, 10)
  txt(g, left, 5, 8, font = 3)
  txt(
    bquote(
      italic(F) * " = " * .(sprintf("%.3f", ref$dispersion_F)) *
        "; adjusted " * italic(p) * " = " *
        .(format(ref$dispersion_p_BH, scientific = FALSE, trim = TRUE))
    ),
    left, 20, 8
  )

  # A fixed 0-2.1 scale keeps the dispersion comparable across the six panels.
  xmax <- 2.1
  xmap <- function(x) left + x / xmax * (right - left)

  ys <- seq(top + 5, bottom - 5, length.out = length(species))
  set.seed(7100L + match(g, lon_order))

  for (i in seq_along(species)) {
    s <- species[i]
    v <- panel_data$distance[panel_data$species == s]
    y <- ys[i]

    label <- gsub("_", " ", s)
    # Only the genus is abbreviated; the panel header gives the full name.
    label <- sub(paste0("^", g, " "), paste0(substr(g, 1, 1), ". "), label)

    # Single-line species and sample-size label, as in the original figure.
    txt(
      bquote(italic(.(label)) * " (" * italic(n) * "=" * .(length(v)) * ")"),
      left - 4, y, 8,
      just = "right"
    )

    # Half-violin, small box and a separate point cloud, as in the original
    # Figure S8. The density is mirrored at zero because these distances are
    # non-negative, not circular angles.
    bw <- max(stats::bw.nrd0(v), .015)
    dx <- seq(0, 2.1, length.out = 256)
    dens <- rowMeans(
      dnorm(outer(dx, v, "-"), sd = bw) + dnorm(outer(dx, -v, "-"), sd = bw)
    )
    grid.polygon(
      unit(c(xmap(dx), rev(xmap(dx))), "native"),
      unit(c(y - 6 - 5.5 * dens / max(dens), rep(y - 6, length(dx))), "native"),
      gp = gpar(
        fill = adjustcolor(spcols[s], .30), col = spcols[s], lwd = .65
      )
    )

    draw_box(v, y - 1, xmap, spcols[s], height = 4)
    pt(
      xmap(v), y + 4 + runif(length(v), -1.7, 1.7),
      adjustcolor(spcols[s], .72),
      size = 2.5
    )
  }

  # Axis and ticks.
  ln(rep(left, 2), c(top - 9, baseline))
  ln(c(left, right), rep(baseline, 2))
  for (t in c(0, .5, 1, 1.5, 2)) {
    ln(rep(xmap(t), 2), baseline + c(0, 3))
    txt(format(t, trim = TRUE), xmap(t), baseline + 11, 8, just = "centre")
  }
  txt(
    "Distance to species spatial median",
    mean(c(left, right)), baseline + 25, 8,
    just = "centre"
  )
}

# ------------------------------------------------------------------------------
# Page layout and device
#
# Three rows: one full-width panel followed by two rows of two panels.
# ------------------------------------------------------------------------------

lon_art <- function() {
  for (i in seq_along(lon_order)) {
    row <- floor((i - 1) / 2)
    yy <- c(0, 225, 400)[row + 1]
    hh <- c(225, 175, 175)[row + 1]
    with_panel(
      ((i - 1) %% 2) * 264.1182, yy, 246.1182, hh,
      function() disp_panel(lon_order[i], LETTERS[i], tall = row == 0)
    )
  }
}

draw_device <- function(path, art, height, full = TRUE, png = FALSE) {
  w <- if (full) 612 else 510.2363858
  h <- if (full) 792 else height

  if (png) {
    ragg::agg_png(
      path,
      width = w / 72, height = h / 72,
      units = "in", res = 300, background = "white"
    )
  } else {
    grDevices::cairo_pdf(
      path,
      width = w / 72, height = h / 72,
      family = "Arial", onefile = FALSE
    )
  }

  grid.newpage()
  if (full) {
    pushViewport(viewport(
      x = unit(51.0234, "pt"), y = unit(792 - 51.0234, "pt"),
      width = unit(510.2363858, "pt"), height = unit(height, "pt"),
      just = c("left", "top"),
      xscale = c(0, 510.2363858), yscale = c(height, 0)
    ))
  } else {
    pushViewport(viewport(xscale = c(0, w), yscale = c(h, 0)))
  }
  art()
  popViewport()
  dev.off()
}

draw_device(
  file.path(MODULE, "outputs/Longitude_dispersion.pdf"), lon_art, 575
)
