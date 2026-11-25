# ------------------------------------------------------------------------------
# 02_AMOVA_heatmaps.R
#
# Three single-page heatmaps of pairwise PhiST for the original geographic
# hierarchy, annotated with pairwise BH significance, retaining the global
# Primary and All results. Completed AMOVA tables are read; nothing is tested
# again here.
#
# Run from the module directory:
#   Rscript scripts/02_AMOVA_heatmaps.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(grid)
})

SCRIPT_ARGUMENTS <- commandArgs(FALSE)
SCRIPT_PATH <- sub(
  "^--file=", "", SCRIPT_ARGUMENTS[grepl("^--file=", SCRIPT_ARGUMENTS)][1]
)
ROOT <- dirname(dirname(normalizePath(SCRIPT_PATH, winslash = "/")))
OUTPUT_DIRECTORY <- Sys.getenv("BATCH_FIGURE_STAGE", file.path(ROOT, "outputs"))

PAGE_LEFT <- 51.0234
PAGE_WIDTH <- 510.2364

# Low values are dark and high values light. Negative estimates are kept in the
# source table but the display colour scale starts at zero.
COLOR_PALETTE <- colorRampPalette(c("#174A70", "#6B9BB9", "#F1F5F7"))(101)

read_table <- function(filename) {
  read.csv(
    file.path(ROOT, "tables", filename),
    fileEncoding = "UTF-8-BOM", check.names = FALSE
  )
}

read_inputs <- function() {
  list(
    overall = read_table("08_AMOVA_overall.csv"),
    pairs = read_table("09_AMOVA_pairs.csv"),
    counts = read_table("03_location_counts.csv"),
    summary = read_table("01_input_summary.csv")
  )
}

# ------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------

significance_symbol <- function(value) {
  if (!is.finite(value)) {
    return("NA")
  }
  if (value < .001) {
    return("***")
  }
  if (value < .01) {
    return("**")
  }
  if (value < .05) {
    return("*")
  }
  "ns"
}

location_name <- function(location) {
  # Province and state names go straight onto the axis, and national records
  # keep the country name, so the reader does not need an abbreviation table.
  sub("^.*:", "", location)
}

heatmap_color <- function(value) {
  if (!is.finite(value)) {
    return("#E6E6E6")
  }
  COLOR_PALETTE[1 + round(max(0, min(1, value)) * 100)]
}

text_label <- function(label, x, y, size = 8, just = "left", font = 1,
                       rotation = 0, color = "#111111") {
  grid.text(
    label, unit(x, "native"), unit(y, "native"),
    just = just, rot = rotation,
    gp = gpar(fontfamily = "Arial", fontsize = size, fontface = font, col = color)
  )
}

rectangle <- function(x, y, width, height, color, border = "white") {
  grid.rect(
    unit(x, "native"), unit(y, "native"),
    unit(width, "bigpts"), unit(height, "bigpts"),
    just = c("left", "top"),
    gp = gpar(fill = color, col = border, lwd = .35)
  )
}

panel <- function(x, y, width, height, draw) {
  pushViewport(viewport(
    x = unit(x, "native"), y = unit(y, "native"),
    width = unit(width, "bigpts"), height = unit(height, "bigpts"),
    just = c("left", "top"),
    xscale = c(0, width), yscale = c(height, 0)
  ))
  draw()
  popViewport()
}

# ------------------------------------------------------------------------------
# Heatmap body
# ------------------------------------------------------------------------------

draw_cells <- function(locations, pairs, left, top, cell) {
  for (row in seq_along(locations)) {
    for (column in seq_len(row)) {
      x <- left + (column - 1) * cell
      y <- top + (row - 1) * cell

      if (row == column) {
        rectangle(x, y, cell, cell, "#E9ECEE")
        next
      }

      pair <- pairs[
        (pairs$location1 == locations[row] & pairs$location2 == locations[column]) |
          (pairs$location2 == locations[row] & pairs$location1 == locations[column]),
      ]
      stopifnot(nrow(pair) == 1)

      rectangle(x, y, cell, cell, heatmap_color(pair$PhiST))

      if (!is.finite(pair$PhiST)) {
        text_label("NA", x + cell / 2, y + cell / 2, 5.5, "centre")
      } else {
        color <- if (pair$PhiST < .4) "white" else "#111111"
        text_label(
          significance_symbol(pair$p_BH_within_species),
          x + cell / 2, y + cell / 2,
          min(8, max(5, cell * .55)), "centre",
          color = color
        )
      }
    }
  }
}

draw_axes <- function(locations, left, top, cell, font_size) {
  for (index in seq_along(locations)) {
    label <- location_name(locations[index])
    text_label(label, left - 4, top + (index - .5) * cell, font_size, "right")
    text_label(
      label, left + (index - .5) * cell, top + length(locations) * cell + 4,
      font_size, "right",
      rotation = 90
    )
  }
}

# Global estimates are printed once above (large panels) or inside (small panels).
draw_overall <- function(rows, x, y, font_size) {
  for (index in seq_len(nrow(rows))) {
    row <- rows[index, ]
    label <- if (row$analysis == "primary") "Primary" else "All"
    statistic <- if (row$status == "estimated") {
      sprintf("%.3f", row$PhiST)
    } else {
      "NA"
    }
    label_expr <- bquote(
      .(label) * " (n = " * .(row$n) * "): " * plain("\u03a6")[plain(ST)] * " = " *
        .(statistic) * .(significance_symbol(row$p_value))
    )
    text_label(label_expr, x, y + (index - 1) * 12, font_size)
  }
}

draw_heatmap <- function(inputs, species, tag, width, height, kind) {
  locations <- sort(
    inputs$counts$location[
      inputs$counts$species == species & inputs$counts$n >= 2
    ]
  )
  pairs <- inputs$pairs[inputs$pairs$species == species, ]

  left <- if (kind == "large") 90 else if (kind == "medium") 77 else 58
  text_label(paste0("(", tag, ")"), left - 23, 5, 10)
  text_label(sub("^([^_])[^_]+_", "\\1. ", species), left, 5, 8, font = 3)

  top <- if (kind == "large") 54 else if (kind == "medium") 48 else 52
  side <- if (kind == "large") 410 else if (kind == "medium") 104 else 45
  cell <- side / length(locations)
  font_size <- if (kind == "large") 6.5 else if (kind == "medium") 6.5 else 5.8

  draw_cells(locations, pairs, left, top, cell)
  draw_axes(locations, left, top, cell, font_size)

  overall <- inputs$overall[inputs$overall$species == species, ]
  overall <- overall[match(c("primary", "all_locations"), overall$analysis), ]

  if (kind == "large") {
    draw_overall(overall, 300, top + 5, 8)
  } else {
    draw_overall(overall, 0, 22, if (kind == "small") 7.2 else 8)
  }
}

draw_legend <- function() {
  text_label(expression(plain("\u03a6")[plain(ST)]), PAGE_LEFT, 703, 8)
  for (index in 0:99) {
    rectangle(89 + index * 1.1, 699, 1.2, 6, COLOR_PALETTE[index + 1], border = NA)
  }
  text_label("0", 84, 703, 7, "right")
  text_label("1", 204, 703, 7)
  text_label(
    "Cell stars: BH q; overall stars: P. * < 0.05; ** < 0.01; *** < 0.001; ns: >= 0.05.",
    PAGE_LEFT, 720, 7
  )
}

# ------------------------------------------------------------------------------
# One page per size class
# ------------------------------------------------------------------------------

draw_page <- function(inputs, species, part) {
  filename <- paste0("Figure_S_star_AMOVA_", part, ".pdf")

  cairo_pdf(
    file.path(OUTPUT_DIRECTORY, filename),
    width = 8.5, height = 11, family = "Arial", pointsize = 8
  )
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 612), yscale = c(792, 0)))

  rows <- list()
  for (index in seq_along(species)) {
    name <- species[index]
    tag <- LETTERS[index]

    # Small panels: 3 x 2 grid.
    if (part == 1) {
      panel(
        PAGE_LEFT + ((index - 1) %% 3) * 175.1182,
        54 + floor((index - 1) / 3) * 158, 160, 149,
        function() draw_heatmap(inputs, name, tag, 160, 149, "small")
      )
    }

    # Medium panels: 2 columns.
    if (part == 2) {
      panel(
        PAGE_LEFT + ((index - 1) %% 2) * 264.1182,
        54 + floor((index - 1) / 2) * 209, 246, 199,
        function() draw_heatmap(inputs, name, tag, 246, 199, "medium")
      )
    }

    # One large panel on its own page.
    if (part == 3) {
      panel(
        PAGE_LEFT, 54, PAGE_WIDTH, 630,
        function() draw_heatmap(inputs, name, tag, PAGE_WIDTH, 630, "large")
      )
    }

    rows[[index]] <- data.frame(species = name, file = filename, panel = tag)
  }

  draw_legend()
  popViewport()
  dev.off()

  cat("[write] ", filename, "\n", sep = "")

  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------
# Main
#
# Species are grouped by the number of locations with at least two samples:
# <= 6 locations on the small page, 7-16 on the medium page, > 16 on the large.
# ------------------------------------------------------------------------------

main <- function() {
  dir.create(OUTPUT_DIRECTORY, recursive = TRUE, showWarnings = FALSE)

  inputs <- read_inputs()

  eligible <- vapply(
    inputs$summary$species,
    function(species) {
      sum(inputs$counts$species == species & inputs$counts$n >= 2)
    },
    0L
  )

  parts <- list(
    inputs$summary$species[eligible <= 6],
    inputs$summary$species[eligible > 6 & eligible <= 16],
    inputs$summary$species[eligible > 16]
  )

  rows <- lapply(
    seq_along(parts),
    function(part) draw_page(inputs, parts[[part]], part)
  )

  write.csv(
    do.call(rbind, rows),
    file.path(OUTPUT_DIRECTORY, "16_AMOVA_figure_panels.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
}

if (sys.nframe() == 0L) main()
