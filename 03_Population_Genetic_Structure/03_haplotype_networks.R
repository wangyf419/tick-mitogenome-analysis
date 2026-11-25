# ------------------------------------------------------------------------------
# 03_haplotype_networks.R
#
# Median-joining haplotype networks, one page per species. Completed PopART
# node/edge tables and sample memberships are read; only the graphical layout
# is recomputed, with fixed seeds. No median-joining inference is run and only
# species with existing network inputs are drawn.
#
# Run from the module directory:
#   Rscript scripts/03_haplotype_networks.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(grid)
  library(igraph)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

SCRIPT_ARGUMENTS <- commandArgs(FALSE)
SCRIPT_PATH <- sub(
  "^--file=", "", SCRIPT_ARGUMENTS[grepl("^--file=", SCRIPT_ARGUMENTS)][1]
)
ROOT <- dirname(dirname(normalizePath(SCRIPT_PATH, winslash = "/")))

read_table <- function(path) {
  read.csv(path, fileEncoding = "UTF-8-BOM")
}

palette <- read_table(file.path(ROOT, "tables", "04_clade_palette.csv"))
summary <- read_table(file.path(ROOT, "tables", "01_input_summary.csv"))

outdir <- file.path(
  Sys.getenv("BATCH_FIGURE_STAGE", unset = file.path(ROOT, "outputs")),
  "networks"
)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Drawing primitives in page coordinates (points, origin at the top left)
# ------------------------------------------------------------------------------

txt <- function(z, x, y, size = 8, just = "left", font = 1, col = "black") {
  grid.text(
    z, unit(x, "native"), unit(y, "native"),
    just = just,
    gp = gpar(fontfamily = "Arial", fontsize = size, fontface = font, col = col)
  )
}

ln <- function(x, y, col = "#90999E", lwd = .5) {
  grid.lines(
    unit(x, "native"), unit(y, "native"),
    gp = gpar(col = col, lwd = lwd)
  )
}

pie <- function(x, y, r, v, cols) {
  start <- pi / 2
  for (i in seq_along(v)) {
    if (v[i] > 0) {
      end <- start + 2 * pi * v[i] / sum(v)
      a <- seq(start, end, length.out = max(5, ceiling((end - start) * 20)))
      grid.polygon(
        unit(c(x, x + r * cos(a)), "native"),
        unit(c(y, y + r * sin(a)), "native"),
        gp = gpar(fill = cols[i], col = NA)
      )
      start <- end
    }
  }
  grid.circle(
    unit(x, "native"), unit(y, "native"), unit(r, "bigpts"),
    gp = gpar(fill = NA, col = "#50595E", lwd = .45)
  )
}

# ------------------------------------------------------------------------------
# One network page per species
# ------------------------------------------------------------------------------

stats <- list()

for (si in seq_len(nrow(summary))) {
  sp <- summary$species[si]
  species_dir <- file.path(ROOT, "species", sp)
  summary_file <- file.path(species_dir, "data", "network_summary.json")

  if (!file.exists(summary_file)) {
    cat(sp, "network not available; inputs retained\n")
    next
  }
  if (nzchar(Sys.getenv("BATCH_NETWORK_SPECIES")) &&
    sp != Sys.getenv("BATCH_NETWORK_SPECIES")) {
    next
  }

  nodes <- read_table(file.path(species_dir, "data", "popart_MJ_nodes.csv"))
  edges <- read_table(file.path(species_dir, "data", "popart_MJ_edges.csv"))
  meta <- read_table(file.path(species_dir, "tables", "samples.csv"))

  # Sample counts per haplotype; inferred median nodes stay at zero.
  freq <- table(meta$haplotype)
  nodes$n <- as.integer(freq[match(nodes$node, names(freq))])
  nodes$n[is.na(nodes$n)] <- 0L

  g <- graph_from_data_frame(
    edges,
    directed = FALSE, vertices = data.frame(name = nodes$node)
  )
  stopifnot(components(g)$no == 1, vcount(g) == nrow(nodes))
  set.seed(20261109L + si)

  # ---------------------------------------------------------------------------
  # Layout
  #
  # Fruchterman-Reingold layout on the complete network, then a light
  # spreading pass so that symbols of different nodes do not touch.
  # ---------------------------------------------------------------------------

  xy <- layout_with_fr(
    g,
    weights = 1 / sqrt(edges$mutations), niter = 5000, grid = "nogrid"
  )
  if (diff(range(xy[, 1])) < diff(range(xy[, 2]))) {
    xy <- xy[, 2:1]
  }

  scale_factor <- min(
    470 / max(1, diff(range(xy[, 1]))),
    525 / max(1, diff(range(xy[, 2])))
  )
  xm <- (xy[, 1] - mean(range(xy[, 1]))) * scale_factor + 306
  ym <- (xy[, 2] - mean(range(xy[, 2]))) * scale_factor + 360

  radius_coefficient <- min(2, 10 / sqrt(max(nodes$n)))
  radius <- ifelse(nodes$n > 0, radius_coefficient * sqrt(nodes$n), .65)

  minimum <- outer(radius, radius, "+") + .5
  diag(minimum) <- 0

  for (iteration in seq_len(350)) {
    dx <- outer(xm, xm, "-")
    dy <- outer(ym, ym, "-")
    distance <- sqrt(dx^2 + dy^2)
    overlap <- pmax(minimum - distance, 0)
    diag(overlap) <- 0
    if (max(overlap) < .1) break

    denom <- pmax(distance, 1e-8)
    xm <- xm + .52 * rowSums(overlap * dx / denom)
    ym <- ym + .52 * rowSums(overlap * dy / denom)

    xm <- pmax(62 + radius, pmin(550 - radius, xm))
    ym <- pmax(85 + radius, pmin(640 - radius, ym))
  }

  nodes$x_pdf_pt <- xm
  nodes$y_pdf_pt <- ym
  nodes$radius_pdf_pt <- radius

  write.csv(
    nodes, file.path(outdir, paste0(sp, "_display_nodes.csv")),
    row.names = FALSE, fileEncoding = "UTF-8"
  )

  # ---------------------------------------------------------------------------
  # Page set-up
  # ---------------------------------------------------------------------------

  # Network inputs are read-only; the layout is used for rendering only.
  pending <- summary$tree_status[si] != "completed"
  name <- paste0(
    sp, if (pending) "_MJ_tree_pending.pdf" else "_median_joining.pdf"
  )

  cairo_pdf(
    file.path(outdir, name),
    width = 8.5, height = 11, family = "Arial", pointsize = 8
  )

  # ---------------------------------------------------------------------------
  # Edges
  #
  # Ticks perpendicular to the edge show single mutations; substitution counts
  # above ten are written as labels afterwards.
  # ---------------------------------------------------------------------------

  grid.newpage()
  pushViewport(viewport(xscale = c(0, 612), yscale = c(792, 0)))
  txt(gsub("_", " ", sp), 51.0234, 56, 8, font = 3)

  np <- setNames(seq_len(nrow(nodes)), nodes$node)
  labelboxes <- list()

  for (i in seq_len(nrow(edges))) {
    a <- np[edges$from[i]]
    b <- np[edges$to[i]]
    ln(xm[c(a, b)], ym[c(a, b)])

    vv <- c(xm[b] - xm[a], ym[b] - ym[a])
    len <- sqrt(sum(vv^2))
    perp <- c(-vv[2], vv[1]) / max(len, 1e-9)
    mut <- edges$mutations[i]

    if (mut >= 1 && mut <= 10) {
      for (t in if (mut == 1) .5 else seq(.35, .65, length.out = mut)) {
        p0 <- c(xm[a], ym[a]) + t * vv
        ln(
          c(p0[1] - perp[1] * 1.2, p0[1] + perp[1] * 1.2),
          c(p0[2] - perp[2] * 1.2, p0[2] + perp[2] * 1.2),
          "#707D84", .45
        )
      }
    }
  }

  # Label long edges after all lines are drawn, largest counts first.
  # Only colliding annotations are suppressed; every edge and its exact weight
  # is retained.
  annotated <- rep(FALSE, nrow(edges))

  for (i in order(edges$mutations, decreasing = TRUE)) {
    if (edges$mutations[i] > 10) {
      a <- np[edges$from[i]]
      b <- np[edges$to[i]]
      mut <- edges$mutations[i]

      vv <- c(xm[b] - xm[a], ym[b] - ym[a])
      perp <- c(-vv[2], vv[1]) / max(sqrt(sum(vv^2)), 1e-9)

      choices <- expand.grid(
        t = c(.5, .35, .65, .25, .75),
        offset = c(5, -5, 10, -10, 15, -15, 22, -22, 30, -30)
      )
      bestp <- NULL
      bestanchor <- NULL

      width <- convertWidth(
        stringWidth(as.character(mut)), "bigpts",
        valueOnly = TRUE
      ) + 3

      for (k in seq_len(nrow(choices))) {
        anchor <- c(xm[a], ym[a]) + choices$t[k] * vv
        p0 <- anchor + choices$offset[k] * perp

        hits <- any(abs(xm - p0[1]) < radius + width / 2 + 1 &
          abs(ym - p0[2]) < radius + 5)
        overlaps <- any(vapply(
          labelboxes,
          function(z) abs(z[1] - p0[1]) < width + 2 && abs(z[2] - p0[2]) < 11,
          TRUE
        ))

        if (!hits && !overlaps &&
          p0[1] > 62 + width / 2 && p0[1] < 550 - width / 2 &&
          p0[2] > 85 && p0[2] < 638) {
          bestp <- p0
          bestanchor <- anchor
          break
        }
      }

      if (!is.null(bestp)) {
        if (sqrt(sum((bestp - bestanchor)^2)) > 8) {
          ln(
            c(bestanchor[1], bestp[1]), c(bestanchor[2], bestp[2]),
            "#A4ADB2", .4
          )
        }
        grid.rect(
          unit(bestp[1], "native"), unit(bestp[2], "native"),
          unit(width, "bigpts"), unit(9, "bigpts"),
          gp = gpar(fill = "white", col = NA)
        )
        txt(
          as.character(mut), bestp[1], bestp[2], 8,
          just = "centre",
          col = "#52616A"
        )
        labelboxes[[length(labelboxes) + 1]] <- bestp
        annotated[i] <- TRUE
      }
    }
  }

  write.csv(
    data.frame(edges, numeric_label_drawn = annotated),
    file.path(outdir, paste0(sp, "_edge_annotations.csv")),
    row.names = FALSE, fileEncoding = "UTF-8"
  )

  # ---------------------------------------------------------------------------
  # Nodes
  #
  # Inferred median vectors are small grey dots; observed haplotypes are pies
  # whose slices give the clade composition and whose radius gives the number
  # of samples.
  # ---------------------------------------------------------------------------

  for (i in which(nodes$n == 0)) {
    grid.circle(
      unit(xm[i], "native"), unit(ym[i], "native"), unit(.65, "bigpts"),
      gp = gpar(fill = "#747E83", col = NA)
    )
  }

  pp <- palette[palette$species == sp, ]
  co <- pp$color
  if (pending) co <- summary$species_hue[si]

  for (i in order(nodes$n, decreasing = TRUE)) {
    if (nodes$n[i] > 0) {
      v <- table(factor(
        meta$clade[meta$haplotype == nodes$node[i]],
        levels = pp$clade
      ))
      pie(xm[i], ym[i], radius[i], as.numeric(v), co)
    }
  }

  # All haplotypes are retained; labels are reserved for frequent nodes to
  # avoid clutter.
  labelidx <- which(nodes$n >= 5)
  if (sum(nodes$n > 0) <= 35) labelidx <- which(nodes$n > 0)

  for (i in labelidx) {
    lab <- nodes$node[i]
    width <- convertWidth(stringWidth(lab), "bigpts", valueOnly = TRUE)
    best <- Inf
    loc <- NULL

    for (ang in seq(0, 2 * pi, length.out = 17)[-17]) {
      p0 <- c(xm[i], ym[i]) + c(cos(ang), sin(ang)) * (radius[i] + 8)
      score <- sum(abs(xm - p0[1]) < radius + width / 2 &
        abs(ym - p0[2]) < radius + 4) * 100 +
        sum(vapply(
          labelboxes,
          function(z) abs(z[1] - p0[1]) < width && abs(z[2] - p0[2]) < 9,
          TRUE
        )) * 100
      if (score < best) {
        best <- score
        loc <- p0
      }
      if (score == 0) break
    }

    if (!is.null(loc) && best == 0 &&
      loc[1] > 70 && loc[1] < 540 && loc[2] > 83 && loc[2] < 640) {
      txt(lab, loc[1], loc[2], 8, just = "centre")
      labelboxes[[length(labelboxes) + 1]] <- loc
    }
  }

  # ---------------------------------------------------------------------------
  # Legend
  # ---------------------------------------------------------------------------

  y <- 663
  for (i in seq_len(nrow(pp))) {
    x <- 51.0234 + (i - 1) * 80
    label <- if (pending) {
      "Tree pending"
    } else if (pp$clade[i] == "Unresolved") {
      "Unassigned"
    } else {
      sub("Clade_0?", "Group ", pp$clade[i])
    }
    pie(x + 4, y, 3, 1, co[i])
    txt(label, x + 12, y, 8)
  }

  txt("Samples", 51.0234, 687, 8)

  size_breaks <- if (max(nodes$n) <= 5) {
    sort(unique(c(1, max(nodes$n))))
  } else {
    c(1, 5, max(nodes$n))
  }

  for (i in seq_along(size_breaks)) {
    n <- size_breaks[i]
    x <- 110 + (i - 1) * 65
    pie(x, 687, radius_coefficient * sqrt(n), 1, "#B8C1C7")
    txt(n, x + radius_coefficient * sqrt(n) + 4, 687, 8)
  }

  grid.circle(
    unit(328, "native"), unit(687, "native"), unit(.8, "bigpts"),
    gp = gpar(fill = "#747E83", col = NA)
  )
  txt("Inferred median", 337, 687, 8)

  txt(
    bquote(
      italic(n) * " = " * .(nrow(meta)) * "; " * .(sum(nodes$n > 0)) *
        " observed haplotypes; " * .(sum(nodes$n == 0)) * " inferred medians."
    ),
    51.0234, 709, 8
  )

  popViewport()
  dev.off()
  cat(sp, "network plotted\n")

  stats[[length(stats) + 1]] <- data.frame(
    species = sp,
    n = nrow(meta),
    observed_haplotypes = sum(nodes$n > 0),
    inferred_medians = sum(nodes$n == 0),
    edges = nrow(edges),
    tree_pending = pending
  )
}

# ------------------------------------------------------------------------------
# Summary table
#
# When only one species is redrawn, the existing summary is kept for the others.
# ------------------------------------------------------------------------------

newstats <- do.call(rbind, stats)
dest <- file.path(outdir, "17_network_figure_summary.csv")

if (file.exists(dest) && nzchar(Sys.getenv("BATCH_NETWORK_SPECIES"))) {
  old <- read_table(dest)
  newstats <- rbind(old[!old$species %in% newstats$species, ], newstats)
}

newstats <- newstats[match(summary$species, newstats$species), ]
newstats <- newstats[!is.na(newstats$species), ]

write.csv(newstats, dest, row.names = FALSE, fileEncoding = "UTF-8")
