# ------------------------------------------------------------------------------
# longitude_plot_functions.R
#
# Shared panel builder for the within-genus longitude rainclouds.

# Small helper for P value labels.

fmt <- function(p) {
  ifelse(p < .001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p))
}

palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#AA4499", "#CCBB44", "#56B4E9", "#332288"
)

build_plot <- function(g, exact = FALSE) {
  sp <- rev(selected[[g]])
  k <- length(sp)
  cols <- setNames(palette[seq_len(k)], sp)

  dd <- locations[locations$genus == g, ]
  zz <- species_summary[species_summary$genus == g, ]

  # ---------------------------------------------------------------------------
  # Display geometry
  #
  # Periodic Gaussian density with a fixed bandwidth of 10 degrees and the
  # same smoothing for all species; longitude is wrapped at +/-360 so that
  # localities close to the antimeridian are not split.
  # ---------------------------------------------------------------------------

  dens <- boxes <- list()
  gridx <- seq(-180, 180, length.out = 721)

  for (s in sp) {
    v <- dd$longitude[dd$species == s]
    y <- match(s, sp)

    if (length(v) >= 3) {
      delta <- outer(gridx, v, "-")
      den <- rowMeans(
        dnorm(delta, sd = 10) + dnorm(delta - 360, sd = 10) +
          dnorm(delta + 360, sd = 10)
      )
      dens[[s]] <- data.frame(
        species = s, x = gridx, base = y + .24,
        top = y + .24 + .27 * den / max(den)
      )
    }

    if (length(v) >= 2) {
      b <- boxplot.stats(v)$stats
      boxes[[s]] <- data.frame(
        species = s, y = y + .055,
        lo = b[1], q1 = b[2], median = b[3], q3 = b[4], hi = b[5]
      )
    }
  }

  dens <- bind_rows(dens)
  boxes <- bind_rows(boxes)

  set.seed(1109)
  dd$y <- match(dd$species, sp) - .16 + runif(nrow(dd), -.07, .07)

  # Only significant pairwise comparisons carry a bracket.
  ann <- pairs[pairs$genus == g & is.finite(pairs$p_BH) & pairs$p_BH < .05, ]
  ann$y1 <- match(ann$group1, sp) + .055
  ann$y2 <- match(ann$group2, sp) + .055
  ann <- ann[order(abs(ann$y1 - ann$y2), ann$y1), ]
  ann$x <- 199 + (seq_len(nrow(ann)) - 1) * if (k >= 6) 11 else 17
  ann$label <- if (exact) {
    ifelse(is.na(ann$p_BH), "NT", fmt(ann$p_BH))
  } else {
    ann$symbol
  }

  maxx <- if (nrow(ann) > 0) max(ann$x) + 20 else 190
  ov <- overall[overall$genus == g, ]
  lab <- paste0(sp, "  (n = ", zz$n_localities[match(sp, zz$species)], ")")

  p <- ggplot() +
    geom_ribbon(
      data = dens,
      aes(x = x, ymin = base, ymax = top, fill = species, colour = species),
      alpha = .30, linewidth = .5
    ) +
    geom_segment(
      data = boxes,
      aes(x = lo, xend = hi, y = y, yend = y, colour = species),
      linewidth = .6
    ) +
    geom_rect(
      data = boxes,
      aes(
        xmin = q1, xmax = q3, ymin = y - .06, ymax = y + .06,
        fill = species, colour = species
      ),
      alpha = .45, linewidth = .65
    ) +
    geom_segment(
      data = boxes,
      aes(
        x = median, xend = median, y = y - .06, yend = y + .06,
        colour = species
      ),
      linewidth = .8
    ) +
    geom_point(
      data = dd,
      aes(x = longitude, y = y, colour = species),
      size = 1.55, alpha = .72
    ) +
    geom_segment(
      data = ann, aes(x = x, xend = x, y = y1, yend = y2), linewidth = .32
    ) +
    geom_segment(
      data = ann, aes(x = x - 3, xend = x, y = y1, yend = y1), linewidth = .32
    ) +
    geom_segment(
      data = ann, aes(x = x - 3, xend = x, y = y2, yend = y2), linewidth = .32
    ) +
    geom_text(
      data = ann, aes(x = x + 3, y = (y1 + y2) / 2, label = label),
      angle = 90, size = if (exact) 2.25 else 2.9
    ) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(
      limits = c(-183, maxx), breaks = seq(-180, 180, 60),
      labels = c("180 W", "120 W", "60 W", "0", "60 E", "120 E", "180 E"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(.4, k + .65), breaks = seq_len(k), labels = lab,
      expand = c(0, 0)
    ) +
    annotate("segment", x = -180, xend = 180, y = .4, yend = .4, linewidth = .45) +
    labs(
      x = "Longitude (degrees)", y = NULL, title = g,
      subtitle = paste0(
        "Circular-coordinate permutation test: R2 = ",
        sprintf("%.3f", ov$R2), ", BH-adjusted P = ", fmt(ov$p_BH)
      ),
      caption = paste0(
        "Shared-site permutations; BH correction across ",
        sum(is.finite(pairs$p_BH)), " pairwise comparisons.\n",
        if (any(is.na(ann$p_BH))) "NT: fewer than 3 localities. " else "",
        if (any(dd$coordinate_review)) {
          "Coordinate check pending for SRR23908119.\n"
        } else {
          ""
        },
        if (exact) {
          "Labels: adjusted P values."
        } else {
          "Significant comparisons only: * < 0.05; ** < 0.01; *** < 0.001."
        }
      )
    ) +
    theme_classic(base_size = 11, base_family = "Arial") +
    theme(
      axis.text.y = element_text(face = "italic", colour = "black", size = 10),
      axis.text.x = element_text(colour = "black", size = 9),
      axis.line.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(hjust = .32),
      plot.title = element_text(face = "bold.italic"),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 8, hjust = 0),
      plot.margin = margin(12, 15, 8, 10)
    )

  # Localities whose coordinates are still under review are marked separately.
  if (any(dd$coordinate_review)) {
    p <- p + geom_point(
      data = dd[dd$coordinate_review, ], aes(x = longitude, y = y),
      inherit.aes = FALSE, shape = 4, colour = "red", size = 3, stroke = 1
    )
  }

  p
}
