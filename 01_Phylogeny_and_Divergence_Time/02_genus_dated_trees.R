# ------------------------------------------------------------------------------
# 02_genus_dated_trees.R
#
# One-page A4 composite of the eight genus-level dated trees, with the numeric
# 95% HPD interval written at every internal node.
#
# Run from the module directory:
#   Rscript scripts/02_genus_dated_trees.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(patchwork)
  library(deeptime)
})

MODULE <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(MODULE, "scripts"))) {
  stop("Run this script from the module directory (the folder containing data/ and scripts/).")
}
dir.create(file.path(MODULE, "outputs"), showWarnings = FALSE, recursive = TRUE)

tree_dir <- file.path(MODULE, "data/genus_trees")
out_pdf <- file.path(MODULE, "outputs/All_Genus_Trees_A4_95HPD.pdf")

# The input trees store times in units of 100 Myr.
time_multiplier <- 100

tree_files <- list.files(path = tree_dir, pattern = "\\.tre$", full.names = TRUE)

# ------------------------------------------------------------------------------
# HPD extraction
#
# Extract the 95% HPD interval of every internal node from a treeio object. The
# column name differs between the source trees, so it is located by pattern
# rather than by position.
# ------------------------------------------------------------------------------

get_hpd_data <- function(tree) {
  dat <- tree@data

  hpd_col <- grep(
    "95.*HPD|HPD.*95",
    colnames(dat),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(hpd_col) == 0) {
    warning("No 95% HPD field found in tree@data.")
    return(data.frame(
      node = integer(),
      HPD_lower = numeric(),
      HPD_upper = numeric()
    ))
  }

  hpd_col <- hpd_col[1]

  # Annotations are stored as numeric, list or character values.
  parse_hpd <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(c(NA_real_, NA_real_))
    }

    if (is.numeric(x) && length(x) >= 2) {
      return(as.numeric(x[1:2]))
    }

    txt <- paste(x, collapse = ",")
    nums <- regmatches(
      txt,
      gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", txt)
    )[[1]]
    nums <- suppressWarnings(as.numeric(nums))

    if (length(nums) >= 2) {
      return(nums[1:2])
    }

    c(NA_real_, NA_real_)
  }

  hpd_matrix <- t(vapply(dat[[hpd_col]], parse_hpd, numeric(2)))

  data.frame(
    node = dat$node,
    HPD_lower = hpd_matrix[, 1] * time_multiplier,
    HPD_upper = hpd_matrix[, 2] * time_multiplier
  )
}

# ------------------------------------------------------------------------------
# Read and draw one genus tree
# ------------------------------------------------------------------------------

plot_genus_tree <- function(file) {
  name <- tools::file_path_sans_ext(basename(file))
  # The supplied FigTree.tre contains this genus.
  if (name == "FigTree") name <- "Haemaphysalis"

  lines <- readLines(file, warn = FALSE)

  # Ixodes.tre is written as UTREE; treeio handles TREE blocks more reliably.
  lines <- sub(
    "^([[:space:]]*)UTREE([[:space:]])",
    "\\1TREE\\2",
    lines,
    ignore.case = TRUE
  )

  # Make sure the file ends with a newline before parsing it.
  lines <- c(lines, "")

  tmp_file <- tempfile(fileext = ".tre")
  writeLines(lines, tmp_file)

  tree <- tryCatch(
    treeio::read.mcmctree(tmp_file),
    error = function(e) {
      message("[skip] ", name, ": could not be read (", e$message, ")")
      NULL
    }
  )

  unlink(tmp_file)

  if (is.null(tree)) {
    return(NULL)
  }

  n_tips <- length(tree@phylo$tip.label)

  # Species names are displayed with a space instead of an underscore.
  tree@phylo$tip.label <- gsub("_", " ", tree@phylo$tip.label)

  hpd_data <- get_hpd_data(tree)

  if (nrow(hpd_data) == 0) {
    message("[skip] ", name, ": no 95% HPD interval extracted.")
    return(NULL)
  }

  # Convert the time unit of the input tree to Myr.
  tree@phylo$edge.length <- tree@phylo$edge.length * time_multiplier

  p <- ggtree(tree, linewidth = 0.4) %>% revts()

  # Match the HPD rows to the ggtree node numbers.
  idx <- match(p$data$node, hpd_data$node)

  p$data$HPD_lower <- hpd_data$HPD_lower[idx]
  p$data$HPD_upper <- hpd_data$HPD_upper[idx]

  # Label internal nodes only, e.g. 78.4-88.3.
  p$data$HPD_label <- ifelse(
    !p$data$isTip &
      !is.na(p$data$HPD_lower) &
      !is.na(p$data$HPD_upper),
    sprintf("%.1f\u2013%.1f", p$data$HPD_lower, p$data$HPD_upper),
    ""
  )

  p <- p +
    # Species names at the tips.
    geom_tiplab(
      size = 2.2,
      align = TRUE,
      linetype = "dashed",
      linesize = 0.2,
      fontface = "italic"
    ) +
    # 95% HPD interval at the internal nodes.
    geom_text2(
      aes(
        subset = !isTip & HPD_label != "",
        label = HPD_label
      ),
      size = 1.45,
      color = "#1B9E77",
      vjust = -0.5,
      hjust = 1.05,
      fontface = "bold"
    ) +
    # Two-tier geological time axis.
    coord_geo(
      dat = list("epochs", "periods"),
      abbrv = list(TRUE, FALSE),
      size = list(1.2, 1.5),
      height = list(
        unit(0.5, "lines"),
        unit(0.5, "lines")
      ),
      bord = c("white", "white"),
      neg = TRUE,
      pos = "bottom"
    ) +
    scale_x_continuous(labels = abs) +
    scale_y_continuous(
      expand = expansion(add = c(1.5, 0.5))
    ) +
    theme_tree2() +
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_text(size = 5),
      plot.margin = margin(t = 2, r = 10, b = 2, l = 2)
    ) +
    hexpand(0.35)

  message(
    "[ok] ", name,
    ": tips = ", n_tips,
    ", internal nodes with 95% HPD = ", sum(p$data$HPD_label != "")
  )

  list(
    plot = p,
    tips = n_tips,
    name = name
  )
}

# ------------------------------------------------------------------------------
# Batch processing
# ------------------------------------------------------------------------------

message("Reading all genus trees and extracting 95% HPD intervals...")

results <- lapply(tree_files, plot_genus_tree)
results <- Filter(Negate(is.null), results)

plot_list <- lapply(results, `[[`, "plot")
tip_counts <- sapply(results, `[[`, "tips")
names_list <- sapply(results, `[[`, "name")

# Panels are ordered by species number, largest first.
sort_idx <- order(tip_counts, decreasing = TRUE)

plot_list <- plot_list[sort_idx]
tip_counts <- tip_counts[sort_idx]
names_ordered <- names_list[sort_idx]

message("Panel order: ", paste(names_ordered, collapse = " -> "))

# Two panels per row; each row is as tall as its largest tree.
ncol <- 2
n_plots <- length(plot_list)
n_rows <- ceiling(n_plots / ncol)

row_heights <- numeric(n_rows)

for (i in seq_len(n_rows)) {
  idx_start <- (i - 1) * ncol + 1
  idx_end <- min(i * ncol, n_plots)

  row_heights[i] <- max(
    tip_counts[idx_start:idx_end]
  )
}

if (n_plots > 0) {
  combined_plot <- wrap_plots(
    plot_list,
    ncol = ncol
  ) +
    plot_layout(
      heights = row_heights
    ) +
    plot_annotation(
      tag_levels = "A"
    ) &
    theme(
      plot.tag = element_text(
        size = 14,
        face = "bold"
      )
    )

  ggsave(
    filename = out_pdf,
    plot = combined_plot,
    width = 8.27,
    height = 11.69,
    units = "in",
    dpi = 300
  )

  message("A4 composite of the dated trees written to: ", out_pdf)
} else {
  message("[error] None of the tree files could be read.")
}
