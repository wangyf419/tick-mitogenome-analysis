# ------------------------------------------------------------------------------
# 01_admixture_bars.R
#
# STRUCTURE/admixture bar plots for every species folder that has prepared Q
# matrices, sample IDs and group labels. Samples are ordered within each group
# using the cluster assignment of the largest K, and that order is applied to
# all K. No mixture model is fitted here.


suppressPackageStartupMessages({
  library(pophelper)
  library(dplyr)
})

MODULE <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(MODULE, "scripts"))) {
  stop("Run this script from the module directory (the folder containing data/ and scripts/).")
}
dir.create(file.path(MODULE, "outputs"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

q_dir <- file.path(MODULE, "data/admixture/Q")
group_dir <- file.path(MODULE, "data/admixture/groups")
sample_dir <- file.path(MODULE, "data/admixture/samples")
output_dir <- file.path(MODULE, "outputs/admixture")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Order samples by their cluster membership inside each group
# ------------------------------------------------------------------------------

sort_by_clusters <- function(qmat, groups) {
  df <- as.data.frame(qmat)
  df$Group <- groups
  df$Sample <- rownames(qmat)

  cluster_cols <- colnames(qmat)

  arrange_args <- lapply(
    cluster_cols,
    function(x) desc(df[[x]])
  )

  df_sorted <- df %>%
    group_by(Group) %>%
    arrange(!!!arrange_args, .by_group = TRUE) %>%
    ungroup()

  df_sorted
}

# ------------------------------------------------------------------------------
# One species folder at a time
# ------------------------------------------------------------------------------

q_folders <- list.dirs(
  q_dir,
  recursive = FALSE,
  full.names = TRUE
)

for (folder in q_folders) {
  species_name <- basename(folder)

  message("Processing: ", species_name)

  # Sample and group label files.
  sample_file <- file.path(
    sample_dir,
    paste0(species_name, "_sample.txt")
  )

  group_file <- file.path(
    group_dir,
    paste0(species_name, "_group.txt")
  )

  if (!file.exists(sample_file) | !file.exists(group_file)) {
    message("[skip] ", species_name, ": sample or group file is missing.")
    next
  }

  sample_labels <- readLines(sample_file)
  group_labels <- readLines(group_file)

  # Q matrices of this species, ordered by file name.
  q_files <- sort(
    list.files(
      folder,
      pattern = "\\.Q$",
      full.names = TRUE
    )
  )

  if (length(q_files) == 0) {
    message("[skip] ", species_name, ": no Q matrix found.")
    next
  }

  slist <- readQ(
    q_files,
    indlabfromfile = FALSE
  )

  # Attach the sample IDs and check that they match the Q matrices.
  slist <- lapply(slist, function(x) {
    stopifnot(
      nrow(x) == length(sample_labels),
      length(group_labels) == length(sample_labels),
      !anyDuplicated(sample_labels)
    )
    rownames(x) <- sample_labels
    x
  })

  # K increases from the first to the last entry.
  slist <- slist[order(vapply(slist, ncol, integer(1)))]

  # ---------------------------------------------------------------------------
  # The largest K defines the display order of all K
  # ---------------------------------------------------------------------------

  maxK_index <- which.max(
    sapply(slist, ncol)
  )

  sort_df <- sort_by_clusters(
    slist[[maxK_index]],
    group_labels
  )

  order_samples <- sort_df$Sample
  group_labels_order <- sort_df$Group

  message(
    "Group order: ",
    paste(unique(group_labels_order), collapse = " | ")
  )

  slist <- lapply(
    slist,
    function(x) {
      x[order_samples, , drop = FALSE]
    }
  )

  # ---------------------------------------------------------------------------
  # Bar plots
  # ---------------------------------------------------------------------------

  plotQ(
    qlist = slist,
    imgoutput = "join",

    # The side labels show the K value only.
    showsp = TRUE,
    splab = paste0("K=", vapply(slist, ncol, integer(1))),
    clustercol = c(
      "#8cc7c5",
      "#f0b880",
      "#f67868",
      "#83aac3",
      "#b6b0d4",
      "#4b658a"
    ),
    showindlab = FALSE,
    showlegend = TRUE,
    showdiv = TRUE,
    divcol = "white",
    divtype = 2,
    divsize = 0.5,
    grplab = data.frame(
      Group = group_labels_order
    ),
    selgrp = "Group",
    ordergrp = FALSE,
    showtitle = FALSE,
    barbordersize = 0,
    barbordercolour = "gray95",
    grplabangle = 45,
    imgtype = "pdf",
    outputfilename = species_name,
    exportpath = output_dir,
    height = 2.5,
    width = 8
  )

  message("[ok] ", species_name, ".pdf")
}
