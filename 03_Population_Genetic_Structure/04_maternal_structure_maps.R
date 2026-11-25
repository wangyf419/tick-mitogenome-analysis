# ------------------------------------------------------------------------------
# 04_maternal_structure_maps.R
#
# Maternal-group maps of pairwise PhiST. The original geography shows the
# lowest quantile subset of significant positive values, whereas the
# country-level summary shows every significant positive value; the curves
# carry no direction. Final tables are read, so no value is re-estimated.


suppressPackageStartupMessages({
  library(sf)
  library(grid)
  library(jsonlite)
})

SCRIPT_ARGUMENTS <- commandArgs(FALSE)
SCRIPT_PATH <- sub(
  "^--file=", "", SCRIPT_ARGUMENTS[grepl("^--file=", SCRIPT_ARGUMENTS)][1]
)
ROOT <- dirname(dirname(normalizePath(SCRIPT_PATH, winslash = "/")))

TABLE_DIRECTORY <- Sys.getenv(
  "BATCH_MAP_TABLE_DIRECTORY", file.path(ROOT, "tables")
)
OUTPUT_DIRECTORY <- Sys.getenv("BATCH_FIGURE_STAGE", file.path(ROOT, "outputs"))
MAP_TARGETS <- strsplit(
  Sys.getenv("BATCH_MAP_TARGETS", "Figure_3,Figure_S35,Figure_S36"),
  ",",
  fixed = TRUE
)[[1]]

# Page frame in points.
PAGE_LEFT <- 51.0234
PAGE_WIDTH <- 510.2364

# Symbol and legend geometry.
PIE_RADIUS_AT_211 <- 12
PHIST_DISPLAY_MULTIPLIER <- 100
MINIMUM_PIE_RADIUS <- 2.8
MINIMUM_PADDING_DEGREES <- 1.2
RELATIVE_PADDING <- .08
WHITE_FRACTION <- .75
LEGEND_COLUMN_WIDTH <- 52
LEGEND_BAR_HEIGHT <- 42
LEGEND_TITLE_SIZE <- 7.2
LEGEND_SPECIES_SIZE <- 5.5
LEGEND_RESERVED_HEIGHT <- 71
LEGEND_RIGHT_INSET <- 2
LEGEND_BOTTOM_INSET <- 4

# External basemap and map-review identifier.
BASEMAP_FILE <- Sys.getenv(
  "BATCH_BASEMAP_FILE",
  file.path(ROOT, "data/basemap/world_center160.shp")
)
MAP_APPROVAL_NUMBER <- Sys.getenv("BATCH_MAP_APPROVAL_NUMBER", "GS(2016)1665")

CENTRAL_MERIDIAN <- 160
EARTH_RADIUS <- 6378137
MAP_CRS <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=160 +datum=WGS84 +units=m +no_defs"

SUPPLEMENTARY_MAP_FILES <- c(
  Figure_S35 = Sys.getenv("BATCH_MAP_DERMACENTOR_FILE", "FigureS24.pdf"),
  Figure_S36 = Sys.getenv("BATCH_MAP_HAEMAPHYSALIS_FILE", "FigureS25.pdf")
)

OCEAN_COLOR <- "#E5F2F9"
LAND_COLOR <- "#F0F0EE"
SAMPLED_COLOR <- "#F5F9F2"

# ------------------------------------------------------------------------------
# Input
#
# Groups belonging to excluded species or to the older tree are not merged in.
# Two table names are kept in their archived form so that the script reads the
# supplied tables unchanged.
# ------------------------------------------------------------------------------

read_table <- function(filename) {
  read.csv(
    file.path(TABLE_DIRECTORY, filename),
    fileEncoding = "UTF-8-BOM", check.names = FALSE
  )
}

read_inputs <- function() {
  inputs <- list(
    samples = read_table("02_all_samples.csv"),
    counts = read_table("03_location_counts.csv"),
    palette = read_table("04_clade_palette.csv"),
    summary = read_table("01_input_summary.csv"),
    pairs = read_table("24_PhiST地图筛选.csv"),
    geography = readRDS(file.path(ROOT, "data", "geography.rds"))
  )

  boolean_columns <- c("selected_for_map", "significant_positive")
  for (name in boolean_columns) {
    inputs$pairs[[name]] <- as.logical(inputs$pairs[[name]])
  }

  # Country-level rows are only drawn for the two species that have them; their
  # selection is rebuilt from the stored statistic, not from the flag.
  country_rows <- inputs$pairs$analysis_level == "country" &
    inputs$pairs$species %in% c(
      "Rhipicephalus_microplus", "Haemaphysalis_longicornis"
    )
  inputs$pairs$selected_for_map[country_rows] <- with(
    inputs$pairs[country_rows, ],
    status == "estimated" & PhiST > 0 & p_BH_within_species < .05
  )

  inputs
}

# ------------------------------------------------------------------------------
# Basemap
# ------------------------------------------------------------------------------

prepare_basemap <- function(inputs) {
  countries <- st_read(BASEMAP_FILE, quiet = TRUE)

  if (!all(c("NAME", "GMI_CNTRY") %in% names(countries))) {
    stop("The basemap has no country name or country code field.")
  }

  country_lookup <- unique(
    st_drop_geometry(inputs$geography$provinces)[, c("adm0_a3", "admin")]
  )
  countries$admin <- country_lookup$admin[
    match(countries$GMI_CNTRY, country_lookup$adm0_a3)
  ]

  missing <- is.na(countries$admin)
  countries$admin[missing] <- countries$NAME[missing]

  # The original map is already cut along the 160 degree central meridian, so
  # its relative degrees are converted straight to equirectangular metres.
  # A general reprojection would cut the coastlines again and draw spurious
  # segments across the whole map.
  geometry <- st_set_crs(st_geometry(countries), NA)
  st_geometry(countries) <- st_set_crs(geometry * (EARTH_RADIUS * pi / 180), MAP_CRS)

  inputs$geography$countries_center <- st_make_valid(countries)
  inputs$geography$provinces_center <- suppressWarnings(
    st_transform(inputs$geography$provinces, MAP_CRS)
  )
  inputs$geography$dash_center <- st_transform(inputs$geography$dash, MAP_CRS)

  inputs
}

add_country_positions <- function(inputs) {
  needed <- unique(sub(
    ":.*", "",
    inputs$samples$location[
      inputs$samples$species %in% c(
        "Rhipicephalus_microplus", "Haemaphysalis_longicornis"
      )
    ]
  ))

  for (country in setdiff(needed, inputs$geography$positions$location)) {
    admin <- if (country == "USA") "United States of America" else country

    features <- inputs$geography$provinces[
      inputs$geography$provinces$admin == admin,
    ]
    geometry <- st_union(st_make_valid(st_geometry(features)))

    # An equal-area projection centred on the country keeps the centroid on the
    # main landmass.
    centre_projection <- if (country == "China") {
      "+proj=laea +lat_0=35 +lon_0=105 +datum=WGS84 +units=m"
    } else {
      "+proj=laea +lat_0=38 +lon_0=-98 +datum=WGS84 +units=m"
    }

    parts <- suppressWarnings(
      st_cast(st_transform(geometry, centre_projection), "POLYGON")
    )
    shape <- parts[which.max(as.numeric(st_area(parts)))]
    point <- suppressWarnings(st_centroid(shape))

    if (!length(st_within(point, shape)[[1]])) {
      point <- st_point_on_surface(shape)
    }

    coordinates <- st_coordinates(st_transform(point, 4326))[1, ]

    row <- data.frame(
      location = country,
      longitude = coordinates[1],
      latitude = coordinates[2],
      position_rule = "country_largest_landmass_centroid",
      source_admin = admin,
      source_name = "country total",
      source_features = paste(features$adm1_code, collapse = ";")
    )

    inputs$geography$positions <- rbind(inputs$geography$positions, row)
  }

  inputs
}

# ------------------------------------------------------------------------------
# Panel subsets
# ------------------------------------------------------------------------------

make_panel_inputs <- function(inputs, species, country_species, scope, region) {
  view <- inputs
  view$summary <- inputs$summary[inputs$summary$species %in% species, ]
  view$samples <- inputs$samples[inputs$samples$species %in% species, ]

  # Country-level panels collapse every location to its country.
  aggregate_rows <- view$samples$species %in% country_species
  view$samples$location[aggregate_rows] <- sub(
    ":.*", "", view$samples$location[aggregate_rows]
  )

  if (scope == "USA_states") {
    view$samples <- view$samples[grepl("^USA:", view$samples$location), ]
  }
  if (!is.null(region)) {
    positions <- inputs$geography$positions
    keep <- positions$location[
      positions$longitude >= 70 & positions$longitude <= 150 &
        positions$latitude >= 0 & positions$latitude <= 60
    ]
    view$samples <- view$samples[view$samples$location %in% keep, ]
  }

  view$counts <- aggregate(
    rep(1, nrow(view$samples)),
    list(species = view$samples$species, location = view$samples$location),
    sum
  )
  names(view$counts)[3] <- "n"

  view$pairs <- inputs$pairs[inputs$pairs$species %in% species, ]
  target <- ifelse(
    view$pairs$species %in% country_species,
    "country",
    if (scope == "USA_states") "USA_states" else "original_geography"
  )
  view$pairs <- view$pairs[view$pairs$analysis_level == target, ]

  keys <- paste(view$counts$species, view$counts$location)
  view$pairs <- view$pairs[
    paste(view$pairs$species, view$pairs$location1) %in% keys &
      paste(view$pairs$species, view$pairs$location2) %in% keys,
  ]

  view
}

legend_reserved_box <- function(species, scales, bounds,
                                right_inset = LEGEND_RIGHT_INSET) {
  count <- sum(vapply(species, function(name) scales[[name]]$n > 0, TRUE))
  if (!count) {
    return(NULL)
  }

  width <- min(LEGEND_COLUMN_WIDTH * count, bounds[3] - bounds[1] - 22)
  c(
    bounds[3] - right_inset - width, bounds[4] - LEGEND_RESERVED_HEIGHT,
    bounds[3] - right_inset, bounds[4]
  )
}

# ------------------------------------------------------------------------------
# Colour scales
#
# The display range, the colour ramp, the curve assignment and the symbol
# spacing are all derived here.
# ------------------------------------------------------------------------------

color_endpoints <- function(species_hue) {
  base <- as.numeric(col2rgb(species_hue)) / 255
  white_fraction <- WHITE_FRACTION

  # The light end is produced directly from the whiten fraction; it is no
  # longer limited by the contrast of the basemap.
  list(
    dark = base,
    light = base * (1 - white_fraction) + white_fraction
  )
}

integer_display_limits <- function(values) {
  if (!length(values)) {
    return(c(0, 1))
  }

  # Limits are widened to whole display units so that the rounded end labels
  # do not suggest a more precise range than is shown.
  limits <- c(
    floor(min(values) * PHIST_DISPLAY_MULTIPLIER + 1e-9),
    ceiling(max(values) * PHIST_DISPLAY_MULTIPLIER - 1e-9)
  ) / PHIST_DISPLAY_MULTIPLIER

  if (diff(limits) < 1e-12) {
    limits <- c(max(0, limits[1] - .01), min(1, limits[2] + .01))
  }

  limits
}

prepare_color_scales <- function(inputs) {
  scales <- list()

  for (index in seq_len(nrow(inputs$summary))) {
    species <- inputs$summary$species[index]
    values <- inputs$pairs$PhiST[
      inputs$pairs$species == species & inputs$pairs$selected_for_map
    ]
    scales[[species]] <- c(
      color_endpoints(inputs$summary$species_hue[index]),
      list(limits = integer_display_limits(values), n = length(values))
    )
  }

  scales
}

phist_color <- function(value, scale) {
  fraction <- max(0, min(1, (value - scale$limits[1]) / diff(scale$limits)))
  rgb_values <- scale$dark + (scale$light - scale$dark) * fraction
  rgb(rgb_values[1], rgb_values[2], rgb_values[3])
}

# ------------------------------------------------------------------------------
# Map window
# ------------------------------------------------------------------------------

species_locations <- function(inputs, species, box = NULL) {
  counts <- inputs$counts[inputs$counts$species %in% species, ]
  locations <- merge(counts, inputs$geography$positions, by = "location")

  if (!is.null(box)) {
    locations <- locations[
      locations$longitude >= box[1] & locations$longitude <= box[3] &
        locations$latitude >= box[2] & locations$latitude <= box[4],
    ]
  }

  locations
}

close_box <- function(inputs, species, regional = FALSE) {
  locations <- species_locations(inputs, species)

  if (regional) {
    locations <- locations[
      locations$longitude >= 70 & locations$longitude <= 140 &
        locations$latitude >= 0 & locations$latitude <= 50,
    ]
  }

  longitude <- range(locations$longitude)
  latitude <- range(locations$latitude)
  padding <- c(
    max(MINIMUM_PADDING_DEGREES, diff(longitude) * RELATIVE_PADDING),
    max(MINIMUM_PADDING_DEGREES, diff(latitude) * RELATIVE_PADDING)
  )

  c(
    longitude[1] - padding[1], latitude[1] - padding[2],
    longitude[2] + padding[1], latitude[2] + padding[2]
  )
}

projected_window <- function(box, width, height) {
  longitude <- (box[c(1, 3)] - CENTRAL_MERIDIAN + 180) %% 360 - 180
  if (box[3] - box[1] >= 359) longitude <- c(-180, 180)

  coordinates <- cbind(longitude, box[c(2, 4)]) * EARTH_RADIUS * pi / 180
  metres_per_point <- max(
    diff(range(coordinates[, 1])) / width,
    diff(range(coordinates[, 2])) / height
  )
  centre <- colMeans(coordinates[, 1:2, drop = FALSE])
  limits <- c(
    centre[1] - metres_per_point * width / 2,
    centre[2] - metres_per_point * height / 2,
    centre[1] + metres_per_point * width / 2,
    centre[2] + metres_per_point * height / 2
  )

  list(limits = limits, metres_per_point = metres_per_point)
}

crop_context <- function(geography, window) {
  limits <- window$limits
  crop_box <- st_bbox(
    setNames(limits, c("xmin", "ymin", "xmax", "ymax")),
    crs = st_crs(MAP_CRS)
  )

  crop_layer <- function(layer) {
    visible <- lengths(st_intersects(layer, st_as_sfc(crop_box))) > 0
    layer <- st_make_valid(layer[visible, ])
    suppressWarnings(st_crop(layer, crop_box))
  }

  # Everything is cropped in the same projection, so a global map may cross the
  # original 180 degree meridian.
  list(
    features = crop_layer(geography$provinces_center),
    countries = crop_layer(geography$countries_center),
    dashes = crop_layer(geography$dash_center)
  )
}

prepare_context <- function(inputs, box, width, height, cache) {
  key <- paste(c(box, width, height), collapse = "_")
  if (exists(key, cache, inherits = FALSE)) {
    return(get(key, cache))
  }

  window <- projected_window(box, width, height)
  cropped <- crop_context(inputs$geography, window)

  # Simplification is applied after the countries are dissolved, so that the
  # slivers between independently simplified provinces are not drawn as
  # internal national borders.
  cropped$countries <- st_simplify(
    cropped$countries,
    dTolerance = 500, preserveTopology = TRUE
  )
  cropped$features <- st_simplify(
    cropped$features,
    dTolerance = 500, preserveTopology = TRUE
  )

  context <- c(window, cropped)
  assign(key, context, cache)

  context
}

# ------------------------------------------------------------------------------
# Symbol positions
# ------------------------------------------------------------------------------

position_symbols <- function(locations, context, left, top, width, height) {
  projected <- st_coordinates(st_transform(
    st_as_sf(locations, coords = c("longitude", "latitude"), crs = 4326),
    MAP_CRS
  ))

  locations$x <- left + (projected[, 1] - context$limits[1]) /
    context$metres_per_point
  locations$y <- top + height - (projected[, 2] - context$limits[2]) /
    context$metres_per_point

  locations$raw_x <- locations$x
  locations$raw_y <- locations$y
  locations$radius <- pmax(
    MINIMUM_PIE_RADIUS,
    PIE_RADIUS_AT_211 * sqrt(locations$n / 211)
  )

  locations
}

separate_symbol_pair <- function(symbols, first, second) {
  displacement <- c(
    symbols$x[first] - symbols$x[second],
    symbols$y[first] - symbols$y[second]
  )
  distance <- sqrt(sum(displacement^2))
  required <- symbols$radius[first] + symbols$radius[second] + 1

  if (distance >= required) {
    return(symbols)
  }
  if (distance < .01) {
    angle <- (first + second) * 2.39996
    displacement <- c(cos(angle), sin(angle))
    distance <- 1
  }

  shift <- displacement / distance * ((required - distance) / 2 + .05)
  symbols$x[c(first, second)] <- symbols$x[c(first, second)] + c(shift[1], -shift[1])
  symbols$y[c(first, second)] <- symbols$y[c(first, second)] + c(shift[2], -shift[2])

  symbols
}

resolve_symbol_overlap <- function(symbols, left, top, width, height) {
  if (nrow(symbols) < 2) {
    return(symbols)
  }

  for (iteration in seq_len(150)) {
    previous <- symbols[, c("x", "y")]

    for (first in 2:nrow(symbols)) {
      for (second in seq_len(first - 1)) {
        symbols <- separate_symbol_pair(symbols, first, second)
      }
    }

    symbols$x <- pmax(
      left + symbols$radius + 1,
      pmin(left + width - symbols$radius - 1, symbols$x)
    )
    symbols$y <- pmax(
      top + symbols$radius + 1,
      pmin(top + height - symbols$radius - 1, symbols$y)
    )

    if (max(abs(as.matrix(previous) - as.matrix(symbols[, c("x", "y")]))) < .02) {
      break
    }
  }

  symbols
}

# ------------------------------------------------------------------------------
# Curves
# ------------------------------------------------------------------------------

panel_pairs <- function(inputs, species, scope) {
  pairs <- inputs$pairs[
    inputs$pairs$species %in% species & inputs$pairs$selected_for_map,
  ]

  if (scope == "none") {
    return(pairs[FALSE, ])
  }
  if (scope == "rhip_overview") {
    pairs <- pairs[
      pairs$species %in%
        c("Rhipicephalus_linnaei", "Rhipicephalus_microplus"),
    ]
  }

  pairs
}

# Quadratic Bezier curve between two symbols; the control point that stays
# inside the panel is chosen.
curve_coordinates <- function(first, second, bounds) {
  vector <- second - first
  distance <- sqrt(sum(vector^2))
  perpendicular <- c(-vector[2], vector[1]) / max(distance, 1e-9)
  midpoint <- (first + second) / 2
  controls <- rbind(
    midpoint + perpendicular * distance * .19,
    midpoint - perpendicular * distance * .19
  )

  penalty <- abs(controls[, 2] - mean(bounds[c(2, 4)])) +
    100 * pmax(bounds[1] - controls[, 1], 0) +
    100 * pmax(controls[, 1] - bounds[3], 0) +
    100 * pmax(bounds[2] - controls[, 2], 0) +
    100 * pmax(controls[, 2] - bounds[4], 0)
  control <- controls[which.min(penalty), ]

  fraction <- seq(0, 1, length.out = 80)

  cbind(
    (1 - fraction)^2 * first[1] + 2 * (1 - fraction) * fraction * control[1] +
      fraction^2 * second[1],
    (1 - fraction)^2 * first[2] + 2 * (1 - fraction) * fraction * control[2] +
      fraction^2 * second[2]
  )
}

# ------------------------------------------------------------------------------
# Vector output
#
# Every page is a vector PDF, and the per-item data used for it are written
# separately as audit tables.
# ------------------------------------------------------------------------------

draw_text <- function(label, x, y, size = 8, just = "left", italic = FALSE,
                      color = "black") {
  grid.text(
    label, unit(x, "native"), unit(y, "native"),
    just = just,
    gp = gpar(
      fontfamily = "Arial", fontsize = size,
      fontface = if (italic) 3 else 1, col = color
    )
  )
}

draw_line <- function(x, y, color = "#8D989F", width = .6) {
  grid.lines(
    unit(x, "native"), unit(y, "native"),
    gp = gpar(col = color, lwd = width, lineend = "round")
  )
}

draw_rectangle <- function(x, y, width, height, fill, border = NA) {
  grid.rect(
    unit(x, "native"), unit(y, "native"),
    unit(width, "bigpts"), unit(height, "bigpts"),
    just = c("left", "top"),
    gp = gpar(fill = fill, col = border, lwd = .5)
  )
}

draw_polygons <- function(features, mapping, fill_colors, borders = NA) {
  for (index in seq_len(nrow(features))) {
    geometry <- st_geometry(features)[[index]]
    polygons <- if (inherits(geometry, "MULTIPOLYGON")) {
      geometry
    } else if (inherits(geometry, "POLYGON")) {
      list(geometry)
    } else {
      list()
    }

    for (polygon in polygons) {
      ring_lengths <- vapply(polygon, nrow, 1L)
      coordinates <- do.call(rbind, polygon)
      grid.path(
        unit(mapping$x(coordinates[, 1]), "native"),
        unit(mapping$y(coordinates[, 2]), "native"),
        id = rep(seq_along(ring_lengths), ring_lengths),
        rule = "evenodd",
        gp = gpar(fill = fill_colors[index], col = borders, lwd = .32)
      )
    }
  }
}

draw_boundaries <- function(geometries, mapping, color, width = .6) {
  for (geometry in st_geometry(geometries)) {
    parts <- if (inherits(geometry, "MULTILINESTRING")) {
      geometry
    } else if (inherits(geometry, "LINESTRING")) {
      list(geometry)
    } else {
      list()
    }

    for (part in parts) {
      draw_line(mapping$x(part[, 1]), mapping$y(part[, 2]), color, width)
    }
  }
}

draw_basemap <- function(context, locations, mapping, country_only) {
  countries <- context$countries

  fills <- if (country_only) {
    ifelse(countries$admin %in% locations$source_admin, SAMPLED_COLOR, LAND_COLOR)
  } else {
    rep(LAND_COLOR, nrow(countries))
  }
  draw_polygons(countries, mapping, fills)

  if (!country_only) {
    subnational <- unique(
      locations$source_admin[grepl(":", locations$location, fixed = TRUE)]
    )
    features <- context$features[context$features$admin %in% subnational, ]

    # The new basemap carries national borders only; province borders come from
    # the verified administrative layer and are clipped to the national land.
    features <- suppressWarnings(
      st_intersection(features, st_union(st_geometry(countries)))
    )

    sampled_features <- unique(
      unlist(strsplit(locations$source_features, ";", fixed = TRUE))
    )
    draw_polygons(
      features, mapping,
      ifelse(features$adm1_code %in% sampled_features, SAMPLED_COLOR, LAND_COLOR)
    )

    draw_boundaries(
      st_boundary(st_geometry(features)), mapping, "#D0D5D7", .4
    )
  }

  draw_boundaries(st_boundary(st_geometry(countries)), mapping, "#ADB9C1", .65)
  draw_boundaries(context$dashes, mapping, "#9AA8B0", .8)
}

draw_pie <- function(x, y, radius, counts, colors) {
  angle_start <- pi / 2
  for (index in seq_along(counts)) {
    if (counts[index] > 0) {
      angle_end <- angle_start + 2 * pi * counts[index] / sum(counts)
      angles <- seq(
        angle_start, angle_end,
        length.out = max(5, ceiling((angle_end - angle_start) * 20))
      )
      grid.polygon(
        unit(c(x, x + radius * cos(angles)), "native"),
        unit(c(y, y + radius * sin(angles)), "native"),
        gp = gpar(fill = colors[index], col = NA)
      )
      angle_start <- angle_end
    }
  }

  grid.circle(
    unit(x, "native"), unit(y, "native"), unit(radius, "bigpts"),
    gp = gpar(fill = NA, col = "#49545B", lwd = .5)
  )
}

draw_symbols <- function(inputs, symbols, panel_id, audit) {
  for (index in order(symbols$n, decreasing = TRUE)) {
    symbol <- symbols[index, ]
    palette <- inputs$palette[inputs$palette$species == symbol$species, ]
    records <- inputs$samples[
      inputs$samples$species == symbol$species &
        inputs$samples$location == symbol$location,
    ]
    slices <- table(factor(records$clade, levels = palette$clade))

    draw_pie(symbol$x, symbol$y, symbol$radius, as.numeric(slices), palette$color)

    audit$samples[[length(audit$samples) + 1]] <- data.frame(
      panel = panel_id, species = symbol$species,
      location = symbol$location, n = symbol$n
    )
    audit$offsets[[length(audit$offsets) + 1]] <- data.frame(
      panel = panel_id, species = symbol$species, location = symbol$location,
      longitude = symbol$longitude, latitude = symbol$latitude,
      offset_x_pdf_pt = symbol$x - symbol$raw_x,
      offset_y_pdf_pt = symbol$y - symbol$raw_y
    )
  }
}

draw_phist_curves <- function(symbols, pairs, scales, bounds, panel_id, audit) {
  # Smallest values are drawn first so that the strongest comparison stays on top.
  for (index in order(pairs$PhiST, decreasing = FALSE)) {
    pair <- pairs[index, ]
    first <- which(
      symbols$species == pair$species & symbols$location == pair$location1
    )
    second <- which(
      symbols$species == pair$species & symbols$location == pair$location2
    )
    if (length(first) != 1 || length(second) != 1) next

    coordinates <- curve_coordinates(
      c(symbols$x[first], symbols$y[first]),
      c(symbols$x[second], symbols$y[second]),
      bounds
    )
    color <- phist_color(pair$PhiST, scales[[pair$species]])
    draw_line(coordinates[, 1], coordinates[, 2], color, 1.5)

    audit$curves[[length(audit$curves) + 1]] <- data.frame(
      panel = panel_id, species = pair$species,
      location1 = pair$location1, location2 = pair$location2,
      PhiST = pair$PhiST, p_BH = pair$p_BH_within_species,
      analysis_level = pair$analysis_level, line_color = color
    )
  }
}

# ------------------------------------------------------------------------------
# Location labels
# ------------------------------------------------------------------------------

location_label <- function(location) {
  aliases <- c(
    "Inner Mongolia" = "Inner Mong.", "Heilongjiang" = "Heilongjiang",
    "New Jersey" = "NJ", "North Carolina" = "NC", "South Carolina" = "SC",
    "West Virginia" = "WV", "Pennsylvania" = "PA", "Minnesota" = "MN",
    "Michigan" = "MI", "Wisconsin" = "WI", "Virginia" = "VA", "Texas" = "TX",
    "Rhode Island" = "RI", "Maine" = "ME", "Alabama" = "AL",
    "Tennessee" = "TN", "Florida" = "FL"
  )

  name <- sub("^.*:", "", location)
  if (name %in% names(aliases)) unname(aliases[name]) else name
}

# Candidate positions around the symbol, keeping the label clear of other
# symbols and of already placed labels.
label_position <- function(label, symbol, symbols, boxes, bounds, size) {
  text_width <- nchar(label) * size * .49
  choices <- expand.grid(
    angle = seq(0, 7) * pi / 4,
    distance = c(symbol$radius + 5, symbol$radius + 13, 24, 36, 52, 70)
  )

  best_score <- Inf
  for (index in seq_len(nrow(choices))) {
    x <- symbol$x + cos(choices$angle[index]) * choices$distance[index]
    y <- symbol$y + sin(choices$angle[index]) * choices$distance[index]

    x <- max(
      bounds[1] + text_width / 2 + 2,
      min(bounds[3] - text_width / 2 - 2, x)
    )
    y <- max(bounds[2] + size, min(bounds[4] - 28, y))

    box <- c(
      x - text_width / 2 - 1, y - size / 2 - 1,
      x + text_width / 2 + 1, y + size / 2 + 1
    )
    hits <- sum(
      symbols$x + symbols$radius > box[1] & symbols$x - symbols$radius < box[3] &
        symbols$y + symbols$radius > box[2] & symbols$y - symbols$radius < box[4]
    )
    overlaps <- sum(vapply(
      boxes,
      function(old) {
        box[1] < old[3] & box[3] > old[1] & box[2] < old[4] & box[4] > old[2]
      },
      TRUE
    ))

    score <- 100 * (hits + overlaps) + choices$distance[index]
    if (score < best_score) {
      best_score <- score
      best <- list(x = x, y = y, box = box)
    }
    if (!hits && !overlaps) break
  }

  best
}

draw_location_labels <- function(symbols, bounds, overview, reserved = NULL) {
  indices <- order(symbols$n, decreasing = TRUE)

  # Overview panels label countries, detail panels label the sampled localities.
  keys <- if (overview) sub(":.*", "", symbols$location) else symbols$location
  indices <- indices[!duplicated(keys[indices])]

  boxes <- if (is.null(reserved)) list() else list(reserved)
  size <- if (bounds[3] - bounds[1] < 200) 6.5 else 7

  for (index in indices) {
    symbol <- symbols[index, ]
    label <- if (overview) keys[index] else location_label(symbol$location)
    position <- label_position(label, symbol, symbols, boxes, bounds, size)

    draw_text(label, position$x, position$y, size, "centre", color = "#49565E")
    boxes[[length(boxes) + 1]] <- position$box
  }
}

# ------------------------------------------------------------------------------
# Scale bar and north arrow
# ------------------------------------------------------------------------------

draw_scale_and_north <- function(context, left, top, width, height) {
  scale_y <- top + height - 15

  # Scale bar length is corrected for the latitude of its own position.
  latitude <- (context$limits[4] - (scale_y - top) * context$metres_per_point) /
    6378137
  projection_factor <- sqrt(1 - .00669437999014 * sin(latitude)^2) /
    cos(latitude)

  distance <- if (width < 200) 500000 else 2000000
  while (distance * projection_factor / context$metres_per_point > width * .32) {
    distance <- distance / 2
  }
  length <- distance * projection_factor / context$metres_per_point

  scale_x <- left + 7
  draw_rectangle(scale_x, scale_y, length / 2, 2.5, "#52616A", "#52616A")
  draw_rectangle(scale_x + length / 2, scale_y, length / 2, 2.5, "white", "#52616A")
  draw_text("0", scale_x, scale_y - 5, 6, "centre")
  draw_text(
    paste0(format(distance / 1000, trim = TRUE), " km"),
    scale_x + length, scale_y - 5, 6, "centre"
  )

  if (nzchar(MAP_APPROVAL_NUMBER)) {
    draw_text(
      MAP_APPROVAL_NUMBER, scale_x, scale_y + 10, 6,
      color = "#52616A"
    )
  }

  north_x <- left + width - 9
  grid.polygon(
    unit(c(north_x, north_x - 2.3, north_x + 2.3), "native"),
    unit(c(top + 17, top + 30, top + 30), "native"),
    gp = gpar(fill = "#52616A", col = NA)
  )
  draw_text("N", north_x, top + 10, 6.5, "centre")
}

map_viewport <- function(left, top, width, height) {
  viewport(
    x = unit(left, "native"), y = unit(top, "native"),
    width = unit(width, "bigpts"), height = unit(height, "bigpts"),
    just = c("left", "top"),
    xscale = c(left, left + width), yscale = c(top + height, top), clip = "on"
  )
}

# ------------------------------------------------------------------------------
# One map panel
# ------------------------------------------------------------------------------

draw_map <- function(inputs, scales, cache, audit, species, box = NULL, width,
                     height, tag, title, scope = "all", overview = FALSE,
                     country_only = FALSE, panel_id = tag,
                     country_species = character(), region = NULL,
                     legend_species = species) {
  view <- make_panel_inputs(inputs, species, country_species, scope, region)
  if (is.null(box)) box <- close_box(view, species)

  draw_text(paste0("(", tag, ")"), 0, 5, 10)
  draw_text(title, 20, 5, 8, italic = TRUE)

  size <- c(width - 6, height - 28)
  left <- (width - size[1]) / 2
  top <- 17
  map_width <- size[1]
  map_height <- size[2]

  context <- prepare_context(view, box, map_width, map_height, cache)
  locations <- species_locations(view, species)
  symbols <- position_symbols(locations, context, left, top, map_width, map_height)
  symbols <- resolve_symbol_overlap(
    symbols, left, top, map_width, map_height
  )
  pairs <- panel_pairs(view, species, scope)

  # Colour scales are rebuilt from the displayed subset, not taken from the
  # module-level list, because the selected pairs depend on the panel.
  scales <- prepare_color_scales(view)

  bounds <- c(left, top, left + map_width, top + map_height)
  mapping <- list(
    x = function(x) left + (x - context$limits[1]) / context$metres_per_point,
    y = function(y) top + map_height - (y - context$limits[2]) / context$metres_per_point
  )
  legend_inset <- LEGEND_RIGHT_INSET

  draw_rectangle(left, top, map_width, map_height, OCEAN_COLOR, "#C8D2D9")

  pushViewport(map_viewport(left, top, map_width, map_height))
  draw_basemap(context, locations, mapping, country_only)
  draw_phist_curves(symbols, pairs, scales, bounds, panel_id, audit)
  draw_symbols(view, symbols, panel_id, audit)
  draw_location_labels(
    symbols, bounds, overview,
    legend_reserved_box(legend_species, scales, bounds, legend_inset)
  )
  draw_scale_and_north(context, left, top, map_width, map_height)
  draw_phist_legend(legend_species, scales, bounds, panel_id, audit, legend_inset)
  popViewport()

  audit$windows[[length(audit$windows) + 1]] <- data.frame(
    panel = panel_id,
    xmin = box[1], ymin = box[2], xmax = box[3], ymax = box[4],
    actual_ymin = context$limits[2] / 6378137 * 180 / pi,
    actual_ymax = context$limits[4] / 6378137 * 180 / pi,
    range_rule = "fixed_S36_tighter_latitudes_and_previous_other_ranges"
  )
}

draw_panel <- function(x, y, width, height, draw) {
  pushViewport(viewport(
    x = unit(x, "native"), y = unit(y, "native"),
    width = unit(width, "bigpts"), height = unit(height, "bigpts"),
    just = c("left", "top"),
    xscale = c(0, width), yscale = c(height, 0)
  ))
  draw()
  popViewport()
}

open_pdf <- function(filename) {
  cairo_pdf(
    file.path(OUTPUT_DIRECTORY, filename),
    width = 8.5, height = 11, family = "Arial", pointsize = 8
  )
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 612), yscale = c(792, 0)))
}

draw_samples_legend <- function(y) {
  draw_text("Samples", PAGE_LEFT, y, 8)
  for (index in seq_along(c(1, 20, 80, 211))) {
    count <- c(1, 20, 80, 211)[index]
    x <- 111 + (index - 1) * 41
    radius <- max(MINIMUM_PIE_RADIUS, PIE_RADIUS_AT_211 * sqrt(count / 211))
    draw_pie(x, y, radius, 1, "#C3CBD0")
    draw_text(count, x + radius + 3, y, 7)
  }
}

draw_phist_legend <- function(species, scales, bounds, panel_id, audit,
                              right_inset = LEGEND_RIGHT_INSET) {
  species <- species[vapply(species, function(name) scales[[name]]$n > 0, TRUE)]
  if (!length(species)) {
    return(invisible(NULL))
  }

  column_width <- min(
    LEGEND_COLUMN_WIDTH, (bounds[3] - bounds[1] - 22) / length(species)
  )
  start <- bounds[3] - length(species) * column_width - right_inset
  bar_top <- bounds[4] - LEGEND_BOTTOM_INSET - LEGEND_BAR_HEIGHT

  # The shared title sits on the first line and the species names on the second,
  # so the repeated titles do not eat into the map.
  draw_text(
    expression(plain("\u03a6")[plain(ST)] * " \u00d7 100"),
    start + length(species) * column_width / 2, bar_top - 24,
    LEGEND_TITLE_SIZE, "centre"
  )

  for (index in seq_along(species)) {
    name <- species[index]
    x <- start + (index - .5) * column_width
    scale <- scales[[name]]

    species_label <- sub("^([^_])[^_]+_", "\\1. ", name)
    draw_text(
      species_label, x, bar_top - 10, LEGEND_SPECIES_SIZE, "centre",
      italic = TRUE
    )

    for (step in 0:79) {
      draw_rectangle(
        x - 7, bar_top + step * LEGEND_BAR_HEIGHT / 80, 6,
        LEGEND_BAR_HEIGHT / 80 + .02,
        phist_color(scale$limits[2] - diff(scale$limits) * step / 79, scale)
      )
    }

    draw_text(
      round(scale$limits[2] * PHIST_DISPLAY_MULTIPLIER), x + 2, bar_top + 2.8, 5.5
    )
    draw_text(
      round(scale$limits[1] * PHIST_DISPLAY_MULTIPLIER), x + 2,
      bar_top + LEGEND_BAR_HEIGHT - 2.8, 5.5
    )

    audit$legends[[length(audit$legends) + 1]] <- data.frame(
      panel = panel_id, species = name,
      minimum = scale$limits[1], maximum = scale$limits[2],
      selected_pairs = scale$n, display_multiplier = PHIST_DISPLAY_MULTIPLIER,
      block_right_inset = right_inset, bar_bottom_inset = LEGEND_BOTTOM_INSET,
      species_font = "Arial Italic"
    )
  }
}

# ------------------------------------------------------------------------------
# Figure 3: Ixodes and Rhipicephalus
# ------------------------------------------------------------------------------

draw_figure3 <- function(inputs, scales, cache, audit) {
  ixodes <- c("Ixodes_scapularis", "Ixodes_ricinus", "Ixodes_persulcatus")
  rhipicephalus <- c(
    "Rhipicephalus_microplus", "Rhipicephalus_sanguineus",
    "Rhipicephalus_haemaphysaloides", "Rhipicephalus_linnaei"
  )
  column_width <- (PAGE_WIDTH - 18) / 3

  open_pdf("Figure_3.pdf")

  # Ixodes overview and one detail panel per species.
  draw_panel(
    PAGE_LEFT, 54, PAGE_WIDTH, 115,
    function() {
      draw_map(
        inputs, scales, cache, audit, ixodes,
        box = c(-20, 22, 340, 70),
        width = PAGE_WIDTH, height = 115, tag = "A", title = "Ixodes",
        scope = "none", overview = TRUE, country_only = TRUE,
        panel_id = "Figure3_A", legend_species = character()
      )
    }
  )

  for (index in 1:3) {
    species <- ixodes[index]
    draw_panel(
      PAGE_LEFT + (index - 1) * (column_width + 9), 179, column_width, 170,
      function() {
        draw_map(
          inputs, scales, cache, audit, species,
          width = column_width, height = 170, tag = LETTERS[index + 1],
          title = sub("^([^_])[^_]+_", "\\1. ", species),
          country_only = species == "Ixodes_ricinus",
          panel_id = paste0("Figure3_", LETTERS[index + 1])
        )
      }
    )
  }

  # Rhipicephalus overview on the country level, then Asian detail panels.
  draw_panel(
    PAGE_LEFT, 359, PAGE_WIDTH, 175,
    function() {
      draw_map(
        inputs, scales, cache, audit, rhipicephalus,
        box = c(-20, -35, 340, 54),
        width = PAGE_WIDTH, height = 175, tag = "E", title = "Rhipicephalus",
        scope = "rhip_overview", overview = TRUE, country_only = TRUE,
        panel_id = "Figure3_E", country_species = "Rhipicephalus_microplus",
        legend_species = c("Rhipicephalus_microplus", "Rhipicephalus_linnaei")
      )
    }
  )

  for (index in 1:3) {
    species <- rhipicephalus[index]
    draw_panel(
      PAGE_LEFT + (index - 1) * (column_width + 9), 544, column_width, 160,
      function() {
        draw_map(
          inputs, scales, cache, audit, species,
          width = column_width, height = 160, tag = LETTERS[index + 5],
          title = sub("^([^_])[^_]+_", "\\1. ", species),
          region = "asia", panel_id = paste0("Figure3_", LETTERS[index + 5])
        )
      }
    )
  }

  draw_samples_legend(724)
  popViewport()
  dev.off()
}

# ------------------------------------------------------------------------------
# Supplementary maps: Dermacentor and Haemaphysalis
# ------------------------------------------------------------------------------

draw_supplementary <- function(inputs, scales, cache, audit) {
  dermacentor <- inputs$summary$species[
    grepl("^Dermacentor_", inputs$summary$species)
  ]
  haemaphysalis <- inputs$summary$species[
    grepl("^Haemaphysalis_", inputs$summary$species)
  ]

  if ("Figure_S35" %in% MAP_TARGETS) {
    open_pdf(SUPPLEMENTARY_MAP_FILES[["Figure_S35"]])
    draw_panel(
      PAGE_LEFT, 54, PAGE_WIDTH, 550,
      function() {
        draw_map(
          inputs, scales, cache, audit, dermacentor,
          width = PAGE_WIDTH, height = 550, tag = "A", title = "Dermacentor",
          panel_id = "FigureS35_A"
        )
      }
    )
    draw_samples_legend(724)
    popViewport()
    dev.off()
  }

  if ("Figure_S36" %in% MAP_TARGETS) {
    open_pdf(SUPPLEMENTARY_MAP_FILES[["Figure_S36"]])

    # (A) country level, (B) Asian detail, (C) USA states of H. longicornis.
    draw_panel(
      PAGE_LEFT, 54, PAGE_WIDTH, 190,
      function() {
        draw_map(
          inputs, scales, cache, audit, haemaphysalis,
          box = c(-20, -45, 340, 50),
          width = PAGE_WIDTH, height = 190, tag = "A", title = "Haemaphysalis",
          overview = TRUE, country_only = TRUE, panel_id = "FigureS36_A",
          country_species = "Haemaphysalis_longicornis",
          legend_species = "Haemaphysalis_longicornis"
        )
      }
    )

    draw_panel(
      PAGE_LEFT, 254, PAGE_WIDTH, 255,
      function() {
        draw_map(
          inputs, scales, cache, audit, haemaphysalis,
          box = c(75, 18, 148, 51),
          width = PAGE_WIDTH, height = 255, tag = "B", title = "Haemaphysalis",
          region = "asia", panel_id = "FigureS36_B"
        )
      }
    )

    draw_panel(
      PAGE_LEFT, 519, PAGE_WIDTH, 175,
      function() {
        draw_map(
          inputs, scales, cache, audit, "Haemaphysalis_longicornis",
          box = c(-95, 33, -70, 45), width = PAGE_WIDTH, height = 175,
          tag = "C", title = "H. longicornis", scope = "USA_states",
          panel_id = "FigureS36_C"
        )
      }
    )

    draw_samples_legend(724)
    popViewport()
    dev.off()
  }
}

# ------------------------------------------------------------------------------
# Audit tables
#
# The last three names are kept in their archived form.
# ------------------------------------------------------------------------------

write_audits <- function(audit, scales) {
  # Archived table names are kept unchanged so that the stored files are read
  # as supplied; only the first two are newly named.
  files <- c(
    samples = "12_map_sample_coverage.csv",
    offsets = "13_map_display_offsets.csv",
    curves = "27_地图实际显示PhiST连线.csv",
    legends = "28_地图PhiST图例范围.csv",
    windows = "29_地图裁切范围.csv"
  )

  for (name in names(files)) {
    write.csv(
      do.call(rbind, audit[[name]]), file.path(OUTPUT_DIRECTORY, files[[name]]),
      row.names = FALSE, fileEncoding = "UTF-8"
    )
  }
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

main <- function() {
  sf_use_s2(FALSE)
  dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)

  inputs <- read_inputs()
  scales <- prepare_color_scales(inputs)
  cache <- new.env(parent = emptyenv())
  audit <- new.env(parent = emptyenv())

  inputs <- add_country_positions(inputs)
  inputs <- prepare_basemap(inputs)

  audit$samples <- list()
  audit$offsets <- list()
  audit$curves <- list()
  audit$legends <- list()
  audit$windows <- list()

  if ("Figure_3" %in% MAP_TARGETS) draw_figure3(inputs, scales, cache, audit)
  draw_supplementary(inputs, scales, cache, audit)
  write_audits(audit, scales)

  cat(
    "[write] ", length(MAP_TARGETS),
    " PhiST map PDFs, curve assignments and legend ranges written.\n",
    sep = ""
  )
}

if (sys.nframe() == 0L) main()
