# 05_main_tables_figures.R  --  STEP 5: main tables and figures
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# Assembles the two tables intended for the main text and the paper figures,
# from the fitted objects of the two analyses.
#
#   1. main table 1: models within two AICc units, both analyses
#   2. main table 2: coefficients of the best model and of the model average
#   3. marginal response curves, model-averaged over the top set
#   4. predicted abundance per station on a satellite basemap
#
# Marginal curves vary one covariate over its observed range while the others
# are held at their mean, which is zero on the standardised scale used for
# fitting; the horizontal axis is returned to the original units.
#
# The basemap is downloaded once and cached, so the figures regenerate offline.
# If it cannot be downloaded the map is drawn without it.
#
# Input : workspace/01_import.rds, 03_dm_selection.rds, 04_rn.rds
# Output: tables/05_*.csv, figures/05_*
# ===========================================================================

#install.packages(c("patchwork", "ggplot2"))

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_DATA <- "C:/Users/Dell/Desktop/codici tesi/file GIS"
ROOT_OUT  <- "C:/Users/Dell/Desktop/codici articolo"

DIR_TAB <- file.path(ROOT_OUT, "tables")
DIR_FIG <- file.path(ROOT_OUT, "figures")
DIR_WS  <- file.path(ROOT_OUT, "workspace")
for (d in c(DIR_TAB, DIR_FIG, DIR_WS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

N_GRID       <- 100        # points along a marginal curve
MAP_PAD      <- 0.02       # degrees added around the station envelope
MAP_ZOOM     <- 12         # basemap tile zoom
MAP_PROVIDER <- "Esri World Imagery"
MAP_URL      <- "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/%d/%d/%d"
TILE_CACHE   <- file.path(DIR_WS, "05_basemap.rds")
FIG_DPI      <- 600

# One colour per analysis, kept identical in every figure. The two hues are a
# green and a purple from a diverging palette, so they remain distinguishable
# in greyscale and under the common forms of colour vision deficiency.
PAL_COUNT <- list(line = "#1B7837", band = "#C7E9C0", point = "#1B7837")
PAL_RN    <- list(line = "#762A83", band = "#E3D3E8", point = "#762A83")

COV_LABELS <- c(
  elevation            = "Elevation (m a.s.l.)",
  tri                  = "Terrain ruggedness index",
  tpi                  = "Topographic position index",
  dist_water_points    = "Distance to water point (m)",
  dist_fish_river      = "Distance to Fish River (m)",
  dist_konkiep_river   = "Distance to Konkiep River (m)",
  dist_tourist_road    = "Distance to tourist road (m)",
  dist_staff_buildings = "Distance to staff buildings (m)")

suppressPackageStartupMessages({
  library(unmarked)
  library(ggplot2)
  library(sf)
})

imp <- readRDS(file.path(DIR_WS, "01_import.rds"))
dm  <- readRDS(file.path(DIR_WS, "03_dm_selection.rds"))
rn  <- readRDS(file.path(DIR_WS, "04_rn.rds"))
set.seed(imp$seed)

# The count-based analysis is the Dail-Madsen model or, when the population
# was treated as closed in STEP 2, the closed N-mixture model.
METHODS <- c(dm$MODEL_LABEL, "Royle-Nichols")
PAL <- setNames(list(PAL_COUNT, PAL_RN), METHODS)

cat("\n=========================================================\n")
cat(" STEP 5 -- main tables and figures\n")
cat("=========================================================\n\n")

# --- 1. MAIN TABLE 1: MODELS WITHIN TWO AICc UNITS -------------------------

top_rows <- function(obj, method, lam_col = "lambda", det_col) {
  r <- obj$rank_tab[obj$rank_tab$in_top_set, ]
  data.frame(
    method     = method,
    model      = r$model,
    hypothesis = r$hypothesis,
    abundance  = r[[lam_col]],
    detection  = r[[det_col]],
    nPars      = r$nPars,
    logLik     = round(r$logLik, 2),
    AICc       = round(r$AICc, 2),
    delta_AICc = round(r$delta_AICc, 2),
    weight     = round(r$weight, 3),
    stringsAsFactors = FALSE)
}

main_models <- rbind(
  top_rows(dm, METHODS[1], det_col = "p"),
  top_rows(rn, METHODS[2], det_col = "r_formula"))

cat("-- main table 1: models within two AICc units --\n")
print(main_models, row.names = FALSE)
cat("\n")

# --- 2. MAIN TABLE 2: COEFFICIENTS -----------------------------------------
# Terms are labelled by submodel and by covariate, so that the two analyses
# can be read side by side.

tidy_coef <- function(obj, method) {
  d <- obj$coef_tab
  sub <- sub("\\(.*$", "", d$term)
  par <- sub("^[^(]*\\(", "", sub("\\)$", "", d$term))
  data.frame(
    method       = method,
    submodel     = sub,
    parameter    = par,
    best_estimate = round(d$best_estimate, 4),
    best_SE       = round(d$best_SE, 4),
    best_lower    = round(d$best_lower, 4),
    best_upper    = round(d$best_upper, 4),
    avg_estimate  = round(d$avg_estimate, 4),
    avg_uncond_SE = round(d$avg_uncond_SE, 4),
    avg_lower     = round(d$avg_lower, 4),
    avg_upper     = round(d$avg_upper, 4),
    weight_present = round(d$weight_present, 3),
    stringsAsFactors = FALSE)
}

main_coef <- rbind(tidy_coef(dm, METHODS[1]), tidy_coef(rn, METHODS[2]))

cat("-- main table 2: coefficients --\n")
print(main_coef[, c("method", "submodel", "parameter", "best_estimate",
                    "best_SE", "avg_estimate", "avg_uncond_SE")],
      row.names = FALSE)
cat("\n")

write.csv(main_models, file.path(DIR_TAB, "05_main_top_models.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(main_coef, file.path(DIR_TAB, "05_main_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# --- 3. MARGINAL RESPONSE CURVES -------------------------------------------

SITE <- imp$site[match(imp$stations, imp$site$Station), imp$COVS, drop = FALSE]
rownames(SITE) <- NULL
center <- imp$center
scale_ <- imp$scale

# Covariates entering the top set of an analysis, in the abundance submodel.
covs_in_top <- function(obj, prefix) {
  tm <- unlist(lapply(obj$top_set, function(m) names(coef(obj$fits[[m]]))))
  tm <- unique(tm[grepl(paste0("^", prefix, "\\("), tm)])
  tm <- sub("\\)$", "", sub(paste0("^", prefix, "\\("), "", tm))
  intersect(names(COV_LABELS), tm)
}

# Model-averaged prediction along one covariate, others held at the mean.
marginal <- function(obj, cov, type, prefix) {
  z   <- seq(min(SITE[[cov]]), max(SITE[[cov]]), length.out = N_GRID)
  nd  <- SITE[rep(1, N_GRID), , drop = FALSE]
  for (v in names(nd)) nd[[v]] <- 0
  nd[[cov]] <- z
  rownames(nd) <- NULL

  w <- obj$weights_top
  preds <- lapply(obj$top_set, function(m) {
    p <- try(predict(obj$fits[[m]], type = type, newdata = nd, appendData = FALSE),
             silent = TRUE)
    if (inherits(p, "try-error")) NULL else
      data.frame(fit = as.numeric(p$Predicted), se = as.numeric(p$SE))
  })
  names(preds) <- obj$top_set
  keep <- !vapply(preds, is.null, logical(1))
  if (!any(keep)) return(NULL)
  ww <- w[keep] / sum(w[keep])

  fm <- vapply(preds[keep], function(x) x$fit, numeric(N_GRID))
  sm <- vapply(preds[keep], function(x) x$se,  numeric(N_GRID))
  if (is.null(dim(fm))) { fm <- matrix(fm, N_GRID, 1); sm <- matrix(sm, N_GRID, 1) }

  est <- as.numeric(fm %*% ww)
  unc <- as.numeric(sqrt(sm^2 + (fm - est)^2) %*% ww)

  data.frame(covariate = cov,
             x_raw = z * scale_[[cov]] + center[[cov]],
             fit = est,
             lower = pmax(est - 1.959964 * unc, 0),
             upper = est + 1.959964 * unc,
             stringsAsFactors = FALSE)
}

curves <- list()
rugs   <- list()
for (spec in list(list(obj = dm, method = METHODS[1], type = dm$LAM_TYPE, prefix = "lam"),
                  list(obj = rn, method = METHODS[2], type = "state",  prefix = "lam"))) {
  cv <- covs_in_top(spec$obj, spec$prefix)
  cat("-- marginal curves, ", spec$method, " : ",
      if (length(cv)) paste(cv, collapse = ", ") else "no abundance covariate in the top set",
      "\n", sep = "")
  for (cc in cv) {
    m <- marginal(spec$obj, cc, spec$type, spec$prefix)
    if (is.null(m)) { cat("   prediction failed for ", cc, "\n", sep = ""); next }
    m$method <- spec$method
    curves[[length(curves) + 1]] <- m
    rugs[[length(rugs) + 1]] <- data.frame(
      covariate = cc, method = spec$method,
      x_raw = SITE[[cc]] * scale_[[cc]] + center[[cc]], stringsAsFactors = FALSE)
  }
}
cat("\n")

if (length(curves)) {
  curve_df <- do.call(rbind, curves)
  rug_df   <- do.call(rbind, rugs)

  curve_theme <- theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey93"),
          panel.border = element_rect(colour = "grey35", linewidth = 0.45),
          strip.background = element_rect(fill = "grey96", colour = "grey35",
                                          linewidth = 0.45),
          strip.text = element_text(size = 9.5, margin = margin(4, 4, 4, 4)),
          axis.text = element_text(size = 8.5, colour = "grey20"),
          axis.title = element_text(size = 10),
          plot.margin = margin(4, 8, 4, 4))

  # One panel row per analysis, so that each keeps its own colour and its own
  # set of covariates, and the two are assembled afterwards.
  one_method <- function(method) {
    cd <- curve_df[curve_df$method == method, , drop = FALSE]
    if (!nrow(cd)) return(NULL)
    rd <- rug_df[rug_df$method == method, , drop = FALSE]
    lev <- unname(COV_LABELS[intersect(names(COV_LABELS), unique(cd$covariate))])
    cd$label <- factor(COV_LABELS[cd$covariate], levels = lev)
    rd$label <- factor(COV_LABELS[rd$covariate], levels = lev)

    ggplot(cd, aes(x = x_raw, y = fit)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  fill = PAL[[method]]$band, colour = NA) +
      geom_line(colour = PAL[[method]]$line, linewidth = 1) +
      geom_rug(data = rd, aes(x = x_raw), inherit.aes = FALSE, sides = "b",
               length = unit(0.028, "npc"), linewidth = 0.3,
               colour = PAL[[method]]$line, alpha = 0.7) +
      facet_wrap(~ label, scales = "free", nrow = 1) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
      scale_x_continuous(expand = expansion(mult = 0.02)) +
      labs(x = NULL,
           y = expression(hat(lambda)~"(individuals per station)")) +
      curve_theme
  }

  plots <- lapply(METHODS, one_method)
  names(plots) <- METHODS
  n_panel <- vapply(METHODS, function(m)
    length(unique(curve_df$covariate[curve_df$method == m])), integer(1))

  for (m in METHODS) {
    if (is.null(plots[[m]])) next
    stem <- file.path(DIR_FIG,
                      paste0("05_marginal_curves_", tolower(gsub("-", "", m))))
    w <- max(3.4, 2.6 * n_panel[[m]])
    ggsave(paste0(stem, ".png"), plots[[m]], width = w, height = 3.1,
           dpi = FIG_DPI, limitsize = FALSE)
    ggsave(paste0(stem, ".pdf"), plots[[m]], width = w, height = 3.1,
           limitsize = FALSE)
    cat("  ", basename(stem), ".(png|pdf)\n", sep = "")
  }

  # Combined figure, one analysis per row, tagged (a) and (b).
  keep_p <- plots[!vapply(plots, is.null, logical(1))]
  w_all  <- max(3.4, 2.6 * max(n_panel))
  h_all  <- 3.1 * length(keep_p)
  combined <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- Reduce(`/`, keep_p) +
      patchwork::plot_annotation(tag_levels = "a",
                                 tag_prefix = "(", tag_suffix = ")")
    combined <- combined & theme(plot.tag = element_text(size = 11, face = "bold"),
                                 plot.tag.position = c(0.005, 0.99))
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    combined <- cowplot::plot_grid(plotlist = keep_p, ncol = 1,
                                   labels = paste0("(", letters[seq_along(keep_p)], ")"),
                                   label_size = 11, hjust = 0)
  } else {
    cat("  install patchwork or cowplot to obtain the combined figure\n")
  }
  if (!is.null(combined)) {
    ggsave(file.path(DIR_FIG, "05_marginal_curves.png"), combined,
           width = w_all, height = h_all, dpi = FIG_DPI, limitsize = FALSE)
    ggsave(file.path(DIR_FIG, "05_marginal_curves.pdf"), combined,
           width = w_all, height = h_all, limitsize = FALSE)
    cat("  05_marginal_curves.(png|pdf)\n")
  }

  write.csv(curve_df, file.path(DIR_TAB, "05_marginal_curves.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n")
} else {
  cat("  no marginal curve produced: no abundance covariate in either top set\n\n")
}

# --- 4. PREDICTED ABUNDANCE PER STATION ON A SATELLITE BASEMAP -------------

coords <- read.csv(file.path(ROOT_DATA, "hyena_analysis_complete_dataset.csv"),
                   stringsAsFactors = FALSE)
coords <- coords[match(imp$stations, coords$Station), ]
coords$X <- coords$X.N.19.15      # longitude, WGS84
coords$Y <- coords$Y.N.19.15      # latitude, WGS84
stopifnot(!anyNA(coords$X), !anyNA(coords$Y))

pts <- data.frame(
  Station = rep(imp$stations, 2),
  X = rep(coords$X, 2), Y = rep(coords$Y, 2),
  method = factor(rep(METHODS, each = length(imp$stations)), levels = METHODS),
  lambda = c(dm$site_pred$lambda_avg, rn$site_pred$lambda_avg),
  lambda_SE = c(dm$site_pred$lambda_avg_SE, rn$site_pred$lambda_avg_SE),
  stringsAsFactors = FALSE)

pts_sf  <- st_as_sf(pts, coords = c("X", "Y"), crs = 4326)
pts_3857 <- st_transform(pts_sf, 3857)
xy <- st_coordinates(pts_3857)
pts$x <- xy[, 1]; pts$y <- xy[, 2]

bb <- c(min(coords$X) - MAP_PAD, min(coords$Y) - MAP_PAD,
        max(coords$X) + MAP_PAD, max(coords$Y) + MAP_PAD)
cat("-- station envelope (degrees) --\n")
cat(sprintf("  longitude %.5f to %.5f   latitude %.5f to %.5f\n",
            bb[1], bb[3], bb[2], bb[4]))

area <- st_transform(st_as_sfc(st_bbox(c(xmin = bb[1], ymin = bb[2],
                                         xmax = bb[3], ymax = bb[4]),
                                       crs = st_crs(4326))), 3857)
area_bb <- st_bbox(area)

# Web-Mercator tiles covering the station envelope, downloaded one by one from
# the tile service, georeferenced from their indices, mosaicked and cropped to
# the envelope.
fetch_tiles <- function(bb, zoom, crop_bb) {
  n <- 2^zoom
  half <- 20037508.342789244                     # half the width of EPSG:3857
  size <- 2 * half / n                           # tile width in metres
  tx <- function(lon) floor((lon + 180) / 360 * n)
  ty <- function(lat) {
    r <- lat * pi / 180
    floor((1 - log(tan(r) + 1 / cos(r)) / pi) / 2 * n)
  }
  xs <- tx(bb[1]):tx(bb[3])
  ys <- ty(bb[4]):ty(bb[2])
  tiles <- list()
  for (x in xs) for (y in ys) {
    f <- tempfile(fileext = ".jpg")
    utils::download.file(sprintf(MAP_URL, zoom, y, x), f, mode = "wb", quiet = TRUE)
    r <- terra::rast(f)
    terra::ext(r) <- terra::ext(-half + x * size, -half + (x + 1) * size,
                                half - (y + 1) * size, half - y * size)
    terra::crs(r) <- "EPSG:3857"
    tiles[[length(tiles) + 1]] <- r
  }
  mosaic <- if (length(tiles) == 1) tiles[[1]] else do.call(terra::merge, tiles)
  terra::crop(mosaic, terra::ext(crop_bb[["xmin"]], crop_bb[["xmax"]],
                                 crop_bb[["ymin"]], crop_bb[["ymax"]]))
}

tile_df <- NULL
if (file.exists(TILE_CACHE)) {
  tile_df <- readRDS(TILE_CACHE)
  cat("  basemap read from cache\n")
} else if (requireNamespace("terra", quietly = TRUE)) {
  cat("  downloading basemap (", MAP_PROVIDER, ", zoom ", MAP_ZOOM, ") ... ", sep = "")
  utils::flush.console()
  tl <- try(fetch_tiles(bb, MAP_ZOOM, area_bb), silent = TRUE)
  if (inherits(tl, "try-error")) {
    cat("FAILED -- ", trimws(as.character(tl)), "\n", sep = "")
  } else {
    m <- terra::as.data.frame(tl, xy = TRUE)
    chan <- setdiff(names(m), c("x", "y"))
    if (length(chan) >= 3) {
      tile_df <- data.frame(x = m$x, y = m$y,
                            hex = grDevices::rgb(m[[chan[1]]], m[[chan[2]]],
                                                 m[[chan[3]]], maxColorValue = 255),
                            stringsAsFactors = FALSE)
      saveRDS(tile_df, TILE_CACHE)
      cat(nrow(tile_df), " pixels, cached\n", sep = "")
    } else {
      cat("unexpected number of channels (", length(chan), ")\n", sep = "")
    }
  }
} else {
  cat("  package terra is not installed; the map is drawn without a basemap\n")
}

map_theme <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 7),
        strip.background = element_rect(fill = "grey95", colour = "grey60",
                                        linewidth = 0.3),
        legend.key.height = unit(0.9, "lines"))

p_map <- ggplot()
if (!is.null(tile_df))
  p_map <- p_map + geom_raster(data = tile_df, aes(x = x, y = y, fill = hex)) +
    scale_fill_identity()
xlim_map <- if (!is.null(tile_df)) range(tile_df$x) else as.numeric(area_bb[c("xmin", "xmax")])
ylim_map <- if (!is.null(tile_df)) range(tile_df$y) else as.numeric(area_bb[c("ymin", "ymax")])

p_map <- p_map +
  geom_point(data = pts, aes(x = x, y = y, size = lambda, colour = lambda),
             alpha = 0.9) +
  scale_size_continuous(range = c(1.2, 5.5)) +
  scale_colour_viridis_c(option = "C", name = expression(lambda)) +
  guides(size = "none") +
  facet_wrap(~ method, nrow = 1) +
  coord_equal(xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  map_theme

ggsave(file.path(DIR_FIG, "05_station_estimates_map.png"), p_map,
       width = 9, height = 4.8, dpi = FIG_DPI, limitsize = FALSE)
ggsave(file.path(DIR_FIG, "05_station_estimates_map.pdf"), p_map,
       width = 9, height = 4.8, limitsize = FALSE)

# The same map with the estimate printed at each station, for a reader who
# needs the value rather than the pattern. Labels are placed away from the
# points where ggrepel is available.
pts$label <- formatC(pts$lambda, format = "f", digits = 2)

p_map_lab <- p_map
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_map_lab <- p_map_lab +
    ggrepel::geom_text_repel(
      data = pts, aes(x = x, y = y, label = label),
      size = 2.5, colour = "white", fontface = "bold",
      bg.colour = "grey15", bg.r = 0.16,
      segment.colour = "white", segment.size = 0.25, segment.alpha = 0.8,
      min.segment.length = 0.15, box.padding = 0.3, point.padding = 0.2,
      max.overlaps = Inf, seed = imp$seed)
} else {
  cat("  install ggrepel for non-overlapping labels; plain labels used\n")
  p_map_lab <- p_map_lab +
    geom_text(data = pts, aes(x = x, y = y, label = label),
              size = 2.5, colour = "white", fontface = "bold",
              vjust = -1.1)
}

ggsave(file.path(DIR_FIG, "05_station_estimates_map_labelled.png"), p_map_lab,
       width = 11, height = 5.8, dpi = FIG_DPI, limitsize = FALSE)
ggsave(file.path(DIR_FIG, "05_station_estimates_map_labelled.pdf"), p_map_lab,
       width = 11, height = 5.8, limitsize = FALSE)

write.csv(pts[, c("Station", "method", "lambda", "lambda_SE")],
          file.path(DIR_TAB, "05_station_estimates.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# --- 4b. PREDICTED ABUNDANCE PER STATION, WITH CONFIDENCE INTERVALS --------
# The map shows where the stations are; this shows how precise each estimate is
# and how the two analyses order the same stations. Stations are ordered by the
# count-based estimate, so the two series can be read against one another.

est <- data.frame(
  Station = rep(imp$stations, 2),
  method  = factor(rep(METHODS, each = length(imp$stations)), levels = METHODS),
  fit     = c(dm$site_pred$lambda_avg, rn$site_pred$lambda_avg),
  lower   = c(dm$site_pred$lambda_avg_lower, rn$site_pred$lambda_avg_lower),
  upper   = c(dm$site_pred$lambda_avg_upper, rn$site_pred$lambda_avg_upper),
  stringsAsFactors = FALSE)

ord <- dm$site_pred$Station[order(dm$site_pred$lambda_avg)]
est$Station <- factor(est$Station, levels = ord)

det_lab <- data.frame(
  Station = factor(dm$site_pred$Station, levels = ord),
  n = dm$site_pred$individuals, stringsAsFactors = FALSE)
est <- merge(est, det_lab, by = "Station")

p_est <- ggplot(est, aes(x = fit, y = Station, colour = method)) +
  geom_linerange(aes(xmin = lower, xmax = upper),
                 position = position_dodge(width = 0.62), linewidth = 0.55,
                 alpha = 0.85) +
  geom_point(position = position_dodge(width = 0.62), size = 2,
             shape = 21, fill = "white", stroke = 0.9) +
  scale_colour_manual(values = c(PAL[[METHODS[1]]]$line, PAL[[METHODS[2]]]$line),
                      breaks = METHODS) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  labs(x = expression(hat(lambda)~"(individuals per station)"), y = NULL,
       colour = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey95"),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        axis.text.y = element_text(family = "mono", size = 8),
        legend.position = "top", legend.background = element_blank(),
        plot.margin = margin(4, 10, 4, 4))

h_est <- 0.22 * length(imp$stations) + 1.2
ggsave(file.path(DIR_FIG, "05_station_estimates_interval.png"), p_est,
       width = 6.2, height = h_est, dpi = FIG_DPI, limitsize = FALSE)
ggsave(file.path(DIR_FIG, "05_station_estimates_interval.pdf"), p_est,
       width = 6.2, height = h_est, limitsize = FALSE)

# The same estimates plotted against one another, so that agreement between the
# analyses is read directly rather than inferred from two panels.
cmp <- data.frame(
  Station = imp$stations,
  dm = dm$site_pred$lambda_avg, dm_lo = dm$site_pred$lambda_avg_lower,
  dm_hi = dm$site_pred$lambda_avg_upper,
  rn = rn$site_pred$lambda_avg, rn_lo = rn$site_pred$lambda_avg_lower,
  rn_hi = rn$site_pred$lambda_avg_upper,
  individuals = dm$site_pred$individuals, stringsAsFactors = FALSE)
lim <- c(0, max(cmp$dm_hi, cmp$rn_hi))

p_cmp <- ggplot(cmp, aes(x = dm, y = rn)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = "grey45", linewidth = 0.4) +
  geom_linerange(aes(xmin = dm_lo, xmax = dm_hi), colour = PAL[[METHODS[1]]]$line,
                 alpha = 0.35, linewidth = 0.45) +
  geom_linerange(aes(ymin = rn_lo, ymax = rn_hi), colour = PAL[[METHODS[2]]]$line,
                 alpha = 0.35, linewidth = 0.45) +
  geom_point(aes(size = individuals), shape = 21, fill = "grey25",
             colour = "white", stroke = 0.6, alpha = 0.9) +
  scale_size_continuous(range = c(1.4, 5), name = "Individuals") +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(x = bquote(hat(lambda)~.(METHODS[1])),
       y = bquote(hat(lambda)~.(METHODS[2]))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        legend.position = "right", legend.background = element_blank(),
        plot.margin = margin(4, 8, 4, 4))

ggsave(file.path(DIR_FIG, "05_station_estimates_agreement.png"), p_cmp,
       width = 5.4, height = 4.4, dpi = FIG_DPI, limitsize = FALSE)
ggsave(file.path(DIR_FIG, "05_station_estimates_agreement.pdf"), p_cmp,
       width = 5.4, height = 4.4, limitsize = FALSE)

cat("  05_station_estimates_interval.(png|pdf)\n")
cat("  05_station_estimates_agreement.(png|pdf)\n")

# --- 5. SESSION ------------------------------------------------------------

si <- sessionInfo()
pk <- sort(vapply(c(si$otherPkgs, si$loadedOnly),
                  function(p) paste(p$Package, p$Version), character(1)))
session_tab <- data.frame(
  item  = c("R version", "platform", "seed", "run date", "occasion length (days)",
            paste("package", seq_along(pk))),
  value = c(si$R.version$version.string, si$platform, imp$seed,
            format(Sys.time()), dm$OCC_LEN, unname(pk)),
  stringsAsFactors = FALSE)
write.csv(session_tab, file.path(DIR_TAB, "05_session_info.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n-- output --\n")
cat("  tables/05_main_top_models.csv\n")
cat("  tables/05_main_coefficients.csv\n")
cat("  tables/05_marginal_curves.csv\n")
cat("  tables/05_station_estimates.csv\n")
cat("  tables/05_session_info.csv\n")
cat("  figures/05_marginal_curves.(png|pdf) and one file per analysis\n")
cat("  figures/05_station_estimates_map.(png|pdf)\n")
cat("  figures/05_station_estimates_map_labelled.(png|pdf)\n")
cat("  figures/05_station_estimates_interval.(png|pdf)\n")
cat("  figures/05_station_estimates_agreement.(png|pdf)\n\n")
cat("STEP 5 complete.\n\n")
