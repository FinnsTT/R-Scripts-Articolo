# 01_import_vif.R  --  STEP 1: data import and collinearity screening
#
# Brown hyena (Parahyaena brunnea) N-mixture analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
#   1. read the raw data files
#   2. reconcile the station set across files
#   3. reconstruct per-station active days (survey effort) from the metadata
#   4. apply the temporal-independence filter to the detection records
#   5. assemble the records and the daily values behind the observation
#      covariates: human activity and moon illumination
#   6. assemble the eight site covariates over the survey period
#   7. screen the covariates for collinearity with the variance inflation
#      factor, removing covariates stepwise until every VIF is below the gate
#   8. assemble the candidate sets of the study, for detection and for
#      abundance; abundance models containing a removed covariate are dropped
#   9. export tables and save the prepared objects
#
# Input : raw data files
# Output: tables/01_*.csv, workspace/01_import.rds

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_DATA <- "C:/Users/Dell/Desktop/codici tesi/file GIS"
ROOT_OUT  <- "C:/Users/Dell/Desktop/codici articolo"

DIR_TAB <- file.path(ROOT_OUT, "tables")
DIR_WS  <- file.path(ROOT_OUT, "workspace")
for (d in c(DIR_TAB, DIR_WS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEED            <- 20260722
INDEP_MIN       <- 15      # minutes between successive retained records at a station
VIF_THRESHOLD   <- 5       # a priori gate
DEPLOY_SESSIONS <- 10:21   # sampling sessions of the maximum-coverage period
HUMAN_CLASSES   <- c("Homo sapiens", "Vehicles")   # records counted as human activity

# Reference occasion length used to assemble the data object for this step.
# The VIF is a property of the covariate matrix, and none of the eight
# covariates depends on occasion boundaries. The occasion length used for
# modelling is set in STEP 1b.
OCC_REF <- 5L

suppressPackageStartupMessages({
  library(lubridate)
  library(suncalc)
})

set.seed(SEED)

cat("\n=========================================================\n")
cat(" STEP 1 -- data import and collinearity screening\n")
cat(" seed = ", SEED, " | independence threshold = ", INDEP_MIN, " min\n", sep = "")
cat("=========================================================\n\n")

# --- HELPERS ---------------------------------------------------------------

# Dates stored as ISO text; empty cells are missing values.
parse_date <- function(x) {
  x <- trimws(as.character(x))
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  ok  <- !is.na(x) & !x %in% c("", "NA", "/")
  if (any(ok)) out[ok] <- suppressWarnings(lubridate::ymd(x[ok]))
  out
}

# Set of days on which each camera was recording, taken from the deployment
# sessions and with any documented pause removed. A pause whose start and end
# dates coincide records no interruption and is ignored.
build_active_dates <- function(meta, sessions) {
  res <- setNames(vector("list", nrow(meta)), meta$station)
  for (i in seq_along(res)) res[[i]] <- as.Date(character(0))
  used <- integer(0); skipped <- integer(0)
  for (s in sessions) {
    sc <- paste0("S", s, "_start"); ec <- paste0("S", s, "_end")
    if (!all(c(sc, ec) %in% names(meta))) { skipped <- c(skipped, s); next }
    used <- c(used, s)
    st <- parse_date(meta[[sc]]); en <- parse_date(meta[[ec]])
    pause_from <- if (paste0("P", s, "S_from") %in% names(meta))
      parse_date(meta[[paste0("P", s, "S_from")]]) else rep(as.Date(NA), nrow(meta))
    pause_to   <- if (paste0("P", s, "E_to") %in% names(meta))
      parse_date(meta[[paste0("P", s, "E_to")]])   else rep(as.Date(NA), nrow(meta))
    for (i in seq_len(nrow(meta))) {
      if (is.na(st[i]) || is.na(en[i]) || en[i] < st[i]) next
      d <- seq(st[i], en[i], by = "day")
      if (!is.na(pause_from[i]) && !is.na(pause_to[i]) && pause_from[i] < pause_to[i])
        d <- d[!(d >= pause_from[i] & d <= pause_to[i])]
      res[[i]] <- unique(c(res[[i]], d))
    }
  }
  attr(res, "sessions_used") <- used
  attr(res, "sessions_skipped") <- skipped
  res
}

# Two records at the same station are treated as independent when at least
# INDEP_MIN minutes separate them. Each record is compared with the previous
# RETAINED record, and the filter runs over the full ordered series at a
# station, so the independent set does not depend on the occasion grid.
filter_independent <- function(dets, indep_min) {
  d <- dets[order(dets$Station, dets$DateTimeOriginal), , drop = FALSE]
  keep <- logical(nrow(d))
  for (s in unique(d$Station)) {
    idx <- which(d$Station == s)
    tt  <- d$DateTimeOriginal[idx]
    k   <- logical(length(idx))
    last <- NA
    for (j in seq_along(idx)) {
      if (is.na(last) ||
          as.numeric(difftime(tt[j], last, units = "mins")) >= indep_min) {
        k[j] <- TRUE; last <- tt[j]
      }
    }
    keep[idx] <- k
  }
  d[keep, , drop = FALSE]
}

# Occasion grid over a date range.
build_periods <- function(start, end, occ_len) {
  s <- seq(start, end, by = occ_len)
  data.frame(period_id = seq_along(s), start = s,
             end = pmin(s + occ_len - 1, end))
}

# Variance inflation factors computed from the covariate matrix alone, as
# 1 / (1 - R2_j) with R2_j from the regression of covariate j on the others.
vif_from_matrix <- function(X) {
  vapply(names(X), function(j) {
    1 / (1 - summary(stats::lm(X[[j]] ~ ., data = X[setdiff(names(X), j)]))$r.squared)
  }, numeric(1))
}

# --- 1. RAW DATA -----------------------------------------------------------

covs <- read.csv(file.path(ROOT_DATA, "hyena_analysis_complete_dataset.csv"),
                 stringsAsFactors = FALSE)
dets <- read.csv(file.path(ROOT_DATA, "DATA-15min_Sampling1-21_counts_Parahyena_brunnea.csv"),
                 stringsAsFactors = FALSE)
recs <- read.csv(file.path(ROOT_DATA, "DATA-15min_Sampling1-21.csv"),
                 stringsAsFactors = FALSE)
meta <- read.csv(file.path(ROOT_DATA, "sampling_effort.csv"), sep = ";",
                 stringsAsFactors = FALSE, colClasses = "character")
# rows holding station totals are not cameras
meta <- meta[!is.na(meta$station) & nzchar(trimws(meta$station)) &
               !grepl("^TOT", trimws(meta$station), ignore.case = TRUE), ]
meta$station <- trimws(meta$station)

dets$Date <- as.Date(dets$Date)
dets$DateTimeOriginal <- as.POSIXct(dets$DateTimeOriginal, tz = "UTC")
dets$session <- as.integer(sub("^S", "", dets$Sampling))

human <- recs[recs$Species %in% HUMAN_CLASSES, ]
human$Date <- as.Date(human$Date)
human$session <- as.integer(sub("^S", "", human$Sampling))

cat("-- input files --\n")
cat("  hyena_analysis_complete_dataset.csv :", nrow(covs), "stations\n")
cat("  sampling_effort.csv                 :", nrow(meta), "cameras\n")
cat("  detections (all sessions)           :", nrow(dets), "records,",
    sum(dets$Count), "individuals\n")

dets  <- dets[dets$session %in% DEPLOY_SESSIONS, ]
human <- human[human$session %in% DEPLOY_SESSIONS, ]
cat("  detections (sessions ", min(DEPLOY_SESSIONS), "-", max(DEPLOY_SESSIONS),
    ")          : ", nrow(dets), " records, ", sum(dets$Count), " individuals\n",
    sep = "")
cat("  human activity records (same sessions) :", nrow(human), "\n\n")

# --- 2. STATION RECONCILIATION ---------------------------------------------
# A station of the covariate file that is absent from the metadata but present
# there under lettered variants (e.g. FR31A, FR31B) is represented by the
# variant with the most active days; its records are relabelled to the station
# name and the records of the other variants are discarded. A station then
# enters the analysis when it appears in every required file and has at least
# one active day.

active_dates <- build_active_dates(meta, DEPLOY_SESSIONS)

variant_tab <- NULL
for (s in setdiff(covs$Station, meta$station)) {
  var_ids <- meta$station[grepl(paste0("^", s, "[A-Za-z]+$"), meta$station)]
  if (!length(var_ids)) next
  n_days <- vapply(var_ids, function(v) length(active_dates[[v]]), integer(1))
  chosen <- var_ids[which.max(n_days)]
  variant_tab <- rbind(variant_tab, data.frame(
    Station = s, variant = var_ids, active_days = as.integer(n_days),
    retained = var_ids == chosen, stringsAsFactors = FALSE))
  drop <- setdiff(var_ids, chosen)
  meta <- meta[!meta$station %in% drop, ]
  meta$station[meta$station == chosen] <- s
  dets  <- dets[!dets$Station %in% drop, ];   dets$Station[dets$Station == chosen]   <- s
  human <- human[!human$Station %in% drop, ]; human$Station[human$Station == chosen] <- s
}
if (!is.null(variant_tab)) active_dates <- build_active_dates(meta, DEPLOY_SESSIONS)

id_all <- sort(unique(c(covs$Station, meta$station, dets$Station)))
n_active <- vapply(id_all, function(s)
  if (s %in% names(active_dates)) length(active_dates[[s]]) else 0L, integer(1))

stations_tab <- data.frame(
  Station        = id_all,
  in_covariates  = id_all %in% covs$Station,
  in_metadata    = id_all %in% meta$station,
  active_days    = as.integer(n_active),
  n_records      = as.integer(table(factor(dets$Station, levels = id_all))),
  stringsAsFactors = FALSE
)
stations_tab$analysed <- with(stations_tab,
  in_covariates & in_metadata & active_days > 0)

missing_from <- apply(
  stations_tab[c("in_covariates", "in_metadata")], 1,
  function(r) paste(c("covariates", "metadata")[!r], collapse = " + "))
stations_tab$reason_excluded <- ""
sel <- !stations_tab$analysed
stations_tab$reason_excluded[sel] <- trimws(paste(
  ifelse(missing_from[sel] == "", "", paste("absent from:", missing_from[sel])),
  ifelse(stations_tab$active_days[sel] == 0, "no active days", "")))

stations <- stations_tab$Station[stations_tab$analysed]
nS <- length(stations)

cat("-- station reconciliation --\n")
if (!is.null(variant_tab)) {
  cat("  stations represented by a lettered variant:\n")
  for (i in seq_len(nrow(variant_tab)))
    cat("    ", variant_tab$Station[i], " <- ", variant_tab$variant[i], "  (",
        variant_tab$active_days[i], " active days",
        if (variant_tab$retained[i]) ", retained)\n" else ", discarded)\n", sep = "")
}
cat("  IDs present in at least one file :", nrow(stations_tab), "\n")
cat("  stations analysed                :", nS, "\n")
if (any(sel)) {
  cat("  excluded:\n")
  ex <- stations_tab[sel, ]
  for (i in seq_len(nrow(ex)))
    cat("    ", ex$Station[i], " -- ", ex$reason_excluded[i],
        "  (records lost: ", ex$n_records[i], ")\n", sep = "")
}
cat("\n")
stopifnot(nS > 0)

# --- 3. SURVEY PERIOD, EFFORT AND DETECTIONS -------------------------------

FULL_START <- min(dets$Date)
FULL_END   <- max(dets$Date)

cat("-- survey period --\n")
cat("  from", format(FULL_START), "to", format(FULL_END),
    " (", as.integer(FULL_END - FULL_START) + 1, "days )\n")
cat("  deployment sessions used :",
    paste(attr(active_dates, "sessions_used"), collapse = ", "), "\n")
if (length(attr(active_dates, "sessions_skipped")))
  cat("  sessions skipped (no end-date column) :",
      paste(attr(active_dates, "sessions_skipped"), collapse = ", "), "\n")

dets_analysed <- dets[dets$Station %in% stations &
                        dets$Date >= FULL_START & dets$Date <= FULL_END, ]
dets_indep <- filter_independent(dets_analysed, INDEP_MIN)

human <- human[human$Station %in% stations &
                 human$Date >= FULL_START & human$Date <= FULL_END, ]

cat("  detection records at analysed stations :", nrow(dets_analysed), "\n")
cat("  independent detections                 :", nrow(dets_indep),
    sprintf(" (%d removed by the >= %d min filter)\n",
            nrow(dets_analysed) - nrow(dets_indep), INDEP_MIN))
cat("  individuals in independent detections  :", sum(dets_indep$Count), "\n")
cat("  stations with at least one detection    :",
    length(unique(dets_indep$Station)), "of", nS, "\n")
cat("  human activity records in the period   :", nrow(human), "\n\n")

# Effort and counts on the reference occasion grid. The count of an occasion is
# the sum of the individuals recorded in its independent detections. Occasions
# with no active day are structural missing values, not observed zeros, and
# are coded NA.
per <- build_periods(FULL_START, FULL_END, OCC_REF)
nP  <- nrow(per)

eff_days <- matrix(0L, nS, nP, dimnames = list(stations, paste0("p", per$period_id)))
for (i in seq_len(nS)) {
  ad <- active_dates[[stations[i]]]
  if (!length(ad)) next
  for (k in seq_len(nP))
    eff_days[i, k] <- sum(ad >= per$start[k] & ad <= per$end[k])
}

period_of <- rep(NA_integer_, nrow(dets_indep))
for (k in seq_len(nP))
  period_of[dets_indep$Date >= per$start[k] & dets_indep$Date <= per$end[k]] <- k
y <- matrix(0L, nS, nP, dimnames = dimnames(eff_days))
tb <- tapply(dets_indep$Count,
             list(factor(dets_indep$Station, levels = stations),
                  factor(period_of, levels = seq_len(nP))), sum)
tb[is.na(tb)] <- 0
y[] <- as.integer(tb)

n_det_zero_effort <- sum(y[eff_days == 0])
y[eff_days == 0] <- NA

cat("-- count array (reference occasion ", OCC_REF, " d) --\n", sep = "")
cat("  sites x occasions  :", nS, "x", nP, "\n")
cat("  total trap-nights  :", sum(eff_days), "\n")
cat("  zero-effort cells  :", sum(is.na(y)), "of", nS * nP, "-> coded NA\n")
if (n_det_zero_effort > 0)
  cat("  WARNING: ", n_det_zero_effort,
      " individuals fall in zero-effort cells and are dropped.\n",
      "   Deployment metadata and photo timestamps disagree for these records.\n", sep = "")
cat("\n")

# --- 4. OBSERVATION COVARIATES: DAILY MOON ILLUMINATION --------------------
# Fraction of the lunar disk illuminated on each day of the survey period
# (suncalc). The occasion value is the mean over the days of the occasion and
# is assembled with the occasion grid in STEP 1b.

survey_days <- seq(FULL_START, FULL_END, by = "day")
moon_daily <- data.frame(
  Date = survey_days,
  fraction = suncalc::getMoonIllumination(date = survey_days,
                                          keep = "fraction")$fraction)

# --- 5. SITE COVARIATES ----------------------------------------------------

ci <- match(stations, covs$Station)

site_raw <- data.frame(
  Station              = stations,
  elevation            = covs$elevation[ci],
  tri                  = covs$tri[ci],
  tpi                  = covs$tpi[ci],
  dist_water_points    = covs$dist_water_points[ci],
  dist_fish_river      = covs$dist_fish_river[ci],
  dist_konkiep_river   = covs$dist_konkiep_river[ci],
  dist_tourist_road    = covs$dist_tourist_road[ci],
  dist_staff_buildings = covs$dist_staff_buildings[ci],
  stringsAsFactors = FALSE
)
rownames(site_raw) <- stations

COVS <- setdiff(names(site_raw), "Station")
if (anyNA(site_raw[COVS]))
  stop("missing covariate values for: ",
       paste(site_raw$Station[!complete.cases(site_raw[COVS])], collapse = ", "))

# z-scores, with the centring and scaling constants kept so predictions can be
# returned on the original scale
center <- vapply(site_raw[COVS], mean, numeric(1))
scale_ <- vapply(site_raw[COVS], stats::sd, numeric(1))
site <- site_raw
for (v in COVS) site[[v]] <- (site_raw[[v]] - center[[v]]) / scale_[[v]]

# --- 6. COLLINEARITY SCREENING ---------------------------------------------
# Stepwise (Zuur et al. 2010): the covariate with the largest VIF is removed
# and the VIF is recomputed on the remaining ones, until every value is below
# the gate.

vif_first <- vif_from_matrix(site[COVS])
kept <- COVS
removal_order <- character(0)
while (length(kept) >= 2) {
  v <- vif_from_matrix(site[kept])
  if (max(v) < VIF_THRESHOLD) break
  worst <- names(v)[which.max(v)]
  removal_order <- c(removal_order, worst)
  kept <- setdiff(kept, worst)
}
vif_final <- if (length(kept) >= 2) vif_from_matrix(site[kept]) else setNames(1, kept)

vif_tab <- data.frame(
  covariate   = COVS,
  description = c("elevation (DEM, m a.s.l.)",
                  "terrain ruggedness index",
                  "topographic position index",
                  "distance to nearest artificial water point (m)",
                  "distance to Fish River (m)",
                  "distance to Konkiep River (m)",
                  "distance to tourist road (m)",
                  "distance to staff buildings (m)"),
  mean_raw = vapply(site_raw[COVS], mean, numeric(1)),
  sd_raw   = vapply(site_raw[COVS], stats::sd, numeric(1)),
  min_raw  = vapply(site_raw[COVS], min, numeric(1)),
  max_raw  = vapply(site_raw[COVS], max, numeric(1)),
  VIF_all_covariates = as.numeric(vif_first),
  removal_step = match(COVS, removal_order),
  VIF_final = as.numeric(vif_final[COVS]),
  row.names = NULL, stringsAsFactors = FALSE
)
vif_tab$status <- ifelse(vif_tab$covariate %in% kept, "retained", "EXCLUDED")

cat("-- collinearity screening (VIF, gate at ", VIF_THRESHOLD, ", stepwise) --\n", sep = "")
print(vif_tab[, c("covariate", "VIF_all_covariates", "removal_step", "VIF_final", "status")],
      row.names = FALSE, digits = 3)
if (length(removal_order)) {
  cat("\n  removed, in order :", paste(removal_order, collapse = ", "), "\n")
} else {
  cat("\n  all covariates below the threshold; none excluded.\n")
}
cat("\n")

COVS_RETAINED <- kept

# --- 7. CANDIDATE SETS -----------------------------------------------------
# Detection: the five structures compared in the study, settled on the
# abundance null before the abundance hypotheses are examined. Survey effort is
# the proportion of the occasion during which the camera was recording; moon
# illumination is the mean fraction of the lunar disk illuminated over the
# occasion, on its original 0-1 scale; human activity is the number of records
# of people and vehicles at the station during the occasion, standardised.
#
# Abundance: the nineteen models of the study. Model names follow the study.

det_set <- data.frame(
  id = c("null", "effort", "effort_moon", "effort_human", "full"),
  formula = c("~1", "~effort", "~effort + moon", "~effort + human",
              "~effort + moon + human"),
  stringsAsFactors = FALSE)

hyp <- data.frame(
  id = paste0("M", 1:19),
  hypothesis = c(
    "intercept only",
    "abundance varies with elevation, through landscape structure",
    "abundance increases with terrain ruggedness, through shelter and denning sites",
    "abundance varies with topographic position, canyon floor against plateau",
    "abundance increases near artificial water points, through resource predictability",
    "abundance increases near the Fish River, through resource predictability",
    "abundance increases near the Konkiep River, through resource predictability",
    "abundance decreases near the tourist road, through human disturbance",
    "abundance decreases near staff buildings, through human disturbance",
    "elevation and terrain ruggedness",
    "elevation and topographic position",
    "distance to the Fish River and to the Konkiep River",
    "elevation, terrain ruggedness and distance to water points",
    "elevation, terrain ruggedness and distance to the Fish River",
    "elevation, terrain ruggedness and distance to the Konkiep River",
    "elevation, terrain ruggedness and distance to both rivers",
    "elevation, terrain ruggedness, Konkiep River and tourist road",
    "elevation, terrain ruggedness, Konkiep River and staff buildings",
    "the effect of terrain ruggedness depends on distance to water points"),
  lambda = c("~1",
             "~elevation",
             "~tri",
             "~tpi",
             "~dist_water_points",
             "~dist_fish_river",
             "~dist_konkiep_river",
             "~dist_tourist_road",
             "~dist_staff_buildings",
             "~elevation + tri",
             "~elevation + tpi",
             "~dist_fish_river + dist_konkiep_river",
             "~elevation + tri + dist_water_points",
             "~elevation + tri + dist_fish_river",
             "~elevation + tri + dist_konkiep_river",
             "~elevation + tri + dist_fish_river + dist_konkiep_river",
             "~elevation + tri + dist_konkiep_river + dist_tourist_road",
             "~elevation + tri + dist_konkiep_river + dist_staff_buildings",
             "~tri * dist_water_points"),
  prediction = c("--", "+/-", "+", "+/-", "-", "-", "-", "+", "+",
                 "--", "--", "--", "--", "--", "--", "--", "--", "--", "--"),
  stringsAsFactors = FALSE)

hyp$vif_ok <- vapply(hyp$lambda, function(f)
  all(all.vars(stats::as.formula(f)) %in% COVS_RETAINED), logical(1))
hyp$excluded_covariates <- vapply(hyp$lambda, function(f)
  paste(setdiff(all.vars(stats::as.formula(f)), COVS_RETAINED), collapse = ", "),
  character(1))

cat("-- candidate set --\n")
cat("  detection structures :", nrow(det_set), "\n")
cat("  abundance models     :", nrow(hyp), "of which", sum(hyp$vif_ok),
    "contain only retained covariates\n")
if (any(!hyp$vif_ok)) {
  cat("  dropped:\n")
  for (i in which(!hyp$vif_ok))
    cat("    ", hyp$id[i], " ", hyp$lambda[i], "  (", hyp$excluded_covariates[i], ")\n",
        sep = "")
}
cat("\n")

# --- 8. EXPORT -------------------------------------------------------------

effort_tab <- data.frame(
  Station                = stations,
  active_days            = vapply(stations, function(s) length(active_dates[[s]]), integer(1)),
  trap_nights            = as.integer(rowSums(eff_days)),
  independent_detections = as.integer(table(factor(dets_indep$Station, levels = stations))),
  individuals            = as.integer(rowSums(y, na.rm = TRUE)),
  human_records          = as.integer(table(factor(human$Station, levels = stations))),
  row.names = NULL)

write.csv(stations_tab, file.path(DIR_TAB, "01_stations.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
if (!is.null(variant_tab))
  write.csv(variant_tab, file.path(DIR_TAB, "01_station_variants.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
write.csv(vif_tab, file.path(DIR_TAB, "01_vif.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(site_raw, file.path(DIR_TAB, "01_covariates_raw.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(effort_tab, file.path(DIR_TAB, "01_effort_by_station.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hyp, file.path(DIR_TAB, "01_candidate_set.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(
  stations = stations, stations_tab = stations_tab, variant_tab = variant_tab,
  active_dates = active_dates,
  dets_all = dets, dets_indep = dets_indep, human = human,
  moon_daily = moon_daily,
  FULL_START = FULL_START, FULL_END = FULL_END,
  covs_all = covs,
  site = site, site_raw = site_raw, center = center, scale = scale_,
  COVS = COVS, COVS_RETAINED = COVS_RETAINED,
  vif_tab = vif_tab, VIF_THRESHOLD = VIF_THRESHOLD,
  hyp = hyp, det_set = det_set,
  OCC_REF = OCC_REF, INDEP_MIN = INDEP_MIN,
  DEPLOY_SESSIONS = DEPLOY_SESSIONS, HUMAN_CLASSES = HUMAN_CLASSES,
  seed = SEED, run_date = Sys.time(), session = sessionInfo()
), file.path(DIR_WS, "01_import.rds"))

cat("-- output --\n")
cat("  tables/01_stations.csv\n")
if (!is.null(variant_tab)) cat("  tables/01_station_variants.csv\n")
cat("  tables/01_vif.csv\n")
cat("  tables/01_covariates_raw.csv\n")
cat("  tables/01_effort_by_station.csv\n")
cat("  tables/01_candidate_set.csv\n")
cat("  workspace/01_import.rds   <- input for STEP 1b\n\n")
cat("STEP 1 complete.\n\n")
