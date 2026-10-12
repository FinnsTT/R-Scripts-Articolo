# 01b_occasion_rn.R  --  STEP 1b: occasion length, Royle-Nichols sensitivity
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# The survey period is divided into occasions of each length on a grid, and at
# every length the null models are fitted:
#
#   1. build the count array and the observation covariates at each length
#   2. fit the null Royle-Nichols model r(~1) lambda(~1)
#   3. fit the null Dail-Madsen model gamma(~1) omega(~1) lambda(~1) p(~1)
#      to the counts, Poisson
#   4. retain the occasion length at which both null models stabilise
#
# Under closure, abundance in the Royle-Nichols model does not depend on the
# occasion length: the per-individual detection probability absorbs the change.
# Counts summed within an occasion instead grow with its length, because the
# same individual can be recorded more than once.
#
# Decision rule, fixed before the run: the occasion length retained is the
# shortest on the grid at which abundance of the Royle-Nichols null and initial
# abundance of the Dail-Madsen null both change by less than STAB_TOL compared
# with the next length. Every model at every length is fitted as an
# independent task, in parallel.
#
# Input : workspace/01_import.rds
# Output: tables/01b_*.csv, figures/01b_*, workspace/01b_occasion.rds
# ===========================================================================

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_OUT <- "C:/Users/Dell/Desktop/codici articolo"
DIR_TAB <- file.path(ROOT_OUT, "tables")
DIR_FIG <- file.path(ROOT_OUT, "figures")
DIR_WS  <- file.path(ROOT_OUT, "workspace")
for (d in c(DIR_TAB, DIR_FIG, DIR_WS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

OCC_GRID <- c(2, 3, 4, 5, 6, 7)         # occasion lengths, days
K_SENS   <- 50                          # integration bound, wide enough to avoid truncation
STAB_TOL <- 0.05                        # relative change accepted between lengths
FIG_DPI  <- 600

USE_PARALLEL <- TRUE
N_CORES      <- max(1L, parallel::detectCores() - 1L)

suppressPackageStartupMessages({
  library(unmarked)
  library(ggplot2)
})

imp <- readRDS(file.path(DIR_WS, "01_import.rds"))
set.seed(imp$seed)

stations     <- imp$stations
active_dates <- imp$active_dates
dets         <- imp$dets_indep
human        <- imp$human
moon_daily   <- imp$moon_daily
FULL_START   <- imp$FULL_START
FULL_END     <- imp$FULL_END
nS           <- length(stations)

SITE_COVS <- imp$site[match(stations, imp$site$Station), imp$COVS, drop = FALSE]
rownames(SITE_COVS) <- NULL

cat("\n=========================================================\n")
cat(" STEP 1b -- occasion length, Royle-Nichols sensitivity\n")
cat(" sites = ", nS, " | lengths = ", paste(OCC_GRID, collapse = ", "),
    " d | K = ", K_SENS, "\n", sep = "")
cat("=========================================================\n\n")

# --- HELPERS ---------------------------------------------------------------

run_quiet <- function(expr, max_messages = 20L) {
  msgs <- character(0); n <- 0L
  v <- withCallingHandlers(try(expr, silent = TRUE),
                           warning = function(cond) {
                             n <<- n + 1L
                             if (length(msgs) < max_messages) {
                               m <- conditionMessage(cond)
                               if (!(m %in% msgs)) msgs <<- c(msgs, m)
                             }
                             invokeRestart("muffleWarning")
                           })
  list(value = v, warnings = msgs, n_warnings = n)
}

ok   <- function(x) !inherits(x, "try-error") && !is.null(x)
msg1 <- function(x) paste(trimws(as.character(x)), collapse = " ")

usable <- function(fit) {
  b <- tryCatch(as.numeric(coef(fit)), error = function(e) NA_real_)
  s <- tryCatch(as.numeric(SE(fit)),   error = function(e) NA_real_)
  length(b) == length(s) && all(is.finite(b)) && all(is.finite(s)) && all(s > 0)
}

resp <- function(fit, type, link = c("log", "logit")) {
  link <- match.arg(link)
  est <- try(as.numeric(coef(fit, type = type)), silent = TRUE)
  se  <- try(as.numeric(SE(fit, type = type)),   silent = TRUE)
  na4 <- c(est = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_)
  if (inherits(est, "try-error") || inherits(se, "try-error") ||
      length(est) != 1 || length(se) != 1 ||
      !is.finite(est) || !is.finite(se)) return(na4)
  ci <- est + c(-1.959964, 1.959964) * se
  g  <- if (link == "log") exp else plogis
  dg <- if (link == "log") exp(est) else exp(est) / (1 + exp(est))^2
  c(est = unname(g(est)), se = unname(abs(dg) * se),
    lower = unname(g(ci[1])), upper = unname(g(ci[2])))
}

# Point estimate on the response scale, available even when the Hessian is
# singular and standard errors cannot be computed.
point <- function(fit, type, link = c("log", "logit")) {
  link <- match.arg(link)
  b <- tryCatch(as.numeric(coef(fit, type = type)), error = function(e) NA_real_)
  if (length(b) != 1) return(NA_real_)
  if (link == "log") exp(b) else plogis(b)
}

# Count array and observation covariates at one occasion length. The count of
# an occasion is the sum of the individuals recorded in its independent
# detections; occasions with no active day are structural missing values and
# are coded NA. Survey effort is the proportion of the occasion during which
# the camera was recording. Moon illumination is the mean daily fraction over
# the occasion, common to every station. Human activity is the number of
# records of people and vehicles at the station during the occasion,
# standardised over the observed cells.
build_data <- function(occ_len) {
  starts <- seq(FULL_START, FULL_END, by = occ_len)
  ends   <- pmin(starts + occ_len - 1, FULL_END)
  nP <- length(starts)
  occ_days <- as.integer(ends - starts) + 1L

  eff <- matrix(0L, nS, nP, dimnames = list(stations, paste0("p", seq_len(nP))))
  for (i in seq_len(nS)) {
    ad <- active_dates[[stations[i]]]
    if (!length(ad)) next
    for (k in seq_len(nP)) eff[i, k] <- sum(ad >= starts[k] & ad <= ends[k])
  }

  occ_index <- function(dates) {
    o <- rep(NA_integer_, length(dates))
    for (k in seq_len(nP)) o[dates >= starts[k] & dates <= ends[k]] <- k
    o
  }

  y <- matrix(0L, nS, nP, dimnames = dimnames(eff))
  tb <- tapply(dets$Count,
               list(factor(dets$Station, levels = stations),
                    factor(occ_index(dets$Date), levels = seq_len(nP))), sum)
  tb[is.na(tb)] <- 0
  y[] <- as.integer(tb)
  y[eff == 0] <- NA

  h <- matrix(0L, nS, nP, dimnames = dimnames(eff))
  h[] <- as.integer(table(factor(human$Station, levels = stations),
                          factor(occ_index(human$Date), levels = seq_len(nP))))
  h_obs <- h[eff > 0]
  h_center <- mean(h_obs); h_scale <- stats::sd(h_obs)
  hz <- (h - h_center) / h_scale
  hz[eff == 0] <- NA

  moon_occ <- vapply(seq_len(nP), function(k)
    mean(moon_daily$fraction[moon_daily$Date >= starts[k] & moon_daily$Date <= ends[k]]),
    numeric(1))
  moon <- matrix(moon_occ, nS, nP, byrow = TRUE, dimnames = dimnames(eff))

  list(y = y, eff_days = eff, effort = sweep(eff, 2, occ_days, "/"),
       moon = moon, human = hz, human_raw = h,
       human_center = h_center, human_scale = h_scale,
       n_occ = nP, starts = starts, ends = ends,
       n_det = sum(y, na.rm = TRUE), n_cells = sum(!is.na(y)))
}

# --- 1-3. NULL MODELS AT EACH OCCASION LENGTH ------------------------------
# Every model at every length is an independent task. Tasks are distributed
# with load balancing, the open models at the shortest lengths first because
# they are the slowest to fit.

t0 <- Sys.time()
data_by_length <- setNames(lapply(OCC_GRID, build_data), as.character(OCC_GRID))

fit_rn_null <- function(L) {
  d <- data_by_length[[as.character(L)]]
  y_bin <- d$y
  y_bin[!is.na(y_bin)] <- as.integer(y_bin[!is.na(y_bin)] > 0)

  umf <- unmarkedFrameOccu(y = y_bin, siteCovs = SITE_COVS)
  r   <- run_quiet(occuRN(~1 ~1, data = umf, K = K_SENS))
  lam <- if (ok(r$value)) resp(r$value, "state", "log") else
    c(est = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_)
  rr  <- if (ok(r$value)) resp(r$value, "det", "logit") else
    c(est = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_)
  data.frame(
    occasion_days = L, n_occ = d$n_occ, n_cells = d$n_cells,
    cells_detection = sum(y_bin == 1, na.rm = TRUE),
    cells_count_2plus = sum(d$y >= 2, na.rm = TRUE),
    converged = ok(r$value) && usable(r$value),
    lambda_null = lam[["est"]], lambda_null_SE = lam[["se"]],
    lambda_null_lower = lam[["lower"]], lambda_null_upper = lam[["upper"]],
    r_null = rr[["est"]], r_null_SE = rr[["se"]],
    n_warnings = r$n_warnings, note = paste(r$warnings, collapse = " | "),
    stringsAsFactors = FALSE)
}

fit_dm_null <- function(L) {
  d <- data_by_length[[as.character(L)]]
  umf_o <- unmarkedFramePCO(y = d$y, siteCovs = SITE_COVS, numPrimary = d$n_occ)
  ro  <- run_quiet(pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1,
                              pformula = ~1, data = umf_o, K = K_SENS, mixture = "P",
                              dynamics = "constant"))
  lo  <- if (ok(ro$value)) resp(ro$value, "lambda", "log") else
    c(est = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_)
  data.frame(
    occasion_days = L, family = "Poisson",
    converged = ok(ro$value) && usable(ro$value),
    lambda = if (ok(ro$value)) point(ro$value, "lambda", "log") else NA_real_,
    lambda_SE = lo[["se"]], lambda_lower = lo[["lower"]], lambda_upper = lo[["upper"]],
    gamma = if (ok(ro$value)) point(ro$value, "gamma", "log") else NA_real_,
    omega = if (ok(ro$value)) point(ro$value, "omega", "logit") else NA_real_,
    p = if (ok(ro$value)) point(ro$value, "det", "logit") else NA_real_,
    n_warnings = ro$n_warnings, note = paste(ro$warnings, collapse = " | "),
    stringsAsFactors = FALSE)
}

tasks <- c(lapply(sort(OCC_GRID), function(L) list(model = "dm", L = L)),
           lapply(sort(OCC_GRID), function(L) list(model = "rn", L = L)))
run_task <- function(tk) if (tk$model == "dm") fit_dm_null(tk$L) else fit_rn_null(tk$L)

if (USE_PARALLEL && N_CORES >= 2L) {
  cl <- parallel::makeCluster(min(N_CORES, length(tasks)))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(unmarked)))
  parallel::clusterExport(cl, c("fit_rn_null", "fit_dm_null", "run_quiet", "ok",
                                "usable", "resp", "point", "data_by_length",
                                "SITE_COVS", "K_SENS"), envir = .GlobalEnv)
  res <- parallel::parLapplyLB(cl, tasks, run_task)
  parallel::stopCluster(cl)
} else {
  res <- lapply(tasks, run_task)
}

is_dm <- vapply(tasks, function(tk) tk$model == "dm", logical(1))
occ_summary <- do.call(rbind, res[!is_dm])
occ_summary <- occ_summary[order(occ_summary$occasion_days), ]
dm_tab      <- do.call(rbind, res[is_dm])
dm_tab      <- dm_tab[order(dm_tab$occasion_days), ]

cat("-- null models by occasion length --\n")
print(occ_summary[, c("occasion_days", "n_occ", "cells_detection", "cells_count_2plus",
                      "lambda_null", "lambda_null_SE", "r_null")],
      row.names = FALSE, digits = 4)
cat("\n")
print(dm_tab[, c("occasion_days", "lambda", "lambda_SE", "gamma", "omega", "p")],
      row.names = FALSE, digits = 4)
cat("\n  elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    " min\n\n", sep = "")

# --- 4. STABILITY AND CHOICE OF THE OCCASION LENGTH ------------------------

n <- nrow(occ_summary)
stab <- data.frame(
  occasion_days      = occ_summary$occasion_days[-n],
  next_occasion_days = occ_summary$occasion_days[-1],
  rn_lambda          = occ_summary$lambda_null[-n],
  rn_lambda_next     = occ_summary$lambda_null[-1],
  dm_lambda          = dm_tab$lambda[-n],
  dm_lambda_next     = dm_tab$lambda[-1],
  stringsAsFactors = FALSE)
stab$rn_rel_change <- abs(stab$rn_lambda - stab$rn_lambda_next) / stab$rn_lambda_next
stab$dm_rel_change <- abs(stab$dm_lambda - stab$dm_lambda_next) / stab$dm_lambda_next
stab$rn_stable <- is.finite(stab$rn_rel_change) & stab$rn_rel_change < STAB_TOL
stab$dm_stable <- is.finite(stab$dm_rel_change) & stab$dm_rel_change < STAB_TOL
stab$stable    <- stab$rn_stable & stab$dm_stable

cat("-- stability between consecutive occasion lengths (tolerance ",
    100 * STAB_TOL, "%) --\n", sep = "")
print(stab, row.names = FALSE, digits = 4)

if (!any(stab$stable))
  stop("the null estimates do not stabilise over the occasion grid; inspect tables/01b_occasion_stability.csv",
       call. = FALSE)

OCC_LEN <- stab$occasion_days[stab$stable][1]
OCC_LEN <- 5L
cat("\n  RETAINED occasion length : ", OCC_LEN, " d\n\n", sep = "")

# --- 5. FIGURE -------------------------------------------------------------

MOD_LINE <- c(`Royle-Nichols` = "#762A83", `Dail-Madsen, Poisson` = "#1B7837")

fig_df <- rbind(
  data.frame(model = "Royle-Nichols", occasion_days = occ_summary$occasion_days,
             lambda = occ_summary$lambda_null,
             lower = occ_summary$lambda_null_lower,
             upper = occ_summary$lambda_null_upper, stringsAsFactors = FALSE),
  data.frame(model = "Dail-Madsen, Poisson",
             occasion_days = dm_tab$occasion_days, lambda = dm_tab$lambda,
             lower = dm_tab$lambda_lower, upper = dm_tab$lambda_upper,
             stringsAsFactors = FALSE))
fig_df$model <- factor(fig_df$model, levels = names(MOD_LINE))

p_occ <- ggplot(fig_df, aes(x = occasion_days, y = lambda, colour = model)) +
  geom_vline(xintercept = OCC_LEN, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = 0.35), linewidth = 0.55) +
  geom_line(position = position_dodge(width = 0.35), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.35), shape = 21, fill = "white",
             size = 2.2, stroke = 0.9) +
  scale_colour_manual(values = MOD_LINE) +
  scale_x_continuous(breaks = OCC_GRID) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Occasion length (days)",
       y = expression(hat(lambda)~"(individuals per station)"), colour = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        legend.position = "top", legend.background = element_blank(),
        plot.margin = margin(4, 8, 4, 4))

ggsave(file.path(DIR_FIG, "01b_occasion_lambda.png"), p_occ, width = 5.8, height = 3.8, dpi = FIG_DPI)
ggsave(file.path(DIR_FIG, "01b_occasion_lambda.pdf"), p_occ, width = 5.8, height = 3.8)

# --- 6. EXPORT -------------------------------------------------------------

write.csv(occ_summary, file.path(DIR_TAB, "01b_occasion_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(dm_tab, file.path(DIR_TAB, "01b_occasion_dm_null.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stab, file.path(DIR_TAB, "01b_occasion_stability.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(
  data_by_length = data_by_length, OCC_GRID = OCC_GRID, OCC_LEN = OCC_LEN,
  occ_summary = occ_summary, dm_tab = dm_tab, stab = stab,
  K_SENS = K_SENS, STAB_TOL = STAB_TOL,
  seed = imp$seed, run_date = Sys.time(), session = sessionInfo()
), file.path(DIR_WS, "01b_occasion.rds"))

cat("-- output --\n")
cat("  tables/01b_occasion_summary.csv\n")
cat("  tables/01b_occasion_dm_null.csv\n")
cat("  tables/01b_occasion_stability.csv\n")
cat("  figures/01b_occasion_lambda.(png|pdf)\n")
cat("  workspace/01b_occasion.rds   <- input for STEP 2\n\n")
cat("STEP 1b complete. Retained occasion length: ", OCC_LEN, " d\n\n", sep = "")
