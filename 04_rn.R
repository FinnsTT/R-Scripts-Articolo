# 04_rn.R  --  STEP 4: Royle-Nichols model, sensitivity and model selection
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# The Royle-Nichols model relates the probability of detecting the species at a
# station to the per-individual detection rate and to local abundance, and is
# fitted to detection and non-detection data. Repeated records of the same
# individual within an occasion therefore cannot inflate the response, which is
# the reason for fitting it alongside the open population model.
#
#   1. build the detection and non-detection array at the occasion length
#      retained in STEP 1b
#   2. sensitivity of the null model r(~1) lambda(~1) to the integration bound
#   3. choice of the integration bound K
#   4. read the candidate set of the study
#   5. settle the detection submodel on the abundance null, among the five
#      structures of the study, before the hypotheses are examined, and hold
#      it fixed thereafter
#   6. fit every candidate, rank it by AICc, verify the integration bound
#      against the largest abundance predicted at a station by the models of
#      the top set, and average over the models within two AICc units
#   7. predict abundance at every station
#   8. state the population assumption and where it is examined
#
# Latent abundance in occuRN is Poisson; the distribution is therefore not a
# choice at this step, and only the integration bound is selected.
#
# Input : workspace/01_import.rds, workspace/02_dm_null.rds
# Output: tables/04_*.csv, figures/04_*, workspace/04_rn.rds
# ===========================================================================

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_OUT <- "C:/Users/Dell/Desktop/codici articolo"

DIR_TAB <- file.path(ROOT_OUT, "tables")
DIR_FIG <- file.path(ROOT_OUT, "figures")
DIR_WS  <- file.path(ROOT_OUT, "workspace")
for (d in c(DIR_TAB, DIR_FIG, DIR_WS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

K_GRID     <- seq(10, 50, by = 5)
TAIL_MAX   <- 1e-10
STABLE_TOL <- 0.01
DRIFT_MAX  <- 1.10
DELTA_SET  <- 2
FIG_DPI    <- 600

suppressPackageStartupMessages({
  library(unmarked)
  library(ggplot2)
})

imp <- readRDS(file.path(DIR_WS, "01_import.rds"))
nul <- readRDS(file.path(DIR_WS, "02_dm_null.rds"))
set.seed(imp$seed)

d         <- nul$data
SITE_COVS <- nul$SITE_COVS
OCC_LEN   <- nul$OCC_LEN
nS        <- nrow(SITE_COVS)
n_occ     <- d$n_occ

cat("\n=========================================================\n")
cat(" STEP 4 -- Royle-Nichols model\n")
cat(" sites = ", nS, " | occasions = ", n_occ,
    " | occasion = ", OCC_LEN, " d\n", sep = "")
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

# Bind two model tables that do not share every column, keeping the columns
# they have in common.
bind_common <- function(a, b) {
  k <- intersect(names(a), names(b))
  rbind(a[, k, drop = FALSE], b[, k, drop = FALSE])
}

aicc <- function(ll, npar, n) -2 * ll + 2 * npar + 2 * npar * (npar + 1) / (n - npar - 1)

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

# --- 1. DETECTION AND NON-DETECTION ARRAY ----------------------------------
# Counts are reduced to presence and absence. Structural missing values are
# carried over unchanged.

y_bin <- d$y
y_bin[!is.na(y_bin)] <- as.integer(y_bin[!is.na(y_bin)] > 0)

n_collapsed <- sum(d$y >= 2, na.rm = TRUE)

cat("-- detection and non-detection array --\n")
cat("  sites x occasions   :", nS, "x", n_occ, "\n")
cat("  observed cells      :", sum(!is.na(y_bin)), "\n")
cat("  cells with detection:", sum(y_bin == 1, na.rm = TRUE), "\n")
cat("  cells reduced from a count of two or more :", n_collapsed,
    sprintf(" (%.2f%% of observed cells)\n", 100 * n_collapsed / sum(!is.na(y_bin))))
cat("  stations with at least one detection :",
    sum(rowSums(y_bin, na.rm = TRUE) > 0), "of", nS, "\n\n")

umf <- unmarkedFrameOccu(y = y_bin, siteCovs = SITE_COVS,
                         obsCovs = list(effort = d$effort, moon = d$moon,
                                        human = d$human))

# --- 2. SENSITIVITY OF THE NULL MODEL TO K ---------------------------------

cat("-- sensitivity of the null model to the integration bound --\n")
cat("  ", length(K_GRID), " fits\n\n", sep = "")

rows <- list()
for (K in K_GRID) {
  cat(sprintf("  K = %2d ... ", K)); utils::flush.console()
  row <- tryCatch({
    t1 <- Sys.time()
    r  <- run_quiet(occuRN(~1 ~1, data = umf, K = K))
    secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (!ok(r$value)) {
      cat(sprintf("failed (%.1f s) -- %s\n", secs, msg1(r$value)))
      data.frame(K = K, converged = FALSE, lambda = NA_real_, lambda_SE = NA_real_,
                 r_det = NA_real_, r_det_SE = NA_real_, logLik = NA_real_,
                 AIC = NA_real_, CV_lambda = NA_real_, tail_prob = NA_real_,
                 seconds = secs, n_warnings = r$n_warnings,
                 note = msg1(r$value), stringsAsFactors = FALSE)
    } else {
      lam <- resp(r$value, "state", "log")
      rr  <- resp(r$value, "det",   "logit")
      cat(sprintf("lambda = %7.4f   r = %6.4f   %5.1f s\n",
                  lam[["est"]], rr[["est"]], secs))
      data.frame(K = K,
                 converged = is.finite(lam[["est"]]) && is.finite(lam[["se"]]),
                 lambda = lam[["est"]], lambda_SE = lam[["se"]],
                 r_det = rr[["est"]], r_det_SE = rr[["se"]],
                 logLik = as.numeric(logLik(r$value)),
                 AIC = as.numeric(r$value@AIC),
                 CV_lambda = lam[["se"]] / lam[["est"]],
                 tail_prob = stats::ppois(K, lam[["est"]], lower.tail = FALSE),
                 seconds = secs, n_warnings = r$n_warnings,
                 note = paste(r$warnings, collapse = " | "), stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    cat("ERROR -- ", conditionMessage(e), "\n", sep = "")
    data.frame(K = K, converged = FALSE, lambda = NA_real_, lambda_SE = NA_real_,
               r_det = NA_real_, r_det_SE = NA_real_, logLik = NA_real_,
               AIC = NA_real_, CV_lambda = NA_real_, tail_prob = NA_real_,
               seconds = NA_real_, n_warnings = 0L,
               note = paste("error:", conditionMessage(e)), stringsAsFactors = FALSE)
  })
  rows[[length(rows) + 1]] <- row
}
sens <- do.call(rbind, rows)
cat("\n  converged ", sum(sens$converged), "/", nrow(sens), "\n\n", sep = "")

# --- 3. CHOICE OF THE INTEGRATION BOUND ------------------------------------

sk <- sens[sens$converged, ]
if (!nrow(sk)) stop("the null model did not converge at any value of K.", call. = FALSE)
sk <- sk[order(sk$K), ]
lambda_ref   <- sk$lambda[nrow(sk)]
sk$rel_diff  <- abs(sk$lambda - lambda_ref) / lambda_ref
sk$tail_ok   <- is.finite(sk$tail_prob) & sk$tail_prob < TAIL_MAX
sk$stable_ok <- sk$rel_diff < STABLE_TOL
sk$eligible  <- sk$tail_ok & sk$stable_ok

drift <- max(sk$lambda) / min(sk$lambda)

cat("-- K: truncation and stability --\n")
print(sk[, c("K", "lambda", "lambda_SE", "r_det", "tail_prob", "rel_diff",
             "tail_ok", "stable_ok", "eligible")], row.names = FALSE, digits = 4)
cat(sprintf("\n  drift of lambda across K : %.4f (identifiable below %.2f: %s)\n",
            drift, DRIFT_MAX, ifelse(drift < DRIFT_MAX, "yes", "NO")))

if (!any(sk$eligible))
  stop("no value of K satisfies both the tail and the stability rule; extend K_GRID.",
       call. = FALSE)

K_RETAINED <- min(sk$K[sk$eligible])
cat("  RETAINED K on the null model : ", K_RETAINED, "\n\n", sep = "")

# --- 4. CANDIDATE SET ------------------------------------------------------
# The models of the study, identical to those of the count-based analysis so
# that the two can be read against one another. Models containing a covariate
# removed by the VIF screen of STEP 1 are not fitted and are listed with the
# reason.

hyp     <- imp$hyp
det_set <- imp$det_set

fit_rn <- function(lam_txt, p_txt, K) {
  f <- stats::as.formula(paste(p_txt, lam_txt))   # occuRN takes ~detection ~state
  run_quiet(occuRN(f, data = umf, K = K))
}

fit_one <- function(id, hypothesis, lam_txt, p_txt, prediction, K) {
  cat(sprintf("  %-15s lambda %-28s r %-22s ... ", id, lam_txt, p_txt))
  utils::flush.console()
  t1 <- Sys.time()
  r  <- fit_rn(lam_txt, p_txt, K)
  secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  blank <- data.frame(
    model = id, hypothesis = hypothesis, lambda = lam_txt, r_formula = p_txt,
    prediction = prediction, converged = FALSE, nPars = NA_integer_,
    logLik = NA_real_, AICc = NA_real_, seconds = secs,
    n_warnings = r$n_warnings, note = NA_character_, stringsAsFactors = FALSE)
  if (!ok(r$value)) {
    cat(sprintf("failed (%.0f s)\n", secs)); blank$note <- msg1(r$value)
    return(list(row = blank, fit = NULL))
  }
  if (!usable(r$value)) {
    cat(sprintf("non-finite standard errors (%.0f s)\n", secs))
    blank$note   <- "non-finite estimates or standard errors"
    blank$nPars  <- length(coef(r$value))
    blank$logLik <- as.numeric(logLik(r$value))
    return(list(row = blank, fit = NULL))
  }
  np <- length(coef(r$value)); ll <- as.numeric(logLik(r$value))
  ac <- aicc(ll, np, nS)
  cat(sprintf("AICc = %8.2f  (%.0f s)\n", ac, secs))
  list(row = data.frame(
    model = id, hypothesis = hypothesis, lambda = lam_txt, r_formula = p_txt,
    prediction = prediction, converged = TRUE, nPars = np, logLik = ll,
    AICc = ac, seconds = secs, n_warnings = r$n_warnings,
    note = paste(r$warnings, collapse = " | "), stringsAsFactors = FALSE),
    fit = r$value)
}

pred_one <- function(fit, newdata) {
  p <- try(predict(fit, type = "state", newdata = newdata, appendData = FALSE),
           silent = TRUE)
  if (inherits(p, "try-error")) return(NULL)
  data.frame(fit = as.numeric(p$Predicted), se = as.numeric(p$SE))
}

# --- 5. DETECTION STRUCTURE, CANDIDATE SET AND INTEGRATION BOUND -----------
# The detection submodel is settled on the abundance null before the
# hypotheses are examined, as in the count-based analysis. In this model the
# detection parameter is the per-individual rate: the probability of detecting
# the species at a station is one minus the probability that none of the
# individuals present is detected. The bound retained on the null model is
# then verified against the largest abundance predicted at a station by the
# models of the top set, and the analysis is refitted at a larger bound when
# the truncated Poisson tail exceeds the tolerance.

select_at <- function(K) {
  cat("-- detection structure, on the abundance null, K = ", K, " --\n", sep = "")
  rows <- list(); fits <- list()
  for (j in seq_len(nrow(det_set))) {
    res <- fit_one(paste0("null_r_", det_set$id[j]),
                   "detection structure, abundance null", "~1", det_set$formula[j],
                   "--", K)
    rows[[length(rows) + 1]] <- res$row
    if (!is.null(res$fit)) fits[[res$row$model]] <- res$fit
  }
  det_tab <- do.call(rbind, rows)
  det_tab$delta_AICc <- det_tab$AICc - min(det_tab$AICc, na.rm = TRUE)
  det_tab$weight <- ifelse(det_tab$converged, exp(-0.5 * det_tab$delta_AICc), NA_real_)
  det_tab$weight <- det_tab$weight / sum(det_tab$weight, na.rm = TRUE)

  print(det_tab[, c("model", "r_formula", "nPars", "logLik", "AICc", "delta_AICc",
                    "weight")], row.names = FALSE, digits = 4)

  if (!any(det_tab$converged))
    stop("no detection structure produced usable estimates.", call. = FALSE)
  P_FORM <- det_tab$r_formula[det_tab$converged][which.min(det_tab$AICc[det_tab$converged])]

  cat("\n  RETAINED detection submodel: r ", P_FORM, "\n", sep = "")
  cat("  it is held fixed across the candidate set below\n\n")

  cat("-- fitting the candidate set, ", sum(hyp$vif_ok), " models, r ", P_FORM,
      " --\n", sep = "")
  rows <- list()
  t0 <- Sys.time()
  for (i in seq_len(nrow(hyp))) {
    if (!hyp$vif_ok[i]) {
      rows[[length(rows) + 1]] <- data.frame(
        model = hyp$id[i], hypothesis = hyp$hypothesis[i], lambda = hyp$lambda[i],
        r_formula = P_FORM, prediction = hyp$prediction[i], converged = FALSE,
        nPars = NA_integer_, logLik = NA_real_, AICc = NA_real_, seconds = NA_real_,
        n_warnings = 0L,
        note = paste("not fitted: contains a covariate removed by the VIF screen:",
                     hyp$excluded_covariates[i]),
        stringsAsFactors = FALSE)
      next
    }
    res <- fit_one(hyp$id[i], hyp$hypothesis[i], hyp$lambda[i], P_FORM,
                   hyp$prediction[i], K)
    rows[[length(rows) + 1]] <- res$row
    if (!is.null(res$fit)) fits[[hyp$id[i]]] <- res$fit
  }
  cat("  elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
      " min\n\n", sep = "")

  aic_tab <- do.call(rbind, rows)
  rank_tab <- aic_tab[aic_tab$converged, ]
  if (!nrow(rank_tab)) stop("no candidate model produced usable estimates.", call. = FALSE)
  rank_tab <- rank_tab[order(rank_tab$AICc), ]
  rank_tab$delta_AICc <- rank_tab$AICc - min(rank_tab$AICc)
  rank_tab$weight <- exp(-0.5 * rank_tab$delta_AICc)
  rank_tab$weight <- rank_tab$weight / sum(rank_tab$weight)
  rank_tab$in_top_set <- rank_tab$delta_AICc < DELTA_SET

  top <- rank_tab$model[rank_tab$in_top_set]
  lam_max <- vapply(top, function(m) {
    p <- pred_one(fits[[m]], SITE_COVS); if (is.null(p)) NA_real_ else max(p$fit)
  }, numeric(1))
  tails <- stats::ppois(K, lam_max, lower.tail = FALSE)

  list(K = K, det_tab = det_tab, P_FORM = P_FORM, aic_tab = aic_tab,
       rank_tab = rank_tab, fits = fits,
       tail_check = data.frame(K = K, model = top, lambda_max_station = lam_max,
                               tail_prob = tails,
                               tail_ok = is.finite(tails) & tails < TAIL_MAX,
                               stringsAsFactors = FALSE))
}

K <- K_RETAINED
tail_log <- NULL
repeat {
  sel <- select_at(K)
  tail_log <- rbind(tail_log, sel$tail_check)
  cat("-- truncation at the largest station-level abundance of the top set --\n")
  print(sel$tail_check, row.names = FALSE, digits = 4)
  if (all(sel$tail_check$tail_ok)) { cat("\n  K = ", K, " retained\n\n", sep = ""); break }
  K_next <- K_GRID[K_GRID > K &
                     vapply(K_GRID, function(k)
                       all(stats::ppois(k, sel$tail_check$lambda_max_station,
                                        lower.tail = FALSE) < TAIL_MAX), logical(1))]
  if (!length(K_next))
    stop("no value of K on the grid bounds the largest predicted abundance; extend K_GRID.",
         call. = FALSE)
  K <- min(K_next)
  cat("\n  truncation exceeds the tolerance; refitting at K = ", K, "\n\n", sep = "")
}

det_tab  <- sel$det_tab
P_FORM   <- sel$P_FORM
aic_tab  <- sel$aic_tab
rank_tab <- sel$rank_tab
fits     <- sel$fits
eff_tab  <- det_tab          # kept for the supplementary table of every fit
K_USED   <- K

det_choice <- data.frame(
  retained = P_FORM,
  detection = det_tab$r_formula, AICc = det_tab$AICc, delta_AICc = det_tab$delta_AICc,
  weight = det_tab$weight, stringsAsFactors = FALSE)

# --- 6. RANKING AND MODEL AVERAGING ----------------------------------------

failed <- aic_tab[!aic_tab$converged, ]

cat("-- ranking --\n")
print(rank_tab[, c("model", "hypothesis", "lambda", "nPars", "logLik",
                   "AICc", "delta_AICc", "weight")], row.names = FALSE, digits = 4)
if (nrow(failed)) {
  cat("\n  excluded from inference:\n")
  for (i in seq_len(nrow(failed)))
    cat("    ", failed$model[i], " -- ", failed$note[i], "\n", sep = "")
}
cat("\n")

top_set <- rank_tab$model[rank_tab$in_top_set]
w_top   <- rank_tab$weight[rank_tab$in_top_set]
w_top   <- w_top / sum(w_top)
names(w_top) <- top_set

cat("-- top model set (delta AICc < ", DELTA_SET, ") --\n", sep = "")
cat("  ", length(top_set), " models: ", paste(top_set, collapse = ", "), "\n\n", sep = "")

coef_of <- function(fit) {
  b <- coef(fit); s <- SE(fit)
  data.frame(term = names(b), estimate = as.numeric(b), se = as.numeric(s),
             stringsAsFactors = FALSE)
}

best_id   <- rank_tab$model[1]
best_fit  <- fits[[best_id]]
cb        <- coef_of(best_fit)
terms_all <- unique(unlist(lapply(top_set, function(m) coef_of(fits[[m]])$term)))

avg_tab <- do.call(rbind, lapply(terms_all, function(tm) {
  est <- vapply(top_set, function(m) {
    c0 <- coef_of(fits[[m]]); v <- c0$estimate[c0$term == tm]; if (length(v)) v else 0
  }, numeric(1))
  se <- vapply(top_set, function(m) {
    c0 <- coef_of(fits[[m]]); v <- c0$se[c0$term == tm]; if (length(v)) v else 0
  }, numeric(1))
  present <- vapply(top_set, function(m) tm %in% coef_of(fits[[m]])$term, logical(1))
  bar <- sum(w_top * est)
  unc <- sum(w_top * sqrt(se^2 + (est - bar)^2))
  data.frame(term = tm, avg_estimate = bar, avg_uncond_SE = unc,
             avg_lower = bar - 1.959964 * unc, avg_upper = bar + 1.959964 * unc,
             weight_present = sum(w_top[present]), n_models_present = sum(present),
             stringsAsFactors = FALSE)
}))

coef_tab <- merge(data.frame(term = terms_all, stringsAsFactors = FALSE),
                  cb, by = "term", all.x = TRUE)
names(coef_tab)[names(coef_tab) == "estimate"] <- "best_estimate"
names(coef_tab)[names(coef_tab) == "se"]       <- "best_SE"
coef_tab$best_lower <- coef_tab$best_estimate - 1.959964 * coef_tab$best_SE
coef_tab$best_upper <- coef_tab$best_estimate + 1.959964 * coef_tab$best_SE
coef_tab <- merge(coef_tab, avg_tab, by = "term", all = TRUE)
coef_tab$submodel <- sub("\\(.*$", "", coef_tab$term)
coef_tab <- coef_tab[order(coef_tab$submodel, coef_tab$term), ]

cat("-- coefficients: best model (", best_id, ") and average over the top set --\n", sep = "")
print(coef_tab[, c("term", "best_estimate", "best_SE",
                   "avg_estimate", "avg_uncond_SE", "weight_present")],
      row.names = FALSE, digits = 4)
cat("\n")

# --- 7. PREDICTED ABUNDANCE PER STATION ------------------------------------

pred_list <- lapply(top_set, function(m) pred_one(fits[[m]], SITE_COVS))
names(pred_list) <- top_set
keep <- !vapply(pred_list, is.null, logical(1))
if (!any(keep)) stop("prediction failed for every model in the top set.", call. = FALSE)
w_use <- w_top[keep] / sum(w_top[keep])

fit_mat <- vapply(pred_list[keep], function(z) z$fit, numeric(nS))
se_mat  <- vapply(pred_list[keep], function(z) z$se,  numeric(nS))
if (is.null(dim(fit_mat))) { fit_mat <- matrix(fit_mat, nS, 1); se_mat <- matrix(se_mat, nS, 1) }

lam_avg <- as.numeric(fit_mat %*% w_use)
se_avg  <- as.numeric(sqrt(se_mat^2 + (fit_mat - lam_avg)^2) %*% w_use)
pb      <- pred_one(best_fit, SITE_COVS)

site_pred <- data.frame(
  Station = imp$stations,
  lambda_best = pb$fit, lambda_best_SE = pb$se,
  lambda_avg = lam_avg, lambda_avg_SE = se_avg,
  lambda_avg_lower = pmax(lam_avg - 1.959964 * se_avg, 0),
  lambda_avg_upper = lam_avg + 1.959964 * se_avg,
  occasions_with_detection = as.integer(rowSums(y_bin, na.rm = TRUE)),
  trap_nights = as.integer(rowSums(d$eff_days)),
  stringsAsFactors = FALSE)

cat("-- predicted abundance per station --\n")
cat(sprintf("  model-averaged: min %.3f   median %.3f   mean %.3f   max %.3f\n",
            min(lam_avg), stats::median(lam_avg), mean(lam_avg), max(lam_avg)))
cat(sprintf("  summed over the %d stations: %.2f\n\n", nS, sum(lam_avg)))

# --- 8. NOTE ON THE POPULATION ASSUMPTION ----------------------------------
# The Royle-Nichols model assumes that abundance at a station is fixed across
# the occasions. The assumption is examined on the null models of STEP 2
# (tables/02_closure_test.csv). Where closure is rejected, the estimate above is
# not abundance at one point in time but the number of individuals that used
# the station over the whole period.

# --- 9. FIGURES ------------------------------------------------------------

RN_LINE <- "#762A83"
RN_BAND <- "#E3D3E8"

fig_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        legend.background = element_blank(), legend.title = element_blank(),
        plot.margin = margin(4, 8, 4, 4))

p_k <- ggplot(sens[sens$converged, ], aes(x = K, y = lambda)) +
  geom_ribbon(aes(ymin = pmax(lambda - 1.96 * lambda_SE, 0),
                  ymax = lambda + 1.96 * lambda_SE),
              fill = RN_BAND, colour = NA) +
  geom_line(colour = RN_LINE, linewidth = 1) +
  geom_point(colour = RN_LINE, fill = "white", shape = 21, size = 2.2, stroke = 0.9) +
  scale_x_continuous(breaks = K_GRID) +
  labs(x = "Integration bound K",
       y = expression(hat(lambda)~"(individuals per station)")) +
  fig_theme

wt <- rank_tab
wt$label <- factor(sprintf("%-4s %s", wt$model, gsub("^~", "\u03bb ~ ", wt$lambda)),
                   levels = rev(sprintf("%-4s %s", wt$model,
                                        gsub("^~", "\u03bb ~ ", wt$lambda))))
p_w <- ggplot(wt, aes(x = weight, y = label, fill = in_top_set)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.25) +
  geom_text(aes(label = formatC(weight, format = "f", digits = 3)),
            hjust = -0.18, size = 2.9, colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = RN_LINE, `FALSE` = "#D9D9D9")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = expression(paste("AIC"[c], " weight")), y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        axis.text.y = element_text(family = "mono", size = 8),
        legend.position = "none",
        plot.margin = margin(4, 10, 4, 4))

ggsave(file.path(DIR_FIG, "04_rn_K_lambda.png"), p_k, width = 5.8, height = 3.6, dpi = FIG_DPI)
ggsave(file.path(DIR_FIG, "04_rn_K_lambda.pdf"), p_k, width = 5.8, height = 3.6)
ggsave(file.path(DIR_FIG, "04_rn_aicc_weights.png"), p_w,
       width = 7.2, height = 0.3 * nrow(wt) + 1, dpi = FIG_DPI, limitsize = FALSE)
ggsave(file.path(DIR_FIG, "04_rn_aicc_weights.pdf"), p_w,
       width = 7.2, height = 0.3 * nrow(wt) + 1, limitsize = FALSE)

# --- 10. EXPORT ------------------------------------------------------------

write.csv(sens, file.path(DIR_TAB, "04_rn_K_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sk[, c("K", "lambda", "lambda_SE", "tail_prob", "rel_diff",
                 "tail_ok", "stable_ok", "eligible")],
          file.path(DIR_TAB, "04_rn_K_selection.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_common(aic_tab, eff_tab), file.path(DIR_TAB, "04_rn_candidates_all.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(rank_tab, file.path(DIR_TAB, "04_rn_aicc_ranking.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(det_choice, file.path(DIR_TAB, "04_rn_detection_structure.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tail_log, file.path(DIR_TAB, "04_rn_K_truncation.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(coef_tab, file.path(DIR_TAB, "04_rn_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(site_pred, file.path(DIR_TAB, "04_rn_site_predictions.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(
  y_bin = y_bin, umf = umf, n_collapsed = n_collapsed,
  sens = sens, k_tab = sk, K_RETAINED = K_RETAINED, K = K_USED, drift = drift,
  tail_log = tail_log,
  hyp = hyp, det_set = det_set, aic_tab = aic_tab, eff_tab = eff_tab,
  det_tab = det_tab, det_choice = det_choice, P_FORM = P_FORM,
  rank_tab = rank_tab, coef_tab = coef_tab, site_pred = site_pred,
  fits = fits, top_set = top_set, weights_top = w_top, best_id = best_id,
  DELTA_SET = DELTA_SET,
  OCC_LEN = OCC_LEN, n_occ = n_occ, nS = nS, SITE_COVS = SITE_COVS,
  seed = imp$seed, run_date = Sys.time(), session = sessionInfo()
), file.path(DIR_WS, "04_rn.rds"))

cat("-- output --\n")
cat("  tables/04_rn_K_sensitivity.csv\n")
cat("  tables/04_rn_K_selection.csv\n")
cat("  tables/04_rn_candidates_all.csv\n")
cat("  tables/04_rn_aicc_ranking.csv\n")
cat("  tables/04_rn_detection_structure.csv\n")
cat("  tables/04_rn_K_truncation.csv\n")
cat("  tables/04_rn_coefficients.csv\n")
cat("  tables/04_rn_site_predictions.csv\n")
cat("  figures/04_rn_K_lambda.(png|pdf)\n")
cat("  figures/04_rn_aicc_weights.(png|pdf)\n")
cat("  workspace/04_rn.rds   <- input for STEP 5\n\n")
cat("STEP 4 complete. K = ", K_USED, " | best model: ", best_id, "\n\n", sep = "")
