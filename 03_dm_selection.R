# 03_dm_selection.R  --  STEP 3: model selection, count-based model
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# The candidate set is the set of ecological hypotheses of the study, fitted
# with the framework, distribution, integration bound and occasion length
# retained in STEP 2: the Dail-Madsen open population model when recruitment
# and apparent survival are supported, the closed N-mixture model otherwise.
# In the open model recruitment and apparent survival are held constant; the
# hypotheses act on initial abundance.
#
#   1. read the candidate set of the study
#   2. settle the detection submodel on the abundance null, among the five
#      structures of the study, before the hypotheses are examined, and hold
#      it fixed thereafter
#   3. fit every candidate and rank it by AICc
#   4. verify the integration bound against the largest abundance predicted at
#      a station by the models of the top set, and refit at a larger bound
#      when the truncated tail exceeds the tolerance
#   5. average over the models within two AICc units, with unconditional
#      standard errors
#   6. open framework: report recruitment and apparent survival, and test
#      whether the dynamic parameters are supported
#   7. predict abundance at every station
#
# Averaging follows Buckland et al. (1997): estimates are weighted by the
# renormalised AICc weights of the models within two units, a coefficient
# absent from a model contributes zero, and the unconditional standard error
# adds the between-model component.
#
# Input : workspace/01_import.rds, workspace/02_dm_null.rds
# Output: tables/03_*.csv, figures/03_*, workspace/03_dm_selection.rds
# ===========================================================================

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_OUT <- "C:/Users/Dell/Desktop/codici articolo"

DIR_TAB <- file.path(ROOT_OUT, "tables")
DIR_FIG <- file.path(ROOT_OUT, "figures")
DIR_WS  <- file.path(ROOT_OUT, "workspace")
for (d in c(DIR_TAB, DIR_FIG, DIR_WS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

DELTA_SET <- 2       # AICc units defining the top model set
FIG_DPI   <- 600

suppressPackageStartupMessages({
  library(unmarked)
  library(ggplot2)
})

imp <- readRDS(file.path(DIR_WS, "01_import.rds"))
nul <- readRDS(file.path(DIR_WS, "02_dm_null.rds"))
set.seed(imp$seed)

FRAMEWORK   <- nul$FRAMEWORK
MODEL_LABEL <- nul$MODEL_LABEL
PFX         <- nul$PFX
OPEN        <- FRAMEWORK == "open"
umf         <- nul$umf
SITE_COVS   <- nul$SITE_COVS
MIXTURE     <- nul$MIXTURE
FAMILY      <- nul$FAMILY
OCC_LEN     <- nul$OCC_LEN
K_GRID      <- nul$K_GRID
TAIL_MAX    <- nul$TAIL_MAX
nS          <- nrow(SITE_COVS)
n_occ       <- nul$data$n_occ
LAM_TYPE    <- if (OPEN) "lambda" else "state"

hyp     <- imp$hyp
det_set <- imp$det_set

cat("\n=========================================================\n")
cat(" STEP 3 -- ", MODEL_LABEL, " model selection\n", sep = "")
cat(" distribution = ", FAMILY, " | K = ", nul$K_RETAINED, " | occasion = ", OCC_LEN,
    " d\n", sep = "")
cat(" sites = ", nS, " | occasions = ", n_occ, "\n", sep = "")
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

# Estimate, standard error and confidence interval on the response scale, for
# an intercept-only submodel.
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

# A fit is usable for inference when the optimiser returned finite estimates
# and finite standard errors for every parameter.
usable <- function(fit) {
  b <- tryCatch(as.numeric(coef(fit)), error = function(e) NA_real_)
  s <- tryCatch(as.numeric(SE(fit)),   error = function(e) NA_real_)
  length(b) == length(s) && all(is.finite(b)) && all(is.finite(s)) && all(s > 0)
}

# Upper tail of the abundance distribution beyond the integration bound. In
# unmarked the negative binomial has mean lambda and size exp(alpha).
tail_prob <- function(K, lambda, size = NA_real_) {
  if (is.finite(size)) stats::pnbinom(K, size = size, mu = lambda, lower.tail = FALSE)
  else stats::ppois(K, lambda, lower.tail = FALSE)
}

nb_size <- function(fit) {
  if (MIXTURE != "NB") return(NA_real_)
  exp(tryCatch(as.numeric(coef(fit, type = "alpha")), error = function(e) NA_real_))
}

fit_count <- function(lam_txt, p_txt, K) {
  if (OPEN)
    run_quiet(pcountOpen(lambdaformula = stats::as.formula(lam_txt),
                         gammaformula  = ~1,
                         omegaformula  = ~1,
                         pformula      = stats::as.formula(p_txt),
                         data = umf, K = K, mixture = MIXTURE,
                         dynamics = "constant"))
  else
    run_quiet(pcount(stats::as.formula(paste(p_txt, lam_txt)),   # ~detection ~abundance
                     data = umf, K = K, mixture = MIXTURE))
}

fit_one <- function(id, hypothesis, lam_txt, p_txt, prediction, K) {
  cat(sprintf("  %-15s lambda %-28s p %-22s ... ", id, lam_txt, p_txt))
  utils::flush.console()
  t1 <- Sys.time()
  r  <- fit_count(lam_txt, p_txt, K)
  secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  blank <- data.frame(
    model = id, hypothesis = hypothesis, lambda = lam_txt, p = p_txt,
    prediction = prediction, converged = FALSE, nPars = NA_integer_,
    logLik = NA_real_, AICc = NA_real_, seconds = secs,
    n_warnings = r$n_warnings, note = NA_character_, stringsAsFactors = FALSE)
  if (!ok(r$value)) {
    cat(sprintf("failed (%.0f s)\n", secs))
    blank$note <- msg1(r$value)
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
    model = id, hypothesis = hypothesis, lambda = lam_txt, p = p_txt,
    prediction = prediction, converged = TRUE, nPars = np, logLik = ll,
    AICc = ac, seconds = secs, n_warnings = r$n_warnings,
    note = paste(r$warnings, collapse = " | "), stringsAsFactors = FALSE),
    fit = r$value)
}

ALPHA <- nul$ALPHA

# Correlation between the estimators of recruitment and apparent survival, on
# the scale the model is fitted on.
rho_gamma_omega <- function(fit) {
  V <- try(vcov(fit), silent = TRUE)
  if (inherits(V, "try-error") || is.null(V) || is.null(rownames(V))) return(NA_real_)
  nm <- rownames(V)
  ig <- grep("^gam", nm)[1]; io <- grep("^om", nm)[1]
  if (is.na(ig) || is.na(io)) return(NA_real_)
  r <- V[ig, io] / sqrt(V[ig, ig] * V[io, io])
  if (!is.finite(r)) return(NA_real_)
  max(-1, min(1, r))
}

# Tail probability of the likelihood ratio under the chi-bar-square null. The
# constrained set is the quadrant gamma >= 0 and 1 - omega >= 0, so the
# relevant correlation is minus the correlation between the estimators of
# gamma and omega.
p_chibar <- function(chi, rho) {
  if (!is.finite(chi) || chi <= 0) return(1)
  if (!is.finite(rho)) rho <- 0
  w2 <- acos(max(-1, min(1, rho))) / (2 * pi)
  0.5 * stats::pchisq(chi, 1, lower.tail = FALSE) +
    w2 * stats::pchisq(chi, 2, lower.tail = FALSE)
}

pred_one <- function(fit, newdata) {
  p <- try(predict(fit, type = LAM_TYPE, newdata = newdata, appendData = FALSE),
           silent = TRUE)
  if (inherits(p, "try-error")) return(NULL)
  data.frame(fit = as.numeric(p$Predicted), se = as.numeric(p$SE))
}

# --- 1-4. DETECTION STRUCTURE, CANDIDATE SET AND INTEGRATION BOUND ---------
# The detection submodel is settled before the abundance hypotheses are
# examined, on the abundance null, so that every candidate is fitted with the
# same detection structure and the ranking reflects the abundance hypotheses
# alone. The five structures are those compared in the study.
#
# Models of the study containing a covariate removed by the VIF screen of
# STEP 1 are not fitted and are listed with the reason.

select_at <- function(K) {
  cat("-- detection structure, on the abundance null, K = ", K, " --\n", sep = "")
  rows <- list(); fits <- list()
  for (j in seq_len(nrow(det_set))) {
    res <- fit_one(paste0("null_p_", det_set$id[j]),
                   "detection structure, abundance null", "~1", det_set$formula[j],
                   "--", K)
    rows[[length(rows) + 1]] <- res$row
    if (!is.null(res$fit)) fits[[res$row$model]] <- res$fit
  }
  det_tab <- do.call(rbind, rows)
  det_tab$delta_AICc <- det_tab$AICc - min(det_tab$AICc, na.rm = TRUE)
  det_tab$weight <- ifelse(det_tab$converged, exp(-0.5 * det_tab$delta_AICc), NA_real_)
  det_tab$weight <- det_tab$weight / sum(det_tab$weight, na.rm = TRUE)

  print(det_tab[, c("model", "p", "nPars", "logLik", "AICc", "delta_AICc", "weight")],
        row.names = FALSE, digits = 4)

  if (!any(det_tab$converged))
    stop("no detection structure produced usable estimates.", call. = FALSE)
  P_FORM <- det_tab$p[det_tab$converged][which.min(det_tab$AICc[det_tab$converged])]

  cat("\n  RETAINED detection submodel: p ", P_FORM, "\n", sep = "")
  cat("  it is held fixed across the candidate set below\n\n")

  cat("-- fitting the candidate set, ", sum(hyp$vif_ok), " models, p ", P_FORM,
      " --\n", sep = "")
  rows <- list()
  t0 <- Sys.time()
  for (i in seq_len(nrow(hyp))) {
    if (!hyp$vif_ok[i]) {
      rows[[length(rows) + 1]] <- data.frame(
        model = hyp$id[i], hypothesis = hyp$hypothesis[i], lambda = hyp$lambda[i],
        p = P_FORM, prediction = hyp$prediction[i], converged = FALSE,
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

  # truncation at the largest station-level abundance of the top set
  top <- rank_tab$model[rank_tab$in_top_set]
  tails <- vapply(top, function(m) {
    p <- pred_one(fits[[m]], SITE_COVS)
    if (is.null(p)) return(NA_real_)
    tail_prob(K, max(p$fit), nb_size(fits[[m]]))
  }, numeric(1))
  lam_max <- vapply(top, function(m) {
    p <- pred_one(fits[[m]], SITE_COVS); if (is.null(p)) NA_real_ else max(p$fit)
  }, numeric(1))

  list(K = K, det_tab = det_tab, P_FORM = P_FORM, aic_tab = aic_tab,
       rank_tab = rank_tab, fits = fits,
       tail_check = data.frame(K = K, model = top, lambda_max_station = lam_max,
                               nb_size = vapply(top, function(m) nb_size(fits[[m]]), numeric(1)),
                               tail_prob = tails, tail_ok = is.finite(tails) & tails < TAIL_MAX,
                               stringsAsFactors = FALSE))
}

K <- nul$K_RETAINED
tail_log <- NULL
repeat {
  sel <- select_at(K)
  tail_log <- rbind(tail_log, sel$tail_check)
  cat("-- truncation at the largest station-level abundance of the top set --\n")
  print(sel$tail_check, row.names = FALSE, digits = 4)
  if (all(sel$tail_check$tail_ok)) { cat("\n  K = ", K, " retained\n\n", sep = ""); break }
  K_next <- K_GRID[K_GRID > K]
  K_next <- K_next[vapply(K_next, function(k) all(vapply(seq_len(nrow(sel$tail_check)),
    function(i) isTRUE(tail_prob(k, sel$tail_check$lambda_max_station[i],
                                 sel$tail_check$nb_size[i]) < TAIL_MAX), logical(1))),
    logical(1))]
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
  detection = det_tab$p, AICc = det_tab$AICc, delta_AICc = det_tab$delta_AICc,
  weight = det_tab$weight, stringsAsFactors = FALSE)

# --- 5. RANKING AND MODEL AVERAGING ----------------------------------------

failed <- aic_tab[!aic_tab$converged, ]

cat("-- ranking --\n")
print(rank_tab[, c("model", "hypothesis", "lambda", "nPars", "logLik",
                   "AICc", "delta_AICc", "weight")],
      row.names = FALSE, digits = 4)
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
cat("  ", length(top_set), " models: ", paste(top_set, collapse = ", "), "\n", sep = "")
cat("  renormalised weights: ",
    paste(sprintf("%s = %.3f", top_set, w_top), collapse = ", "), "\n\n", sep = "")

# Coefficients of the best model and of the average over the top set. A
# coefficient absent from a model in the set contributes zero.
coef_of <- function(fit) {
  b <- coef(fit); s <- SE(fit)
  data.frame(term = names(b), estimate = as.numeric(b), se = as.numeric(s),
             stringsAsFactors = FALSE)
}

best_id  <- rank_tab$model[1]
best_fit <- fits[[best_id]]
cb <- coef_of(best_fit)

terms_all <- unique(unlist(lapply(top_set, function(m) coef_of(fits[[m]])$term)))

avg_tab <- do.call(rbind, lapply(terms_all, function(tm) {
  est <- vapply(top_set, function(m) {
    c0 <- coef_of(fits[[m]]); v <- c0$estimate[c0$term == tm]
    if (length(v)) v else 0
  }, numeric(1))
  se <- vapply(top_set, function(m) {
    c0 <- coef_of(fits[[m]]); v <- c0$se[c0$term == tm]
    if (length(v)) v else 0
  }, numeric(1))
  present <- vapply(top_set, function(m) tm %in% coef_of(fits[[m]])$term, logical(1))
  bar <- sum(w_top * est)
  unc <- sum(w_top * sqrt(se^2 + (est - bar)^2))
  data.frame(term = tm,
             avg_estimate = bar, avg_uncond_SE = unc,
             avg_lower = bar - 1.959964 * unc,
             avg_upper = bar + 1.959964 * unc,
             weight_present = sum(w_top[present]),
             n_models_present = sum(present),
             stringsAsFactors = FALSE)
}))

coef_tab <- merge(
  data.frame(term = terms_all, stringsAsFactors = FALSE),
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

# --- 6. DYNAMICS, AND SUPPORT FOR THE OPEN STRUCTURE -----------------------
# Open framework only. Recruitment and apparent survival are constant in every
# candidate model, so they are reported for the best model and averaged over
# the top set. Apparent survival is per occasion; retention over the survey is
# omega raised to the number of transitions. Under the closed framework the
# record of the closure decision of STEP 2 is carried forward.

nT <- n_occ - 1L   # transitions between occasions
dyn_tab <- NULL

if (OPEN) {
  dyn_one <- function(fit) {
    c(gamma = resp(fit, "gamma", "log")[["est"]],
      gamma_SE = resp(fit, "gamma", "log")[["se"]],
      gamma_lower = resp(fit, "gamma", "log")[["lower"]],
      gamma_upper = resp(fit, "gamma", "log")[["upper"]],
      omega = resp(fit, "omega", "logit")[["est"]],
      omega_SE = resp(fit, "omega", "logit")[["se"]],
      omega_lower = resp(fit, "omega", "logit")[["lower"]],
      omega_upper = resp(fit, "omega", "logit")[["upper"]])
  }

  dyn_best <- dyn_one(best_fit)

  dyn_avg <- function(type, link) {
    e <- vapply(top_set, function(m) resp(fits[[m]], type, link)[["est"]], numeric(1))
    s <- vapply(top_set, function(m) resp(fits[[m]], type, link)[["se"]],  numeric(1))
    good <- is.finite(e) & is.finite(s)
    if (!any(good)) return(c(est = NA_real_, se = NA_real_,
                             lower = NA_real_, upper = NA_real_))
    ww  <- w_top[good] / sum(w_top[good])
    bar <- sum(ww * e[good])
    unc <- sum(ww * sqrt(s[good]^2 + (e[good] - bar)^2))
    c(est = bar, se = unc,
      lower = bar - 1.959964 * unc, upper = bar + 1.959964 * unc)
  }

  gam_a <- dyn_avg("gamma", "log")
  omg_a <- dyn_avg("omega", "logit")
  omg_a[["lower"]] <- max(omg_a[["lower"]], 0)
  omg_a[["upper"]] <- min(omg_a[["upper"]], 1)

  dyn_tab <- data.frame(
    parameter = c("gamma (recruitment, individuals per site per occasion)",
                  "omega (apparent survival per occasion)",
                  paste0("retention over the survey, omega^", nT)),
    best_estimate = c(dyn_best[["gamma"]], dyn_best[["omega"]],
                      dyn_best[["omega"]]^nT),
    best_SE = c(dyn_best[["gamma_SE"]], dyn_best[["omega_SE"]],
                nT * dyn_best[["omega"]]^(nT - 1) * dyn_best[["omega_SE"]]),
    best_lower = c(dyn_best[["gamma_lower"]], dyn_best[["omega_lower"]],
                   dyn_best[["omega_lower"]]^nT),
    best_upper = c(dyn_best[["gamma_upper"]], dyn_best[["omega_upper"]],
                   dyn_best[["omega_upper"]]^nT),
    avg_estimate = c(gam_a[["est"]], omg_a[["est"]], omg_a[["est"]]^nT),
    avg_uncond_SE = c(gam_a[["se"]], omg_a[["se"]],
                      nT * omg_a[["est"]]^(nT - 1) * omg_a[["se"]]),
    avg_lower = c(gam_a[["lower"]], omg_a[["lower"]], omg_a[["lower"]]^nT),
    avg_upper = c(gam_a[["upper"]], omg_a[["upper"]], omg_a[["upper"]]^nT),
    stringsAsFactors = FALSE)

  cat("-- recruitment and apparent survival --\n")
  print(dyn_tab[, c("parameter", "best_estimate", "best_lower", "best_upper",
                    "avg_estimate", "avg_lower", "avg_upper")],
        row.names = FALSE, digits = 4)
  cat("\n")

  # Do the two dynamic parameters earn their place? The best-supported open
  # model is compared with the closed model of the same abundance and
  # detection structure, fitted to the same data at the same integration
  # bound. Under the closed model omega = 1 and gamma = 0, both on the boundary
  # of the parameter space, so the likelihood ratio is referred to the
  # chi-bar-square distribution (Self and Liang 1987), with weights following
  # the correlation between the estimators of gamma and omega. The chi-square
  # on two degrees of freedom, which is conservative there, and the difference
  # in AICc are reported alongside. The intervals of gamma and omega are
  # reported as descriptions, not as tests: on the log and logit links the null
  # values lie at minus and plus infinity, so the intervals can never contain
  # them.

  cat("-- support for the open structure --\n")
  cat("  fitting the closed counterpart: lambda ", rank_tab$lambda[1], " p ", P_FORM,
      " ... ", sep = "")
  utils::flush.console()
  rc <- run_quiet(pcount(stats::as.formula(paste(P_FORM, rank_tab$lambda[1])),
                         data = nul$umf_closed, K = K_USED, mixture = MIXTURE))

  if (!ok(rc$value) || !usable(rc$value)) {
    cat("failed\n")
    open_test <- data.frame(
      comparison = "open versus closed, same abundance and detection structure",
      logLik_open = as.numeric(logLik(best_fit)), logLik_closed = NA_real_,
      df = 2L, chi_square = NA_real_, rho_gamma_omega = NA_real_,
      p_mixture = NA_real_, p_chisq2 = NA_real_,
      AICc_open = rank_tab$AICc[1], AICc_closed = NA_real_,
      note = "the closed counterpart did not produce usable estimates",
      stringsAsFactors = FALSE)
  } else {
    ll_o <- as.numeric(logLik(best_fit))
    ll_c <- as.numeric(logLik(rc$value))
    chi  <- 2 * (ll_o - ll_c)
    ac_c <- aicc(ll_c, length(coef(rc$value)), nS)
    rr   <- rho_gamma_omega(best_fit)
    cat(sprintf("AICc = %.2f\n", ac_c))
    open_test <- data.frame(
      comparison = "open versus closed, same abundance and detection structure",
      logLik_open = ll_o, logLik_closed = ll_c, df = 2L,
      chi_square = chi,
      rho_gamma_omega = rr,
      p_mixture = p_chibar(chi, -rr),
      p_chisq2 = stats::pchisq(max(chi, 0), df = 2, lower.tail = FALSE),
      AICc_open = rank_tab$AICc[1], AICc_closed = ac_c,
      note = "boundary hypothesis; p_mixture is the reference, p_chisq2 is conservative",
      stringsAsFactors = FALSE)
    fits[["closed_counterpart"]] <- rc$value
  }

  open_test$omega <- dyn_best[["omega"]]
  open_test$omega_lower <- dyn_best[["omega_lower"]]
  open_test$omega_upper <- dyn_best[["omega_upper"]]
  open_test$gamma <- dyn_best[["gamma"]]
  open_test$gamma_lower <- dyn_best[["gamma_lower"]]
  open_test$gamma_upper <- dyn_best[["gamma_upper"]]
  open_test$delta_AICc <- open_test$AICc_closed - open_test$AICc_open
  open_test$verdict <- with(open_test,
    if (!is.finite(p_mixture))
      "not assessable"
    else if (p_mixture < ALPHA)
      "closure rejected"
    else
      "closure not rejected")

  print(t(open_test), quote = FALSE)
  cat("\n")
} else {
  open_test <- nul$closure_tab
  cat("-- closure --\n")
  cat("  the population was treated as closed in STEP 2; recruitment and apparent\n")
  cat("  survival are not estimated\n\n")
}

# --- 7. PREDICTED ABUNDANCE PER STATION ------------------------------------

site_new <- SITE_COVS
pred_list <- lapply(top_set, function(m) pred_one(fits[[m]], site_new))
names(pred_list) <- top_set
keep <- !vapply(pred_list, is.null, logical(1))
if (!any(keep)) stop("prediction failed for every model in the top set.", call. = FALSE)
w_use <- w_top[keep] / sum(w_top[keep])

fit_mat <- vapply(pred_list[keep], function(z) z$fit, numeric(nS))
se_mat  <- vapply(pred_list[keep], function(z) z$se,  numeric(nS))
if (is.null(dim(fit_mat))) { fit_mat <- matrix(fit_mat, nS, 1); se_mat <- matrix(se_mat, nS, 1) }

lam_avg <- as.numeric(fit_mat %*% w_use)
se_avg  <- as.numeric(sqrt(se_mat^2 + (fit_mat - lam_avg)^2) %*% w_use)

pb <- pred_one(best_fit, site_new)

site_pred <- data.frame(
  Station = imp$stations,
  lambda_best = pb$fit, lambda_best_SE = pb$se,
  lambda_avg = lam_avg, lambda_avg_SE = se_avg,
  lambda_avg_lower = pmax(lam_avg - 1.959964 * se_avg, 0),
  lambda_avg_upper = lam_avg + 1.959964 * se_avg,
  individuals = as.integer(rowSums(nul$data$y, na.rm = TRUE)),
  trap_nights = as.integer(rowSums(nul$data$eff_days)),
  stringsAsFactors = FALSE)

cat("-- predicted ", if (OPEN) "initial " else "", "abundance per station --\n", sep = "")
cat(sprintf("  model-averaged: min %.3f   median %.3f   mean %.3f   max %.3f\n",
            min(lam_avg), stats::median(lam_avg), mean(lam_avg), max(lam_avg)))
cat(sprintf("  summed over the %d stations: %.2f\n", nS, sum(lam_avg)))
cat("\n")

# --- 8. FIGURE: AICc WEIGHTS -----------------------------------------------

fig_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        axis.text.y = element_text(family = "mono", size = 8),
        legend.position = "none",
        plot.margin = margin(4, 10, 4, 4))

wt <- rank_tab
wt$label <- factor(sprintf("%-4s %s", wt$model, gsub("^~", "\u03bb ~ ", wt$lambda)),
                   levels = rev(sprintf("%-4s %s", wt$model,
                                        gsub("^~", "\u03bb ~ ", wt$lambda))))

p_w <- ggplot(wt, aes(x = weight, y = label, fill = in_top_set)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.25) +
  geom_text(aes(label = formatC(weight, format = "f", digits = 3)),
            hjust = -0.18, size = 2.9, colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "#1B7837", `FALSE` = "#D9D9D9")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = expression(paste("AIC"[c], " weight")), y = NULL) +
  fig_theme

stem <- file.path(DIR_FIG, paste0("03_", PFX, "_aicc_weights"))
ggsave(paste0(stem, ".png"), p_w,
       width = 7.2, height = 0.3 * nrow(wt) + 1, dpi = FIG_DPI, limitsize = FALSE)
ggsave(paste0(stem, ".pdf"), p_w,
       width = 7.2, height = 0.3 * nrow(wt) + 1, limitsize = FALSE)

# --- 9. EXPORT -------------------------------------------------------------

tab <- function(name) file.path(DIR_TAB, paste0("03_", PFX, "_", name, ".csv"))

write.csv(bind_common(aic_tab, eff_tab), tab("candidates_all"),
          row.names = FALSE, fileEncoding = "UTF-8")
if (OPEN) write.csv(dyn_tab, tab("dynamics"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(open_test, tab("open_support"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(rank_tab, tab("aicc_ranking"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(det_choice, tab("detection_structure"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tail_log, tab("K_truncation"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(coef_tab, tab("coefficients"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(site_pred, tab("site_predictions"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(
  hyp = hyp, det_set = det_set, aic_tab = aic_tab, eff_tab = eff_tab,
  det_tab = det_tab, det_choice = det_choice, P_FORM = P_FORM,
  rank_tab = rank_tab, coef_tab = coef_tab, site_pred = site_pred,
  dyn_tab = dyn_tab, open_test = open_test, n_transitions = nT,
  fits = fits, top_set = top_set, weights_top = w_top,
  best_id = best_id, tail_log = tail_log,
  FRAMEWORK = FRAMEWORK, MODEL_LABEL = MODEL_LABEL, PFX = PFX, LAM_TYPE = LAM_TYPE,
  DELTA_SET = DELTA_SET, FAMILY = FAMILY, MIXTURE = MIXTURE,
  K = K_USED, OCC_LEN = OCC_LEN, n_occ = n_occ, nS = nS,
  SITE_COVS = SITE_COVS,
  seed = imp$seed, run_date = Sys.time(), session = sessionInfo()
), file.path(DIR_WS, "03_dm_selection.rds"))

cat("-- output --\n")
for (n in c("candidates_all", if (OPEN) "dynamics", "open_support", "aicc_ranking",
            "detection_structure", "K_truncation", "coefficients", "site_predictions"))
  cat("  tables/03_", PFX, "_", n, ".csv\n", sep = "")
cat("  figures/03_", PFX, "_aicc_weights.(png|pdf)\n", sep = "")
cat("  workspace/03_dm_selection.rds   <- input for STEP 5\n\n")
cat("STEP 3 complete. ", MODEL_LABEL, " | K = ", K_USED, " | best model: ", best_id,
    " | top set: ", paste(top_set, collapse = ", "), "\n\n", sep = "")
