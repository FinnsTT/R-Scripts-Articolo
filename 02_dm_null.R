# 02_dm_null.R  --  STEP 2: null-model structure and closure of the population
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# The open population model of Dail and Madsen (2011) is fitted to the full
# survey period at the occasion length retained in STEP 1b, together with its
# closed counterpart, the N-mixture model of Royle (2004). This step fixes the
# structural choices on the null models gamma(~1) omega(~1) lambda(~1) p(~1)
# and lambda(~1) p(~1):
#
#   1. read the count array at the retained occasion length
#   2. factorial sensitivity of both null models, abundance distribution x K
#   3. choice of the abundance distribution, from the drift of initial
#      abundance of the open null across the integration bound
#   4. choice of the integration bound K on the open null
#   5. closure: the open null against the closed null at the retained
#      distribution and K
#   6. tables and figures
#
# Decision rules, fixed before the run:
#
#   distribution  retained when abundance is stable across K, with a drift
#                 below DRIFT_MAX. An estimate that grows with K is a property
#                 of the truncation rather than of the data, and a drift that
#                 cannot be evaluated because fits did not converge is not
#                 evidence of stability
#   K             the smallest value of the grid at which the truncated upper
#                 tail P(N > K | lambda-hat) is below TAIL_MAX and the estimate
#                 is within STABLE_TOL of its value at the largest K. The tail
#                 is Poisson or negative binomial, following the distribution
#   closure       tested at the retained distribution and K. Under the closed
#                 model omega = 1 and gamma = 0, both on the boundary of the
#                 parameter space, so the likelihood ratio is referred to the
#                 chi-bar-square distribution (Self and Liang 1987), with
#                 weights following the correlation between the estimators of
#                 gamma and omega; the chi-square on two degrees of freedom,
#                 conservative there, is reported alongside. When closure is
#                 rejected the Dail-Madsen model is carried forward, otherwise
#                 the closed N-mixture model at the same distribution and K.
#                 When the open null yields no usable estimates, no evidence
#                 against closure is obtainable and distribution and K are
#                 chosen on the closed null with the same rules
#
# Input : workspace/01_import.rds, workspace/01b_occasion.rds
# Output: tables/02_*.csv, figures/02_*, workspace/02_dm_null.rds
# ===========================================================================

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_OUT <- "C:/Users/Dell/Desktop/codici articolo"

DIR_TAB <- file.path(ROOT_OUT, "tables")
DIR_FIG <- file.path(ROOT_OUT, "figures")
DIR_WS  <- file.path(ROOT_OUT, "workspace")
for (d in c(DIR_TAB, DIR_FIG, DIR_WS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

K_GRID   <- seq(10, 50, by = 5)      # integration bound
FAMILIES <- c(Poisson = "P", NegBin = "NB")

DRIFT_MAX   <- 1.10    # ratio of largest to smallest estimate accepted across K
TAIL_MAX    <- 1e-10   # largest accepted P(N > K | lambda-hat), Dennis et al. 2015
STABLE_TOL  <- 0.01    # relative distance from the estimate at the largest K
ALPHA       <- 0.05    # level of the open-versus-closed test

FIG_DPI <- 600

# The fits of the factorial are independent of one another and run in parallel.
USE_PARALLEL <- TRUE
N_CORES      <- max(1L, parallel::detectCores() - 1L)

suppressPackageStartupMessages({
  library(unmarked)
  library(ggplot2)
})

imp <- readRDS(file.path(DIR_WS, "01_import.rds"))
occ <- readRDS(file.path(DIR_WS, "01b_occasion.rds"))
set.seed(imp$seed)

stations <- imp$stations
dets     <- imp$dets_indep
nS       <- length(stations)
OCC_LEN  <- occ$OCC_LEN
d        <- occ$data_by_length[[as.character(OCC_LEN)]]

SITE_COVS <- imp$site[match(stations, imp$site$Station), imp$COVS, drop = FALSE]
rownames(SITE_COVS) <- NULL

cat("\n=========================================================\n")
cat(" STEP 2 -- null models: closure, distribution and K\n")
cat(" sites = ", nS, " | independent detections = ", nrow(dets),
    " | occasion = ", OCC_LEN, " d\n", sep = "")
cat("=========================================================\n\n")

# --- HELPERS ---------------------------------------------------------------

# Evaluate, recording warnings rather than suppressing them. Only distinct
# messages are stored, up to a cap, while the total is counted: an optimiser
# can raise the same warning thousands of times.
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

# Estimate, standard error and confidence interval on the response scale, for
# an intercept-only submodel. Names are dropped so that the elements of the
# result can be addressed by position-independent names.
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
  d  <- if (link == "log") exp(est) else exp(est) / (1 + exp(est))^2   # delta method
  c(est = unname(g(est)), se = unname(abs(d) * se),
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

# Upper tail of the abundance distribution beyond the integration bound. In
# unmarked the negative binomial has mean lambda and size exp(alpha).
tail_prob <- function(K, lambda, size = NA_real_) {
  if (!is.finite(lambda)) return(NA_real_)
  if (is.finite(size)) stats::pnbinom(K, size = size, mu = lambda, lower.tail = FALSE)
  else stats::ppois(K, lambda, lower.tail = FALSE)
}

aicc <- function(ll, npar, n) -2 * ll + 2 * npar + 2 * npar * (npar + 1) / (n - npar - 1)

# Correlation between the estimators of recruitment and apparent survival, on
# the scale the model is fitted on, used to weight the chi-bar-square.
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

# Tail probability of the likelihood ratio under the chi-bar-square null that
# applies when both parameters lie on the boundary. The constrained set is the
# quadrant gamma >= 0 and 1 - omega >= 0, so the relevant correlation is minus
# the correlation between the estimators of gamma and omega. With that
# correlation zero the weights reduce to the familiar quarter, half, quarter.
p_chibar <- function(chi, rho) {
  if (!is.finite(chi) || chi <= 0) return(1)
  if (!is.finite(rho)) rho <- 0
  w2 <- acos(max(-1, min(1, rho))) / (2 * pi)
  w1 <- 0.5
  w1 * stats::pchisq(chi, 1, lower.tail = FALSE) +
  w2 * stats::pchisq(chi, 2, lower.tail = FALSE)
}

# --- 1. COUNT ARRAY --------------------------------------------------------

cat("-- count array --\n")
cat("  sites x occasions      :", nS, "x", d$n_occ, "\n")
cat("  survey period          :", format(imp$FULL_START), "to", format(imp$FULL_END),
    " (", as.integer(imp$FULL_END - imp$FULL_START) + 1, "days )\n")
cat("  observed cells         :", d$n_cells, "of", nS * d$n_occ, "\n")
cat("  individuals counted    :", d$n_det, "\n")
cat("  cells with >= 1        :", sum(d$y >= 1, na.rm = TRUE), "\n")
cat("  cells with >= 2        :", sum(d$y >= 2, na.rm = TRUE), "\n")
cat("  largest cell count     :", max(d$y, na.rm = TRUE), "\n\n")

OBS_COVS <- list(effort = d$effort, moon = d$moon, human = d$human)

umf_open <- unmarkedFramePCO(y = d$y, siteCovs = SITE_COVS, obsCovs = OBS_COVS,
                             numPrimary = d$n_occ)
umf_closed <- unmarkedFramePCount(y = d$y, siteCovs = SITE_COVS, obsCovs = OBS_COVS)

# --- 2. FACTORIAL SENSITIVITY OF THE NULL MODELS ---------------------------

n_total <- 2L * length(FAMILIES) * length(K_GRID)
cat("-- factorial sensitivity of the null models --\n")
cat("  2 models x ", length(FAMILIES), " distributions x ", length(K_GRID),
    " values of K = ", n_total, " fits\n", sep = "")
cat("  the open likelihood is O(K^2 * T) with T = ", d$n_occ,
    " primary periods; elapsed time is reported per fit\n\n", sep = "")

blank_row <- function(model, fam, K, note) data.frame(
  model = model, occasion_days = OCC_LEN, n_occ = d$n_occ, n_cells = d$n_cells,
  individuals = d$n_det, family = fam, K = K, converged = FALSE,
  lambda = NA_real_, lambda_SE = NA_real_,
  gamma = NA_real_, gamma_SE = NA_real_,
  omega = NA_real_, omega_SE = NA_real_,
  p = NA_real_, p_SE = NA_real_, nb_size = NA_real_,
  lambda_point = NA_real_, gamma_point = NA_real_, omega_point = NA_real_,
  logLik = NA_real_, AIC = NA_real_, CV_lambda = NA_real_,
  tail_prob = NA_real_, seconds = NA_real_,
  n_warnings = 0L, note = note, stringsAsFactors = FALSE)

fit_null <- function(model, fam, K) {
  if (model == "open")
    pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1,
               pformula = ~1, data = umf_open, K = K, mixture = FAMILIES[[fam]],
               dynamics = "constant")
  else
    pcount(~1 ~1, data = umf_closed, K = K, mixture = FAMILIES[[fam]])
}

one_fit <- function(spec) {
  model <- spec$model; fam <- spec$family; K <- spec$K
  tryCatch({
    t1 <- Sys.time()
    r  <- run_quiet(fit_null(model, fam, K))
    secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

    if (!ok(r$value)) {
      list(row = blank_row(model, fam, K, msg1(r$value)), fit = NULL)
    } else {
      lam_type <- if (model == "open") "lambda" else "state"
      lam <- resp(r$value, lam_type, "log")
      pp  <- resp(r$value, "det", "logit")
      gam <- if (model == "open") resp(r$value, "gamma", "log") else
        c(est = NA_real_, se = NA_real_)
      omg <- if (model == "open") resp(r$value, "omega", "logit") else
        c(est = NA_real_, se = NA_real_)
      size <- if (fam == "NegBin") point(r$value, "alpha", "log") else NA_real_
      g_pt <- if (model == "open") point(r$value, "gamma", "log")  else NA_real_
      o_pt <- if (model == "open") point(r$value, "omega", "logit") else NA_real_
      l_pt <- point(r$value, lam_type, "log")
      ll  <- tryCatch(as.numeric(logLik(r$value)), error = function(e) NA_real_)
      ai  <- tryCatch(as.numeric(r$value@AIC),     error = function(e) NA_real_)
      list(row = data.frame(
        model = model, occasion_days = OCC_LEN, n_occ = d$n_occ, n_cells = d$n_cells,
        individuals = d$n_det, family = fam, K = K,
        converged = is.finite(lam[["est"]]) && is.finite(lam[["se"]]),
        lambda = lam[["est"]], lambda_SE = lam[["se"]],
        gamma = gam[["est"]], gamma_SE = gam[["se"]],
        omega = omg[["est"]], omega_SE = omg[["se"]],
        p = pp[["est"]], p_SE = pp[["se"]], nb_size = size,
        lambda_point = l_pt, gamma_point = g_pt, omega_point = o_pt,
        logLik = ll, AIC = ai,
        CV_lambda = lam[["se"]] / lam[["est"]],
        tail_prob = tail_prob(K, lam[["est"]], size), seconds = secs,
        n_warnings = r$n_warnings,
        note = paste(r$warnings, collapse = " | "),
        stringsAsFactors = FALSE), fit = r$value)
    }
  }, error = function(e)
    list(row = blank_row(model, fam, K, paste("error:", conditionMessage(e))), fit = NULL))
}

# Tasks are distributed with load balancing, the open models at the largest
# bounds first because they are the slowest to fit.
specs <- expand.grid(K = K_GRID, family = names(FAMILIES), model = c("open", "closed"),
                     stringsAsFactors = FALSE)
specs <- specs[order(specs$model != "open", -specs$K), ]
specs <- lapply(seq_len(nrow(specs)), function(i) as.list(specs[i, ]))

t0 <- Sys.time()
if (USE_PARALLEL && N_CORES >= 2L) {
  cl <- parallel::makeCluster(min(N_CORES, length(specs)))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(unmarked)))
  parallel::clusterExport(cl, c("one_fit", "fit_null", "blank_row", "run_quiet", "ok",
                                "msg1", "resp", "point", "tail_prob", "umf_open",
                                "umf_closed", "FAMILIES", "OCC_LEN", "d"),
                          envir = .GlobalEnv)
  res <- parallel::parLapplyLB(cl, specs, one_fit)
  parallel::stopCluster(cl)
} else {
  res <- lapply(specs, one_fit)
}

rows <- lapply(res, `[[`, "row")
fits <- lapply(res, `[[`, "fit")
names(fits) <- vapply(specs, function(z) paste(z$model, z$family, z$K, sep = "_"), character(1))
fits <- fits[!vapply(fits, is.null, logical(1))]

sens_all <- do.call(rbind, rows)
if (!is.null(sens_all))
  sens_all <- sens_all[order(sens_all$model != "open", sens_all$family, sens_all$K), ]
if (is.null(sens_all) || !nrow(sens_all))
  stop("the factorial produced no rows.", call. = FALSE)

print(sens_all[, c("model", "family", "K", "converged", "lambda_point", "gamma_point",
                  "omega_point", "logLik", "seconds")],
      row.names = FALSE, digits = 4)

sens_open   <- sens_all[sens_all$model == "open", ]
sens_closed <- sens_all[sens_all$model == "closed", ]

# written before the decision rules are applied, so that the fits remain
# available for inspection when no setting satisfies them
write.csv(sens_open, file.path(DIR_TAB, "02_dm_null_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sens_closed, file.path(DIR_TAB, "02_nm_null_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n  elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    " min; converged open ", sum(sens_open$converged), "/", nrow(sens_open),
    ", closed ", sum(sens_closed$converged), "/", nrow(sens_closed), "\n\n", sep = "")

# --- 3-4. DISTRIBUTION AND INTEGRATION BOUND ------------------------------
# Drift is the ratio of the largest to the smallest estimate of abundance over
# the values of K at which the model converged. Among the identifiable
# distributions the one with the lower AIC at the largest value of K is
# retained, and K is the smallest value that satisfies the tail and the
# stability rule.

select_structure <- function(sens, model) {
  drift_tab <- do.call(rbind, lapply(names(FAMILIES), function(fam) {
    s <- sens[sens$family == fam & sens$converged, ]
    l <- s$lambda[is.finite(s$lambda)]
    data.frame(model = model, family = fam,
               n_converged = nrow(s), n_fits = length(K_GRID),
               lambda_min = if (length(l)) min(l) else NA_real_,
               lambda_max = if (length(l)) max(l) else NA_real_,
               drift      = if (length(l) >= 2) max(l) / min(l) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  drift_tab$identifiable <- is.finite(drift_tab$drift) &
    drift_tab$drift < DRIFT_MAX & drift_tab$n_converged == drift_tab$n_fits
  out <- list(model = model, sens = sens, drift_tab = drift_tab,
              FAMILY = NA_character_, k_tab = NULL, K = NA_real_, reason = NA_character_)
  if (!any(drift_tab$identifiable)) {
    out$reason <- "no distribution satisfies the identifiability rule"
    return(out)
  }
  cand <- drift_tab$family[drift_tab$identifiable]
  aic_at_max <- vapply(cand, function(f) {
    s <- sens[sens$family == f & sens$K == max(K_GRID), ]
    if (nrow(s) && is.finite(s$AIC[1])) s$AIC[1] else Inf
  }, numeric(1))
  out$FAMILY <- cand[which.min(aic_at_max)]

  sk <- sens[sens$family == out$FAMILY & sens$converged, ]
  sk <- sk[order(sk$K), ]
  lambda_ref   <- sk$lambda[nrow(sk)]
  sk$rel_diff  <- abs(sk$lambda - lambda_ref) / lambda_ref
  sk$tail_ok   <- is.finite(sk$tail_prob) & sk$tail_prob < TAIL_MAX
  sk$stable_ok <- sk$rel_diff < STABLE_TOL
  sk$eligible  <- sk$tail_ok & sk$stable_ok
  out$k_tab <- sk
  if (!any(sk$eligible)) {
    out$reason <- "no value of K satisfies both the tail and the stability rule"
    return(out)
  }
  out$K <- min(sk$K[sk$eligible])
  out
}

report_structure <- function(st) {
  label <- if (st$model == "open") "Dail-Madsen" else "N-mixture"
  cat("-- distribution: drift of abundance across K, ", label, " --\n", sep = "")
  print(st$drift_tab, row.names = FALSE, digits = 4)
  if (!is.null(st$k_tab)) {
    cat("\n-- K: truncation and stability, ", label, ", ", st$FAMILY, " --\n", sep = "")
    print(st$k_tab[, c("K", "lambda", "lambda_SE", "tail_prob", "rel_diff",
                       "tail_ok", "stable_ok", "eligible")],
          row.names = FALSE, digits = 4)
  }
  if (is.na(st$K)) cat("\n  ", st$reason, "\n\n", sep = "")
  else cat("\n  RETAINED: ", st$FAMILY, ", K = ", st$K, "\n\n", sep = "")
}

st_open <- select_structure(sens_open, "open")
report_structure(st_open)

# --- 5. CLOSURE ------------------------------------------------------------

closure_row <- function(FAMILY, K) data.frame(
  comparison = "open versus closed null", family = FAMILY, K = K,
  logLik_open = NA_real_, logLik_closed = NA_real_, df = 2L,
  chi_square = NA_real_, rho_gamma_omega = NA_real_,
  p_mixture = NA_real_, p_chisq2 = NA_real_,
  AICc_open = NA_real_, AICc_closed = NA_real_, delta_AICc = NA_real_,
  gamma = NA_real_, omega = NA_real_,
  verdict = NA_character_, note = NA_character_, stringsAsFactors = FALSE)

if (!is.na(st_open$K)) {
  st <- st_open
  FAMILY <- st$FAMILY; K_RETAINED <- st$K
  closure_tab <- closure_row(FAMILY, K_RETAINED)
  fo <- fits[[paste("open",   FAMILY, K_RETAINED, sep = "_")]]
  fc <- fits[[paste("closed", FAMILY, K_RETAINED, sep = "_")]]
  ro <- sens_open[sens_open$family == FAMILY & sens_open$K == K_RETAINED, ]
  if (is.null(fc) || !sens_closed$converged[sens_closed$family == FAMILY &
                                             sens_closed$K == K_RETAINED]) {
    closure_tab$verdict <- "not assessable"
    closure_tab$note <- "the closed null did not produce usable estimates"
    FRAMEWORK <- "open"
  } else {
    llo <- as.numeric(logLik(fo)); llc <- as.numeric(logLik(fc))
    chi <- 2 * (llo - llc)
    rr  <- rho_gamma_omega(fo)
    closure_tab$logLik_open <- llo; closure_tab$logLik_closed <- llc
    closure_tab$chi_square <- chi
    closure_tab$rho_gamma_omega <- rr
    closure_tab$p_mixture <- p_chibar(chi, -rr)
    closure_tab$p_chisq2 <- stats::pchisq(max(chi, 0), df = 2, lower.tail = FALSE)
    closure_tab$AICc_open <- aicc(llo, length(coef(fo)), nS)
    closure_tab$AICc_closed <- aicc(llc, length(coef(fc)), nS)
    closure_tab$delta_AICc <- closure_tab$AICc_closed - closure_tab$AICc_open
    closure_tab$gamma <- ro$gamma_point; closure_tab$omega <- ro$omega_point
    rejected <- is.finite(closure_tab$p_mixture) && closure_tab$p_mixture < ALPHA
    closure_tab$verdict <- if (rejected) "closure rejected" else "closure not rejected"
    closure_tab$note <- "boundary hypothesis; p_mixture is the reference, p_chisq2 is conservative"
    FRAMEWORK <- if (rejected) "open" else "closed"
  }
  sens <- if (FRAMEWORK == "open") sens_open else sens_closed
  drift_all <- st_open$drift_tab
} else {
  closure_tab <- closure_row(NA_character_, NA_real_)
  closure_tab$verdict <- "not assessable"
  closure_tab$note <- paste("the open null did not produce usable estimates:", st_open$reason)
  FRAMEWORK <- "closed"
  st <- select_structure(sens_closed, "closed")
  report_structure(st)
  drift_all <- rbind(st_open$drift_tab, st$drift_tab)
  if (is.na(st$K))
    stop(st$reason, " on the closed null; inspect tables/02_nm_null_sensitivity.csv",
         call. = FALSE)
  FAMILY <- st$FAMILY; K_RETAINED <- st$K
  closure_tab$family <- FAMILY; closure_tab$K <- K_RETAINED
  sens <- sens_closed
}

MODEL_LABEL <- if (FRAMEWORK == "open") "Dail-Madsen" else "N-mixture"
PFX <- if (FRAMEWORK == "open") "dm" else "nm"
sk  <- st$k_tab

cat("-- closure: open versus closed null, ", FAMILY, ", K = ", K_RETAINED, " --\n", sep = "")
print(t(closure_tab), quote = FALSE)
cat("\n  RETAINED FRAMEWORK : ", FRAMEWORK, " (", MODEL_LABEL, ")\n\n", sep = "")

write.csv(closure_tab, file.path(DIR_TAB, "02_closure_test.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(drift_all, file.path(DIR_TAB, "02_distribution_drift.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
if (!is.null(sk))
  write.csv(sk[, c("model", "family", "K", "lambda", "lambda_SE", "nb_size", "tail_prob",
                   "rel_diff", "tail_ok", "stable_ok", "eligible")],
            file.path(DIR_TAB, "02_K_selection.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

best <- sens[sens$family == FAMILY & sens$K == K_RETAINED, ]

cat("-- null model at the retained settings --\n")
cat(sprintf("  abundance lambda         : %.4f (SE %.4f)\n", best$lambda, best$lambda_SE))
if (FRAMEWORK == "open") {
  cat(sprintf("  recruitment       gamma  : %.4f (SE %.4f)\n", best$gamma,  best$gamma_SE))
  cat(sprintf("  apparent survival omega  : %.4f (SE %.4f)\n", best$omega,  best$omega_SE))
}
cat(sprintf("  detection         p      : %.4f (SE %.4f)\n", best$p,      best$p_SE))
if (FRAMEWORK == "open")
  cat(sprintf("  retention over the survey, omega^(T-1) with T = %d : %.3f\n",
              d$n_occ, best$omega^(d$n_occ - 1)))
cat("\n")

# --- 6. FIGURES ------------------------------------------------------------

FAM_LINE <- c(Poisson = "#0571B0", NegBin = "#CA0020")
FAM_BAND <- c(Poisson = "#CDE3F1", NegBin = "#F7D4D7")

fig_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey93"),
        panel.border = element_rect(colour = "grey35", linewidth = 0.45),
        axis.text = element_text(size = 8.5, colour = "grey20"),
        legend.key.height = unit(0.9, "lines"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        plot.margin = margin(4, 8, 4, 4))

sp <- st$sens[st$sens$converged, ]

p_k <- ggplot(sp, aes(x = K, y = lambda)) +
  geom_ribbon(aes(ymin = pmax(lambda - 1.96 * lambda_SE, 0),
                  ymax = lambda + 1.96 * lambda_SE, fill = family),
              colour = NA, alpha = 0.75) +
  geom_line(aes(colour = family), linewidth = 1) +
  geom_point(aes(colour = family, shape = family), size = 2.2,
             fill = "white", stroke = 0.9) +
  scale_colour_manual(values = FAM_LINE) +
  scale_fill_manual(values = FAM_BAND) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(breaks = K_GRID) +
  labs(x = "Integration bound K",
       y = expression(hat(lambda)~"(individuals per station)")) +
  fig_theme

p_p <- ggplot(sp, aes(x = K, y = p)) +
  geom_ribbon(aes(ymin = pmax(p - 1.96 * p_SE, 0),
                  ymax = pmin(p + 1.96 * p_SE, 1), fill = family),
              colour = NA, alpha = 0.75) +
  geom_line(aes(colour = family), linewidth = 1) +
  geom_point(aes(colour = family, shape = family), size = 2.2,
             fill = "white", stroke = 0.9) +
  scale_colour_manual(values = FAM_LINE) +
  scale_fill_manual(values = FAM_BAND) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(breaks = K_GRID) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Integration bound K", y = expression(hat(p))) +
  fig_theme

stem <- file.path(DIR_FIG, paste0("02_", if (st$model == "open") "dm" else "nm", "_K_"))
ggsave(paste0(stem, "lambda.png"), p_k, width = 5.8, height = 3.6, dpi = FIG_DPI)
ggsave(paste0(stem, "lambda.pdf"), p_k, width = 5.8, height = 3.6)
ggsave(paste0(stem, "detection.png"), p_p, width = 5.8, height = 3.6, dpi = FIG_DPI)
ggsave(paste0(stem, "detection.pdf"), p_p, width = 5.8, height = 3.6)

# --- 7. EXPORT -------------------------------------------------------------

saveRDS(list(
  data = d, umf = if (FRAMEWORK == "open") umf_open else umf_closed,
  umf_open = umf_open, umf_closed = umf_closed, SITE_COVS = SITE_COVS,
  sens = sens, sens_open = sens_open, sens_closed = sens_closed,
  closure_tab = closure_tab, FRAMEWORK = FRAMEWORK, MODEL_LABEL = MODEL_LABEL, PFX = PFX,
  STRUCTURE_FROM = st$model,
  drift_tab = drift_all, k_tab = sk, null_best = best,
  fits_null = fits,
  OCC_LEN = OCC_LEN, K_GRID = K_GRID, FAMILIES = FAMILIES,
  FAMILY = FAMILY, MIXTURE = unname(FAMILIES[[FAMILY]]), K_RETAINED = K_RETAINED,
  DRIFT_MAX = DRIFT_MAX, TAIL_MAX = TAIL_MAX, STABLE_TOL = STABLE_TOL,
  ALPHA = ALPHA,
  seed = imp$seed, run_date = Sys.time(), session = sessionInfo()
), file.path(DIR_WS, "02_dm_null.rds"))

cat("-- output --\n")
cat("  tables/02_dm_null_sensitivity.csv\n")
cat("  tables/02_nm_null_sensitivity.csv\n")
cat("  tables/02_closure_test.csv\n")
cat("  tables/02_distribution_drift.csv\n")
cat("  tables/02_K_selection.csv\n")
cat("  figures/02_", if (st$model == "open") "dm" else "nm", "_K_(lambda|detection).(png|pdf)\n", sep = "")
cat("  workspace/02_dm_null.rds   <- input for STEP 3\n\n")
cat("STEP 2 complete. Framework: ", MODEL_LABEL, ", ", FAMILY, ", K = ", K_RETAINED,
    ", occasion = ", OCC_LEN, " d\n\n", sep = "")
