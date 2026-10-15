# 06_paper_tables.R  --  STEP 6: formatted tables and captions
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# Turns the output of steps 01 to 05 into tables laid out for the
# manuscript: scientific column headings, symbols and units in the headers,
# numbers rounded to a reported precision, and hypotheses written out.
#
#   Table 1  candidate set of ecological hypotheses
#   Table 2  models within two AICc units, both analyses
#   Table 3  coefficients of the best model and of the model average
#   Table 4  recruitment and apparent survival, when the open population
#            model is retained
#   Table S1 every model fitted, both analyses
#   Table S2 sensitivity of the open and closed null models to distribution
#            and to K
#   Table S3 sensitivity of the Royle-Nichols null model to K
#   Table S4 predicted abundance per station
#   Table S5 test of population closure, open against closed model
#
# Tables are written as Word documents when flextable and officer are
# available, and as delimited text with the same headings in every case.
# Captions for the tables are written alongside them.
#
# Input : workspace/02_dm_null.rds, 03_dm_selection.rds, 04_rn.rds
# Output: paper_tables/Table_*.csv, .docx and captions.md
# ===========================================================================

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_OUT <- "C:/Users/Dell/Desktop/codici articolo"

DIR_WS  <- file.path(ROOT_OUT, "workspace")
DIR_PAP <- file.path(ROOT_OUT, "paper_tables")
dir.create(DIR_PAP, recursive = TRUE, showWarnings = FALSE)

# Files left by an earlier run under names no longer produced are removed at
# the end, once the new ones exist, so that an interrupted run never leaves the
# folder empty. The names written are collected as the tables are written.
written <- character(0)

DEC <- "."          # decimal mark in the delimited output
SEP <- ","          # field separator

nul <- readRDS(file.path(DIR_WS, "02_dm_null.rds"))
dm  <- readRDS(file.path(DIR_WS, "03_dm_selection.rds"))
rn  <- readRDS(file.path(DIR_WS, "04_rn.rds"))

DM_LAB <- dm$MODEL_LABEL
RN_LAB <- "Royle-Nichols"

cat("\n=========================================================\n")
cat(" STEP 6 -- formatted tables and captions\n")
cat("=========================================================\n\n")

# --- HELPERS ---------------------------------------------------------------

has_word <- requireNamespace("flextable", quietly = TRUE) &&
            requireNamespace("officer",   quietly = TRUE)

num <- function(x, d = 2) {
  out <- formatC(as.numeric(x), format = "f", digits = d)
  out[!is.finite(as.numeric(x))] <- "--"
  out
}

sci <- function(x, d = 1) {
  out <- formatC(as.numeric(x), format = "e", digits = d)
  out[!is.finite(as.numeric(x))] <- "--"
  out
}

# Estimate with its standard error in one cell.
est_se <- function(e, s, d = 3) {
  out <- paste0(num(e, d), " (", num(s, d), ")")
  out[!is.finite(as.numeric(e))] <- "--"
  out
}

# Confidence interval in one cell.
ci <- function(l, u, d = 3) {
  out <- paste0(num(l, d), ", ", num(u, d))
  out[!is.finite(as.numeric(l))] <- "--"
  out
}

# Covariate names as they should read in a table.
pretty_formula <- function(f) {
  f <- gsub("^~", "", f)
  f <- gsub("dist_water_points",    "distance to water points",    f)
  f <- gsub("dist_fish_river",      "distance to Fish River",      f)
  f <- gsub("dist_konkiep_river",   "distance to Konkiep River",   f)
  f <- gsub("dist_tourist_road",    "distance to tourist road",    f)
  f <- gsub("dist_staff_buildings", "distance to staff buildings", f)
  f <- gsub("\\btri\\b", "TRI", f)
  f <- gsub("\\btpi\\b", "TPI", f)
  f <- gsub("\\bmoon\\b",  "moon illumination", f)
  f <- gsub("\\bhuman\\b", "human activity",    f)
  f <- gsub(" \\* ", " \u00d7 ", f)
  f <- gsub(":", " \u00d7 ", f)
  f <- gsub("^1$", "constant", f)
  f
}

# Parameter names as they should read in a table.
pretty_term <- function(term, method) {
  sub_ <- sub("\\(.*$", "", term)
  par  <- sub("^[^(]*\\(", "", sub("\\)$", "", term))
  submodel <- c(lam = "\u03bb", lambda = "\u03bb",
                gam = "\u03b3", gamConst = "\u03b3",
                om = "\u03c9", omega = "\u03c9",
                alpha = "\u03b1",
                p = if (method == RN_LAB) "r" else "p")[sub_]
  submodel[is.na(submodel)] <- sub_[is.na(submodel)]
  par[sub_ == "alpha"] <- "Int"      # log of the negative binomial size parameter
  par <- gsub("Int", "intercept", par)
  par <- pretty_formula(par)
  data.frame(submodel = unname(submodel), parameter = par, stringsAsFactors = FALSE)
}

# Page geometry. Wide tables are laid out on a landscape page so that no column
# is wrapped; the narrow ones stay portrait.
page_section <- function(landscape) {
  w <- if (landscape) 11.69 else 8.27
  h <- if (landscape) 8.27  else 11.69
  officer::prop_section(
    page_size = officer::page_size(
      width = w, height = h,
      orient = if (landscape) "landscape" else "portrait"),
    page_margins = officer::page_mar(top = 0.8, bottom = 0.8,
                                     left = 0.8, right = 0.8,
                                     header = 0.4, footer = 0.4, gutter = 0))
}

# Orientation follows the width the table actually needs, not the number of
# columns: a table of four columns one of which holds a sentence is wider than
# one of six columns holding numbers. The width of a column is the longer of
# its heading and its widest entry, and a table wider than the characters that
# fit across a portrait page is set on a landscape one.
LANDSCAPE_MIN_CHARS <- 90

table_width <- function(df) {
  w <- vapply(seq_along(df), function(j)
    max(nchar(names(df)[j]), max(nchar(as.character(df[[j]])), 0L)), integer(1))
  sum(pmin(w, 40L)) + 2L * ncol(df)
}

write_table <- function(df, stem, caption,
                        landscape = table_width(df) > LANDSCAPE_MIN_CHARS) {
  # UTF-8 with a byte-order mark, so that spreadsheet software reads the
  # symbols in the headings as symbols
  con <- file(file.path(DIR_PAP, paste0(stem, ".csv")), open = "w",
              encoding = "UTF-8")
  cat("\ufeff", file = con)
  utils::write.table(df, con, sep = SEP, dec = DEC, row.names = FALSE,
                     qmethod = "double")
  close(con)
  written <<- c(written, paste0(stem, ".csv"), paste0(stem, ".docx"))
  if (has_word) {
    ft <- flextable::flextable(df)
    ft <- flextable::set_caption(ft, caption)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::align(ft, align = "center", part = "header")
    ft <- flextable::align(ft, align = "left",   part = "body", j = 1)
    ft <- flextable::padding(ft, padding.top = 2, padding.bottom = 2, part = "all")
    ft <- flextable::autofit(ft)
    ft <- flextable::fit_to_width(ft, max_width = if (landscape) 10 else 6.6)
    doc <- officer::read_docx()
    doc <- flextable::body_add_flextable(doc, ft)
    doc <- officer::body_set_default_section(doc, page_section(landscape))
    print(doc, target = file.path(DIR_PAP, paste0(stem, ".docx")))
  }
  cat(sprintf("  %-38s %3d cols, %3d chars  %s\n", stem, ncol(df),
              table_width(df),
              if (!has_word) "csv" else
                paste0("csv + docx, ", if (landscape) "landscape" else "portrait")))
}

if (!has_word)
  cat("NOTE: flextable and officer are not installed; only delimited text is written.\n",
      "      install.packages(c(\"flextable\", \"officer\")) to obtain Word tables.\n\n",
      sep = "")

cat("-- writing --\n")

# --- TABLE 1: CANDIDATE SET ------------------------------------------------

h <- dm$hyp
tab1 <- data.frame(
  id        = h$id,
  hypo      = h$hypothesis,
  submodel  = paste0("\u03bb ~ ", pretty_formula(h$lambda)),
  sign      = h$prediction,
  fitted    = ifelse(h$vif_ok, "yes",
                     paste0("no, VIF screen removed ",
                            pretty_formula(h$excluded_covariates))),
  check.names = FALSE, stringsAsFactors = FALSE)

names(tab1) <- c("Model", "Ecological hypothesis",
                 "Abundance submodel", "Predicted sign", "Fitted")

write_table(tab1, "Table_1_candidate_set",
  "Table 1. Candidate set of ecological hypotheses for brown hyena abundance.")

# --- TABLE 2: MODELS WITHIN TWO AICc UNITS ---------------------------------

top_block <- function(obj, method, det_col) {
  r <- obj$rank_tab[obj$rank_tab$in_top_set, , drop = FALSE]
  det_lab <- if (method == RN_LAB) "r ~ " else "p ~ "
  out <- data.frame(
    Analysis   = method,
    Model      = r$model,
    Abundance  = paste0("\u03bb ~ ", pretty_formula(r$lambda)),
    Detection  = paste0(det_lab, pretty_formula(r[[det_col]])),
    np         = r$nPars,
    logL       = num(r$logLik, 2),
    AICc       = num(r$AICc, 2),
    dAICc      = num(r$delta_AICc, 2),
    wi         = num(r$weight, 3),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c("Analysis", "Model", "Abundance submodel",
                  "Detection submodel", "np", "logL", "AICc",
                  "\u0394AICc", "wi")
  out
}

tab2 <- rbind(top_block(dm, DM_LAB, "p"),
              top_block(rn, RN_LAB, "r_formula"))
tab2$Analysis[duplicated(tab2$Analysis)] <- ""

write_table(tab2, "Table_2_top_models",
  paste0("Table 2. Models within two AICc units of the best-supported model, ",
         "for the ", DM_LAB, " model and the Royle-Nichols model. np, number ",
         "of estimated parameters; wi, Akaike weight."))

# --- TABLE 3: COEFFICIENTS -------------------------------------------------

coef_block <- function(obj, method) {
  d <- obj$coef_tab
  pt <- pretty_term(d$term, method)
  data.frame(
    `Model`       = method,
    `Submodel`    = pt$submodel,
    `Parameter`   = pt$parameter,
    `Best model: estimate (SE)`   = est_se(d$best_estimate, d$best_SE, 3),
    `Best model: 95% CI`          = ci(d$best_lower, d$best_upper, 3),
    `Averaged: estimate (uncond. SE)` = est_se(d$avg_estimate, d$avg_uncond_SE, 3),
    `Averaged: 95% CI`            = ci(d$avg_lower, d$avg_upper, 3),
    `wi of models containing the term` = num(d$weight_present, 3),
    check.names = FALSE, stringsAsFactors = FALSE)
}

tab3 <- rbind(coef_block(dm, DM_LAB), coef_block(rn, RN_LAB))
tab3$Model[duplicated(tab3$Model)] <- ""
tab3$Submodel[c(FALSE, tab3$Submodel[-1] == tab3$Submodel[-nrow(tab3)] &
                       tab3$Model[-1] == "")] <- ""

write_table(tab3, "Table_3_coefficients",
  paste0("Table 3. Coefficients on the link scale for the best-supported model ",
         "and averaged over the models within two AICc units, with unconditional ",
         "standard errors. Covariates were standardised before fitting."))

# --- TABLE 4: RECRUITMENT AND APPARENT SURVIVAL ----------------------------

# Parameters carry their symbol. The exponent of the retention row is written
# with superscript digits, so the row reads as it would in the text.
sup_num <- function(n) {
  glyph <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3",
             "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077",
             "8" = "\u2078", "9" = "\u2079")
  paste(glyph[strsplit(as.character(n), "")[[1]]], collapse = "")
}

if (!is.null(dm$dyn_tab)) {
  d4 <- dm$dyn_tab
  tab4a <- data.frame(
    quantity = c(
      "\u03b3 (recruitment, individuals per station per occasion)",
      "\u03c9 (apparent survival per occasion)",
      paste0("\u03c9", sup_num(dm$n_transitions),
             " (retention over the survey)")),
    best = paste0(num(d4$best_estimate, 4), " (", num(d4$best_lower, 4), ", ",
                  num(d4$best_upper, 4), ")"),
    avg  = paste0(num(d4$avg_estimate, 4), " (", num(d4$avg_lower, 4), ", ",
                  num(d4$avg_upper, 4), ")"),
    stringsAsFactors = FALSE)
  names(tab4a) <- c("Quantity", "Best model: estimate (95% CI)",
                    "Averaged: estimate (95% CI)")

  write_table(tab4a, "Table_4_dynamics",
    paste0("Table 4. Recruitment (\u03b3) and apparent survival (\u03c9) ",
           "estimated by the Dail-Madsen model, for the best-supported model ",
           "and averaged over the models within two AICc units, with the ",
           "retention implied over the ", dm$n_transitions,
           " transitions of the survey."))
}

# --- TABLE S1: FULL CANDIDATE SET ------------------------------------------

full_block <- function(obj, method, det_col) {
  k <- intersect(names(obj$aic_tab), names(obj$eff_tab))
  a <- rbind(obj$aic_tab[, k, drop = FALSE], obj$eff_tab[, k, drop = FALSE])
  rk <- obj$rank_tab[, c("model", "delta_AICc", "weight")]
  a <- merge(a, rk, by = "model", all.x = TRUE, sort = FALSE)
  a <- a[order(is.na(a$AICc), a$AICc), , drop = FALSE]
  det_lab <- if (method == RN_LAB) "r ~ " else "p ~ "
  out <- data.frame(
    method     = method,
    model      = a$model,
    hypothesis = a$hypothesis,
    lam        = paste0("\u03bb ~ ", pretty_formula(a$lambda)),
    det        = paste0(det_lab, pretty_formula(a[[det_col]])),
    np         = ifelse(is.na(a$nPars), "--", as.character(a$nPars)),
    logL       = num(a$logLik, 2),
    AICc       = num(a$AICc, 2),
    dAICc      = num(a$delta_AICc, 2),
    wi         = num(a$weight, 3),
    status     = ifelse(a$converged, "retained",
                        ifelse(is.na(a$note) | a$note == "",
                               "not estimable", a$note)),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c("Analysis", "Model", "Ecological hypothesis",
                  "Abundance submodel", "Detection submodel", "np", "logL",
                  "AICc", "\u0394AICc", "wi", "Status")
  out
}

tabS1 <- rbind(full_block(dm, DM_LAB, "p"),
               full_block(rn, RN_LAB, "r_formula"))
tabS1$Analysis[duplicated(tabS1$Analysis)] <- ""

write_table(tabS1, "Table_S1_full_candidate_set",
  paste0("Table S1. Every model fitted, for both analyses, ordered by AICc. ",
         "The models named null_p_* and null_r_* are the fits on the abundance ",
         "null used to settle the detection submodel before the hypotheses were ",
         "examined; they are not part of the ranked candidate set. Models that ",
         "were not fitted or did not produce usable estimates are listed with ",
         "the reason and were excluded from inference. np, number of estimated ",
         "parameters; wi, Akaike weight."))

# --- TABLE S2: NULL MODEL SENSITIVITY -------------------------------------
# Both null models are reported. Open fits on the closed boundary have a
# singular Hessian and no standard errors; their point estimates of
# recruitment and apparent survival are given.

s <- rbind(nul$sens_open, nul$sens_closed)
tabS2 <- data.frame(
  model        = ifelse(s$model == "open", "Dail-Madsen", "N-mixture"),
  distribution = s$family,
  K            = s$K,
  converged    = ifelse(s$converged, "yes", "no"),
  lambda       = est_se(s$lambda, s$lambda_SE, 3),
  gamma        = est_se(s$gamma,  s$gamma_SE,  4),
  omega        = est_se(s$omega,  s$omega_SE,  4),
  gamma_pt     = sci(s$gamma_point, 1),
  omega_pt     = num(s$omega_point, 5),
  p            = est_se(s$p,      s$p_SE,      4),
  cv           = num(s$CV_lambda, 3),
  tail         = sci(s$tail_prob, 1),
  AIC          = num(s$AIC, 2),
  check.names = FALSE, stringsAsFactors = FALSE)
names(tabS2) <- c("Model", "Abundance distribution", "K", "Converged",
                  "\u03bb (SE)", "\u03b3 (SE)", "\u03c9 (SE)",
                  "\u03b3 point estimate", "\u03c9 point estimate", "p (SE)",
                  "CV of \u03bb", "P(N > K)", "AIC")
tabS2$Model[duplicated(tabS2$Model)] <- ""

write_table(tabS2, "Table_S2_null_sensitivity",
  paste0("Table S2. Sensitivity of the null Dail-Madsen and N-mixture models to ",
         "the abundance distribution and to the integration bound K, at a ",
         nul$OCC_LEN, "-day sampling occasion. \u03bb, abundance (initial ",
         "abundance in the Dail-Madsen model); \u03b3, recruitment; \u03c9, ",
         "apparent survival; p, detection probability."))

# --- TABLE S3: ROYLE-NICHOLS NULL MODEL SENSITIVITY ------------------------

s <- rn$sens
tabS3 <- data.frame(
  K         = s$K,
  converged = ifelse(s$converged, "yes", "no"),
  lambda    = est_se(s$lambda, s$lambda_SE, 3),
  r         = est_se(s$r_det,  s$r_det_SE,  4),
  cv        = num(s$CV_lambda, 3),
  tail      = sci(s$tail_prob, 1),
  AIC       = num(s$AIC, 2),
  check.names = FALSE, stringsAsFactors = FALSE)
names(tabS3) <- c("K", "Converged", "\u03bb (SE)", "r (SE)",
                  "CV of \u03bb", "P(N > K)", "AIC")

write_table(tabS3, "Table_S3_royle_nichols_sensitivity",
  paste0("Table S3. Sensitivity of the null Royle-Nichols model to the ",
         "integration bound K, at a ", rn$OCC_LEN, "-day sampling occasion. \u03bb, ",
         "abundance; r, per-individual detection probability."))

# --- TABLE S4: PREDICTED ABUNDANCE PER STATION -----------------------------

tabS4 <- data.frame(
  station    = dm$site_pred$Station,
  nights     = dm$site_pred$trap_nights,
  individuals = dm$site_pred$individuals,
  occasions  = rn$site_pred$occasions_with_detection,
  dm_est     = est_se(dm$site_pred$lambda_avg, dm$site_pred$lambda_avg_SE, 3),
  dm_ci      = ci(dm$site_pred$lambda_avg_lower, dm$site_pred$lambda_avg_upper, 3),
  rn_est     = est_se(rn$site_pred$lambda_avg, rn$site_pred$lambda_avg_SE, 3),
  rn_ci      = ci(rn$site_pred$lambda_avg_lower, rn$site_pred$lambda_avg_upper, 3),
  check.names = FALSE, stringsAsFactors = FALSE)
names(tabS4) <- c("Station", "Trap-nights", "Individuals counted",
                  "Occasions with a detection",
                  paste0(DM_LAB, ": \u03bb (SE)"), paste0(DM_LAB, ": 95% CI"),
                  paste0(RN_LAB, ": \u03bb (SE)"), paste0(RN_LAB, ": 95% CI"))

write_table(tabS4, "Table_S4_station_estimates",
  paste0("Table S4. Model-averaged abundance predicted at each camera-trap ",
         "station by the two analyses, with unconditional standard errors."))

# --- TABLE S5: TEST OF POPULATION CLOSURE ----------------------------------
# One row for the null models, which decides the framework, and, when the open
# model is retained, one row for the best-supported model against its closed
# counterpart.

ct <- nul$closure_tab
closure_rows <- data.frame(
  comparison = "Null models",
  distribution = ct$family, K = ct$K,
  logLik_open = ct$logLik_open, logLik_closed = ct$logLik_closed,
  chi = ct$chi_square, rho = ct$rho_gamma_omega,
  p_mix = ct$p_mixture, p_chi2 = ct$p_chisq2, dAICc = ct$delta_AICc,
  gamma = ct$gamma, omega = ct$omega, verdict = ct$verdict,
  stringsAsFactors = FALSE)
if (dm$FRAMEWORK == "open") {
  ot <- dm$open_test
  closure_rows <- rbind(closure_rows, data.frame(
    comparison = paste0("Best-supported model (", dm$best_id, ")"),
    distribution = dm$FAMILY, K = dm$K,
    logLik_open = ot$logLik_open, logLik_closed = ot$logLik_closed,
    chi = ot$chi_square, rho = ot$rho_gamma_omega,
    p_mix = ot$p_mixture, p_chi2 = ot$p_chisq2, dAICc = ot$delta_AICc,
    gamma = ot$gamma, omega = ot$omega, verdict = ot$verdict,
    stringsAsFactors = FALSE))
}

tabS5 <- data.frame(
  comparison = closure_rows$comparison,
  distribution = ifelse(closure_rows$distribution == "NegBin", "negative binomial",
                        closure_rows$distribution),
  K     = closure_rows$K,
  llo   = num(closure_rows$logLik_open, 2),
  llc   = num(closure_rows$logLik_closed, 2),
  lr    = num(closure_rows$chi, 3),
  rho   = num(closure_rows$rho, 3),
  pmix  = num(closure_rows$p_mix, 4),
  pchi  = num(closure_rows$p_chi2, 4),
  daicc = num(closure_rows$dAICc, 2),
  gamma = sci(closure_rows$gamma, 2),
  omega = num(closure_rows$omega, 5),
  verdict = closure_rows$verdict,
  check.names = FALSE, stringsAsFactors = FALSE)
names(tabS5) <- c("Comparison", "Abundance distribution", "K",
                  "logL open", "logL closed", "LR",
                  "\u03c1(\u03b3, \u03c9)",
                  "P, \u03c7\u0304\u00b2", "P, \u03c7\u00b2\u2082",
                  "\u0394AICc", "\u03b3", "\u03c9", "Verdict")

write_table(tabS5, "Table_S5_closure_test",
  paste0("Table S5. Likelihood-ratio test of population closure, Dail-Madsen ",
         "against N-mixture model with the same abundance and detection structure, ",
         "distribution and integration bound K."))

# --- CAPTIONS --------------------------------------------------------------

cap <- c(
"# Table captions",
"",
"The choice of the occasion length and the choice of the detection submodel are",
"not tabulated: they are reported in the text. The numbers are in",
paste0("tables/01b_occasion_stability.csv, tables/03_", dm$PFX,
       "_detection_structure.csv and tables/04_rn_detection_structure.csv."),
"",
"Generated from the fitted models.",
"",
paste0("**Table 1.** Candidate set of ecological hypotheses for brown ",
       "hyena (*Parahyaena brunnea*) abundance in the Fish River Canyon, Namibia, ",
       "following the hypotheses and model names of the study. Every model was fitted ",
       "with the detection submodel retained beforehand on the abundance null: ",
       "p ", dm$P_FORM, " for the ", DM_LAB, " model and r ", rn$P_FORM,
       " for the Royle-Nichols model."),
"",
paste0("**Table 2.** Models within two AICc units of the best-supported model, for ",
       "the ", DM_LAB, " model and the Royle-Nichols model, fitted ",
       "to ", dm$nS, " camera-trap stations over ", dm$n_occ, " sampling occasions ",
       "of ", dm$OCC_LEN, " days. np, number of estimated parameters; logL, ",
       "log-likelihood; wi, Akaike weight."),
"",
"**Table 3.** Coefficients on the link scale for the best-supported model and",
"averaged over the models within two AICc units, with unconditional standard",
"errors following Buckland et al. (1997). Continuous covariates were centred and",
"scaled before fitting, so coefficients are comparable in magnitude.",
"",
if (!is.null(dm$dyn_tab))
  paste0("**Table 4.** Recruitment (\u03b3) and apparent survival (\u03c9) ",
         "estimated by the Dail-Madsen model, per sampling occasion of ",
         dm$OCC_LEN, " days, with the retention implied over the ",
         dm$n_transitions, " transitions of the survey."),
if (!is.null(dm$dyn_tab)) "",
paste0("**Table S1.** Every model fitted, for both analyses, ordered by AICc, ",
       "including the fits on the abundance null used to settle the detection ",
       "submodel among the five structures of the study, and the models that ",
       "were not fitted or did not produce usable estimates."),
"",
paste0("**Table S2.** Sensitivity of the null Dail-Madsen model, ",
       "\u03bb(~1) \u03b3(~1) \u03c9(~1) p(~1), and of the null N-mixture model, ",
       "\u03bb(~1) p(~1), to the abundance distribution (Poisson, negative ",
       "binomial) and to the integration bound K. The ", DM_LAB, " model was ",
       "retained, with the ", nul$FAMILY, " distribution; the bound retained on the ",
       "null model was K = ", nul$K_RETAINED, " and the bound used for the ",
       "candidate set was K = ", dm$K, "."),
"",
paste0("**Table S3.** Sensitivity of the null Royle-Nichols model, ",
       "\u03bb(~1) r(~1), to the integration bound K. The bound retained on the ",
       "null model was K = ", rn$K_RETAINED, " and the bound used for the ",
       "candidate set was K = ", rn$K, "."),
"",
"**Table S4.** Model-averaged abundance predicted at each camera-trap station by",
"the two analyses, with survey effort and observed counts for reference.",
"",
paste0("**Table S5.** Likelihood-ratio test of population closure. The ",
       "Dail-Madsen model is compared with the N-mixture model of the same ",
       "abundance and detection structure, fitted to the same counts at the same ",
       "distribution and integration bound K; the first row compares the null ",
       "models and decides the framework carried forward, the second, when present, ",
       "the best-supported open model with its closed counterpart. Under closure ",
       "\u03c9 = 1 and \u03b3 = 0 lie on the boundary of the parameter space, so ",
       "the likelihood ratio (LR) is referred to the chi-bar-square distribution ",
       "of Self and Liang (1987), with weights set by the correlation ",
       "\u03c1(\u03b3, \u03c9) between the estimators of recruitment and apparent ",
       "survival; P, \u03c7\u0304\u00b2 is the reference p-value at \u03b1 = ",
       nul$ALPHA, ". P, \u03c7\u00b2\u2082 is the chi-square test on two degrees ",
       "of freedom, conservative on the boundary, reported for comparison. ",
       "\u0394AICc, AICc of the closed minus the open model; \u03b3, recruitment ",
       "(individuals per station per occasion); \u03c9, apparent survival per ",
       "occasion; logL, log-likelihood."),
"",
"")

writeLines(cap, file.path(DIR_PAP, "captions.md"), useBytes = FALSE)
written <- c(written, "captions.md")
cat("  captions.md\n\n")

# --- REMOVE WHAT THIS RUN DID NOT PRODUCE ----------------------------------

# Only the files this step owns are considered: the tables and their captions.
# Anything written by another step is left alone.
present <- list.files(DIR_PAP, pattern = "^(Table_[1-4]_|Table_S[1-5]_|captions\\.md$)")
stale <- setdiff(present, written)
if (length(stale)) {
  file.remove(file.path(DIR_PAP, stale))
  cat("-- removed from an earlier run --\n")
  for (f in stale) cat("  ", f, "\n", sep = "")
  cat("\n")
}

cat("-- output directory --\n")
cat("  ", DIR_PAP, "\n\n", sep = "")
cat("STEP 6 complete.\n\n")
