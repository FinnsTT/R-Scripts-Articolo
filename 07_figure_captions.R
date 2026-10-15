# 07_figure_captions.R  --  STEP 7: figure captions
#
# Brown hyena (Parahyaena brunnea) abundance analysis -- Fish River Canyon,
# Canyon Nature Park, Namibia
#
# Writes the caption of every figure, numbered in the order the figures appear
# in the manuscript and in the supplementary material, with the counts and the
# settings read from the fitted objects rather than written by hand. Figures
# carry no text of their own, so the caption is the only place a reader is told
# what is plotted.
#
# Input : workspace/01_import.rds, 01b_occasion.rds, 02_dm_null.rds,
#         03_dm_selection.rds, 04_rn.rds
# Output: figures/figure_captions.md, figure_captions.docx
# ===========================================================================

rm(list = ls())

# --- CONFIGURATION ---------------------------------------------------------

ROOT_OUT <- "C:/Users/Dell/Desktop/codici articolo"

DIR_WS  <- file.path(ROOT_OUT, "workspace")
DIR_FIG <- file.path(ROOT_OUT, "figures")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# The map in the main text is available with and without the estimate printed
# at each station. Setting this to TRUE assumes the printed version is used.
MAP_WITH_VALUES <- FALSE

imp <- readRDS(file.path(DIR_WS, "01_import.rds"))
occ <- readRDS(file.path(DIR_WS, "01b_occasion.rds"))
nul <- readRDS(file.path(DIR_WS, "02_dm_null.rds"))
dm  <- readRDS(file.path(DIR_WS, "03_dm_selection.rds"))
rn  <- readRDS(file.path(DIR_WS, "04_rn.rds"))

has_basemap <- file.exists(file.path(DIR_WS, "05_basemap.rds"))

cat("\n=========================================================\n")
cat(" STEP 7 -- figure captions\n")
cat("=========================================================\n\n")

# --- QUANTITIES THE CAPTIONS QUOTE -----------------------------------------

nS        <- dm$nS
n_occ     <- dm$n_occ
occ_len   <- dm$OCC_LEN
occ_grid  <- occ$OCC_GRID
stab_pct  <- 100 * occ$STAB_TOL
fam       <- nul$FAMILY
k_null    <- nul$K_RETAINED
k_count   <- dm$K
k_rn_null <- rn$K_RETAINED
k_rn      <- rn$K
count_lab <- dm$MODEL_LABEL
count_pfx <- dm$PFX
open_fw   <- dm$FRAMEWORK == "open"

# Model on whose null factorial the distribution and the bound were chosen.
struct_open <- identical(nul$STRUCTURE_FROM, "open")
struct_lab  <- if (struct_open) "Dail-Madsen" else "N-mixture"
struct_pfx  <- if (struct_open) "dm" else "nm"

fam_word <- function(f) if (identical(f, "NegBin")) "negative binomial" else f

drift_nb <- nul$drift_tab[nul$drift_tab$model == nul$STRUCTURE_FROM &
                            nul$drift_tab$family == "NegBin", , drop = FALSE]

det_words <- function(f) {
  f <- gsub("^~", "", f)
  if (f == "1") return("constant")
  f <- gsub("effort", "survey effort", f)
  f <- gsub("\\bmoon\\b", "moon illumination", f)
  f <- gsub("\\bhuman\\b", "human activity", f)
  gsub(" \\+ ", " and ", f)
}

# --- CAPTIONS --------------------------------------------------------------
# Each entry is the figure label, the file it refers to, and the caption.

fig <- list()
add <- function(label, file, text) fig[[length(fig) + 1]] <<-
  list(label = label, file = file, text = text)

## main text ----------------------------------------------------------------

add("Figure 1", "not produced by this analysis",
    "Study area, camera-trap stations and survey design.")

add("Figure 2", "05_marginal_curves.(png|pdf)",
    paste0("Marginal effect of each covariate retained in the best-supported ",
           "models on predicted abundance: (a) the ", count_lab, " model, ",
           "(b) the Royle-Nichols model. Curves are averaged over the models ",
           "within two AICc units and shaded bands are 95% confidence ",
           "intervals from the unconditional standard errors. One covariate ",
           "is varied over its observed range while the others are held at ",
           "their mean, and the horizontal axis is on the original scale. ",
           "Tick marks show the observed station values."))

add("Figure 3",
    if (MAP_WITH_VALUES) "05_station_estimates_map_labelled.(png|pdf)"
    else "05_station_estimates_map.(png|pdf)",
    paste0("Abundance predicted at each of the ", nS, " camera-trap stations ",
           "by the ", count_lab, " model (left) and the Royle-Nichols model ",
           "(right), averaged over the models within two AICc units. Point ",
           "size and colour are proportional to the predicted value",
           if (MAP_WITH_VALUES) ", which is printed beside each station" else "",
           ".", if (has_basemap) " Background imagery, Esri World Imagery." else ""))

## supplementary -------------------------------------------------------------

sup <- list()
adds <- function(label, file, text) sup[[length(sup) + 1]] <<-
  list(label = label, file = file, text = text)

adds("Figure S1", "01b_occasion_lambda.(png|pdf)",
     paste0("Abundance estimated by the null Royle-Nichols model and initial ",
            "abundance estimated by the null Dail-Madsen model with Poisson ",
            "abundance, as a function of the sampling occasion length (",
            paste(occ_grid, collapse = ", "), " days), with 95% confidence ",
            "intervals where the standard errors could be computed. The dashed ",
            "line marks the occasion length retained, ", occ_len, " days: the ",
            "shortest at which both estimates change by less than ", stab_pct,
            "% at the next length."))

adds("Figure S2", paste0("02_", struct_pfx, "_K_lambda.(png|pdf)"),
     paste0(if (struct_open) "Initial abundance" else "Abundance",
            " estimated by the null ", struct_lab, " model as a function of ",
            "the integration bound K, for the Poisson and the negative ",
            "binomial abundance distributions, with 95% confidence bands. ",
            if (!struct_open)
              paste0("The null Dail-Madsen model did not yield usable ",
                     "estimates at any bound, so the distribution and the ",
                     "bound were chosen on the null N-mixture model. ")
            else "",
            "The ", fam_word(fam), " distribution was retained because its ",
            "estimate does not move with the bound",
            if (nrow(drift_nb) && !isTRUE(drift_nb$identifiable) &&
                is.finite(drift_nb$drift))
              paste0("; the negative binomial estimate tracks the bound, ",
                     "changing by a factor of ", sprintf("%.2f", drift_nb$drift),
                     " across the grid, and was not carried further")
            else "",
            "."))

adds("Figure S3", paste0("02_", struct_pfx, "_K_detection.(png|pdf)"),
     paste0("Detection probability estimated by the null ", struct_lab,
            " model as a function of the integration bound K, for the two ",
            "abundance distributions, with 95% confidence bands."))

adds("Figure S4", "04_rn_K_lambda.(png|pdf)",
     paste0("Abundance estimated by the null Royle-Nichols model as a function ",
            "of the integration bound K, with the 95% confidence band. The ",
            "bound retained on the null model was K = ", k_rn_null,
            if (k_rn != k_rn_null)
              paste0("; the candidate set was fitted at K = ", k_rn,
                     ", the smallest bound at which the truncated tail is ",
                     "negligible for the largest abundance predicted at a station")
            else "",
            "."))

adds("Figure S5", paste0("03_", count_pfx, "_aicc_weights.(png|pdf)"),
     paste0("Akaike weights of the candidate set for the ", count_lab,
            " model. Filled bars are the models within two AICc units."))

adds("Figure S6", "04_rn_aicc_weights.(png|pdf)",
     paste0("Akaike weights of the candidate set for the Royle-Nichols model. ",
            "Filled bars are the models within two AICc units."))

adds("Figure S7", "05_station_estimates_interval.(png|pdf)",
     paste0("Abundance predicted at each station by the two analyses, with ",
            "95% confidence intervals from the unconditional standard ",
            "errors. Stations are ordered by the estimate of the ", count_lab,
            " model."))

adds("Figure S8", "05_station_estimates_agreement.(png|pdf)",
     paste0("Abundance predicted at each station by the two analyses, plotted ",
            "against one another. Bars are 95% confidence intervals on each ",
            "axis, point size is the number of individuals counted at the ",
            "station, and the dashed line is the line of equality."))

# --- WRITE -----------------------------------------------------------------

fmt <- function(x) c(paste0("**", x$label, ".** ", x$text),
                     paste0("    ", x$file), "")

md <- c("# Figure captions",
        "",
        paste0("Brown hyena (*Parahyaena brunnea*), Fish River Canyon. ", nS,
               " camera-trap stations, ", n_occ, " sampling occasions of ",
               occ_len, " days."),
        paste0(count_lab, " model: ", fam_word(fam), ", K = ", k_count,
               ", detection ", det_words(dm$P_FORM),
               ". Royle-Nichols model: K = ", k_rn,
               ", per-individual detection ", det_words(rn$P_FORM), "."),
        "",
        "## Main text", "",
        unlist(lapply(fig, fmt)),
        "## Supplementary material", "",
        unlist(lapply(sup, fmt)))

writeLines(md, file.path(DIR_FIG, "figure_captions.md"), useBytes = FALSE)

has_word <- requireNamespace("officer", quietly = TRUE)

if (has_word) {
  doc <- officer::read_docx()
  para <- function(txt, bold = FALSE, size = 10, after = 6, italic = FALSE) {
    fp <- officer::fp_text(font.family = "Calibri", font.size = size,
                           bold = bold, italic = italic)
    officer::fpar(officer::ftext(txt, fp),
                  fp_p = officer::fp_par(padding.bottom = after))
  }
  doc <- officer::body_add_fpar(doc, para("Figure captions", TRUE, 16, 10))
  doc <- officer::body_add_fpar(doc, para(gsub("\\*", "", md[3]), FALSE, 9, 2, TRUE))
  doc <- officer::body_add_fpar(doc, para(md[4], FALSE, 9, 12, TRUE))
  for (section in list(list("Main text", fig),
                       list("Supplementary material", sup))) {
    doc <- officer::body_add_fpar(doc, para(section[[1]], TRUE, 12, 8))
    for (x in section[[2]]) {
      doc <- officer::body_add_fpar(
        doc, officer::fpar(
          officer::ftext(paste0(x$label, ". "),
                         officer::fp_text(font.family = "Calibri",
                                          font.size = 10, bold = TRUE)),
          officer::ftext(x$text,
                         officer::fp_text(font.family = "Calibri",
                                          font.size = 10)),
          fp_p = officer::fp_par(padding.bottom = 2)))
      doc <- officer::body_add_fpar(doc, para(x$file, FALSE, 8, 10))
    }
  }
  print(doc, target = file.path(DIR_FIG, "figure_captions.docx"))
}

cat("-- captions written --\n")
cat("  main text     : ", length(fig), " figures\n", sep = "")
cat("  supplementary : ", length(sup), " figures\n", sep = "")
for (x in c(fig, sup)) cat("  ", x$label, "  ", x$file, "\n", sep = "")
cat("\n")
cat("  figures/figure_captions.md\n")
if (has_word) cat("  figures/figure_captions.docx\n")
if (!has_word)
  cat("  install officer for the Word version\n")
cat("\nSTEP 7 complete.\n\n")
