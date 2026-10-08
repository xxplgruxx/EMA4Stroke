# =============================================================================
# ORPheoS / EMA4Stroke -- psychometric analysis pipeline
# Manuscript: "Beyond the Clinic: Towards Post-Stroke Cognition Monitoring with
# the Oldenburg Test Battery for Remote Digital Phenotyping" (Rotermund et al.)
# =============================================================================
#
# This script reproduces every analysis, table and figure reported in the
# manuscript and its supplement from the curated baseline data set
# (ORPheoS_Baseline_Curated.feather), which is created by
# EMA4Stroke_ORPheoS_Preprocessing.R from the study data.
#
# Usage: set data_dir in Section 0.1 to the folder containing the curated data
# set and run the script from top to bottom. All output is written to the
# subfolder "output" of data_dir.
#
# Sections follow the order of the manuscript:
#
#   0  Configuration, packages and helper functions
#   1  Data import (curated data set)
#   2  Task scoring (Methods 2.5; Results 3.2)
#   3  Analysed sample (Results 3.1; Supplementary Figure 1; Suppl. Tables 1-3)
#   4  Item-bank properties (Suppl. Tables 5-6; Suppl. Figures 3-4)
#   5  Differential item functioning (Suppl. Tables 7-8)
#   6  Reliability and standard error of measurement (Suppl. Table 2)
#   7  Known-groups validity and robustness (Figure 3; Suppl. Tables 9-15)
#   8  Convergent and divergent validity (Figure 2; Suppl. Table 16)
#   9  Exploratory criterion validity (Suppl. Table 17)
#   10 Structural validity: S-1 bifactor model (Figure 4; Suppl. Table 18)
#   11 Export of all reported tables (HTML overview)
#
# Terminology used in object names:
#   Patient / Control      group membership
#   Phase1 / Phase2        control recruitment phase (2023 / 2026)
#   Fam / Recoll           DRT dimensions: familiarity-based recognition /
#                          recollection-based discrimination
#   theta_*                ability score per domain (EAP for GRT, DRT and PAL;
#                          WLE-based residual for CI; mixed-model based for
#                          SART-ED)
#   *_A / *_B              split-half parcels (odd/even or interleaved items)
#
# Output files are named after the manuscript element they belong to
# (e.g. SupplTable09_KnownGroups.csv, Figure4_Bifactor_PathDiagram.svg).

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
# --- 0.1 Data location -------------------------------------------------------
# The curated data set contains one row per analysed participant (234
# patients, 280 controls) with all variables used here: group and control
# recruitment phase, demographics, clinical variables, item responses and the
# SART-ED trials in wide format. Participant exclusions and harmonisation are
# done in EMA4Stroke_ORPheoS_Preprocessing.R.
data_dir <- "path/to/data"   # folder containing ORPheoS_Baseline_Curated.feather

curated_file <- file.path(data_dir, "ORPheoS_Baseline_Curated.feather")
if (!file.exists(curated_file)) stop(sprintf("Curated data set not found at:\n  %s\nSet data_dir in Section 0.1.", curated_file))
output_dir <- file.path(data_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 0.2 Run-time switches ----------------------------------------------------
# The three analyses below dominate the run time (repeated IRT model fits).
# Setting a switch to FALSE skips the analysis and its supplementary table;
# all other results are unaffected.
RUN_DIF_ANALYSIS           <- TRUE   # Supplementary Tables 7-8
RUN_LATENT_MEAN_ROBUSTNESS <- TRUE   # Supplementary Table 10
RUN_INTERCHANGEABILITY     <- TRUE   # Supplementary Table 6 (random item halves)

# --- 0.3 Packages -------------------------------------------------------------
# DiagrammeR and DiagrammeRsvg are only needed to render the participant flow
# diagram; without them the diagram is written as DOT source.
packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble",
  "arrow", "ggplot2", "GGally",
  "mirt", "TAM", "psych", "GPArotation", "Gifi", "lavaan",
  "lme4", "effsize", "rlang", "rsvg"
)
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(sprintf("Missing package(s): %s\nInstall with: install.packages(c(%s))",
               paste(missing_packages, collapse = ", "),
               paste0('"', missing_packages, '"', collapse = ", ")))
}
invisible(lapply(packages, library, character.only = TRUE))

set.seed(20260821)

# --- 0.4 Plot settings and helper functions ----------------------------------
# assign_group(), assign_control_phase()
#                       Group and control recruitment phase per study ID,
#                       looked up from the curated data set (defined in
#                       Section 1, after loading).
# rasch_theta_eap()     EAP ability estimates from a unidimensional Rasch model
#                       (used for split-half parcels).
# iterative_mnsq_trim() Iterative item-fit trimming: the item with the largest
#                       deviation of Infit or Outfit MNSQ from 1 outside
#                       0.5-1.5 (Wright & Linacre, 1994) is removed and the
#                       model refitted until all items fit. One item per step,
#                       because removing an item changes the fit of the others.
# compare_groups()      Welch's t-test and Hedges' g with 95% CI.
# spearman_brown()      Split-half reliability: Pearson r between parcels,
#                       Spearman-Brown corrected to full length.
ctrl_color <- "#004E9F"
pat_color  <- "#E06A4E"

theme_corr <- theme_minimal(base_family = "sans") +
  theme(
    text               = element_text(size = 12),
    plot.title         = element_text(size = 14, face = "bold", hjust = 0),
    panel.grid         = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 11, color = "black"),
    plot.background    = element_rect(fill = "transparent", color = NA),
    panel.background   = element_rect(fill = "transparent", color = NA),
    strip.background   = element_rect(fill = "transparent", color = NA),
    strip.text         = element_text(size = 10, face = "bold", color = "black"),
    legend.background  = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

fmt_p   <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "< .001", paste0("= ", sub("^0\\.", ".", sprintf("%.3f", p)))))

rasch_theta_eap <- function(item_mat) {
  item_mat <- as.matrix(item_mat)
  has_data <- rowSums(!is.na(item_mat)) > 0
  theta <- rep(NA_real_, nrow(item_mat))
  if (sum(has_data) < 10) return(theta)
  ok_col <- apply(item_mat[has_data, , drop = FALSE], 2, function(x) length(unique(na.omit(x))) >= 2)
  if (sum(ok_col) < 2) return(theta)
  sub <- item_mat[has_data, ok_col, drop = FALSE]
  mod <- tryCatch(
    mirt::mirt(sub, 1, itemtype = "Rasch", verbose = FALSE, SE = FALSE, technical = list(NCYCLES = 1000)),
    error = function(e) NULL
  )
  if (is.null(mod)) return(theta)
  fs <- mirt::fscores(mod, method = "EAP", full.scores = TRUE)
  theta[has_data] <- as.numeric(fs[, 1])
  theta
}

iterative_mnsq_trim <- function(item_mat, mnsq_lower = 0.5, mnsq_upper = 1.5, min_items = 3,
                                control = list(snodes = 1000, qmc = TRUE)) {
  extract_item_table <- function(mod, item_names) {
    fit <- as.data.frame(TAM::msq.itemfit(mod)$itemfit)
    infit_col  <- grep("^Infit$",  names(fit), value = TRUE)[1]
    outfit_col <- grep("^Outfit$", names(fit), value = TRUE)[1]
    item_col   <- grep("^item$|^parameter$", names(fit), value = TRUE)[1]
    xsi <- mod$xsi
    tibble::tibble(
      item       = fit[[item_col]],
      difficulty = xsi$xsi[match(fit[[item_col]], item_names)],
      se         = xsi$se.xsi[match(fit[[item_col]], item_names)],
      infit      = fit[[infit_col]],
      outfit     = fit[[outfit_col]]
    )
  }
  
  current_items <- colnames(item_mat)
  dropped <- character(0)
  fit_before <- NULL
  fit_after  <- NULL
  step <- 0
  repeat {
    step <- step + 1
    if (length(current_items) < min_items) {
      warning("Iterative MNSQ trimming stopped early (< min_items remaining).")
      break
    }
    sub <- item_mat[, current_items, drop = FALSE]
    mod <- TAM::tam.mml(as.matrix(sub), irtmodel = "1PL", control = control, verbose = FALSE)
    item_table <- extract_item_table(mod, current_items)
    if (step == 1) fit_before <- item_table
    fit_after <- item_table
    
    misfits <- item_table[item_table$infit > mnsq_upper | item_table$infit < mnsq_lower |
                            item_table$outfit > mnsq_upper | item_table$outfit < mnsq_lower, ]
    if (nrow(misfits) == 0) break
    misfits$dev <- pmax(abs(misfits$infit - 1), abs(misfits$outfit - 1))
    worst <- misfits$item[which.max(misfits$dev)]
    dropped <- c(dropped, worst)
    current_items <- setdiff(current_items, worst)
  }
  list(keep = current_items, dropped = dropped, fit_before = fit_before, fit_after = fit_after)
}

compare_groups <- function(theta, group, group_levels = c("Control", "Patient")) {
  d <- tibble::tibble(theta = theta, Group = factor(group, levels = group_levels)) %>%
    dplyr::filter(is.finite(theta), !is.na(Group))
  if (dplyr::n_distinct(d$Group) < 2) return(NULL)
  desc <- d %>% dplyr::group_by(Group) %>% dplyr::summarise(N = dplyr::n(), M = mean(theta), SD = sd(theta), .groups = "drop")
  
  tt <- tryCatch(t.test(theta ~ Group, data = d), error = function(e) NULL)
  gg <- tryCatch(effsize::cohen.d(theta ~ Group, data = d, hedges.correction = TRUE), error = function(e) NULL)
  if (is.null(tt)) {
    wt <- tryCatch(wilcox.test(theta ~ Group, data = d), error = function(e) NULL)
    return(list(desc = desc, t = NA_real_, df = NA_real_,
                p = if (!is.null(wt)) wt$p.value else NA_real_,
                g = if (!is.null(gg)) unname(gg$estimate) else NA_real_,
                g_ci = if (!is.null(gg)) gg$conf.int else c(NA_real_, NA_real_),
                status = "t-test failed (near-zero variance in at least one group, likely a ceiling/floor effect) -- p is from a Wilcoxon rank-sum test instead; g may be NA or unreliable for the same reason"))
  }
  list(desc = desc, t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
       g = if (!is.null(gg)) unname(gg$estimate) else NA_real_,
       g_ci = if (!is.null(gg)) gg$conf.int else c(NA_real_, NA_real_),
       status = "ok")
}

spearman_brown <- function(x, y) {
  r <- suppressWarnings(cor(x, y, use = "complete.obs"))
  (2 * r) / (1 + r)
}

# =============================================================================
# 1. DATA IMPORT (CURATED DATA SET)
# =============================================================================
# One row per analysed participant. Group and control recruitment phase are
# stored explicitly; the functions below look them up by study ID.
cat(sprintf("[Data] Reading: %s\n", curated_file))
df <- arrow::read_feather(curated_file, mmap = FALSE)
cat(sprintf("[Data] %d participants (%d patients, %d controls: %d phase 1, %d phase 2), %d variables.\n",
            nrow(df), sum(df$Group == "Patient"), sum(df$Group == "Control"),
            sum(df$Phase %in% "Phase1"), sum(df$Phase %in% "Phase2"), ncol(df)))
stopifnot(!anyDuplicated(df$SERIAL), all(df$Group %in% c("Patient", "Control")))
group_lookup <- stats::setNames(df$Group, df$SERIAL)
phase_lookup <- stats::setNames(df$Phase, df$SERIAL)
assign_group         <- function(serial) unname(group_lookup[as.character(serial)])
assign_control_phase <- function(serial, group = NULL) unname(phase_lookup[as.character(serial)])

# Clinical variables of patients with an entry in the aetiology file.
stroke_aetiology <- if ("has_aetiology" %in% names(df)) {
  df %>% dplyr::filter(has_aetiology %in% TRUE) %>%
    dplyr::select(SERIAL, NIHSS_score, Hemisphere_clean, haemorrhagic, Loc_Anterior, Loc_Posterior, Loc_Media)
} else NULL

get_union_items <- function(d1, d2, cols) {
  i1 <- cols[colMeans(!is.na(d1[, cols, drop = FALSE])) > 0]
  i2 <- cols[colMeans(!is.na(d2[, cols, drop = FALSE])) > 0]
  union(i1, i2)
}

item_fit_report <- list()
add_item_fit <- function(task_label, dimension_label, fit_before, fit_after, dropped) {
  item_fit_report[[length(item_fit_report) + 1]] <<- dplyr::bind_rows(
    if (!is.null(fit_before)) dplyr::mutate(fit_before, Task = task_label, Dimension = dimension_label, stage = "before", dropped = item %in% dropped),
    if (!is.null(fit_after))  dplyr::mutate(fit_after,  Task = task_label, Dimension = dimension_label, stage = "after",  dropped = FALSE)
  )
}

# =============================================================================
# 2. TASK SCORING
# =============================================================================
# Each task is scaled separately on the pooled sample of patients and controls
# so that all participants share one metric (joint calibration). Under the
# planned-missingness design, every participant answered only part of each
# item bank; marginal maximum likelihood estimation uses all available
# responses, and items not administered are treated as missing at random.
# For every score, split-half parcels (A/B) are derived for the reliability
# estimates (Section 6) and the bifactor model (Section 10).
# =============================================================================
# 2.1 GRT -- GEOMETRIC REASONING TASK (reasoning)
# =============================================================================
# Dimensionality: principal component analysis of the tetrachoric correlation
# matrix. Four items with substantial secondary loadings are removed, then
# items are trimmed iteratively by Infit/Outfit MNSQ (0.5-1.5).
# Model choice: the Rasch model is retained because equal discrimination is
# the basis of interchangeable item subsets; the comparison with the 2PL model
# (AIC/BIC) is reported in Supplementary Table 5.
# Scores: EAP estimates. Parcels: odd vs. even item numbers.

grt_cols_raw <- grep("^GRT", names(df), value = TRUE)
df[grt_cols_raw] <- lapply(df[grt_cols_raw], function(v) suppressWarnings(as.numeric(as.character(v))))

is_binary_col <- function(x) { ux <- unique(na.omit(x)); length(ux) > 0 && all(ux %in% c(0, 1)) }
grt_item_cols <- grt_cols_raw[vapply(df[grt_cols_raw], is_binary_col, logical(1))]

grt_df       <- dplyr::filter(df, !is.na(Group))
grt_patients <- dplyr::filter(grt_df, Group == "Patient")
grt_controls  <- dplyr::filter(grt_df, Group == "Control")
grt_common   <- get_union_items(grt_patients, grt_controls, grt_item_cols)

grt_combined <- dplyr::bind_rows(
  dplyr::select(grt_patients, SERIAL, dplyr::all_of(grt_common)),
  dplyr::select(grt_controls,  SERIAL, dplyr::all_of(grt_common))
) %>% dplyr::arrange(SERIAL)

grt_item_var <- sapply(grt_combined[, grt_common, drop = FALSE], var, na.rm = TRUE)
grt_keep_var <- names(grt_item_var[grt_item_var > 0])
grt_X <- dplyr::select(grt_combined, SERIAL, dplyr::all_of(grt_keep_var))
grt_Y <- dplyr::select(grt_X, -SERIAL)

grt_tet <- suppressWarnings(psych::tetrachoric(grt_Y, na.rm = TRUE))$rho
grt_pca <- psych::pca(grt_tet, nfactors = 2, rotate = "none")
cat(sprintf("[GRT] PCA on %d items: eigenvalues %.2f and %.2f; first component explains %.1f%% of the variance.\n",
            ncol(grt_Y), grt_pca$values[1], grt_pca$values[2], 100 * grt_pca$Vaccounted["Proportion Var", 1]))

grt_manual_drop <- c("GRT30", "GRT10", "GRT24", "GRT13")
grt_stage2 <- setdiff(colnames(grt_Y), grt_manual_drop)
grt_Y2 <- grt_Y[, grt_stage2, drop = FALSE]

grt_rows_any <- rowSums(!is.na(grt_Y2)) > 0
grt_Y2_fit   <- grt_Y2[grt_rows_any, , drop = FALSE]
grt_trim <- iterative_mnsq_trim(grt_Y2_fit)
grt_final_items <- grt_trim$keep
add_item_fit("GRT", NA_character_, grt_trim$fit_before, grt_trim$fit_after, grt_trim$dropped)
grt_Y_final <- grt_Y2_fit[, grt_final_items, drop = FALSE]

cat(sprintf("[GRT] %d items entered screening -> %d after manual exclusion -> %d after MNSQ trim (%d dropped).\n",
            ncol(grt_Y), ncol(grt_Y2), length(grt_final_items), length(grt_trim$dropped)))

grt_mod_1pl_cmp <- mirt::mirt(grt_Y_final, 1, itemtype = "Rasch", verbose = FALSE)
grt_mod_2pl     <- mirt::mirt(grt_Y_final, 1, itemtype = "2PL", verbose = FALSE)
invisible(utils::capture.output(grt_comp <- anova(grt_mod_1pl_cmp, grt_mod_2pl, verbose = FALSE)))
grt_delta_aic <- grt_comp$AIC[2] - grt_comp$AIC[1]
grt_delta_bic <- grt_comp$BIC[2] - grt_comp$BIC[1]

grt_mod_1pl <- mirt::mirt(grt_Y_final, 1, itemtype = "Rasch", SE = TRUE, verbose = FALSE)
grt_rel <- mirt::marginal_rxx(grt_mod_1pl)

grt_theta_all <- rep(NA_real_, nrow(grt_X))
grt_theta_all[grt_rows_any] <- mirt::fscores(grt_mod_1pl, method = "EAP", full.scores = TRUE)[, 1]

grt_scores <- tibble::tibble(SERIAL = grt_X$SERIAL, theta = grt_theta_all) %>%
  dplyr::left_join(df %>% dplyr::select(SERIAL, Group) %>% dplyr::distinct(), by = "SERIAL") %>%
  dplyr::filter(!is.na(theta))

grt_group_comp <- compare_groups(grt_scores$theta, grt_scores$Group)

cat(sprintf("[GRT] Marginal reliability rxx = %.2f. Group comparison: t(%.1f) = %.2f, p %s, Hedges g = %.2f [%.2f, %.2f].\n",
            grt_rel, grt_group_comp$df, grt_group_comp$t, fmt_p(grt_group_comp$p),
            grt_group_comp$g, grt_group_comp$g_ci[1], grt_group_comp$g_ci[2]))

grt_nums <- as.integer(gsub("[^0-9]", "", grt_final_items))
grt_A_items <- grt_final_items[grt_nums %% 2 == 1]
grt_B_items <- grt_final_items[grt_nums %% 2 == 0]

grt_parcel_theta_A <- rasch_theta_eap(grt_X[grt_A_items])
grt_parcel_theta_B <- rasch_theta_eap(grt_X[grt_B_items])
grt_parcels <- tibble::tibble(SERIAL = grt_X$SERIAL, GRT_A = grt_parcel_theta_A, GRT_B = grt_parcel_theta_B)
grt_reliability_sb <- spearman_brown(grt_parcels$GRT_A, grt_parcels$GRT_B)

# =============================================================================
# 2.2 DRT -- DELAYED RECOGNITION TASK (two dimensions)
# =============================================================================
# Responses are dichotomised (1 = correct; codes 0 and 2 = incorrect). The
# DRT is two-dimensional: trials with a previously seen target measure
# familiarity-based recognition (Fam); trials with novel objects only,
# requiring the 'object unknown' response, measure recollection-based
# discrimination (Recoll).
# Item assignment follows the task design. An exploratory two-dimensional 2PL
# model (oblimin rotation) screens item salience (largest absolute loading
# >= .30) and identifies items whose empirical loading contradicts their
# designed dimension; these are excluded. Three items are removed after
# inspection of the item biplot. Each dimension is then trimmed by MNSQ and
# scaled with a unidimensional Rasch model (EAP scores).
# Parcels: odd vs. even item numbers within each dimension.

drt_cols_raw <- grep("^DRT", names(df), value = TRUE)
df[drt_cols_raw] <- lapply(df[drt_cols_raw], function(v) suppressWarnings(as.numeric(v)))
is_012 <- function(x) { ux <- unique(na.omit(x)); length(ux) > 0 && all(ux %in% c(0, 1, 2)) }
drt_cols <- drt_cols_raw[vapply(df[drt_cols_raw], is_012, logical(1))]

drt_patients <- dplyr::filter(df, Group == "Patient")
drt_controls  <- dplyr::filter(df, Group == "Control")
drt_common   <- get_union_items(drt_patients, drt_controls, drt_cols)
drt_combined <- dplyr::bind_rows(
  dplyr::select(drt_patients, SERIAL, dplyr::all_of(drt_common)),
  dplyr::select(drt_controls,  SERIAL, dplyr::all_of(drt_common))
) %>% dplyr::arrange(SERIAL)

drt_item_var <- sapply(drt_combined[, drt_common, drop = FALSE], var, na.rm = TRUE)
drt_keep_var <- names(drt_item_var[drt_item_var > 0])
drt_X <- dplyr::select(drt_combined, SERIAL, dplyr::all_of(drt_keep_var))

drt_Y_bin <- dplyr::select(drt_X, -SERIAL)
drt_Y_bin[drt_Y_bin == 2] <- 0

drt_pca_items <- colnames(drt_Y_bin)[vapply(drt_Y_bin, function(x) sum(!is.na(x)) >= 10 && length(unique(stats::na.omit(x))) >= 2, logical(1))]
drt_R_all <- suppressWarnings(psych::tetrachoric(as.matrix(drt_Y_bin[, drt_pca_items, drop = FALSE])))$rho
drt_R_all[!is.finite(drt_R_all)] <- 0; diag(drt_R_all) <- 1
drt_pca_all <- psych::principal(drt_R_all, nfactors = 2, rotate = "none")
cat(sprintf("[DRT] Item biplot PCA on all %d items: eigenvalues %.2f and %.2f.\n", length(drt_pca_items), drt_pca_all$values[1], drt_pca_all$values[2]))

drt_manual_drop <- intersect(c("DRT19", "DRT26", "DRT29"), colnames(drt_Y_bin))
drt_items_stage2 <- setdiff(colnames(drt_Y_bin), drt_manual_drop)
drt_Y_bin_rm <- drt_Y_bin[, drt_items_stage2, drop = FALSE]

drt_fit_idx <- rowSums(!is.na(drt_Y_bin_rm)) > 0
drt_Y_fit   <- drt_Y_bin_rm[drt_fit_idx, , drop = FALSE]
drt_mod_2d  <- mirt::mirt(drt_Y_fit, 2, itemtype = "2PL", SE = FALSE, verbose = FALSE)
drt_rot     <- summary(drt_mod_2d, rotate = "oblimin", verbose = FALSE)
drt_loadmat <- as.matrix(drt_rot$rotF)
rownames(drt_loadmat) <- colnames(drt_Y_fit)

drt_maxdim <- apply(abs(drt_loadmat), 1, which.max)
drt_maxval <- apply(abs(drt_loadmat), 1, max)

drt_true_dimension <- c(
  DRT1=2, DRT2=1, DRT3=1, DRT4=2, DRT5=2, DRT6=2, DRT7=2, DRT8=1, DRT9=2, DRT10=1,
  DRT11=2, DRT12=1, DRT13=2, DRT14=1, DRT15=2, DRT16=1, DRT17=2, DRT18=2, DRT19=1, DRT20=1,
  DRT21=1, DRT22=2, DRT23=2, DRT24=1, DRT25=2, DRT26=1, DRT27=2, DRT28=2, DRT29=1, DRT30=1,
  DRT31=2, DRT32=1, DRT33=2, DRT34=1, DRT35=2, DRT36=2, DRT37=1, DRT38=2, DRT39=1, DRT40=2
)

drt_design_dim <- function(items) unname(drt_true_dimension[paste0("DRT", as.integer(gsub("\\D", "", items)))])

drt_salient <- rownames(drt_loadmat)[drt_maxval >= 0.30]
drt_Fam_items <- drt_salient[drt_design_dim(drt_salient) %in% 1]
drt_Recoll_items <- drt_salient[drt_design_dim(drt_salient) %in% 2]

drt_Fam_items_emp <- rownames(drt_loadmat)[drt_maxdim == 2 & drt_maxval >= 0.30]
drt_Recoll_items_emp <- rownames(drt_loadmat)[drt_maxdim == 1 & drt_maxval >= 0.30]

drt_conflict_items <- drt_salient[(drt_salient %in% drt_Fam_items_emp & drt_design_dim(drt_salient) %in% 2) |
                                    (drt_salient %in% drt_Recoll_items_emp & drt_design_dim(drt_salient) %in% 1)]
if (length(drt_conflict_items) > 0) {
  cat(sprintf("[DRT] %d item(s) load contrary to their designed dimension and are excluded: %s\n",
              length(drt_conflict_items), paste(drt_conflict_items, collapse = ", ")))
  drt_Fam_items <- setdiff(drt_Fam_items, drt_conflict_items)
  drt_Recoll_items <- setdiff(drt_Recoll_items, drt_conflict_items)
}
cat(sprintf("[DRT] Design-based assignment of %d salient items: Fam = %d, Recoll = %d (%d item(s) below the .30 salience threshold excluded).\n",
            length(drt_salient), length(drt_Fam_items), length(drt_Recoll_items), nrow(drt_loadmat) - length(drt_salient)))

drt_Ybin_Fam <- drt_Y_bin_rm[, intersect(drt_Fam_items, colnames(drt_Y_bin_rm)), drop = FALSE]
drt_Ybin_Recoll <- drt_Y_bin_rm[, intersect(drt_Recoll_items, colnames(drt_Y_bin_rm)), drop = FALSE]
drt_sel_Fam <- iterative_mnsq_trim(drt_Ybin_Fam, control = list())
drt_sel_Recoll <- iterative_mnsq_trim(drt_Ybin_Recoll, control = list())
drt_final_Fam <- drt_sel_Fam$keep
drt_final_Recoll <- drt_sel_Recoll$keep
add_item_fit("DRT", "Fam", drt_sel_Fam$fit_before, drt_sel_Fam$fit_after, drt_sel_Fam$dropped)
add_item_fit("DRT", "Recoll", drt_sel_Recoll$fit_before, drt_sel_Recoll$fit_after, drt_sel_Recoll$dropped)

drt_dim_comparison <- tibble::tibble(item = drt_salient) %>%
  dplyr::mutate(
    design_dim      = drt_design_dim(item),
    exploratory_dim = dplyr::case_when(item %in% drt_Fam_items_emp ~ 1L, item %in% drt_Recoll_items_emp ~ 2L, TRUE ~ NA_integer_),
    loading_dim1    = drt_loadmat[item, 1],
    loading_dim2    = drt_loadmat[item, 2],
    excluded_design_conflict = item %in% drt_conflict_items,
    in_final_pool   = item %in% c(drt_final_Fam, drt_final_Recoll),
    agree           = exploratory_dim == design_dim
  )
readr::write_csv(drt_dim_comparison, file.path(output_dir, "Results_DRT_DimensionAssignment.csv"))
drt_dim_summary <- drt_dim_comparison %>%
  dplyr::summarise(n_salient_items = dplyr::n(), n_agree = sum(agree, na.rm = TRUE),
                   n_disagree = sum(!agree, na.rm = TRUE), pct_agree = round(100 * n_agree / n_salient_items, 1))
cat("[DRT] Design-based vs. exploratory dimension assignment -- agreement:\n")
print(drt_dim_summary)
cat("[DRT] Items where the exploratory model disagrees with the design (scored by design regardless):\n")
print(drt_dim_comparison %>% dplyr::filter(!agree) %>% dplyr::select(item, design_dim, exploratory_dim, loading_dim1, loading_dim2, in_final_pool))

cat(sprintf("[DRT] %d items entered dimensionality screening -> %d after manual exclusion. Fam final: %d items (%d dropped). Recoll final: %d items (%d dropped).\n",
            ncol(drt_Y_bin), ncol(drt_Y_bin_rm),
            length(drt_final_Fam), length(drt_sel_Fam$dropped),
            length(drt_final_Recoll), length(drt_sel_Recoll$dropped)))

fit_dim_theta <- function(Ybin_full, keep_items) {
  Ysub <- Ybin_full[, keep_items, drop = FALSE]
  idx <- rowSums(!is.na(Ysub)) > 0
  mod_1pl <- mirt::mirt(Ysub[idx, , drop = FALSE], 1, itemtype = "Rasch", SE = TRUE, verbose = FALSE)
  mod_2pl <- mirt::mirt(Ysub[idx, , drop = FALSE], 1, itemtype = "2PL", SE = FALSE, verbose = FALSE)
  invisible(utils::capture.output(cmp <- anova(mod_1pl, mod_2pl, verbose = FALSE)))
  
  daic <- cmp$AIC[2] - cmp$AIC[1]
  dbic <- cmp$BIC[2] - cmp$BIC[1]
  rel <- mirt::marginal_rxx(mod_1pl)
  theta <- rep(NA_real_, nrow(Ybin_full))
  theta[idx] <- mirt::fscores(mod_1pl, method = "EAP", full.scores = TRUE)[, 1]
  list(mod = mod_1pl, rel = rel, daic = daic, dbic = dbic, theta = theta)
}

drt_fit1 <- fit_dim_theta(drt_Y_bin_rm, drt_final_Fam)
drt_fit2 <- fit_dim_theta(drt_Y_bin_rm, drt_final_Recoll)

drt_scores <- tibble::tibble(
  SERIAL = drt_X$SERIAL,
  Group  = assign_group(drt_X$SERIAL),
  theta_DRT_Fam = drt_fit1$theta,
  theta_DRT_Recoll = drt_fit2$theta
)

drt_gc_Fam <- compare_groups(drt_scores$theta_DRT_Fam, drt_scores$Group)
drt_gc_Recoll <- compare_groups(drt_scores$theta_DRT_Recoll, drt_scores$Group)

cat(sprintf("[DRT Fam] rxx = %.2f | Hedges g = %.2f [%.2f, %.2f]\n",
            drt_fit1$rel, drt_gc_Fam$g, drt_gc_Fam$g_ci[1], drt_gc_Fam$g_ci[2]))
cat(sprintf("[DRT Recoll] rxx = %.2f | Hedges g = %.2f [%.2f, %.2f]\n",
            drt_fit2$rel, drt_gc_Recoll$g, drt_gc_Recoll$g_ci[1], drt_gc_Recoll$g_ci[2]))

drt_Fam_nums <- as.integer(gsub("[^0-9]", "", drt_final_Fam))
drt_A_items <- drt_final_Fam[drt_Fam_nums %% 2 == 1]
drt_B_items <- drt_final_Fam[drt_Fam_nums %% 2 == 0]

drt_parcel_theta_A <- rasch_theta_eap(drt_Y_bin_rm[drt_A_items])
drt_parcel_theta_B <- rasch_theta_eap(drt_Y_bin_rm[drt_B_items])
drt_parcels <- tibble::tibble(SERIAL = drt_X$SERIAL, DRT_Fam_A = drt_parcel_theta_A, DRT_Fam_B = drt_parcel_theta_B)
drt_reliability_sb <- spearman_brown(drt_parcels$DRT_Fam_A, drt_parcels$DRT_Fam_B)

drt_Recoll_nums   <- as.integer(gsub("[^0-9]", "", drt_final_Recoll))
drt_Recoll_A_items <- drt_final_Recoll[drt_Recoll_nums %% 2 == 1]
drt_Recoll_B_items <- drt_final_Recoll[drt_Recoll_nums %% 2 == 0]

drt_Recoll_parcel_theta_A <- rasch_theta_eap(drt_Y_bin_rm[drt_Recoll_A_items])
drt_Recoll_parcel_theta_B <- rasch_theta_eap(drt_Y_bin_rm[drt_Recoll_B_items])
drt_Recoll_parcels <- tibble::tibble(SERIAL = drt_X$SERIAL, DRT_Recoll_A = drt_Recoll_parcel_theta_A, DRT_Recoll_B = drt_Recoll_parcel_theta_B)
drt_Recoll_reliability_sb <- spearman_brown(drt_Recoll_parcels$DRT_Recoll_A, drt_Recoll_parcels$DRT_Recoll_B)
cat(sprintf("[DRT Recoll] Split-half reliability (Spearman-Brown corrected) = %.3f\n", drt_Recoll_reliability_sb))

# =============================================================================
# 2.3 PAL -- PAIRED ASSOCIATION LEARNING (associative memory)
# =============================================================================
# Performance per difficulty level (2, 3, 4, 6, 8 pairs) is the percentage of
# correct placements, recoded to ordered categories 0-10. A failed level is
# propagated to all higher levels, as the adaptive task ends there.
# Dimensionality: nonlinear principal component analysis for ordinal data
# (Gifi::princals). Model: partial credit model (Rasch family; generalised
# partial credit model with slopes constrained equal), compared with the graded
# response model by AIC/BIC. The five levels are not trimmed by item fit,
# because all levels are needed to retain the adaptive task structure.
# Scores: EAP estimates. Parcels: interleaved levels (1, 3, 5 vs. 2, 4),
# scored pass/fail, to balance difficulty between halves.

pal_level_names <- paste0("L", 1:5)
get_pal_level <- function(data, level_num) {
  cols <- paste0("PALV", 1:4, "L", level_num)
  cols_present <- intersect(cols, names(data))
  if (length(cols_present) == 0) return(rep(NA_real_, nrow(data)))
  purrr::reduce(data[cols_present], dplyr::coalesce)
}

pal_X <- df %>% dplyr::select(SERIAL)
for (i in 1:5) pal_X[[paste0("L", i)]] <- get_pal_level(df, i)

pal_scale_max <- suppressWarnings(max(as.matrix(pal_X[pal_level_names]), na.rm = TRUE))
if (is.finite(pal_scale_max) && pal_scale_max > 1.5) pal_X[pal_level_names] <- pal_X[pal_level_names] / 100

propagate_failure <- function(vec) {
  fail_idx <- which(!is.na(vec) & vec == 0)
  if (length(fail_idx) > 0) {
    first_fail <- min(fail_idx)
    if (first_fail < length(vec)) vec[(first_fail + 1):length(vec)] <- 0
  }
  vec
}

pal_Y_raw <- pal_X %>% dplyr::select(dplyr::all_of(pal_level_names)) %>% apply(1, propagate_failure) %>% t() %>% as.data.frame()

pal_Y_ord <- pal_Y_raw %>% dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(round(.x * 10, 0))))
pal_keep_idx <- rowSums(!is.na(pal_Y_ord)) > 0
pal_Y_ord <- pal_Y_ord[pal_keep_idx, ]
pal_serial_kept <- pal_X$SERIAL[pal_keep_idx]
pal_group_kept <- assign_group(pal_serial_kept)

pal_unique_counts <- apply(pal_Y_ord, 2, function(x) length(unique(na.omit(x))))
pal_princals_cols <- names(pal_Y_ord)[pal_unique_counts > 1]
pal_Y_princals <- pal_Y_ord[, pal_princals_cols, drop = FALSE]
pal_princals_fit <- Gifi::princals(as.data.frame(lapply(pal_Y_princals, as.ordered)), ndim = 2, ordinal = TRUE, ties = "s")
cat(sprintf("[PAL] Nonlinear PCA on %d levels: first dimension explains %.1f%% of the variance.\n",
            ncol(pal_Y_princals), pal_princals_fit$evals[1] / ncol(pal_Y_princals) * 100))

pal_spec_pcm <- mirt::mirt.model(paste0("F = 1-", ncol(pal_Y_ord), "\nCONSTRAIN = (1-", ncol(pal_Y_ord), ", a1)"))
pal_mod_pcm  <- mirt::mirt(pal_Y_ord, model = pal_spec_pcm, itemtype = "gpcm", SE = TRUE, verbose = FALSE)
pal_mod_grm  <- mirt::mirt(pal_Y_ord, 1, itemtype = "graded", SE = FALSE, verbose = FALSE)
invisible(utils::capture.output(pal_anova <- anova(pal_mod_pcm, pal_mod_grm, verbose = FALSE)))
pal_delta_aic <- pal_anova$AIC[2] - pal_anova$AIC[1]
pal_delta_bic <- pal_anova$BIC[2] - pal_anova$BIC[1]
pal_rel <- mirt::marginal_rxx(pal_mod_pcm)

pal_itemfit_raw <- mirt::itemfit(pal_mod_pcm, fit_stats = "infit")
pal_coefs <- mirt::coef(pal_mod_pcm, IRTpars = TRUE, simplify = TRUE)$items
pal_b_cols <- grep("^b", colnames(pal_coefs), value = TRUE)

pal_se_vals <- rep(NA_real_, length(pal_itemfit_raw$item))
tryCatch({
  pal_coefs_se_list <- mirt::coef(pal_mod_pcm, IRTpars = TRUE, printSE = TRUE)
  for (i in seq_along(pal_itemfit_raw$item)) {
    item_mat <- pal_coefs_se_list[[pal_itemfit_raw$item[i]]]
    if (!is.null(item_mat) && "SE" %in% rownames(item_mat)) {
      se_row <- item_mat["SE", ]
      se_b_cols <- grep("^b", names(se_row), value = TRUE)
      if (length(se_b_cols) > 0) pal_se_vals[i] <- mean(se_row[se_b_cols], na.rm = TRUE)
    }
  }
}, error = function(e) cat(sprintf("[PAL] SE extraction failed, leaving SEs as NA: %s\n", conditionMessage(e))))

pal_item_table <- tibble::tibble(
  item       = pal_itemfit_raw$item,
  difficulty = rowMeans(pal_coefs[pal_itemfit_raw$item, pal_b_cols, drop = FALSE], na.rm = TRUE),
  se         = pal_se_vals,
  infit      = pal_itemfit_raw$infit,
  outfit     = pal_itemfit_raw$outfit
)
add_item_fit("PAL", NA_character_, pal_item_table, pal_item_table, character(0))

pal_fscores <- mirt::fscores(pal_mod_pcm, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)
pal_scores <- tibble::tibble(SERIAL = pal_serial_kept, theta = pal_fscores[, "F"], Group = pal_group_kept)

pal_group_comp <- compare_groups(pal_scores$theta, pal_scores$Group)

cat(sprintf("[PAL] Partial credit model retained (graded response model minus PCM: delta AIC = %.1f, delta BIC = %.1f). rxx = %.2f. Hedges g = %.2f [%.2f, %.2f]\n",
            pal_delta_aic, pal_delta_bic, pal_rel, pal_group_comp$g, pal_group_comp$g_ci[1], pal_group_comp$g_ci[2]))

pal_bin <- pal_Y_ord
pal_bin[pal_level_names] <- lapply(pal_bin[pal_level_names], function(x) as.integer(!is.na(x) & x > 0))

pal_A_cols <- c("L1", "L3", "L5")
pal_B_cols <- c("L2", "L4")
pal_A_mat <- t(apply(as.matrix(pal_bin[pal_A_cols]), 1, propagate_failure)); colnames(pal_A_mat) <- pal_A_cols
pal_B_mat <- t(apply(as.matrix(pal_bin[pal_B_cols]), 1, propagate_failure)); colnames(pal_B_mat) <- pal_B_cols

irt_or_sum_theta <- function(mat) {
  theta <- rasch_theta_eap(mat)
  if (sum(!is.na(theta)) < 10) {
    raw_sum <- rowSums(mat, na.rm = FALSE)
    raw_sum[rowSums(!is.na(mat)) == 0] <- NA_real_
    theta <- as.numeric(scale(raw_sum))
  }
  theta
}

pal_parcels <- tibble::tibble(
  SERIAL = pal_serial_kept,
  PAL_A = irt_or_sum_theta(pal_A_mat),
  PAL_B = irt_or_sum_theta(pal_B_mat)
)
pal_reliability_sb <- spearman_brown(pal_parcels$PAL_A, pal_parcels$PAL_B)

# =============================================================================
# 2.4 CI -- COGNITIVE INTERFERENCE TASK (internal interference control)
# =============================================================================
# Item format: nine five-item testlets; baseline recall (A-C pairs) and recall
# under interference (A-D pairs) are modelled as two correlated dimensions in a
# two-dimensional Rasch model (TAM).
# Item fit: as for GRT and DRT, the misfitting item (Infit or Outfit MNSQ
# outside 0.5-1.5) with the largest deviation from 1 is removed and the model
# refitted, one item at a time, until all retained items fit.
# Score: interference recall residualised on baseline recall via the latent
# regression slope (latent covariance / latent baseline variance). Weighted
# likelihood estimates (WLE) are used because, unlike EAP estimates, they are
# not shrunk towards the mean in proportion to the number of items answered,
# which is a prerequisite for residualising one score on another under planned
# missingness.
# Model comparison: 2D vs. 1D model (AIC/BIC; Supplementary Table 5).
# Parcels: matched baseline/interference item pairs are split into odd and even
# pairs; each half is scaled and residualised exactly like the full score.

ci_cols <- grep("^CI_", names(df), value = TRUE)
ci_long <- df %>%
  dplyr::select(SERIAL, Group, dplyr::all_of(ci_cols)) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("CI_"), names_to = "ItemFull", values_to = "raw") %>%
  dplyr::filter(!is.na(raw), raw != "") %>%
  dplyr::mutate(
    raw     = as.numeric(raw),
    ver     = stringr::str_extract(ItemFull, "(?<=CI_)\\d+"),
    section = stringr::str_extract(ItemFull, "(?<=_)[A-Z]{2}(?=_)"),
    id      = stringr::str_extract(ItemFull, "[^_]+$"),
    item    = paste0(section, "_V", ver, "_I", id)
  )

ci_item_uses_one <- ci_long %>% dplyr::group_by(item) %>% dplyr::summarise(any_one = any(raw == 1, na.rm = TRUE), .groups = "drop")
ci_long <- ci_long %>%
  dplyr::left_join(ci_item_uses_one, by = "item") %>%
  dplyr::mutate(GR = dplyr::case_when(any_one ~ as.numeric(raw == 1), TRUE ~ as.numeric(raw == 0)))

ci_wide <- ci_long %>%
  dplyr::filter(section %in% c("AC", "AD")) %>%
  dplyr::group_by(SERIAL, item) %>%
  dplyr::summarise(item_score = round(mean(GR, na.rm = TRUE)), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = item, values_from = item_score) %>%
  dplyr::filter(!is.na(SERIAL))

ci_mat <- as.matrix(ci_wide[, -1]); rownames(ci_mat) <- ci_wide$SERIAL
ci_mat <- ci_mat[, apply(ci_mat, 2, function(x) length(unique(na.omit(x))) >= 2), drop = FALSE]

ci_Q <- matrix(0, nrow = ncol(ci_mat), ncol = 2, dimnames = list(colnames(ci_mat), c("Dim_AC", "Dim_AD")))
ci_Q[grepl("^AC", rownames(ci_Q)), 1] <- 1
ci_Q[grepl("^AD", rownames(ci_Q)), 2] <- 1

ci_mod_init <- TAM::tam.mml(resp = ci_mat, Q = ci_Q, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
ci_fit_init <- TAM::tam.fit(ci_mod_init, progress = FALSE)

ci_pair_key <- function(items) sub("^(AC|AD)_", "", items)

# Iterative item-fit trimming as for GRT and DRT: the misfitting item with the
# largest deviation of Infit or Outfit MNSQ from 1 is removed and the
# two-dimensional model refitted, until all retained items fit.
ci_mat_trim   <- ci_mat
ci_drop_items <- character(0)
ci_fit_it     <- ci_fit_init
repeat {
  mis <- ci_fit_it$itemfit %>% dplyr::filter(Outfit > 1.5 | Outfit < 0.5 | Infit > 1.5 | Infit < 0.5)
  if (nrow(mis) == 0) break
  worst <- as.character(mis$parameter[which.max(pmax(abs(mis$Infit - 1), abs(mis$Outfit - 1)))])
  ci_drop_items <- c(ci_drop_items, worst)
  ci_mat_trim <- ci_mat_trim[, colnames(ci_mat_trim) != worst, drop = FALSE]
  if (min(sum(grepl("^AC", colnames(ci_mat_trim))), sum(grepl("^AD", colnames(ci_mat_trim)))) < 3) {
    warning("[CI] Item-fit trimming stopped: fewer than 3 items left in one dimension.")
    break
  }
  mod_it <- TAM::tam.mml(resp = ci_mat_trim, Q = ci_Q[colnames(ci_mat_trim), , drop = FALSE], irtmodel = "1PL",
                         control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
  ci_fit_it <- TAM::tam.fit(mod_it, progress = FALSE)
}
cat(sprintf("[CI] Item-fit trimming: %d of %d items retained; removed in order: %s\n",
            ncol(ci_mat_trim), ncol(ci_mat), if (length(ci_drop_items)) paste(ci_drop_items, collapse = ", ") else "none"))
ci_ac_all <- grep("^AC", colnames(ci_mat_trim), value = TRUE)
ci_ad_all <- grep("^AD", colnames(ci_mat_trim), value = TRUE)
ci_paired_keys <- intersect(ci_pair_key(ci_ac_all), ci_pair_key(ci_ad_all))
ci_keep_items <- intersect(colnames(ci_mat_trim), c(paste0("AC_", ci_paired_keys), paste0("AD_", ci_paired_keys)))

ci_mat_final <- ci_mat_trim[, ci_keep_items, drop = FALSE]

ci_Q_trim    <- ci_Q[rownames(ci_Q) %in% colnames(ci_mat_trim), , drop = FALSE]
ci_mod_final <- TAM::tam.mml(resp = ci_mat_trim, Q = ci_Q_trim, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)

cat(sprintf("[CI] %d items entered joint calibration -> %d retained after item-fit trimming (%d dropped); %d of these form complete AC/AD pairs for the difference-score parcels.\n",
            ncol(ci_mat), ncol(ci_mat_trim), ncol(ci_mat) - ncol(ci_mat_trim), ncol(ci_mat_final)))

ci_fit_final <- TAM::tam.fit(ci_mod_final, progress = FALSE)
extract_ci_item_table <- function(mod, fit) {
  xsi <- mod$xsi
  tibble::tibble(
    item       = fit$itemfit$parameter,
    difficulty = xsi$xsi[match(fit$itemfit$parameter, rownames(xsi))],
    se         = xsi$se.xsi[match(fit$itemfit$parameter, rownames(xsi))],
    infit      = fit$itemfit$Infit,
    outfit     = fit$itemfit$Outfit
  )
}
ci_dropped_items <- setdiff(colnames(ci_mat), colnames(ci_mat_trim))
add_item_fit("CI", NA_character_,
             extract_ci_item_table(ci_mod_init, ci_fit_init),
             extract_ci_item_table(ci_mod_final, ci_fit_final),
             ci_dropped_items)

ci_mod_1d <- TAM::tam.mml(resp = ci_mat_trim, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
ci_aic_1d <- ci_mod_1d$ic$AIC
ci_aic_2d <- ci_mod_final$ic$AIC
ci_bic_1d <- ci_mod_1d$ic$BIC
ci_bic_2d <- ci_mod_final$ic$BIC

ci_wle <- TAM::tam.wle(ci_mod_final, progress = FALSE)
ci_rel_ac <- TAM::WLErel(ci_wle$theta.Dim01, ci_wle$error.Dim01)
ci_rel_ad <- TAM::WLErel(ci_wle$theta.Dim02, ci_wle$error.Dim02)

ci_latent_cov <- ci_mod_final$variance
ci_beta <- ci_latent_cov[1, 2] / ci_latent_cov[1, 1]

ci_theta <- tibble::tibble(
  SERIAL = rownames(ci_mat_trim),
  Theta_AC = ci_wle$theta.Dim01,
  Theta_AD = ci_wle$theta.Dim02
) %>%
  dplyr::left_join(df %>% dplyr::select(SERIAL, Group) %>% dplyr::distinct(), by = "SERIAL") %>%
  dplyr::filter(!is.na(Group))

ci_intercept <- mean(ci_theta$Theta_AD, na.rm = TRUE) - ci_beta * mean(ci_theta$Theta_AC, na.rm = TRUE)
ci_theta$CI_Resilience <- ci_theta$Theta_AD - (ci_intercept + ci_beta * ci_theta$Theta_AC)

ci_group_comp <- compare_groups(ci_theta$CI_Resilience, ci_theta$Group)

cat(sprintf("[CI] Latent AC-AD covariance = %.3f, beta = %.3f. Reliabilities: AC = %.2f, AD = %.2f. Hedges g = %.2f [%.2f, %.2f]\n",
            ci_latent_cov[1, 2], ci_beta, ci_rel_ac, ci_rel_ad,
            ci_group_comp$g, ci_group_comp$g_ci[1], ci_group_comp$g_ci[2]))

# Split-half parcels of the residualised score: odd vs. even matched item
# pairs, each half scaled and residualised like the full score.
ci_pair_keys_final <- sort(unique(ci_pair_key(colnames(ci_mat_final))))
ci_half_keys <- list(A = ci_pair_keys_final[seq(1, length(ci_pair_keys_final), by = 2)],
                     B = ci_pair_keys_final[seq(2, length(ci_pair_keys_final), by = 2)])
ci_half_items <- lapply(ci_half_keys, function(k) intersect(colnames(ci_mat_final), c(paste0("AC_", k), paste0("AD_", k))))
ci_residual_score <- function(mat) {
  is_ac <- grepl("^AC", colnames(mat)); is_ad <- grepl("^AD", colnames(mat))
  has_both <- rowSums(!is.na(mat[, is_ac, drop = FALSE])) > 0 & rowSums(!is.na(mat[, is_ad, drop = FALSE])) > 0
  out <- rep(NA_real_, nrow(mat))
  if (sum(has_both) < 10) return(out)
  mod <- TAM::tam.mml(resp = mat[has_both, , drop = FALSE], Q = ci_Q[colnames(mat), , drop = FALSE],
                      irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
  wle <- TAM::tam.wle(mod, progress = FALSE)
  beta  <- mod$variance[1, 2] / mod$variance[1, 1]
  alpha <- mean(wle$theta.Dim02, na.rm = TRUE) - beta * mean(wle$theta.Dim01, na.rm = TRUE)
  out[has_both] <- wle$theta.Dim02 - (alpha + beta * wle$theta.Dim01)
  out
}
ci_parcels <- tibble::tibble(
  SERIAL = rownames(ci_mat_final),
  CI_A   = ci_residual_score(ci_mat_final[, ci_half_items$A, drop = FALSE]),
  CI_B   = ci_residual_score(ci_mat_final[, ci_half_items$B, drop = FALSE])
)
ci_reliability_sb <- spearman_brown(ci_parcels$CI_A, ci_parcels$CI_B)
cat(sprintf("[CI] Split-half reliability of the residualised score (Spearman-Brown corrected) = %.3f (%d item pairs per half: %d / %d)\n",
            ci_reliability_sb, length(ci_half_keys$A), length(ci_half_keys$A), length(ci_half_keys$B)))

# =============================================================================
# 2.5 SART-ED -- SUSTAINED ATTENTION TO RESPONSE TASK WITH EXTERNAL DISTRACTION
# =============================================================================
# Trial-level data (200 trials: 140 Go, 60 No-Go; distractor faces in 60
# trials) are read from the curated data set (wide format, one set of columns
# per trial) and reshaped to long format. Go reaction times outside 150-1500 ms
# are treated as anticipations or lapses and excluded.
# Two scores are derived:
#   External interference control (2.5.2): two-step mixed-effects approach.
#     Linear (log RT, correct Go trials) and logistic (Go accuracy) mixed
#     models regress performance on condition (distractor vs. none), centred
#     baseline performance and their interaction, with random intercepts and
#     random condition slopes per participant. The empirical Bayes condition
#     slopes are each person's distraction effect adjusted for baseline; they
#     are sign-aligned (higher = more resilient) and combined by minimum
#     residual factor analysis. Mixed models are used instead of item response
#     models because consecutive trials are not locally independent.
#   Mental speed (2.5.3): mean log RT on correct Go trials without distractor,
#     after iterative intra-individual winsorising at +/-3.5 SD, z-standardised
#     so that higher values indicate faster responding.
# Phase-1 controls are not included: their SART-ED data were stored as
# aggregate scores without the trial-level information both scores require.

# --- 2.5.1 Trial data and helper functions -----------------------------------

sart_RT_MIN <- 150; sart_RT_MAX <- 1500
sart_MIN_TRIALS_FULL <- 10
sart_MIN_TRIALS_PARCEL <- 5

as01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.character(x)) return(as.integer(x %in% c("1", "TRUE", "true", "T")))
  as.integer(x)
}
safe_log <- function(x) log(pmax(x, 1e-6))

extract_eb_slope <- function(fit, slope_term = "condition") {
  fe <- lme4::fixef(fit); re <- lme4::ranef(fit)[["SERIAL"]]
  tibble::tibble(SERIAL = rownames(re), eb_slope = as.numeric(fe[[slope_term]] + re[[slope_term]]))
}

latent_score_1f <- function(df_ind) {
  X <- df_ind %>% dplyr::select(-SERIAL) %>% as.data.frame()
  ok_var <- sapply(X, function(v) stats::sd(v, na.rm = TRUE)) > 0
  X <- X[, ok_var, drop = FALSE]
  cc <- stats::complete.cases(X)
  if (ncol(X) < 2 || sum(cc) < 30) return(tibble::tibble(SERIAL = df_ind$SERIAL, theta = NA_real_))
  Xcc <- X[cc, , drop = FALSE]
  fit <- psych::fa(Xcc, nfactors = 1, fm = "minres", rotate = "none")
  scores <- as.numeric(fit$scores[, 1])
  
  # Orient the factor so that higher scores mean higher values on the first indicator
  ref_cor <- suppressWarnings(cor(scores, Xcc[[1]], use = "complete.obs"))
  reflected <- !is.na(ref_cor) && ref_cor < 0
  if (reflected) scores <- scores * -1
  theta <- rep(NA_real_, nrow(df_ind)); theta[cc] <- scores
  tibble::tibble(SERIAL = df_ind$SERIAL, theta = theta)
}

sart_cols <- grep("^SART_T\\d{3}_", names(df), value = TRUE)
sart_trials <- df %>%
  dplyr::select(SERIAL, dplyr::all_of(sart_cols)) %>%
  tidyr::pivot_longer(-SERIAL, names_to = c("Trial", ".value"), names_pattern = "^SART_T(\\d{3})_(.*)$") %>%
  dplyr::filter(!(is.na(GoNoGo) & is.na(Face) & is.na(Correct) & is.na(RT))) %>%
  dplyr::mutate(
    Trial     = as.integer(Trial),
    GoNoGo    = as01(GoNoGo),
    Correct   = as01(Correct),
    condition = dplyr::if_else(Face == 0, 0L, 1L),
    parcel    = dplyr::if_else(Trial %% 2 == 1, "A", "B"),
    Group     = assign_group(SERIAL)
  )
cat(sprintf("[SART-ED] %d trials from %d participants.\n", nrow(sart_trials), dplyr::n_distinct(sart_trials$SERIAL)))

sart_go        <- dplyr::filter(sart_trials, GoNoGo == 1L)
sart_nogo      <- dplyr::filter(sart_trials, GoNoGo == 0L)
sart_go_rt     <- sart_go %>% dplyr::filter(Correct == 1L, !is.na(RT), RT >= sart_RT_MIN, RT <= sart_RT_MAX) %>% dplyr::mutate(logRT = safe_log(RT))
sart_go_acc    <- sart_go %>% dplyr::filter(!is.na(Correct)) %>% dplyr::mutate(correct_go = as01(Correct))

sart_baseline <- sart_trials %>%
  dplyr::filter(condition == 0L) %>%
  dplyr::group_by(SERIAL) %>%
  dplyr::summarise(
    baseline_logRT = mean(safe_log(RT[GoNoGo == 1L & Correct == 1L & !is.na(RT) & RT >= sart_RT_MIN & RT <= sart_RT_MAX]), na.rm = TRUE),
    baseline_acc   = mean(Correct[GoNoGo == 1L], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(baseline_logRT_c = as.numeric(scale(baseline_logRT, scale = FALSE)),
                baseline_acc_c   = as.numeric(scale(baseline_acc,   scale = FALSE)))

make_cov <- function(dat, label) {
  dat %>% dplyr::count(SERIAL, condition, name = "n") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
    dplyr::rename(!!paste0("cond0_n_", label) := `0`, !!paste0("cond1_n_", label) := `1`)
}

# --- 2.5.2 External interference control (incl. Suppl. Table 11) ------------

sart_elig <- make_cov(sart_go_rt, "go_rt") %>%
  dplyr::full_join(make_cov(sart_go_acc, "go_acc"), by = "SERIAL") %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("cond"), ~ dplyr::coalesce(.x, 0L))) %>%
  dplyr::mutate(elig_go_rt  = cond0_n_go_rt  >= sart_MIN_TRIALS_FULL & cond1_n_go_rt  >= sart_MIN_TRIALS_FULL,
                elig_go_acc = cond0_n_go_acc >= sart_MIN_TRIALS_FULL & cond1_n_go_acc >= sart_MIN_TRIALS_FULL)

dat_rt <- sart_go_rt %>%
  dplyr::inner_join(dplyr::select(sart_elig, SERIAL, elig_go_rt), by = "SERIAL") %>%
  dplyr::filter(elig_go_rt) %>%
  dplyr::left_join(dplyr::select(sart_baseline, SERIAL, baseline_logRT_c), by = "SERIAL") %>%
  dplyr::filter(!is.na(baseline_logRT_c))
fit_go_rt <- lme4::lmer(logRT ~ condition * baseline_logRT_c + (condition | SERIAL), data = dat_rt, REML = TRUE,
                        control = lme4::lmerControl(optimizer = "bobyqa"))
eb_rt <- extract_eb_slope(fit_go_rt) %>% dplyr::rename(eb_slope_rt = eb_slope)

dat_acc <- sart_go_acc %>%
  dplyr::inner_join(dplyr::select(sart_elig, SERIAL, elig_go_acc), by = "SERIAL") %>%
  dplyr::filter(elig_go_acc) %>%
  dplyr::left_join(dplyr::select(sart_baseline, SERIAL, baseline_acc_c), by = "SERIAL") %>%
  dplyr::filter(!is.na(baseline_acc_c))
fit_go_acc <- lme4::glmer(correct_go ~ condition * baseline_acc_c + (condition | SERIAL), data = dat_acc, family = binomial(),
                          control = lme4::glmerControl(optimizer = "bobyqa"))
eb_acc <- extract_eb_slope(fit_go_acc) %>% dplyr::rename(eb_slope_acc = eb_slope)

cat(sprintf("[SART-ED] Distraction effect: RT beta = %.3f (SE %.3f), Wald t = %.1f; Accuracy beta = %.3f (SE %.3f), Wald z = %.2f\n",
            lme4::fixef(fit_go_rt)["condition"], sqrt(diag(vcov(fit_go_rt)))["condition"], summary(fit_go_rt)$coefficients["condition", "t value"],
            lme4::fixef(fit_go_acc)["condition"], sqrt(diag(vcov(fit_go_acc)))["condition"], summary(fit_go_acc)$coefficients["condition", "z value"]))

sart_indicators <- tibble::tibble(SERIAL = unique(sart_trials$SERIAL)) %>%
  dplyr::left_join(eb_rt, by = "SERIAL") %>%
  dplyr::left_join(eb_acc, by = "SERIAL") %>%
  dplyr::mutate(ind_rt_resilience = -1 * eb_slope_rt, ind_acc_resilience = eb_slope_acc)

ei_fit <- latent_score_1f(dplyr::select(sart_indicators, SERIAL, ind_rt_resilience, ind_acc_resilience))
sart_ei_theta <- ei_fit %>% dplyr::rename(SART_ExternalInterference_Theta = theta)

sart_interaction_fit <- function(fit, label) {
  co <- summary(fit)$coefficients
  rows <- intersect(c("condition", "GroupPatient", "condition:GroupPatient"), rownames(co))
  tibble::tibble(outcome = label, term = rows, estimate = co[rows, 1], se = co[rows, 2],
                 wald = co[rows, 3], p = 2 * stats::pnorm(-abs(co[rows, 3])))
}
sart_group_interaction <- tryCatch({
  dat_rt_g  <- dat_rt  %>% dplyr::mutate(Group = factor(assign_group(SERIAL), levels = c("Control", "Patient"))) %>% dplyr::filter(!is.na(Group))
  dat_acc_g <- dat_acc %>% dplyr::mutate(Group = factor(assign_group(SERIAL), levels = c("Control", "Patient"))) %>% dplyr::filter(!is.na(Group))
  fit_rt_g  <- lme4::lmer(logRT ~ condition * (baseline_logRT_c + Group) + (condition | SERIAL), data = dat_rt_g, REML = TRUE,
                          control = lme4::lmerControl(optimizer = "bobyqa"))
  fit_acc_g <- lme4::glmer(correct_go ~ condition * (baseline_acc_c + Group) + (condition | SERIAL), data = dat_acc_g, family = binomial(),
                           control = lme4::glmerControl(optimizer = "bobyqa"))
  dplyr::bind_rows(sart_interaction_fit(fit_rt_g, "Go log-RT"), sart_interaction_fit(fit_acc_g, "Go accuracy (logit)"))
}, error = function(e) { cat(sprintf("[SART-ED] Group-interaction models failed: %s\n", conditionMessage(e))); NULL })
if (!is.null(sart_group_interaction)) {
  cat("[SART-ED] Robustness -- condition x Group interaction in the baseline-adjusted mixed models:\n")
  print(sart_group_interaction)
  readr::write_csv(sart_group_interaction, file.path(output_dir, "SupplTable11_SART_GroupInteraction.csv"))
}

# --- 2.5.3 Mental speed ------------------------------------------------------
winsorize_iterative <- function(x, k = 3.5, max_iter = 10) {
  x_clean <- x
  for (i in seq_len(max_iter)) {
    m <- mean(x_clean, na.rm = TRUE)
    s <- stats::sd(x_clean, na.rm = TRUE)
    if (!is.finite(s) || s == 0) break
    lower <- m - k * s
    upper <- m + k * s
    n_capped <- sum(x_clean < lower | x_clean > upper, na.rm = TRUE)
    if (n_capped == 0) break
    x_clean <- pmin(pmax(x_clean, lower), upper)
  }
  x_clean
}

sart_ms_baseline_trials <- sart_go_rt %>%
  dplyr::filter(condition == 0L) %>%
  dplyr::group_by(SERIAL) %>%
  dplyr::mutate(logRT_winsorized = winsorize_iterative(logRT, k = 3.5), n_capped = sum(logRT_winsorized != logRT)) %>%
  dplyr::ungroup()

sart_ms_ind_winsorized <- sart_ms_baseline_trials %>%
  dplyr::group_by(SERIAL) %>%
  dplyr::summarise(ms_mean_logRT = mean(logRT_winsorized, na.rm = TRUE), .groups = "drop")

ms_clean_winsorized <- dplyr::filter(sart_ms_ind_winsorized, is.finite(ms_mean_logRT))
sart_ms_theta <- ms_clean_winsorized %>%
  dplyr::transmute(SERIAL, SART_MentalSpeed_Theta = -1 * as.numeric(scale(ms_mean_logRT)))

cat(sprintf("[SART Mental Speed] Intra-individual winsorizing (k=3.5 SD, iterative): %d of %d people had at least one trial capped (total %d trials capped across %d baseline trials).\n",
            sum(sart_ms_baseline_trials %>% dplyr::group_by(SERIAL) %>% dplyr::summarise(any_capped = dplyr::first(n_capped) > 0, .groups = "drop") %>% dplyr::pull(any_capped)),
            dplyr::n_distinct(sart_ms_baseline_trials$SERIAL),
            sum(sart_ms_baseline_trials$logRT_winsorized != sart_ms_baseline_trials$logRT),
            nrow(sart_ms_baseline_trials)))

sart_gc_ms <- compare_groups(sart_ms_theta$SART_MentalSpeed_Theta, assign_group(sart_ms_theta$SERIAL))

sart_scores <- tibble::tibble(SERIAL = unique(sart_trials$SERIAL), Group = assign_group(unique(sart_trials$SERIAL))) %>%
  dplyr::left_join(sart_ei_theta, by = "SERIAL") %>%
  dplyr::left_join(sart_ms_theta, by = "SERIAL")

sart_gc_ei <- compare_groups(sart_scores$SART_ExternalInterference_Theta, sart_scores$Group)

cat(sprintf("[SART-ED External Interference] Hedges g = %.2f [%.2f, %.2f]\n", sart_gc_ei$g, sart_gc_ei$g_ci[1], sart_gc_ei$g_ci[2]))
cat(sprintf("[SART-ED Mental Speed] Hedges g = %.2f [%.2f, %.2f]\n", sart_gc_ms$g, sart_gc_ms$g_ci[1], sart_gc_ms$g_ci[2]))

# --- 2.5.4 Split-half parcels (independent models for odd and even trials) ---

sart_parcel_indicators <- tibble::tibble(SERIAL = unique(sart_trials$SERIAL))

for (p in c("A", "B")) {
  dat_rt_p <- sart_go_rt %>% dplyr::filter(parcel == p) %>%
    dplyr::left_join(dplyr::select(sart_baseline, SERIAL, baseline_logRT_c), by = "SERIAL") %>%
    dplyr::filter(!is.na(baseline_logRT_c))
  cov_p <- make_cov(dat_rt_p, "go_rt")
  elig_p <- cov_p$SERIAL[cov_p$cond0_n_go_rt >= sart_MIN_TRIALS_PARCEL & cov_p$cond1_n_go_rt >= sart_MIN_TRIALS_PARCEL]
  fit_rt_p <- lme4::lmer(logRT ~ condition * baseline_logRT_c + (condition | SERIAL),
                         data = dplyr::filter(dat_rt_p, SERIAL %in% elig_p), REML = TRUE,
                         control = lme4::lmerControl(optimizer = "bobyqa"))
  eb_rt_p <- extract_eb_slope(fit_rt_p) %>% dplyr::rename(!!paste0("eb_rt_", p) := eb_slope)
  
  dat_acc_p <- sart_go_acc %>% dplyr::filter(parcel == p) %>%
    dplyr::left_join(dplyr::select(sart_baseline, SERIAL, baseline_acc_c), by = "SERIAL") %>%
    dplyr::filter(!is.na(baseline_acc_c))
  cov_acc_p <- make_cov(dat_acc_p, "go_acc")
  elig_acc_p <- cov_acc_p$SERIAL[cov_acc_p$cond0_n_go_acc >= sart_MIN_TRIALS_PARCEL & cov_acc_p$cond1_n_go_acc >= sart_MIN_TRIALS_PARCEL]
  fit_acc_p <- lme4::glmer(correct_go ~ condition * baseline_acc_c + (condition | SERIAL),
                           data = dplyr::filter(dat_acc_p, SERIAL %in% elig_acc_p), family = binomial(),
                           control = lme4::glmerControl(optimizer = "bobyqa"))
  eb_acc_p <- extract_eb_slope(fit_acc_p) %>% dplyr::rename(!!paste0("eb_acc_", p) := eb_slope)
  
  sart_parcel_indicators <- sart_parcel_indicators %>%
    dplyr::left_join(eb_rt_p, by = "SERIAL") %>% dplyr::left_join(eb_acc_p, by = "SERIAL") %>%
    dplyr::mutate(!!paste0("SART_", p) := rowMeans(cbind(-1 * .data[[paste0("eb_rt_", p)]], .data[[paste0("eb_acc_", p)]]), na.rm = FALSE))
}

compute_ms_half <- function(p) {
  half <- sart_go_rt %>% dplyr::filter(condition == 0L, parcel == p) %>%
    dplyr::group_by(SERIAL) %>%
    dplyr::mutate(logRT_w = winsorize_iterative(logRT, k = 3.5)) %>%
    dplyr::summarise(h_mean_logRT = mean(logRT_w, na.rm = TRUE), h_n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(h_n >= 3, is.finite(h_mean_logRT))
  half %>% dplyr::transmute(SERIAL, theta = -1 * as.numeric(scale(h_mean_logRT)))
}
ms_A <- compute_ms_half("A") %>% dplyr::rename(MS_A = theta)
ms_B <- compute_ms_half("B") %>% dplyr::rename(MS_B = theta)

sart_parcels <- sart_parcel_indicators %>%
  dplyr::select(SERIAL, SART_A, SART_B) %>%
  dplyr::left_join(ms_A, by = "SERIAL") %>%
  dplyr::left_join(ms_B, by = "SERIAL")

sart_reliability_sb <- spearman_brown(sart_parcels$SART_A, sart_parcels$SART_B)
ms_reliability_sb   <- spearman_brown(sart_parcels$MS_A,   sart_parcels$MS_B)

# =============================================================================
# 2.6 MoCA AND MASTER SCORE TABLE
# =============================================================================
# MoCA total scores (patients only) are converted to percent of maximum; the
# impairment criterion is < 26/30 (Nasreddine et al., 2005). All task scores
# are merged into one table (master_theta), which is the basis of all
# validity analyses.

moca_items <- grep("^MC01_\\d{2}$", names(df), value = TRUE)
moca_df <- df %>%
  dplyr::select(SERIAL, dplyr::all_of(moca_items)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(moca_items), ~ suppressWarnings(as.numeric(.x)))) %>%
  dplyr::mutate(
    valid      = rowSums(!is.na(dplyr::across(dplyr::all_of(moca_items)))) > 0,
    MoCA_total = rowSums(dplyr::across(dplyr::all_of(moca_items)), na.rm = TRUE),
    MoCA       = dplyr::if_else(valid, pmin(pmax(MoCA_total, 0), 30) / 30 * 100, NA_real_)
  ) %>%
  dplyr::filter(is.finite(MoCA), MoCA > 0, MoCA <= 100) %>%
  dplyr::select(SERIAL, MoCA)

master_theta <- df %>%
  dplyr::select(SERIAL) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Group = assign_group(SERIAL)) %>%
  dplyr::left_join(dplyr::select(grt_scores, SERIAL, theta_GRT = theta), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(drt_scores, SERIAL, theta_DRT_Fam, theta_DRT_Recoll), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(pal_scores, SERIAL, theta_PAL = theta), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(ci_theta, SERIAL, theta_CI_resilience = CI_Resilience), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(sart_scores, SERIAL, theta_SART_EI = SART_ExternalInterference_Theta, theta_SART_MS = SART_MentalSpeed_Theta), by = "SERIAL") %>%
  dplyr::left_join(moca_df, by = "SERIAL")

readr::write_csv(master_theta, file.path(output_dir, "Master_Theta_Scores.csv"))

moca_impairment_cutoff <- 26

# =============================================================================
# 3. ANALYSED SAMPLE (Results 3.1; Supplementary Figure 1; Suppl. Tables 1-3)
# =============================================================================
# 3.1 Analysed sample, item allocation, days since admission, demographics.
# All participants in the curated data set are analysed (at least one task
# response). Days from hospital admission to testing and the harmonised
# demographics are taken from the curated data set.
analysed_ids <- df$SERIAL
cat(sprintf("[Sample] Analysed sample: %d patients, %d controls.\n",
            sum(assign_group(analysed_ids) == "Patient"), sum(assign_group(analysed_ids) == "Control")))

compact_ranges <- function(x) {
  x <- sort(unique(as.integer(x))); if (!length(x)) return("")
  br <- c(0, which(diff(x) != 1), length(x))
  paste(vapply(seq_len(length(br) - 1), function(i) { a <- x[br[i] + 1]; b <- x[br[i + 1]]
  if (a == b) as.character(a) else paste0(a, "\u2013", b) }, character(1)), collapse = ", ")
}
alloc_patients <- df$SERIAL %in% analysed_ids & assign_group(df$SERIAL) == "Patient"
answered <- function(cols) {
  if (!length(cols)) return(matrix(FALSE, sum(alloc_patients), 0))
  m <- vapply(df[alloc_patients, cols, drop = FALSE], function(x) !is.na(x) & trimws(as.character(x)) != "", logical(sum(alloc_patients)))
  if (is.null(dim(m))) m <- matrix(m, nrow = sum(alloc_patients)); colnames(m) <- cols; m
}
alloc_spec <- list(
  GRT = list(cols = grep("^GRT\\d+$", names(df), value = TRUE), unit = function(c) sub("^GRT0*", "", c)),
  DRT = list(cols = grep("^DRT\\d+$", names(df), value = TRUE), unit = function(c) sub("^DRT", "", c)),
  PAL = list(cols = grep("^PALV\\d+L\\d+$", names(df), value = TRUE), unit = function(c) sub("^PALV(\\d+)L.*$", "\\1", c)),
  CI  = list(cols = grep("^CI_\\d+_", names(df), value = TRUE), unit = function(c) sub("^CI_0*(\\d+)_.*$", "\\1", c))
)
item_allocation_table <- purrr::map_dfr(names(alloc_spec), function(task) {
  sp <- alloc_spec[[task]]; m <- answered(sp$cols)
  if (!ncol(m)) return(tibble::tibble())
  units <- sp$unit(colnames(m))
  per_patient <- apply(m, 1, function(r) compact_ranges(unique(units[r])))
  n_items <- rowSums(m)
  keep <- n_items > 0
  tibble::tibble(Task = task, unit_set = per_patient[keep], n_items = n_items[keep]) %>%
    dplyr::group_by(Task, unit_set) %>%
    dplyr::summarise(n_patients = dplyr::n(), items_answered_median = stats::median(n_items), items_answered_range = paste0(min(n_items), "\u2013", max(n_items)), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_patients))
})
item_allocation_table <- dplyr::mutate(item_allocation_table,
                                       unit = dplyr::recode(Task, GRT = "items", DRT = "items", PAL = "versions", CI = "testlets"))
readr::write_csv(item_allocation_table, file.path(output_dir, "Table1_ItemAllocation.csv"))
cat("[Allocation] Item sets answered by patients (most frequent first; PAL = versions, CI = testlets):\n")
print(item_allocation_table %>% dplyr::group_by(Task) %>% dplyr::slice_head(n = 8) %>% dplyr::ungroup(), n = Inf, width = Inf)

df_days_admission <- if ("days_admission_to_test" %in% names(df)) {
  df %>% dplyr::filter(Group == "Patient") %>% dplyr::select(SERIAL, days_admission_to_test)
} else NULL
df_days <- if (!is.null(df_days_admission)) {
  df_days_admission %>% dplyr::transmute(SERIAL, days_post_stroke = days_admission_to_test, over_15_days = days_admission_to_test > 15)
} else NULL

# Demographic variables of the analysed sample: raw value, recoded label and
# group per participant.
harmonise_demo <- function(var, recode_fun) {
  df %>%
    dplyr::filter(SERIAL %in% analysed_ids) %>%
    dplyr::transmute(SERIAL, raw = as.character(.data[[var]]), label = recode_fun(raw), Group = assign_group(SERIAL))
}

df_edu <- harmonise_demo("school", function(school) dplyr::case_when(
  school %in% c("1", "No degree")                                ~ "No degree",
  school %in% c("4", "Abitur", "University entrance (Abitur)")   ~ "Abitur",
  school %in% c("2", "Realschule", "Intermediate (Realschule)")  ~ "Realschule",
  school %in% c("3", "Hauptschule", "Secondary (Hauptschule)")   ~ "Hauptschule",
  TRUE ~ NA_character_
)) %>% dplyr::mutate(D_Abitur = dplyr::case_when(label == "Abitur" ~ 1L, label %in% c("No degree", "Hauptschule", "Realschule") ~ 0L, TRUE ~ NA_integer_))

df_sex <- harmonise_demo("sex", function(sex) {
  s <- tolower(trimws(sex))
  dplyr::case_when(
    s %in% c("m", "male", "mannlich", "m\u00e4nnlich")         ~ "Male",
    s %in% c("f", "w", "female", "weiblich")                   ~ "Female",
    s %in% c("d", "diverse", "divers", "other", "andere", "3") ~ "Diverse",
    s == "1" ~ "Female",
    s == "2" ~ "Male",
    TRUE ~ NA_character_
  )
}) %>% dplyr::mutate(D_Male = dplyr::case_when(label == "Male" ~ 1L, label == "Female" ~ 0L, TRUE ~ NA_integer_))

age_col <- "age"
if (!(age_col %in% names(df))) {
  stop(sprintf("Expected an age column named '%s' in the data set but did not find one. Check names(df) and update age_col.", age_col))
}

df_age <- harmonise_demo(age_col, function(a) a) %>%
  dplyr::transmute(SERIAL, age = suppressWarnings(as.numeric(sub(",", ".", label))))

cat(sprintf("[Demographics] Using '%s' as the age column (%d non-missing values, range %s-%s).\n",
            age_col, sum(!is.na(df_age$age)),
            suppressWarnings(min(df_age$age, na.rm = TRUE)), suppressWarnings(max(df_age$age, na.rm = TRUE))))

# --- 3.2 Clinical characteristics (patients) ----------------------------------

if (!is.null(stroke_aetiology)) {
  aetiology_joined <- master_theta %>%
    dplyr::select(SERIAL, Group) %>%
    dplyr::inner_join(stroke_aetiology, by = "SERIAL")
  
  ae_patients <- dplyr::filter(aetiology_joined, Group == "Patient", SERIAL %in% analysed_ids)
  n_total_patients_with_theta <- sum(master_theta$Group == "Patient", na.rm = TRUE)
  
  fmt_mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
  fmt_median_iqr <- function(x) sprintf("%.1f [%.1f, %.1f]", stats::median(x, na.rm = TRUE),
                                        stats::quantile(x, 0.25, na.rm = TRUE), stats::quantile(x, 0.75, na.rm = TRUE))
  fmt_n_pct <- function(n, denom) sprintf("%d (%.1f%%)", n, 100 * n / denom)
  
  n_nihss <- sum(!is.na(ae_patients$NIHSS_score))
  hemis_tab <- table(ae_patients$Hemisphere_clean, useNA = "ifany")
  loc_n <- nrow(ae_patients)
  n_haem <- sum(ae_patients$haemorrhagic, na.rm = TRUE)
  n_isch <- sum(!ae_patients$haemorrhagic, na.rm = TRUE)
  
  days_src_admission <- !is.null(df_days_admission) && any(!is.na(df_days_admission$days_admission_to_test))
  days_label_mean <- if (days_src_admission) "Days from admission to testing, mean (SD)" else "Days since stroke onset, mean (SD)"
  days_label_mdn  <- if (days_src_admission) "Days from admission to testing, median [IQR]" else "Days since stroke onset, median [IQR]"
  days_label_15   <- if (days_src_admission) "Days from admission to testing: >15 days" else "Days since stroke onset: >15 days"
  days_rows <- if (!is.null(df_days) && any(!is.na(df_days$days_post_stroke))) {
    dd <- df_days$days_post_stroke; n_dd <- sum(!is.na(dd))
    tibble::tibble(Variable = c(days_label_mean, days_label_mdn, days_label_15),
                   Statistic = c(fmt_mean_sd(dd), fmt_median_iqr(dd), fmt_n_pct(sum(df_days$over_15_days, na.rm = TRUE), n_dd)),
                   N = c(n_dd, n_dd, sum(df_days$over_15_days, na.rm = TRUE)))
  } else NULL
  stroke_characteristics <- dplyr::bind_rows(
    days_rows,
    tibble::tibble(Variable = "NIHSS score, median [IQR]", Statistic = if (n_nihss > 0) fmt_median_iqr(ae_patients$NIHSS_score) else "NA", N = n_nihss),
    tibble::tibble(Variable = "Stroke type: Haemorrhagic", Statistic = fmt_n_pct(n_haem, loc_n), N = n_haem),
    tibble::tibble(Variable = "Stroke type: Ischemic", Statistic = fmt_n_pct(n_isch, loc_n), N = n_isch),
    tibble::tibble(Variable = paste0("Hemisphere: ", names(hemis_tab)), Statistic = fmt_n_pct(as.integer(hemis_tab), loc_n), N = as.integer(hemis_tab)),
    tibble::tibble(Variable = "Lesion location: Anterior*", Statistic = fmt_n_pct(sum(ae_patients$Loc_Anterior, na.rm = TRUE), loc_n), N = sum(ae_patients$Loc_Anterior, na.rm = TRUE)),
    tibble::tibble(Variable = "Lesion location: Posterior*", Statistic = fmt_n_pct(sum(ae_patients$Loc_Posterior, na.rm = TRUE), loc_n), N = sum(ae_patients$Loc_Posterior, na.rm = TRUE)),
    tibble::tibble(Variable = "Lesion location: Media*", Statistic = fmt_n_pct(sum(ae_patients$Loc_Media, na.rm = TRUE), loc_n), N = sum(ae_patients$Loc_Media, na.rm = TRUE))
  ) %>% dplyr::mutate(
    footnote = "* Lesion location categories are not mutually exclusive (a patient may have more than one); percentages need not sum to 100%.",
    pct_of_analytic_patient_sample = round(100 * loc_n / n_total_patients_with_theta, 1)
  )
  readr::write_csv(stroke_characteristics, file.path(output_dir, "SupplTable01_ClinicalCharacteristics.csv"))
  cat(sprintf("[Stroke aetiology] Clinical characteristics computed for %d of %d analysed patients (%.1f%%):\n",
              loc_n, n_total_patients_with_theta, round(100 * loc_n / n_total_patients_with_theta, 1)))
  print(stroke_characteristics)
} else {
  cat("[Stroke aetiology] No aetiology variables in the curated data set -- clinical characteristics skipped.\n")
}

# --- 3.3 Participant characteristics (Supplementary Table 1) ----------------
# Continuous variables: Welch's t-test with Cohen's d. Categorical variables:
# chi-square test, or Fisher's exact test (simulated p) when any observed cell
# count is below 5. p values are Holm-corrected across all tested variables.

df_handedness <- harmonise_demo("handedness", function(h) dplyr::case_when(
  tolower(h) %in% c("1", "right", "rechts")             ~ "Right",
  tolower(h) %in% c("2", "left", "links")                ~ "Left",
  tolower(h) %in% c("3", "ambidextrous", "beidhandig", "beidh\u00e4ndig") ~ "Ambidextrous",
  TRUE ~ NA_character_
))

parse_decimal <- function(x) suppressWarnings(as.numeric(stringr::str_replace(as.character(x), ",", ".")))
df_weight <- harmonise_demo("weight", function(w) w) %>%
  dplyr::transmute(SERIAL, weight_raw = label, weight_kg = parse_decimal(label))
df_height <- harmonise_demo("height", function(h) h) %>%
  dplyr::transmute(SERIAL, height_raw = label, height_m = parse_decimal(label)) %>%
  dplyr::mutate(height_m = ifelse(!is.na(height_m) & height_m < 3, height_m, height_m / 100))

df_bmi <- df_height %>%
  dplyr::inner_join(df_weight, by = "SERIAL") %>%
  dplyr::mutate(Group = assign_group(SERIAL),
                BMI   = round(weight_kg / (height_m^2), 1))

bmi_bad <- !is.na(df_bmi$BMI) & (df_bmi$BMI < 14 | df_bmi$BMI > 70)
if (any(bmi_bad)) {
  cat(sprintf("[Demographics] %d BMI value(s) outside 14-70 set to NA. Raw entries:\n", sum(bmi_bad)))
  print(df_bmi[bmi_bad, c("SERIAL", "Group", "height_raw", "weight_raw", "height_m", "weight_kg", "BMI")])
  df_bmi$BMI[bmi_bad] <- NA_real_
}
bmi_range <- range(df_bmi$BMI, na.rm = TRUE)
cat(sprintf("[Demographics] BMI computed (height unit detected per row), range %.1f-%.1f across %d valid values.\n",
            bmi_range[1], bmi_range[2], sum(!is.na(df_bmi$BMI))))

comorbidity_cols <- setdiff(grep("^CM_", names(df), value = TRUE), "CM_None")
cat(sprintf("[Demographics] %d comorbidity column(s) found: %s\n",
            length(comorbidity_cols), if (length(comorbidity_cols) > 0) paste(comorbidity_cols, collapse = ", ") else "none"))

compare_demo_continuous <- function(x, group, label, group_levels = c("Patient", "Control")) {
  d <- tibble::tibble(x = x, Group = factor(group, levels = group_levels)) %>% dplyr::filter(is.finite(x), !is.na(Group))
  if (dplyr::n_distinct(d$Group) < 2 || nrow(d) < 4) {
    return(tibble::tibble(Variable = label, Patient = "NA", Control = "NA", Statistic = NA_character_, p_raw = NA_real_, Effect_size = NA_character_))
  }
  means <- d %>% dplyr::group_by(Group) %>% dplyr::summarise(M = mean(x), SD = stats::sd(x), .groups = "drop")
  tt  <- tryCatch(stats::t.test(x ~ Group, data = d), error = function(e) NULL)
  eff <- tryCatch(effsize::cohen.d(x ~ Group, data = d), error = function(e) NULL)
  tibble::tibble(
    Variable = label,
    Patient = sprintf("%.1f \u00b1 %.1f", means$M[means$Group == "Patient"], means$SD[means$Group == "Patient"]),
    Control = sprintf("%.1f \u00b1 %.1f", means$M[means$Group == "Control"], means$SD[means$Group == "Control"]),
    Statistic = if (!is.null(tt)) sprintf("t(%.0f)=%.2f", unname(tt$parameter), unname(tt$statistic)) else NA_character_,
    p_raw = if (!is.null(tt)) tt$p.value else NA_real_,
    Effect_size = if (!is.null(eff)) sprintf("d=%.2f", unname(eff$estimate)) else NA_character_
  )
}

compare_demo_categorical <- function(x, group, label, group_levels = c("Patient", "Control"), test_rule = c("fisher", "chisq_or_fisher")) {
  test_rule <- match.arg(test_rule)
  d <- tibble::tibble(x = as.character(x), Group = factor(group, levels = group_levels)) %>% dplyr::filter(!is.na(x), !is.na(Group))
  empty <- tibble::tibble(Variable = character(0), Patient = character(0), Control = character(0), Statistic = character(0), p_raw = double(0), Effect_size = character(0))
  if (dplyr::n_distinct(d$Group) < 2 || dplyr::n_distinct(d$x) < 2) {
    return(dplyr::bind_rows(tibble::tibble(Variable = label, Patient = NA_character_, Control = NA_character_, Statistic = NA_character_, p_raw = NA_real_, Effect_size = NA_character_), empty))
  }
  tab <- table(d$Group, d$x)
  test <- if (test_rule == "chisq_or_fisher") {
    if (any(tab < 5)) tryCatch({ set.seed(20260821); stats::fisher.test(tab, simulate.p.value = TRUE) }, error = function(e) NULL)
    else tryCatch(stats::chisq.test(tab), error = function(e) NULL)
  } else {
    tryCatch(stats::fisher.test(tab), error = function(e) tryCatch(stats::chisq.test(tab), error = function(e2) NULL))
  }
  stat_label <- NA_character_
  if (!is.null(test)) stat_label <- if (grepl("Fisher", test$method)) "Fisher" else sprintf("\u03c7\u00b2(%d)=%.2f", unname(test$parameter), unname(test$statistic))
  n_group <- rowSums(tab)
  pct_tab <- sweep(tab, 1, n_group, "/") * 100
  header <- tibble::tibble(Variable = label, Patient = NA_character_, Control = NA_character_,
                           Statistic = stat_label, p_raw = if (!is.null(test)) test$p.value else NA_real_, Effect_size = NA_character_)
  categories <- tibble::tibble(
    Variable = paste0("  ", colnames(tab)),
    Patient  = sprintf("%d (%.1f%%)", tab[group_levels[1], ], pct_tab[group_levels[1], ]),
    Control  = sprintf("%d (%.1f%%)", tab[group_levels[2], ], pct_tab[group_levels[2], ]),
    Statistic = NA_character_, p_raw = NA_real_, Effect_size = NA_character_
  )
  dplyr::bind_rows(header, categories)
}

clean_numeric_string <- function(x, max_val = Inf) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "(\\d+)\\s*-\\s*(\\d+)", "\\2")
  x_num <- suppressWarnings(as.numeric(stringr::str_replace(stringr::str_extract(x, "\\d+([\\.,]\\d+)?"), ",", ".")))
  x_final <- floor(x_num + 0.5)
  ifelse(x_final > max_val, NA_real_, x_final)
}

demo_age    <- compare_demo_continuous(df_age$age, assign_group(df_age$SERIAL), "Age (years)")
demo_sex    <- compare_demo_categorical(df_sex$label, df_sex$Group, "Sex", test_rule = "chisq_or_fisher")
demo_hand   <- compare_demo_categorical(df_handedness$label, df_handedness$Group, "Handedness", test_rule = "chisq_or_fisher")
demo_bmi    <- compare_demo_continuous(df_bmi$BMI, assign_group(df_bmi$SERIAL), "BMI (kg/m\u00b2)")
demo_school <- compare_demo_categorical(df_edu$label, df_edu$Group, "School qualification", test_rule = "chisq_or_fisher")

extra_demo_spec <- list(
  list(col = "education_in_years", label = "Education (years)",           clean = function(x) clean_numeric_string(x, max_val = 30)),
  list(col = "occupation",         label = "Occupation",                  clean = NULL),
  list(col = "exercise",           label = "Exercise",                    clean = NULL),
  list(col = "sleep_avg",          label = "Sleep, average (hours)",      clean = function(x) clean_numeric_string(x, max_val = 24)),
  list(col = "sleep_last_night",   label = "Sleep, last night (hours)",   clean = function(x) clean_numeric_string(x, max_val = 24))
)
extra_demo_table <- purrr::map_dfr(extra_demo_spec, function(sp) {
  if (!sp$col %in% names(df)) {
    cat(sprintf("[Demographics] Column '%s' not found -- '%s' omitted from Supplementary Table 1.\n", sp$col, sp$label))
    return(tibble::tibble())
  }
  in_sample <- df$SERIAL %in% analysed_ids
  x <- df[[sp$col]][in_sample]; grp <- assign_group(df$SERIAL[in_sample])
  if (!is.null(sp$clean)) x <- sp$clean(x)
  if (is.numeric(x)) {
    if (is.null(sp$clean) && dplyr::n_distinct(stats::na.omit(x)) <= 10)
      cat(sprintf("[Demographics] NOTE: '%s' is numeric with <= 10 distinct values -- treated as continuous; check whether these are category codes.\n", sp$label))
    compare_demo_continuous(x, grp, sp$label)
  } else {
    compare_demo_categorical(x, grp, sp$label, test_rule = "chisq_or_fisher")
  }
})

comorbidity_table <- tibble::tibble()
if (length(comorbidity_cols) > 0) {
  comorbidity_data <- df %>% dplyr::filter(SERIAL %in% analysed_ids) %>% dplyr::select(SERIAL, dplyr::all_of(comorbidity_cols)) %>%
    dplyr::mutate(Group = assign_group(SERIAL))
  comorbidity_table <- purrr::map_dfr(comorbidity_cols, function(col) {
    res <- compare_demo_categorical(as.character(comorbidity_data[[col]]), comorbidity_data$Group,
                                    paste0("Comorbidity: ", sub("^CM_", "", col)), test_rule = "chisq_or_fisher")
    hdr <- dplyr::filter(res, !is.na(Statistic))
    present <- res[trimws(res$Variable) %in% c("1", "TRUE", "Yes", "yes", "ja", "Ja"), ]
    if (nrow(hdr) == 1 && nrow(present) == 1) { hdr$Patient <- present$Patient; hdr$Control <- present$Control }
    hdr
  })
}

tested_table <- dplyr::bind_rows(demo_age, demo_sex, demo_hand, demo_bmi, demo_school, extra_demo_table, comorbidity_table)
test_rows <- !is.na(tested_table$p_raw)
tested_table$p_holm <- NA_real_
tested_table$p_holm[test_rows] <- stats::p.adjust(tested_table$p_raw[test_rows], method = "holm")

stroke_count_table <- tibble::tibble()
if ("BQ30" %in% names(df)) {
  sc_levels <- c("1 (present only)", "2", "3", "4", "5 or more")
  sc <- df %>% dplyr::filter(assign_group(SERIAL) == "Patient", SERIAL %in% analysed_ids) %>%
    dplyr::transmute(stroke_count = dplyr::case_when(
      as.character(BQ30) == "1" ~ sc_levels[1], as.character(BQ30) == "2" ~ sc_levels[2],
      as.character(BQ30) == "3" ~ sc_levels[3], as.character(BQ30) == "4" ~ sc_levels[4],
      as.character(BQ30) == "5" ~ sc_levels[5], TRUE ~ NA_character_)) %>%
    dplyr::filter(!is.na(stroke_count))
  if (nrow(sc) > 0) {
    sc_tab <- table(factor(sc$stroke_count, levels = sc_levels))
    stroke_count_table <- dplyr::bind_rows(
      tibble::tibble(Variable = "Number of strokes", Patient = NA_character_),
      tibble::tibble(Variable = paste0("  ", names(sc_tab)), Patient = sprintf("%d (%.1f%%)", as.integer(sc_tab), 100 * as.integer(sc_tab) / sum(sc_tab)))
    ) %>% dplyr::mutate(Control = "N/A", Statistic = NA_character_, p_raw = NA_real_, Effect_size = NA_character_, p_holm = NA_real_)
  }
} else {
  cat("[Demographics] Column 'BQ30' not found -- number of strokes omitted from Supplementary Table 1.\n")
}

participant_characteristics <- dplyr::bind_rows(
  tested_table,
  stroke_count_table,
  if (exists("stroke_characteristics")) stroke_characteristics %>% dplyr::transmute(Variable, Patient = Statistic, Control = "N/A", Statistic = NA_character_, p_raw = NA_real_, Effect_size = NA_character_, p_holm = NA_real_) else NULL
)
readr::write_csv(participant_characteristics, file.path(output_dir, "SupplTable01_ParticipantCharacteristics.csv"))
cat(sprintf("[Demographics] Supplementary Table 1 written (%d tested variables, Holm-corrected together):\n", sum(test_rows)))
print(participant_characteristics, n = Inf)

# --- 3.4 Score summary tables (Suppl. Tables 2, 4, 5, 9) ---------------------
# Analysed sample sizes per task and group, item fit before and after
# trimming, IRT model comparisons, known-groups effect sizes (Holm-corrected
# across the seven scores), and split-half reliabilities with percentile
# bootstrap confidence intervals (2,000 resamples of participants).

sample_size_table <- dplyr::bind_rows(
  grt_scores %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "GRT", Score = "theta"),
  drt_scores %>% dplyr::filter(!is.na(theta_DRT_Fam)) %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "DRT", Score = "Fam"),
  drt_scores %>% dplyr::filter(!is.na(theta_DRT_Recoll)) %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "DRT", Score = "Recoll"),
  pal_scores %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "PAL", Score = "theta"),
  ci_theta   %>% dplyr::filter(is.finite(CI_Resilience)) %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "CI",  Score = "internal_interference"),
  sart_scores %>% dplyr::filter(!is.na(SART_ExternalInterference_Theta)) %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "SART", Score = "external_interference"),
  sart_scores %>% dplyr::filter(!is.na(SART_MentalSpeed_Theta)) %>% dplyr::count(Group, name = "N") %>% dplyr::mutate(Task = "SART", Score = "mental_speed")
) %>% dplyr::select(Task, Score, Group, N) %>% dplyr::arrange(Task, Score, Group)
readr::write_csv(sample_size_table, file.path(output_dir, "SupplTable02_SampleSizes.csv"))
cat("\n[Sample] Analysed sample sizes by task and group:\n"); print(sample_size_table)

model_comparison_table <- tibble::tibble(
  Task       = c("GRT", "DRT", "DRT", "CI", "PAL"),
  Dimension  = c(NA, "Fam", "Recoll", NA, NA),
  Comparison = c("1PL vs 2PL", "1PL vs 2PL", "1PL vs 2PL", "1D vs 2D (joint)", "PCM vs GRM"),
  delta_AIC  = c(grt_delta_aic, drt_fit1$daic, drt_fit2$daic, ci_aic_1d - ci_aic_2d, pal_delta_aic),
  delta_BIC  = c(grt_delta_bic, drt_fit1$dbic, drt_fit2$dbic, ci_bic_1d - ci_bic_2d, pal_delta_bic)
)
readr::write_csv(model_comparison_table, file.path(output_dir, "SupplTable05_ModelComparisons.csv"))

item_fit_full <- dplyr::bind_rows(item_fit_report)
readr::write_csv(item_fit_full, file.path(output_dir, "SupplTable04_ItemFit.csv"))

effect_size_table <- dplyr::bind_rows(
  with(grt_group_comp, tibble::tibble(Task = "GRT", Score = "theta", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(drt_gc_Fam, tibble::tibble(Task = "DRT", Score = "Fam", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(drt_gc_Recoll, tibble::tibble(Task = "DRT", Score = "Recoll", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(pal_group_comp, tibble::tibble(Task = "PAL", Score = "theta", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(ci_group_comp, tibble::tibble(Task = "CI", Score = "internal_interference", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(sart_gc_ei, tibble::tibble(Task = "SART", Score = "external_interference", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(sart_gc_ms, tibble::tibble(Task = "SART", Score = "mental_speed", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2]))
) %>% dplyr::mutate(
  magnitude = dplyr::case_when(abs(g) < 0.2 ~ "negligible", abs(g) < 0.5 ~ "small", abs(g) < 0.8 ~ "medium", TRUE ~ "large"),
  p_holm = stats::p.adjust(p, method = "holm")
)

es_desc_cols <- c(`GRT|theta` = "theta_GRT", `DRT|Fam` = "theta_DRT_Fam", `DRT|Recoll` = "theta_DRT_Recoll", `PAL|theta` = "theta_PAL",
                  `CI|internal_interference` = "theta_CI_resilience", `SART|external_interference` = "theta_SART_EI", `SART|mental_speed` = "theta_SART_MS")
es_desc <- purrr::map_dfr(names(es_desc_cols), function(k) {
  z <- as.numeric(scale(master_theta[[es_desc_cols[[k]]]]))
  gc <- compare_groups(z, master_theta$Group)
  tibble::tibble(Task = sub("\\|.*", "", k), Score = sub(".*\\|", "", k),
                 z_M_ctrl = gc$desc$M[gc$desc$Group == "Control"], z_SD_ctrl = gc$desc$SD[gc$desc$Group == "Control"],
                 z_M_pat = gc$desc$M[gc$desc$Group == "Patient"], z_SD_pat = gc$desc$SD[gc$desc$Group == "Patient"])
})
effect_size_table <- dplyr::left_join(effect_size_table, es_desc, by = c("Task", "Score"))
readr::write_csv(effect_size_table, file.path(output_dir, "SupplTable09_KnownGroups.csv"))

bootstrap_sb_ci <- function(x, y, n_boot = 2000) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  boot_r <- replicate(n_boot, {
    idx <- sample.int(n, n, replace = TRUE)
    spearman_brown(x[idx], y[idx])
  })
  stats::quantile(boot_r, c(.025, .975), na.rm = TRUE)
}
reliability_table <- tibble::tibble(
  Task = c("GRT", "DRT (Fam parcel)", "DRT (Recoll parcel)", "PAL", "CI", "SART (external interference)", "SART (mental speed)"),
  reliability_sb = c(grt_reliability_sb, drt_reliability_sb, drt_Recoll_reliability_sb, pal_reliability_sb, ci_reliability_sb, sart_reliability_sb, ms_reliability_sb)
)
reliability_table$boot_ci_lo <- NA_real_
reliability_table$boot_ci_hi <- NA_real_
boot_pairs <- list(
  list(grt_parcels$GRT_A, grt_parcels$GRT_B),
  list(drt_parcels$DRT_Fam_A, drt_parcels$DRT_Fam_B),
  list(drt_Recoll_parcels$DRT_Recoll_A, drt_Recoll_parcels$DRT_Recoll_B),
  list(pal_parcels$PAL_A, pal_parcels$PAL_B),
  list(ci_parcels$CI_A, ci_parcels$CI_B),
  list(sart_parcels$SART_A, sart_parcels$SART_B),
  list(sart_parcels$MS_A, sart_parcels$MS_B)
)
for (i in seq_along(boot_pairs)) {
  ci_vals <- bootstrap_sb_ci(boot_pairs[[i]][[1]], boot_pairs[[i]][[2]])
  reliability_table$boot_ci_lo[i] <- ci_vals[1]
  reliability_table$boot_ci_hi[i] <- ci_vals[2]
}
readr::write_csv(reliability_table, file.path(output_dir, "SupplTable02_Reliability.csv"))

# --- 3.5 Attrition: completers vs. non-completers (Supplementary Table 3) ---
# Completers finished all five tasks (the CI is administered last). Because
# missing data arise both by design and by early discontinuation, completers
# and non-completers among the analysed patients are compared on MoCA, age,
# sex, education and performance on the first task (GRT). Differences on
# observed variables are compatible with the missing-at-random assumption of
# the estimation methods.

all_patients <- intersect(unique(df$SERIAL[assign_group(df$SERIAL) == "Patient"]), analysed_ids)

completer_status <- tibble::tibble(SERIAL = all_patients) %>%
  dplyr::mutate(completer = SERIAL %in% ci_theta$SERIAL[is.finite(ci_theta$CI_Resilience)])

cat(sprintf("[Completers] %d of %d patients (%.1f%%) completed through CI; %d (%.1f%%) did not.\n",
            sum(completer_status$completer), nrow(completer_status), 100 * mean(completer_status$completer),
            sum(!completer_status$completer), 100 * (1 - mean(completer_status$completer))))

completer_compare_data <- completer_status %>%
  dplyr::left_join(moca_df, by = "SERIAL") %>%
  dplyr::left_join(df_age, by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(df_sex, SERIAL, D_Male), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(df_edu, SERIAL, D_Abitur), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(grt_scores, SERIAL, theta_GRT = theta), by = "SERIAL")

compare_completers <- function(x, completer, label, type = c("continuous", "binary")) {
  type <- match.arg(type)
  d <- tibble::tibble(x = x, completer = completer) %>% dplyr::filter(!is.na(x), !is.na(completer))
  if (dplyr::n_distinct(d$completer) < 2 || nrow(d) < 10) {
    return(tibble::tibble(variable = label, N_completer = sum(d$completer, na.rm = TRUE), N_noncompleter = sum(!d$completer, na.rm = TRUE), status = "insufficient data"))
  }
  if (type == "binary") {
    test <- tryCatch(stats::chisq.test(table(d$completer, d$x), correct = FALSE), error = function(e) NULL)
    tibble::tibble(variable = label, N_completer = sum(d$completer), N_noncompleter = sum(!d$completer),
                   M_completer = mean(d$x[d$completer]) * 100, M_noncompleter = mean(d$x[!d$completer]) * 100,
                   statistic = if (!is.null(test)) unname(test$statistic) else NA_real_,
                   p = if (!is.null(test)) test$p.value else NA_real_,
                   effect = NA_real_, status = "ok (% shown; chi-square test)")
  } else {
    tt  <- tryCatch(t.test(x ~ completer, data = d), error = function(e) NULL)
    eff <- tryCatch(effsize::cohen.d(x ~ completer, data = d, hedges.correction = TRUE), error = function(e) NULL)
    tibble::tibble(variable = label, N_completer = sum(d$completer), N_noncompleter = sum(!d$completer),
                   M_completer = mean(d$x[d$completer]), M_noncompleter = mean(d$x[!d$completer]),
                   statistic = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
                   p = if (!is.null(tt)) tt$p.value else NA_real_,
                   effect = if (!is.null(eff)) unname(eff$estimate) else NA_real_,
                   status = "ok (means shown; t-test + Hedges g)")
  }
}

completer_results <- dplyr::bind_rows(
  compare_completers(completer_compare_data$MoCA / 100 * 30, completer_compare_data$completer, "MoCA (score, 0-30)", "continuous"),
  compare_completers(completer_compare_data$age, completer_compare_data$completer, "Age", "continuous"),
  compare_completers(completer_compare_data$D_Male, completer_compare_data$completer, "Sex (% Male)", "binary"),
  compare_completers(completer_compare_data$D_Abitur, completer_compare_data$completer, "Education (% Abitur)", "binary"),
  compare_completers(completer_compare_data$theta_GRT, completer_compare_data$completer, "GRT theta (first-administered task)", "continuous")
)
readr::write_csv(completer_results, file.path(output_dir, "SupplTable03_Attrition.csv"))
cat("[Completers] Completer vs. non-completer comparison on baseline variables (analysed patients only; t and g are non-completers minus completers):\n")
print(completer_results)

# --- 3.6 Participant flow (Supplementary Figure 1) ---------------------------
# Screening and enrolment numbers come from the recruitment log; exclusions
# before the curated data set are reported by
# EMA4Stroke_ORPheoS_Preprocessing.R (Section 7). All other counts are taken
# from the data of this run.

fc_n_stroke_screened <- 589
fc_n_stroke_enrolled <- 289
fc_stroke_dates      <- "Aug 2024 - Jun 2026"
fc_ctrl_p1_dates     <- "Feb&#8211;Nov 2023"
fc_ctrl_p2_dates     <- "Mar&#8211;Apr 2026"
fc_ctrl_p1_tasks     <- "GRT, DRT (SART-ED without trial-level data, not analysed)"
fc_ctrl_p2_tasks     <- "All tasks"
fc_excl_reasons      <- c("Physical/cognitive inability to comply", "Absent from room during assessment", "Declined to participate")

fc_n_not_started <- 21L   # enrolled, digital assessment not started
fc_n_nodata      <- 30L   # started, but no task responses
fc_review_exclusions <- tibble::tibble(
  reason = c("Hemispatial neglect", "Transient ischaemic attack", "Stroke not confirmed"),
  n      = c(1L, 2L, 1L)
)  # excluded after clinical review
fc_n_review_excluded <- sum(fc_review_exclusions$n)

fc_pkgs <- c("DiagrammeR", "DiagrammeRsvg")
fc_ok <- vapply(fc_pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))

fc_n_stroke_excluded <- fc_n_stroke_screened - fc_n_stroke_enrolled
fc_n_stroke_analysed <- sum(df$Group == "Patient")
fc_n_ctrl_p1 <- sum(df$Phase %in% "Phase1")
fc_n_ctrl_p2 <- sum(df$Phase %in% "Phase2")
if (fc_n_stroke_enrolled - fc_n_not_started - fc_n_nodata - fc_n_review_excluded != fc_n_stroke_analysed)
  cat("[Flowchart] WARNING: enrolled minus exclusions does not equal the analysed patients -- check the fixed values above.\n")
cat(sprintf("[Flowchart] Patients: %d screened -> %d excluded at screening -> %d enrolled -> %d did not start, %d no task data, %d excluded after clinical review -> %d analysed. Controls: %d phase 1, %d phase 2.\n",
            fc_n_stroke_screened, fc_n_stroke_excluded, fc_n_stroke_enrolled, fc_n_not_started, fc_n_nodata, fc_n_review_excluded,
            fc_n_stroke_analysed, fc_n_ctrl_p1, fc_n_ctrl_p2))

fc_get_n <- function(task, score, group) {
  v <- sample_size_table$N[sample_size_table$Task == task & sample_size_table$Score == score & sample_size_table$Group == group]
  if (length(v) == 0) 0L else v
}
fc_tasks <- list(
  list(label = "SART-ED", stroke = max(fc_get_n("SART", "external_interference", "Patient"), fc_get_n("SART", "mental_speed", "Patient")),
       control = max(fc_get_n("SART", "external_interference", "Control"), fc_get_n("SART", "mental_speed", "Control")), note = ""),
  list(label = "CI",  stroke = fc_get_n("CI", "internal_interference", "Patient"), control = fc_get_n("CI", "internal_interference", "Control"), note = ""),
  list(label = "GRT", stroke = fc_get_n("GRT", "theta", "Patient"), control = fc_get_n("GRT", "theta", "Control"), note = ""),
  list(label = "DRT", stroke = max(fc_get_n("DRT", "Fam", "Patient"), fc_get_n("DRT", "Recoll", "Patient")),
       control = max(fc_get_n("DRT", "Fam", "Control"), fc_get_n("DRT", "Recoll", "Control")), note = ""),
  list(label = "PAL", stroke = fc_get_n("PAL", "theta", "Patient"), control = fc_get_n("PAL", "theta", "Control"), note = "")
)

{
  fc_font <- 34
  col_ctrl  <- c(fill = "#6695c5", border = "#3d6694", font = "#ffffff")
  col_str   <- c(fill = "#eca594", border = "#c45a3d", font = "#5a1e0f")
  col_excl  <- c(fill = "#f5e6e2", border = "#c45a3d", font = "#7a2010")
  col_mid   <- c(fill = "#f0ede8", border = "#888680", font = "#2c2c2a")
  col_task  <- c(fill = "#f4f4f2", border = "#888680", font = "#2c2c2a")
  col_edge  <- "#5a5a58"
  
  tr        <- function(content) sprintf('<TR><TD ALIGN="CENTER">%s</TD></TR>', content)
  tr_space  <- function() '<TR><TD> </TD></TR>'
  html_lbl  <- function(rows) sprintf('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">\n    %s\n  </TABLE>>', paste(unlist(rows), collapse = "\n    "))
  row_bold  <- function(txt) tr(sprintf("<B>%s</B>", txt))
  row_plain <- function(txt) tr(txt)
  row_small <- function(txt, pt = fc_font - 4) tr(sprintf('<FONT POINT-SIZE="%d">%s</FONT>', pt, txt))
  ni        <- function(val) sprintf("<I>n</I> = %s", val)
  lbl_top    <- function(rows, total) html_lbl(c(rows, rep(list(tr_space()), total - length(rows))))
  lbl_bottom <- function(rows, total) html_lbl(c(rep(list(tr_space()), total - length(rows)), rows))
  sty <- function(col) sprintf('style="filled,rounded" fillcolor="%s" color="%s" fontcolor="%s"', col["fill"], col["border"], col["font"])
  
  lbl_stroke_screen <- lbl_top(list(row_bold("Patients"), row_plain(sprintf("In-Patient recruitment (%s)", fc_stroke_dates)), row_plain(ni(fc_n_stroke_screened))), 4)
  lbl_ctrl_p1 <- lbl_top(list(row_bold("Controls"), row_plain(sprintf("Recruited in phase 1 (%s)", fc_ctrl_p1_dates)), row_plain(ni(fc_n_ctrl_p1)), row_small(sprintf("Tasks: %s", fc_ctrl_p1_tasks))), 4)
  lbl_ctrl_p2 <- lbl_top(list(row_bold("Controls"), row_plain(sprintf("Recruited in phase 2 (%s)", fc_ctrl_p2_dates)), row_plain(ni(fc_n_ctrl_p2)), row_small(sprintf("Tasks: %s", fc_ctrl_p2_tasks))), 4)
  lbl_excl <- html_lbl(c(list(row_bold(sprintf("Excluded (%s)", ni(fc_n_stroke_excluded)))),
                         lapply(fc_excl_reasons, function(r) row_plain(sprintf("&#8226; %s", r)))))
  
  lbl_stroke_enrol <- html_lbl(list(row_bold("Patients enrolled"), row_plain(ni(fc_n_stroke_enrolled))))
  excl2_rows <- list(row_bold(sprintf("Excluded after enrolment (%s)", ni(fc_n_not_started + fc_n_nodata + fc_n_review_excluded))),
                     row_plain(sprintf("&#8226; No task data (%s):", ni(fc_n_not_started + fc_n_nodata))))
  if (fc_n_not_started > 0) excl2_rows <- c(excl2_rows, list(row_small(sprintf("Did not start the digital assessment (%d)", fc_n_not_started))))
  excl2_rows <- c(excl2_rows, list(row_small(sprintf("Started, but no task responses (%d)", fc_n_nodata))))
  if (fc_n_review_excluded > 0) {
    excl2_rows <- c(excl2_rows, list(row_plain(sprintf("&#8226; After clinical review (%s):", ni(fc_n_review_excluded)))),
                    lapply(seq_len(nrow(fc_review_exclusions)), function(i) row_small(sprintf("%s (%d)", fc_review_exclusions$reason[i], fc_review_exclusions$n[i]))))
  }
  lbl_excl2 <- html_lbl(excl2_rows)
  lbl_stroke_analysed <- lbl_bottom(list(row_bold("Patients analysed"), row_plain(ni(fc_n_stroke_analysed))), 2)
  lbl_ctrl_enrol      <- lbl_bottom(list(row_bold("Controls analysed"), row_plain(ni(fc_n_ctrl_p1 + fc_n_ctrl_p2))), 2)
  lbl_pmiss <- html_lbl(list(
    row_bold("Planned-missingness design"),
    row_plain("Patients: consecutive cohorts (n &#8776; 50) received different item sets of equal size;"),
    row_plain("Controls (phase 1) received only GRT items 1-40 and DRT items 1-18;"),
    row_plain("Controls (phase 2) received all items;"),
    row_plain("some items/sets skipped due to early dropout"),
    row_plain("(fatigue or loss of motivation)")
  ))
  task_nodes <- vapply(seq_along(fc_tasks), function(i) {
    t <- fc_tasks[[i]]
    rows <- list(row_bold(t$label), row_plain(sprintf("Patients: %s", ni(t$stroke))), row_plain(sprintf("Controls: %s", ni(t$control))))
    if (nzchar(t$note)) rows <- c(rows, list(row_small(t$note)))
    sprintf("  T%d [label=%s %s]", i, lbl_bottom(rows, 4), sty(col_task))
  }, character(1))
  
  dot_code <- paste0(
    'digraph flow {
  graph [rankdir=TB nodesep=0.70 ranksep=0.85 fontname="Arial" bgcolor="white" splines=ortho pad="0.6,0.5"]
  node  [shape=box fontname="Arial" fontsize=', fc_font, ' margin="0.50,0.35" width=5.0 penwidth=1.8]
  edge  [color="', col_edge, '" penwidth=1.6 arrowsize=0.9 fontname="Arial" fontsize=', fc_font - 4, ']
  STROKE_SCREEN   [label=', lbl_stroke_screen, ' ', sty(col_str), ']
  CTRL_P1         [label=', lbl_ctrl_p1, ' ', sty(col_ctrl), ']
  CTRL_P2         [label=', lbl_ctrl_p2, ' ', sty(col_ctrl), ']
  EXCL            [label=', lbl_excl, ' ', sty(col_excl), ' width=6.2]
  STROKE_ENROL    [label=', lbl_stroke_enrol, ' ', sty(col_str), ']
  EXCL2           [label=', lbl_excl2, ' ', sty(col_excl), ' width=6.2]
  STROKE_ANALYSED [label=', lbl_stroke_analysed, ' ', sty(col_str), ']
  CTRL_ENROL      [label=', lbl_ctrl_enrol, ' ', sty(col_ctrl), ' width=6.8]
  PMISS           [label=', lbl_pmiss, ' ', sty(col_mid), ' width=13.0]
', paste(task_nodes, collapse = "\n"), '
  { rank=same; STROKE_SCREEN; CTRL_P1; CTRL_P2 }
  { rank=same; STROKE_ENROL; EXCL }
  { rank=same; STROKE_ANALYSED; EXCL2; CTRL_ENROL }
  { rank=same; T1; T2; T3; T4; T5 }
  STROKE_SCREEN -> EXCL [style=dashed arrowhead=open color="', col_excl["border"], '"]
  STROKE_SCREEN -> STROKE_ENROL
  STROKE_ENROL -> EXCL2 [style=dashed arrowhead=open color="', col_excl["border"], '"]
  STROKE_ENROL -> STROKE_ANALYSED
  CTRL_P1 -> CTRL_ENROL
  CTRL_P2 -> CTRL_ENROL
  STROKE_ANALYSED -> PMISS
  CTRL_ENROL -> PMISS
  PMISS -> T1
  PMISS -> T2
  PMISS -> T3
  PMISS -> T4
  PMISS -> T5
}')
  
  writeLines(dot_code, file.path(output_dir, "SupplFigure1_ParticipantFlow.dot"))
  fc_png <- file.path(output_dir, "SupplFigure1_ParticipantFlow.png")
  fc_svg <- file.path(output_dir, "SupplFigure1_ParticipantFlow.svg")
  if (!all(fc_ok)) {
    cat("[Flowchart] DiagrammeR/DiagrammeRsvg not installed -- DOT source written (SupplFigure1_ParticipantFlow.dot); counts above are valid.\n")
  } else tryCatch({
    svg_txt <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot_code))
    writeLines(svg_txt, fc_svg)
    rsvg::rsvg_png(charToRaw(svg_txt), file = fc_png, width = 3800)
    cat(sprintf("[Flowchart] Written: %s (+ .png)\n", fc_svg))
  }, error = function(e) cat(sprintf("[Flowchart] Rendering failed: %s\n", conditionMessage(e))))
}

# =============================================================================
# 4. ITEM-BANK PROPERTIES (Suppl. Tables 5-6; Suppl. Figures 3-4)
# =============================================================================
# 4.1 Practical impact of the Rasch model, targeting, test information.
# Rasch vs. 2PL: ability estimates and group effect sizes are compared to show
# whether assuming equal discrimination has practical consequences.
# Targeting: proportion of persons within the range of item locations and the
# mean person-item offset; conditional standard errors per theta, where
# SE <= 0.55 corresponds to a conditional reliability of about .70.

irt_scales <- list(
  list(label = "GRT",    mod = grt_mod_1pl,  serial = grt_X$SERIAL[grt_rows_any]),
  list(label = "DRT-Fam", mod = drt_fit1$mod, serial = drt_X$SERIAL[rowSums(!is.na(drt_Y_bin_rm[, drt_final_Fam, drop = FALSE])) > 0]),
  list(label = "DRT-Recoll", mod = drt_fit2$mod, serial = drt_X$SERIAL[rowSums(!is.na(drt_Y_bin_rm[, drt_final_Recoll, drop = FALSE])) > 0]),
  list(label = "PAL",    mod = pal_mod_pcm,  serial = pal_serial_kept)
)
for (dimname in c("AC", "AD")) {
  ci_sub <- ci_mat_trim[, grep(paste0("^", dimname), colnames(ci_mat_trim)), drop = FALSE]
  ci_sub <- ci_sub[rowSums(!is.na(ci_sub)) > 0, , drop = FALSE]
  ci_mod_dim <- tryCatch(mirt::mirt(as.data.frame(ci_sub), 1, itemtype = "Rasch", verbose = FALSE), error = function(e) NULL)
  if (!is.null(ci_mod_dim)) irt_scales[[length(irt_scales) + 1]] <- list(label = paste0("CI-", dimname), mod = ci_mod_dim, serial = rownames(ci_sub))
}

theta_grid <- seq(-4, 4, by = 0.05)
info_long <- list(); items_long <- list(); persons_long <- list(); targeting_rows <- list()
for (sc in irt_scales) {
  res <- tryCatch({
    info <- mirt::testinfo(sc$mod, Theta = matrix(theta_grid))
    ip <- mirt::coef(sc$mod, IRTpars = TRUE, simplify = TRUE)$items
    bcols <- grep("^b", colnames(ip))
    loc <- if (length(bcols) == 1) ip[, bcols] else rowMeans(ip[, bcols, drop = FALSE], na.rm = TRUE)
    fs <- mirt::fscores(sc$mod, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)
    persons <- tibble::tibble(task = sc$label, SERIAL = sc$serial, theta = fs[, 1], se = fs[, 2], Group = assign_group(sc$serial)) %>%
      dplyr::filter(!is.na(Group))
    se_curve <- 1 / sqrt(info)
    good <- theta_grid[se_curve <= 0.55]
    list(info = tibble::tibble(task = sc$label, theta = theta_grid, Information = info, SE = se_curve),
         items = tibble::tibble(task = sc$label, item = rownames(ip), location = loc),
         persons = persons,
         target = tibble::tibble(
           task = sc$label, n_items = length(loc), item_loc_mean = mean(loc), item_loc_sd = stats::sd(loc),
           person_theta_mean = mean(persons$theta), person_theta_sd = stats::sd(persons$theta),
           targeting_offset = mean(persons$theta) - mean(loc),
           pct_persons_within_item_range = 100 * mean(persons$theta >= min(loc) & persons$theta <= max(loc)),
           peak_info_theta = theta_grid[which.max(info)], peak_info = max(info),
           theta_range_se_le_055 = if (length(good)) sprintf("%.2f to %.2f", min(good), max(good)) else "none",
           mean_se_control = mean(persons$se[persons$Group == "Control"]), mean_se_patient = mean(persons$se[persons$Group == "Patient"])))
  }, error = function(e) { cat(sprintf("[Targeting] %s failed: %s\n", sc$label, conditionMessage(e))); NULL })
  if (!is.null(res)) {
    info_long[[sc$label]] <- res$info; items_long[[sc$label]] <- res$items
    persons_long[[sc$label]] <- res$persons; targeting_rows[[sc$label]] <- res$target
  }
}
targeting_table <- dplyr::bind_rows(targeting_rows)
readr::write_csv(targeting_table, file.path(output_dir, "SupplTable06_Targeting.csv"))
cat("[Targeting] Test information, person-item targeting and conditional SE by scale:\n")
print(targeting_table, width = Inf)
if (length(info_long) > 0) {
  info_df <- dplyr::bind_rows(info_long) %>% tidyr::pivot_longer(c(Information, SE), names_to = "measure", values_to = "value")
  fig_info <- ggplot(info_df, aes(x = theta, y = value)) + geom_line(linewidth = 0.8) +
    facet_grid(measure ~ task, scales = "free_y") +
    labs(x = expression(theta), y = NULL, title = "Test information and conditional standard error") + theme_corr
  ggsave(file.path(output_dir, "SupplFigure4_TestInformation.png"), fig_info, width = 14, height = 6, dpi = 300, bg = "white")
  fig_wright <- ggplot(dplyr::bind_rows(persons_long), aes(x = theta, fill = Group)) +
    geom_histogram(bins = 30, alpha = 0.55, position = "identity") +
    geom_rug(data = dplyr::bind_rows(items_long), aes(x = location), inherit.aes = FALSE, sides = "b",
             length = grid::unit(0.08, "npc"), linewidth = 0.7) +
    facet_wrap(~ task, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c(Control = ctrl_color, Patient = pat_color)) +
    labs(x = expression(theta), y = "Persons", title = "Person-item map (bars: persons; ticks: item locations)") + theme_corr
  ggsave(file.path(output_dir, "SupplFigure3_PersonItemMap.png"), fig_wright, width = 12, height = 10, dpi = 300, bg = "white")
}

rasch_vs_2pl_rows <- list(); rasch_vs_2pl_items <- list()
for (sc in irt_scales[vapply(irt_scales, function(x) x$label %in% c("GRT", "DRT-Fam", "DRT-Recoll"), logical(1))]) {
  res <- tryCatch({
    Y <- as.data.frame(mirt::extract.mirt(sc$mod, "data"))
    mod2 <- mirt::mirt(Y, 1, itemtype = "2PL", verbose = FALSE)
    th1 <- mirt::fscores(sc$mod, method = "EAP")[, 1]
    th2 <- mirt::fscores(mod2, method = "EAP")[, 1]
    grp <- assign_group(sc$serial)
    g1 <- compare_groups(th1, grp); g2 <- compare_groups(th2, grp)
    a <- mirt::coef(mod2, simplify = TRUE)$items[, "a1"]
    med_a <- stats::median(a)
    flag <- a < med_a / 2 | a > med_a * 2
    list(row = tibble::tibble(
      task = sc$label, n_items = length(a), N = length(th1),
      r_pearson = stats::cor(th1, th2), r_spearman = stats::cor(th1, th2, method = "spearman"),
      mean_abs_diff_sd = mean(abs(as.numeric(scale(th1)) - as.numeric(scale(th2)))),
      g_rasch = g1$g, g_2pl = g2$g, delta_g = g2$g - g1$g,
      slope_min = min(a), slope_median = med_a, slope_max = max(a), slope_cv = stats::sd(a) / mean(a),
      n_items_flagged = sum(flag),
      items_low = paste(names(a)[a < med_a / 2], collapse = "; "),
      items_high = paste(names(a)[a > med_a * 2], collapse = "; ")),
      items = tibble::tibble(task = sc$label, item = names(a), slope_2pl = a, flagged = flag))
  }, error = function(e) { cat(sprintf("[Rasch vs 2PL] %s failed: %s\n", sc$label, conditionMessage(e))); NULL })
  if (!is.null(res)) { rasch_vs_2pl_rows[[sc$label]] <- res$row; rasch_vs_2pl_items[[sc$label]] <- res$items }
}
rasch_vs_2pl_table <- dplyr::bind_rows(rasch_vs_2pl_rows)
rasch_vs_2pl_slopes <- dplyr::bind_rows(rasch_vs_2pl_items)
if (nrow(rasch_vs_2pl_table) > 0) {
  readr::write_csv(rasch_vs_2pl_table, file.path(output_dir, "SupplTable05_RaschVs2PL.csv"))
  readr::write_csv(rasch_vs_2pl_slopes, file.path(output_dir, "SupplTable05_RaschVs2PL_Slopes.csv"))
  cat("[Rasch vs 2PL] Practical impact of the model choice (theta agreement, group effect, slope spread):\n")
  print(rasch_vs_2pl_table, width = Inf)
}

# --- 4.2 Local dependence (Yen's Q3) -----------------------------------------
# Residual correlations between item pairs; a pair is flagged when its Q3
# exceeds the scale's mean Q3 by more than 0.20 (Christensen et al., 2017),
# a criterion that accounts for the negative bias of Q3 in short scales.

ld_summary <- list(); ld_flagged <- list()
for (sc in irt_scales) {
  q3 <- tryCatch(mirt::residuals(sc$mod, type = "Q3", verbose = FALSE),
                 error = function(e) tryCatch(residuals(sc$mod, type = "Q3", verbose = FALSE), error = function(e2) NULL))
  if (is.null(q3)) { cat(sprintf("[Local dependence] Q3 not computable for %s.\n", sc$label)); next }
  ut <- which(upper.tri(q3), arr.ind = TRUE)
  vals <- q3[ut]; q3_bar <- mean(vals, na.rm = TRUE)
  flag <- !is.na(vals) & (vals - q3_bar) > 0.20
  ld_summary[[sc$label]] <- tibble::tibble(task = sc$label, n_pairs = length(vals), mean_Q3 = q3_bar,
                                           max_Q3 = max(vals, na.rm = TRUE), n_flagged_pairs = sum(flag))
  if (any(flag)) ld_flagged[[sc$label]] <- tibble::tibble(task = sc$label, item_1 = rownames(q3)[ut[flag, 1]],
                                                          item_2 = colnames(q3)[ut[flag, 2]], Q3 = vals[flag], Q3_minus_mean = vals[flag] - q3_bar)
}
ld_summary_table <- dplyr::bind_rows(ld_summary); ld_flagged_table <- dplyr::bind_rows(ld_flagged)
readr::write_csv(ld_summary_table, file.path(output_dir, "SupplTable06_LocalDependence.csv"))
readr::write_csv(ld_flagged_table, file.path(output_dir, "SupplTable06_LocalDependence_FlaggedPairs.csv"))
cat("[Local dependence] Yen's Q3 summary (flag: Q3 - mean Q3 > 0.20):\n"); print(ld_summary_table)

# --- 4.3 Interchangeability of item subsets ----------------------------------
# Each item bank is split 100 times into two random, non-overlapping halves.
# Both halves are scored with the calibrated item parameters held fixed, and
# the halves are compared by correlation (also Spearman-Brown corrected) and
# standardised mean difference. Equivalent estimates from disjoint subsets are
# the item-bank analogue of alternate-form equivalence.

interchangeability_check <- function(mod, label, B = 100, seed = 20260821) {
  Y <- as.data.frame(mirt::extract.mirt(mod, "data"))
  items <- colnames(Y)
  set.seed(seed)
  purrr::map_dfr(seq_len(B), function(b) {
    half_a <- sample(items, floor(length(items) / 2)); half_b <- setdiff(items, half_a)
    keep <- rowSums(!is.na(Y[, half_a, drop = FALSE])) > 0 & rowSums(!is.na(Y[, half_b, drop = FALSE])) > 0
    if (sum(keep) < 20) return(tibble::tibble())
    YA <- Y[keep, , drop = FALSE]; YA[, half_b] <- NA
    YB <- Y[keep, , drop = FALSE]; YB[, half_a] <- NA
    tA <- mirt::fscores(mod, method = "EAP", response.pattern = as.matrix(YA))[, "F1"]
    tB <- mirt::fscores(mod, method = "EAP", response.pattern = as.matrix(YB))[, "F1"]
    r <- stats::cor(tA, tB)
    tibble::tibble(task = label, split = b, n = sum(keep), r = r, r_sb = 2 * r / (1 + r),
                   smd = mean(tA - tB) / stats::sd(c(tA, tB)))
  })
}
if (RUN_INTERCHANGEABILITY) {
  ic_splits <- dplyr::bind_rows(lapply(irt_scales[vapply(irt_scales, function(sc) sc$label %in% c("GRT", "DRT-Fam", "DRT-Recoll"), logical(1))],
                                       function(sc) tryCatch(interchangeability_check(sc$mod, sc$label),
                                                             error = function(e) { cat(sprintf("[Interchangeability] %s failed: %s\n", sc$label, conditionMessage(e))); NULL })))
  if (nrow(ic_splits) > 0) {
    ic_summary <- ic_splits %>% dplyr::group_by(task) %>%
      dplyr::summarise(n_splits = dplyr::n(), median_n = stats::median(n),
                       r_median = stats::median(r), r_2.5 = stats::quantile(r, .025), r_97.5 = stats::quantile(r, .975),
                       r_sb_median = stats::median(r_sb),
                       smd_median = stats::median(smd), smd_abs_max = max(abs(smd)), .groups = "drop")
    readr::write_csv(ic_splits, file.path(output_dir, "SupplTable06_Interchangeability_AllSplits.csv"))
    readr::write_csv(ic_summary, file.path(output_dir, "SupplTable06_Interchangeability.csv"))
    cat("[Interchangeability] Random disjoint item halves scored with fixed calibrated parameters:\n"); print(ic_summary, width = Inf)
  }
} else {
  cat("[Interchangeability] Skipped -- RUN_INTERCHANGEABILITY = FALSE.\n")
}

# =============================================================================
# 5. DIFFERENTIAL ITEM FUNCTIONING (Suppl. Tables 7-8)
# =============================================================================
# Grouping variables: patient status, control recruitment phase, sex,
# education (Abitur vs. lower) and age (median split). Multiple-group Rasch
# models are fitted with all other items as anchors; each item is tested by a
# likelihood-ratio test of freeing its difficulty across groups.
# Benjamini-Hochberg correction is applied within each scale and grouping
# variable. Groups with fewer than 30 persons are not tested. The impact of
# flagged items is evaluated by re-scoring each scale without them.

if (RUN_DIF_ANALYSIS) {
  
  dif_min_group_n <- 30
  
  run_dif <- function(item_mat, group, task_label, dimension_label, grouping_label, itemtype = "Rasch") {
    item_mat <- as.data.frame(item_mat)
    
    ok <- !is.na(group) & rowSums(!is.na(item_mat)) > 0
    item_mat <- item_mat[ok, , drop = FALSE]
    group <- factor(as.character(group[ok]))
    
    keep_cols <- vapply(item_mat, function(col) {
      all(vapply(split(col, group), function(x) length(unique(stats::na.omit(x))) >= 2, logical(1)))
    }, logical(1))
    item_mat <- item_mat[, keep_cols, drop = FALSE]
    
    if (nlevels(group) < 2 || any(table(group) < dif_min_group_n) || ncol(item_mat) < 3) {
      return(tibble::tibble(Task = task_label, Dimension = dimension_label, Grouping = grouping_label,
                            item = NA_character_, status = "skipped: insufficient N per group or < 3 usable items",
                            group_n = paste(names(table(group)), table(group), sep = "=", collapse = "; ")))
    }
    
    which_par <- if (itemtype == "gpcm") {
      max_ncat <- max(vapply(item_mat, function(x) length(unique(stats::na.omit(x))), integer(1)))
      paste0("d", seq_len(max(max_ncat - 1, 1)))
    } else {
      "d"
    }
    
    result <- tryCatch({
      mod0 <- mirt::multipleGroup(
        item_mat, model = 1, group = group, itemtype = itemtype,
        invariance = c("free_means", "free_var", colnames(item_mat)),
        SE = FALSE, verbose = FALSE, technical = list(NCYCLES = 500)
      )
      dif_out <- mirt::DIF(mod0, which.par = which_par, scheme = "drop", verbose = FALSE)
      dif_out <- tibble::as_tibble(dif_out, rownames = "item")
      
      if ("X2" %in% names(dif_out)) {
        neg_x2 <- !is.na(dif_out$X2) & dif_out$X2 < 0
        if (any(neg_x2)) {
          dif_out$X2[neg_x2] <- 0
          if ("df" %in% names(dif_out) && "p" %in% names(dif_out)) {
            dif_out$p[neg_x2] <- stats::pchisq(0, df = dif_out$df[neg_x2], lower.tail = FALSE)
          }
        }
      }
      dif_out %>% dplyr::mutate(status = "ok")
    }, error = function(e) {
      tibble::tibble(item = NA_character_, status = paste("error:", conditionMessage(e)))
    })
    
    result %>% dplyr::mutate(Task = task_label, Dimension = dimension_label, Grouping = grouping_label,
                             group_n = paste(names(table(group)), table(group), sep = "=", collapse = "; "))
  }
  
  dif_task_defs <- list(
    list(task = "GRT", dim = NA_character_, item_mat = grt_Y_final,
         serial = grt_X$SERIAL[grt_rows_any], itemtype = "Rasch"),
    list(task = "DRT", dim = "Fam", item_mat = drt_Ybin_Fam[, drt_final_Fam, drop = FALSE],
         serial = drt_X$SERIAL, itemtype = "Rasch"),
    list(task = "DRT", dim = "Recoll", item_mat = drt_Ybin_Recoll[, drt_final_Recoll, drop = FALSE],
         serial = drt_X$SERIAL, itemtype = "Rasch"),
    list(task = "CI", dim = NA_character_, item_mat = ci_mat_trim,
         serial = rownames(ci_mat_trim), itemtype = "Rasch"),
    list(task = "PAL", dim = NA_character_, item_mat = pal_Y_ord,
         serial = pal_serial_kept, itemtype = "Rasch")
  )
  
  dif_groupings <- list(
    "Patient_vs_Control" = function(serial, task) assign_group(serial),
    "Control_recruitment_phase" = function(serial, task) assign_control_phase(serial),
    "Sex" = function(serial, task) {
      d <- df_sex$D_Male[match(serial, df_sex$SERIAL)]
      dplyr::if_else(d == 1, "Male", dplyr::if_else(d == 0, "Female", NA_character_))
    },
    "Education_Abitur" = function(serial, task) {
      d <- df_edu$D_Abitur[match(serial, df_edu$SERIAL)]
      dplyr::if_else(d == 1, "Abitur", dplyr::if_else(d == 0, "BelowAbitur", NA_character_))
    },
    "Age_median_split" = function(serial, task) {
      a <- df_age$age[match(serial, df_age$SERIAL)]
      med <- stats::median(a, na.rm = TRUE)
      dplyr::if_else(a <= med, "Younger", dplyr::if_else(a > med, "Older", NA_character_))
    }
  )
  
  dif_results <- list()
  for (td in dif_task_defs) {
    for (gname in names(dif_groupings)) {
      group_vec <- dif_groupings[[gname]](td$serial, td$task)
      cat(sprintf("[DIF] Running %s (%s) x %s ...\n", td$task, dplyr::coalesce(td$dim, "-"), gname))
      dif_results[[length(dif_results) + 1]] <- run_dif(
        td$item_mat, group_vec, td$task, td$dim, gname, itemtype = td$itemtype
      )
    }
  }
  dif_table <- dplyr::bind_rows(dif_results)
  
  dif_table <- dif_table %>%
    dplyr::group_by(Task, Dimension, Grouping) %>%
    dplyr::mutate(p_adj = if (all(status == "ok")) stats::p.adjust(p, method = "BH") else NA_real_) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dif_flagged = status == "ok" & !is.na(p_adj) & p_adj < .05)
  
  readr::write_csv(dif_table, file.path(output_dir, "SupplTable07_DIF_ItemLevel.csv"))
  
  dif_summary <- dif_table %>%
    dplyr::group_by(Task, Dimension, Grouping) %>%
    dplyr::summarise(
      n_items_tested = sum(status == "ok"),
      n_flagged      = sum(dif_flagged, na.rm = TRUE),
      pct_flagged    = dplyr::if_else(n_items_tested > 0, round(100 * n_flagged / n_items_tested, 1), NA_real_),
      flagged_items  = paste(item[dif_flagged %in% TRUE], collapse = "; "),
      .groups = "drop"
    )
  readr::write_csv(dif_summary, file.path(output_dir, "SupplTable07_DIF.csv"))
  
  n_flagged <- sum(dif_table$dif_flagged, na.rm = TRUE)
  cat(sprintf("\n[DIF] Done. %d item x grouping tests completed; %d flagged after BH correction (p_adj < .05, within each Task x Dimension x Grouping combination). See SupplTable07_DIF_ItemLevel.csv (per item) and SupplTable07_DIF.csv (per combination).\n",
              sum(dif_table$status == "ok"), n_flagged))
  
  dif_flagged_items <- dif_table %>% dplyr::filter(dif_flagged) %>% dplyr::distinct(Task, Dimension, item)
  get_flagged <- function(task, dim = NA_character_) {
    dif_flagged_items %>%
      dplyr::filter(Task == task, (is.na(Dimension) & is.na(dim)) | Dimension == dim) %>%
      dplyr::pull(item)
  }
  
  summarize_dif_sensitivity <- function(task_label, dim_label, flagged_items,
                                        serial_orig, theta_orig, group_orig,
                                        serial_new, theta_new) {
    theta_new_aligned <- theta_new[match(serial_orig, serial_new)]
    gc_orig <- compare_groups(theta_orig, group_orig)
    gc_new  <- compare_groups(theta_new_aligned, group_orig)
    tibble::tibble(
      Task = task_label, Dimension = dim_label,
      n_items_dropped = length(flagged_items), items_dropped = paste(flagged_items, collapse = "; "),
      r_theta_orig_vs_excl = suppressWarnings(cor(theta_orig, theta_new_aligned, use = "complete.obs")),
      g_original      = if (!is.null(gc_orig)) gc_orig$g else NA_real_,
      g_excluding_dif = if (!is.null(gc_new))  gc_new$g  else NA_real_,
      delta_g         = if (!is.null(gc_orig) && !is.null(gc_new)) gc_new$g - gc_orig$g else NA_real_
    )
  }
  
  dif_sensitivity <- list()
  
  grt_flagged <- get_flagged("GRT")
  if (length(grt_flagged) > 0) {
    grt_keep2 <- setdiff(grt_final_items, grt_flagged)
    grt_theta_new <- rep(NA_real_, nrow(grt_X))
    grt_theta_new[grt_rows_any] <- rasch_theta_eap(grt_Y_final[, grt_keep2, drop = FALSE])
    dif_sensitivity[["GRT"]] <- summarize_dif_sensitivity(
      "GRT", NA_character_, grt_flagged,
      grt_X$SERIAL, grt_theta_all, assign_group(grt_X$SERIAL),
      grt_X$SERIAL, grt_theta_new
    )
  }
  
  drt_d1_flagged <- get_flagged("DRT", "Fam")
  if (length(drt_d1_flagged) > 0) {
    keep2 <- setdiff(drt_final_Fam, drt_d1_flagged)
    theta_new <- fit_dim_theta(drt_Y_bin_rm, keep2)$theta
    dif_sensitivity[["DRT_Fam"]] <- summarize_dif_sensitivity(
      "DRT", "Fam", drt_d1_flagged,
      drt_X$SERIAL, drt_scores$theta_DRT_Fam, assign_group(drt_X$SERIAL),
      drt_X$SERIAL, theta_new
    )
  }
  
  drt_recoll_flagged <- get_flagged("DRT", "Recoll")
  if (length(drt_recoll_flagged) > 0) {
    keep2 <- setdiff(drt_final_Recoll, drt_recoll_flagged)
    theta_new <- fit_dim_theta(drt_Y_bin_rm, keep2)$theta
    dif_sensitivity[["DRT_Recoll"]] <- summarize_dif_sensitivity(
      "DRT", "Recoll", drt_recoll_flagged,
      drt_X$SERIAL, drt_scores$theta_DRT_Recoll, assign_group(drt_X$SERIAL),
      drt_X$SERIAL, theta_new
    )
  }
  
  pal_flagged <- get_flagged("PAL")
  if (length(pal_flagged) > 0) {
    pal_keep2 <- setdiff(colnames(pal_Y_ord), pal_flagged)
    pal_spec2 <- mirt::mirt.model(paste0("F = 1-", length(pal_keep2), "\nCONSTRAIN = (1-", length(pal_keep2), ", a1)"))
    pal_mod2  <- mirt::mirt(pal_Y_ord[, pal_keep2, drop = FALSE], model = pal_spec2, itemtype = "gpcm", verbose = FALSE)
    pal_fs2   <- mirt::fscores(pal_mod2, method = "EAP", full.scores = TRUE)
    dif_sensitivity[["PAL"]] <- summarize_dif_sensitivity(
      "PAL", NA_character_, pal_flagged,
      pal_serial_kept, pal_scores$theta, pal_group_kept,
      pal_serial_kept, pal_fs2[, "F"]
    )
  }
  
  ci_flagged <- get_flagged("CI")
  if (length(ci_flagged) > 0) {
    ci_keep2 <- setdiff(colnames(ci_mat_trim), ci_flagged)
    ci_mat2  <- ci_mat_trim[, ci_keep2, drop = FALSE]
    ci_Q2    <- ci_Q_trim[rownames(ci_Q_trim) %in% ci_keep2, , drop = FALSE]
    ci_mod2  <- TAM::tam.mml(resp = ci_mat2, Q = ci_Q2, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
    ci_wle2  <- TAM::tam.wle(ci_mod2, progress = FALSE)
    ci_cov2  <- ci_mod2$variance
    ci_beta2 <- ci_cov2[1, 2] / ci_cov2[1, 1]
    ci_intercept2 <- mean(ci_wle2$theta.Dim02, na.rm = TRUE) - ci_beta2 * mean(ci_wle2$theta.Dim01, na.rm = TRUE)
    ci_resil2   <- ci_wle2$theta.Dim02 - (ci_intercept2 + ci_beta2 * ci_wle2$theta.Dim01)
    ci_serial2  <- rownames(ci_mat2)
    dif_sensitivity[["CI"]] <- summarize_dif_sensitivity(
      "CI", NA_character_, ci_flagged,
      ci_theta$SERIAL, ci_theta$CI_Resilience, ci_theta$Group,
      ci_serial2, ci_resil2
    )
  }
  
  if (length(dif_sensitivity) > 0) {
    dif_sensitivity_table <- dplyr::bind_rows(dif_sensitivity)
    readr::write_csv(dif_sensitivity_table, file.path(output_dir, "SupplTable08_DIF_Sensitivity.csv"))
    cat(sprintf("[DIF] Sensitivity re-score (excluding flagged items) written for %d task(s) with at least one flagged item: %s.\n",
                length(dif_sensitivity), paste(names(dif_sensitivity), collapse = ", ")))
  } else {
    cat("[DIF] No items flagged -- no sensitivity re-score needed.\n")
  }
}

# =============================================================================
# 6. STANDARD ERROR OF MEASUREMENT (Supplementary Table 2)
# =============================================================================
# Cross-sectional SEM = SD x sqrt(1 - reliability), reported in SD units.
# Split-half reliabilities are computed with the parcels in Section 2. The
# smallest detectable change requires retest data and is not estimated.

sem_specs <- list(
  list(label = "Reasoning (GRT)",                    col = "theta_GRT",           rel = grt_reliability_sb),
  list(label = "Familiarity-based recognition",      col = "theta_DRT_Fam",       rel = drt_reliability_sb),
  list(label = "Recollection-based discrimination",  col = "theta_DRT_Recoll",    rel = drt_Recoll_reliability_sb),
  list(label = "Associative memory (PAL)",           col = "theta_PAL",           rel = pal_reliability_sb),
  list(label = "Internal interference (CI)",         col = "theta_CI_resilience", rel = ci_reliability_sb),
  list(label = "External interference (SART-ED)",    col = "theta_SART_EI",       rel = sart_reliability_sb),
  list(label = "Mental speed (SART-ED)",             col = "theta_SART_MS",       rel = ms_reliability_sb)
)
sem_table <- purrr::map_dfr(sem_specs, function(sp) {
  x <- master_theta[[sp$col]]
  sdx <- stats::sd(x, na.rm = TRUE)
  tibble::tibble(score = sp$label, N = sum(is.finite(x)), SD = sdx, reliability_split_half = sp$rel,
                 SEM = sdx * sqrt(1 - sp$rel), SEM_in_SD_units = sqrt(1 - sp$rel))
})
readr::write_csv(sem_table, file.path(output_dir, "SupplTable02_SEM.csv"))
cat("[SEM] Cross-sectional standard error of measurement (SDC not estimable without retest data):\n"); print(sem_table)

# =============================================================================
# 7. KNOWN-GROUPS VALIDITY AND ROBUSTNESS (Figure 3; Suppl. Tables 9-15)
# =============================================================================
# Main comparison: Welch's t-tests and Hedges' g, Holm-corrected across the
# seven scores (Section 3.4; Supplementary Table 9). The condition x group
# interaction of the SART-ED (Supplementary Table 11) is computed in 2.5.2.
# 7.1 Figure 3: raincloud plots of the ability scores by group.

domain_labels <- c(
  theta_SART_EI       = "External interference control",
  theta_CI_resilience = "Internal interference control",
  theta_SART_MS       = "Mental speed",
  theta_GRT           = "Reasoning",
  theta_DRT_Fam       = "Familiarity-based recognition",
  theta_DRT_Recoll    = "Recollection-based discrimination",
  theta_PAL           = "Associative memory"
)
domain_order <- unname(domain_labels)

fig_kg_data <- master_theta %>%
  dplyr::select(SERIAL, Group, dplyr::all_of(names(domain_labels))) %>%
  tidyr::pivot_longer(cols = -c(SERIAL, Group), names_to = "raw_test", values_to = "score") %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::mutate(
    Test  = factor(domain_labels[raw_test], levels = domain_order),
    Group = factor(dplyr::if_else(Group == "Control", "Controls", "Patients"), levels = c("Controls", "Patients"))
  )
levels(fig_kg_data$Test) <- stringr::str_wrap(levels(fig_kg_data$Test), width = 18)

fig_kg_sig <- fig_kg_data %>%
  dplyr::group_by(Test) %>%
  dplyr::summarise(p_val = tryCatch(t.test(score ~ Group)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  dplyr::mutate(p_val = stats::p.adjust(p_val, method = "holm")) %>%
  dplyr::mutate(
    asterisks = dplyr::case_when(is.na(p_val) ~ "", p_val < .001 ~ "***", p_val < .01 ~ "**", p_val < .05 ~ "*", TRUE ~ "ns"),
    y_pos  = max(fig_kg_data$score, na.rm = TRUE) + 0.10 * diff(range(fig_kg_data$score, na.rm = TRUE)),
    y_text = max(fig_kg_data$score, na.rm = TRUE) + 0.16 * diff(range(fig_kg_data$score, na.rm = TRUE))
  )

fig_kg_data <- fig_kg_data %>%
  dplyr::mutate(group_pos = as.numeric(factor(Group, levels = c("Controls", "Patients"))))

compute_half_density <- function(data, width_max = 0.38) {
  data %>%
    dplyr::group_by(Test, Group, group_pos) %>%
    dplyr::group_modify(function(d, key) {
      if (nrow(d) < 3) return(tibble::tibble())
      dens <- stats::density(d$score, n = 256)
      tibble::tibble(y = dens$x, dens_scaled = dens$y / max(dens$y) * width_max)
    }) %>%
    dplyr::ungroup()
}
fig_kg_density <- compute_half_density(fig_kg_data)

set.seed(20260821)
fig_kg_data <- fig_kg_data %>%
  dplyr::mutate(x_jitter = group_pos - stats::runif(dplyr::n(), 0.04, 0.34))

fig_known_groups <- ggplot() +
  geom_ribbon(data = fig_kg_density, aes(y = y, xmin = group_pos, xmax = group_pos + dens_scaled, fill = Group, group = interaction(Group, Test)),
              orientation = "y", alpha = 0.6, color = NA) +
  geom_point(data = fig_kg_data, aes(x = x_jitter, y = score, color = Group), alpha = 0.5, size = 1) +
  geom_boxplot(data = fig_kg_data, aes(x = group_pos, y = score, group = group_pos), width = 0.1, outlier.shape = NA,
               fill = "white", alpha = 0.8, position = position_nudge(x = 0.03)) +
  scale_fill_manual(values = c(Controls = ctrl_color, Patients = pat_color)) +
  scale_color_manual(values = c(Controls = ctrl_color, Patients = pat_color)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Controls", "Patients"), expand = expansion(mult = 0.15)) +
  facet_grid(Test ~ ., scales = "free_x", switch = "y") +
  coord_flip(clip = "off") +
  labs(x = "", y = expression("Score (" * theta * ")")) +
  theme_corr +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 14, color = "black", lineheight = 1.1),
    strip.placement   = "outside",
    axis.text.y       = element_text(size = 12, color = "black"),
    panel.spacing     = unit(1.2, "lines"),
    legend.position   = "bottom", legend.title = element_blank(),
    plot.margin       = margin(t = 5, r = 40, b = 5, l = 5)
  ) +
  geom_segment(data = fig_kg_sig, aes(x = 1, xend = 2, y = y_pos, yend = y_pos), inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = fig_kg_sig, aes(x = 1.5, y = y_text, label = asterisks), inherit.aes = FALSE, size = 5, vjust = 0.5, fontface = "bold")

ggsave(file.path(output_dir, "Figure3_GroupComparison_Raincloud.png"), fig_known_groups, width = 12, height = 14, dpi = 300, bg = "transparent")

# --- 7.2 Latent mean differences (Suppl. Table 10) and control recruitment
#         phase (Suppl. Table 14) -----------------------------------------------
# EAP scores shrink towards the mean when few items are answered. Group
# differences are therefore also estimated directly in multiple-group IRT
# models (item parameters equal across groups; patient latent mean and
# variance free). Comparisons with each control phase separately, and DRT
# re-scored on the items both phases received (items 1-18), test whether the
# recognition results depend on the reference sample.

latent_mean_robustness <- function(item_mat, serial, task, score, itemtype = "Rasch") {
  grp <- assign_group(serial)
  Y <- as.data.frame(item_mat)
  keep <- !is.na(grp) & rowSums(!is.na(Y)) > 0
  Y <- Y[keep, , drop = FALSE]; g <- factor(grp[keep], levels = c("Control", "Patient"))
  
  ok <- vapply(Y, function(x) all(tapply(!is.na(x), g, any)) && length(unique(stats::na.omit(x))) >= 2, logical(1))
  Y <- Y[, ok, drop = FALSE]
  mg <- mirt::multipleGroup(Y, 1, group = g, itemtype = itemtype, SE = TRUE, verbose = FALSE,
                            invariance = c("slopes", "intercepts", "free_means", "free_var"))
  cf <- mirt::coef(mg, printSE = TRUE)
  gpH <- cf$Control$GroupPars; gpP <- cf$Patient$GroupPars
  n <- table(g)
  var_H <- gpH["par", "COV_11"]; var_P <- gpP["par", "COV_11"]
  mu_P <- gpP["par", "MEAN_1"]; se_mu <- gpP["SE", "MEAN_1"]
  pooled_sd <- sqrt(((n[["Control"]] - 1) * var_H + (n[["Patient"]] - 1) * var_P) / (sum(n) - 2))
  d <- (0 - mu_P) / pooled_sd
  tibble::tibble(Task = task, Score = score, N_ctrl = n[["Control"]], N_pat = n[["Patient"]], n_items = ncol(Y),
                 d_latent = d, d_ci_lo = d - 1.96 * se_mu / pooled_sd, d_ci_hi = d + 1.96 * se_mu / pooled_sd,
                 p_latent = 2 * stats::pnorm(-abs(mu_P / se_mu)), latent_var_ratio_pat_ctrl = var_P / var_H)
}
if (RUN_LATENT_MEAN_ROBUSTNESS) {
  lm_specs <- list(
    list(item_mat = grt_Y_final, serial = grt_X$SERIAL[grt_rows_any], task = "GRT", score = "theta", itemtype = "Rasch"),
    list(item_mat = drt_Y_bin_rm[, drt_final_Fam, drop = FALSE], serial = drt_X$SERIAL, task = "DRT", score = "Fam", itemtype = "Rasch"),
    list(item_mat = drt_Y_bin_rm[, drt_final_Recoll, drop = FALSE], serial = drt_X$SERIAL, task = "DRT", score = "Recoll", itemtype = "Rasch"),
    list(item_mat = pal_Y_ord, serial = pal_serial_kept, task = "PAL", score = "theta", itemtype = "Rasch")
  )
  latent_mean_table <- purrr::map_dfr(lm_specs, function(sp) {
    cat(sprintf("[Latent-mean check] Fitting multiple-group model: %s %s ...\n", sp$task, sp$score))
    tryCatch(latent_mean_robustness(sp$item_mat, sp$serial, sp$task, sp$score, sp$itemtype),
             error = function(e) { cat(sprintf("  failed: %s\n", conditionMessage(e))); tibble::tibble(Task = sp$task, Score = sp$score) })
  }) %>%
    dplyr::left_join(dplyr::select(effect_size_table, Task, Score, g_EAP = g, g_EAP_ci_lo = g_ci_lo, g_EAP_ci_hi = g_ci_hi), by = c("Task", "Score"))
  readr::write_csv(latent_mean_table, file.path(output_dir, "SupplTable10_LatentMeans.csv"))
  cat("[Latent-mean check] Latent standardized group difference (multiple-group IRT) vs. EAP-based Hedges g:\n")
  print(latent_mean_table, width = Inf)
  cat("[Latent-mean check] |d_latent| > |g_EAP| is the pattern expected from EAP shrinkage attenuating group differences. CI (residualised 2D score) and SART-ED (condition x group interaction, Section 2.5.2) are not included here.\n")
} else {
  cat("[Latent-mean check] Skipped -- RUN_LATENT_MEAN_ROBUSTNESS = FALSE.\n")
}

phase_scores <- c(`GRT (reasoning)` = "theta_GRT",
                  `DRT-Fam (familiarity-based recognition)` = "theta_DRT_Fam",
                  `DRT-Recoll (recollection-based discrimination)` = "theta_DRT_Recoll",
                  `PAL (associative memory)` = "theta_PAL",
                  `CI (internal interference control)` = "theta_CI_resilience",
                  `SART-ED (external interference control)` = "theta_SART_EI",
                  `SART-ED (mental speed)` = "theta_SART_MS")
phase_cohort <- dplyr::if_else(master_theta$Group == "Patient", "Patient", assign_control_phase(master_theta$SERIAL, master_theta$Group))
phase_comparisons <- list(
  c(first = "Phase1", second = "Patient", label = "Phase-1 controls vs. patients"),
  c(first = "Phase2", second = "Patient", label = "Phase-2 controls vs. patients"),
  c(first = "Phase1", second = "Phase2",  label = "Phase-1 vs. phase-2 controls")
)
phase_table <- purrr::map_dfr(names(phase_scores), function(sc) {
  x <- master_theta[[phase_scores[[sc]]]]
  n_old <- sum(is.finite(x) & phase_cohort %in% "Phase1")
  if (n_old < 10) return(tibble::tibble())
  purrr::map_dfr(phase_comparisons, function(cmp) {
    keep <- phase_cohort %in% cmp[c("first", "second")]
    gc <- compare_groups(x[keep], phase_cohort[keep], group_levels = unname(cmp[c("first", "second")]))
    if (is.null(gc)) return(tibble::tibble())
    d1 <- gc$desc[gc$desc$Group == cmp[["first"]], ]; d2 <- gc$desc[gc$desc$Group == cmp[["second"]], ]
    tibble::tibble(Score = sc, Comparison = cmp[["label"]],
                   N_first = d1$N, M_first = d1$M, SD_first = d1$SD,
                   N_second = d2$N, M_second = d2$M, SD_second = d2$SD,
                   t = gc$t, df = gc$df, p = gc$p, g = gc$g, g_ci_lo = gc$g_ci[1], g_ci_hi = gc$g_ci[2])
  })
})
if (nrow(phase_table) > 0) {
  phase_table$p_holm <- stats::p.adjust(phase_table$p, method = "holm")
  readr::write_csv(phase_table, file.path(output_dir, "SupplTable14_ControlPhase.csv"))
  cat("[Phase] Known-groups comparisons by control recruitment phase (positive g = first-named group higher):\n")
  print(phase_table, n = Inf, width = Inf)
} else {
  cat("[Phase] No score has >= 10 phase-1 controls -- comparison by phase skipped.\n")
}

drt_common_items_max <- 18
drt_common_rows <- list()
for (dim_spec in list(list(label = "DRT-Fam (familiarity-based recognition)", mod = drt_fit1$mod, items = drt_final_Fam),
                      list(label = "DRT-Recoll (recollection-based discrimination)", mod = drt_fit2$mod, items = drt_final_Recoll))) {
  res <- tryCatch({
    Y <- as.data.frame(mirt::extract.mirt(dim_spec$mod, "data"))
    serial <- drt_X$SERIAL[rowSums(!is.na(drt_Y_bin_rm[, dim_spec$items, drop = FALSE])) > 0]
    common <- colnames(Y)[as.integer(gsub("\\D", "", colnames(Y))) <= drt_common_items_max]
    if (length(common) < 3) stop(sprintf("only %d common items", length(common)))
    Yc <- Y; Yc[, setdiff(colnames(Y), common)] <- NA
    keep <- rowSums(!is.na(Yc[, common, drop = FALSE])) > 0
    th <- rep(NA_real_, nrow(Y))
    th[keep] <- mirt::fscores(dim_spec$mod, method = "EAP", response.pattern = as.matrix(Yc[keep, , drop = FALSE]))[, "F1"]
    coh <- dplyr::if_else(assign_group(serial) == "Patient", "Patient", assign_control_phase(serial, assign_group(serial)))
    purrr::map_dfr(phase_comparisons, function(cmp) {
      sel <- coh %in% cmp[c("first", "second")] & is.finite(th)
      gc <- compare_groups(th[sel], coh[sel], group_levels = unname(cmp[c("first", "second")]))
      if (is.null(gc)) return(tibble::tibble())
      d1 <- gc$desc[gc$desc$Group == cmp[["first"]], ]; d2 <- gc$desc[gc$desc$Group == cmp[["second"]], ]
      tibble::tibble(Score = dim_spec$label, n_common_items = length(common), common_items = paste(common, collapse = "; "),
                     Comparison = cmp[["label"]], N_first = d1$N, M_first = d1$M, SD_first = d1$SD, N_second = d2$N, M_second = d2$M, SD_second = d2$SD,
                     t = gc$t, df = gc$df, p = gc$p, g = gc$g, g_ci_lo = gc$g_ci[1], g_ci_hi = gc$g_ci[2])
    })
  }, error = function(e) { cat(sprintf("[DRT common items] %s failed: %s\n", dim_spec$label, conditionMessage(e))); NULL })
  if (!is.null(res)) drt_common_rows[[dim_spec$label]] <- res
}
drt_common_items_table <- dplyr::bind_rows(drt_common_rows)
if (nrow(drt_common_items_table) > 0) {
  drt_common_items_table$p_holm <- stats::p.adjust(drt_common_items_table$p, method = "holm")
  readr::write_csv(drt_common_items_table, file.path(output_dir, "SupplTable14_ControlPhase_CommonItems.csv"))
  cat("[DRT common items] Phase comparisons re-scored on items 1-18 only (fixed calibrated parameters; positive g = first-named group higher):\n")
  print(dplyr::select(drt_common_items_table, -common_items), n = Inf, width = Inf)
}

# --- 7.3 Joint covariate model (Supplementary Table 12) ----------------------
# Each score is regressed simultaneously on group, age, sex and education;
# recruitment phase is added for tasks completed by both control phases
# (GRT, DRT). Patients are coded 0 on the phase indicator, so the group
# coefficient compares patients with phase-1 controls. For tasks without
# phase-1 data the phase indicator is collinear with group and dropped by lm().
covariate_tasks <- list(
  list(theta = master_theta$theta_SART_EI,       label = "SART-ED (external interference)"),
  list(theta = master_theta$theta_CI_resilience, label = "CI (internal interference)"),
  list(theta = master_theta$theta_SART_MS,       label = "SART-ED (mental speed)"),
  list(theta = master_theta$theta_GRT,           label = "GRT (reasoning)"),
  list(theta = master_theta$theta_DRT_Fam,       label = "DRT-Fam (familiarity-based recognition)"),
  list(theta = master_theta$theta_DRT_Recoll,    label = "DRT-Recoll (recollection-based discrimination)"),
  list(theta = master_theta$theta_PAL,           label = "PAL (associative memory)")
)

joint_covariate_model <- function(theta, serial, group, task_label, score_label) {
  control_phase <- assign_control_phase(serial, group)
  is_phase2 <- dplyr::case_when(
    group == "Patient" ~ 0L,
    control_phase == "Phase2" ~ 1L,
    control_phase == "Phase1" ~ 0L,
    TRUE ~ NA_integer_
  )
  d <- tibble::tibble(
    theta     = theta,
    Group     = factor(group, levels = c("Control", "Patient")),
    is_phase2 = is_phase2,
    age       = df_age$age[match(serial, df_age$SERIAL)],
    D_Male    = df_sex$D_Male[match(serial, df_sex$SERIAL)],
    D_Abitur  = df_edu$D_Abitur[match(serial, df_edu$SERIAL)]
  ) %>% dplyr::filter(!is.na(theta), !is.na(Group))
  
  fit <- tryCatch(lm(theta ~ Group + is_phase2 + age + D_Male + D_Abitur, data = d), error = function(e) NULL)
  if (is.null(fit)) return(tibble::tibble(Task = task_label, Score = score_label, term = NA_character_, status = "model failed to fit"))
  
  co <- summary(fit)$coefficients
  tibble::tibble(
    Task = task_label, Score = score_label, term = rownames(co),
    estimate = co[, "Estimate"], se = co[, "Std. Error"], t = co[, "t value"], p = co[, "Pr(>|t|)"],
    N = stats::nobs(fit), status = "ok"
  )
}

joint_covariate_table <- purrr::map_dfr(covariate_tasks, ~ joint_covariate_model(
  .x$theta, master_theta$SERIAL, master_theta$Group, .x$label, "theta"))
readr::write_csv(joint_covariate_table, file.path(output_dir, "SupplTable12_JointCovariateModel.csv"))
cat("[Covariates] Joint covariate model (group, recruitment phase, age, sex, education) written (Supplementary Table 12).\n")

# --- 7.4 Setting-robust checks (Supplementary Table 13) ----------------------
# Patients were tested supervised on study devices, controls remotely on own
# devices. Two analyses do not depend on this confound:
#  (a) Severity gradient within patients (identical setting and device):
#      Spearman correlations with MoCA and NIHSS, and comparisons of controls
#      with MoCA-unimpaired and MoCA-impaired patients (and NIHSS 0-3 vs. >= 4).
#  (b) Setting-independent SART-ED measures: NoGo accuracy (trials without
#      distractor) and post-error slowing, i.e. mean log RT of the first three
#      correct Go responses after a commission error minus mean log RT of
#      correct Go responses after correct trials (participants with >= 2
#      errors). Post-error slowing is a within-person RT difference, so device
#      latency and general response speed cancel out.
setting_checks <- list()
holm_tbl <- function(tb) { tb$p_holm <- stats::p.adjust(tb$p, "holm"); tb }
cmp_row <- function(x, g, lv, label, comparison) {
  r <- compare_groups(x, g, group_levels = lv)
  if (is.null(r)) return(NULL)
  d1 <- r$desc[r$desc$Group == lv[1], ]; d2 <- r$desc[r$desc$Group == lv[2], ]
  tibble::tibble(measure = label, comparison = comparison, n_first = d1$N, M_first = d1$M, SD_first = d1$SD,
                 n_second = d2$N, M_second = d2$M, SD_second = d2$SD, t = r$t, df = r$df, p = r$p,
                 g = r$g, g_ci_lo = r$g_ci[1], g_ci_hi = r$g_ci[2])
}

sev <- master_theta %>%
  dplyr::select(SERIAL, Group, theta_SART_MS, theta_GRT, MoCA) %>%
  dplyr::left_join(dplyr::select(stroke_aetiology, SERIAL, NIHSS_score), by = "SERIAL") %>%
  dplyr::mutate(moca30 = MoCA / 100 * 30,
                moca_grp = dplyr::case_when(Group == "Control" ~ "Controls",
                                            Group == "Patient" & moca30 >= moca_impairment_cutoff ~ "Patients, MoCA unimpaired",
                                            Group == "Patient" & moca30 <  moca_impairment_cutoff ~ "Patients, MoCA impaired"),
                nihss_grp = dplyr::case_when(Group == "Control" ~ "Controls",
                                             Group == "Patient" & NIHSS_score <= 3 ~ "Patients, NIHSS 0-3",
                                             Group == "Patient" & NIHSS_score >= 4 ~ "Patients, NIHSS >= 4"))
pt <- sev %>% dplyr::filter(Group == "Patient")
sp <- function(x, y) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 10) return(c(rho = NA, p = NA, n = sum(ok)))
ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman")); c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok)) }
sev_cor <- dplyr::bind_rows(
  tibble::as_tibble(as.list(sp(pt$theta_SART_MS, pt$MoCA)))         %>% dplyr::mutate(measure = "Mental speed", correlate = "MoCA (%)"),
  tibble::as_tibble(as.list(sp(pt$theta_SART_MS, -pt$NIHSS_score))) %>% dplyr::mutate(measure = "Mental speed", correlate = "NIHSS (reversed: higher = milder)"),
  tibble::as_tibble(as.list(sp(pt$theta_GRT, pt$MoCA)))             %>% dplyr::mutate(measure = "Reasoning", correlate = "MoCA (%)"),
  tibble::as_tibble(as.list(sp(pt$theta_GRT, -pt$NIHSS_score)))     %>% dplyr::mutate(measure = "Reasoning", correlate = "NIHSS (reversed: higher = milder)")
) %>% dplyr::select(measure, correlate, n, rho, p)
setting_checks$severity_correlations <- sev_cor

sev_groups <- dplyr::bind_rows(lapply(list(c("theta_SART_MS", "Mental speed"), c("theta_GRT", "Reasoning")), function(m) {
  x <- sev[[m[1]]]
  holm_tbl(dplyr::bind_rows(
    cmp_row(x, sev$moca_grp,  c("Controls", "Patients, MoCA unimpaired"),                m[2], "Controls vs. patients with unimpaired MoCA"),
    cmp_row(x, sev$moca_grp,  c("Controls", "Patients, MoCA impaired"),                  m[2], "Controls vs. patients with impaired MoCA"),
    cmp_row(x, sev$moca_grp,  c("Patients, MoCA unimpaired", "Patients, MoCA impaired"), m[2], "Unimpaired vs. impaired MoCA (patients only)"),
    cmp_row(x, sev$nihss_grp, c("Controls", "Patients, NIHSS 0-3"),                      m[2], "Controls vs. patients with NIHSS 0-3"),
    cmp_row(x, sev$nihss_grp, c("Controls", "Patients, NIHSS >= 4"),                     m[2], "Controls vs. patients with NIHSS >= 4"),
    cmp_row(x, sev$nihss_grp, c("Patients, NIHSS 0-3", "Patients, NIHSS >= 4"),          m[2], "NIHSS 0-3 vs. >= 4 (patients only)")))
}))
setting_checks$severity_groups <- sev_groups

nogo_acc <- sart_nogo %>% dplyr::filter(condition == 0L, !is.na(Correct)) %>%
  dplyr::group_by(SERIAL) %>% dplyr::summarise(n_nogo = dplyr::n(), nogo_accuracy = mean(Correct), .groups = "drop") %>%
  dplyr::filter(n_nogo >= 10) %>% dplyr::mutate(Group = assign_group(SERIAL))

pes_one <- function(d) {
  d <- d[order(d$Trial), ]
  err <- d$GoNoGo == 0L & d$Correct == 0L & !is.na(d$Correct)
  vgo <- d$GoNoGo == 1L & d$Correct == 1L & !is.na(d$RT) & d$RT >= sart_RT_MIN & d$RT <= sart_RT_MAX
  idx <- seq_len(nrow(d)); post_err <- integer(0)
  for (i in which(err)) {
    nxt <- idx[vgo & idx > i]; nxt_err <- idx[err & idx > i][1]
    if (!is.na(nxt_err)) nxt <- nxt[nxt < nxt_err]
    post_err <- c(post_err, utils::head(nxt, 3))
  }
  prev_correct <- c(FALSE, (d$Correct == 1L & !is.na(d$Correct))[-nrow(d)])
  post_corr <- setdiff(idx[vgo & prev_correct], post_err)
  if (sum(err) < 2 || length(post_err) < 2 || length(post_corr) < 10) return(NULL)
  tibble::tibble(n_errors = sum(err), n_post_error = length(post_err),
                 pes = mean(safe_log(d$RT[post_err])) - mean(safe_log(d$RT[post_corr])))
}
pes_scores <- sart_trials %>% dplyr::group_by(SERIAL) %>% dplyr::group_modify(~ { r <- pes_one(.x); if (is.null(r)) tibble::tibble() else r }) %>%
  dplyr::ungroup() %>% dplyr::mutate(Group = assign_group(SERIAL))
pes_within <- pes_scores %>% dplyr::group_by(Group) %>%
  dplyr::summarise(n = dplyr::n(), M_pes = mean(pes), SD_pes = stats::sd(pes),
                   t = if (dplyr::n() > 2) unname(stats::t.test(pes)$statistic) else NA_real_,
                   p = if (dplyr::n() > 2) stats::t.test(pes)$p.value else NA_real_, .groups = "drop")
setting_checks$pes_within <- pes_within

setting_checks$group_tests <- holm_tbl(dplyr::bind_rows(
  cmp_row(nogo_acc$nogo_accuracy, nogo_acc$Group, c("Control", "Patient"), "NoGo accuracy (trials without distractor)", "Controls vs. patients"),
  cmp_row(pes_scores$pes, pes_scores$Group, c("Control", "Patient"), "Post-error slowing (log RT difference)", "Controls vs. patients")))

for (nm in names(setting_checks)) {
  cat(sprintf("\n[Setting checks] %s:\n", nm)); print(setting_checks[[nm]], n = Inf, width = Inf)
  readr::write_csv(setting_checks[[nm]], file.path(output_dir, sprintf("SupplTable13_SettingChecks_%s.csv", nm)))
}
sart_severity_table <- setting_checks$severity_groups
sart_setting_table  <- setting_checks$group_tests

# --- 7.5 Signal detection analysis of the DRT (Supplementary Table 15) -------
# Sensitivity (d') and response criterion (c) separate discrimination from
# response bias. Hits: correct responses on target trials; false alarms:
# incorrect responses on lure trials. Extreme rates are corrected with the
# log-linear rule (add 0.5 to counts, 1 to trial numbers; Hautus, 1995).

drt_hits_Fam <- as.matrix(drt_Y_bin_rm[, drt_final_Fam, drop = FALSE])
drt_hits_Recoll <- as.matrix(drt_Y_bin_rm[, drt_final_Recoll, drop = FALSE])

drt_sdt <- tibble::tibble(
  SERIAL         = drt_X$SERIAL,
  n_Fam          = rowSums(!is.na(drt_hits_Fam)),
  n_hits         = rowSums(!is.na(drt_hits_Fam) & drt_hits_Fam == 1),
  n_Recoll       = rowSums(!is.na(drt_hits_Recoll)),
  n_false_alarms = rowSums(!is.na(drt_hits_Recoll) & drt_hits_Recoll == 0)
) %>%
  dplyr::filter(n_Fam > 0, n_Recoll > 0) %>%
  dplyr::mutate(
    HR_corrected  = (n_hits + 0.5) / (n_Fam + 1),
    FAR_corrected = (n_false_alarms + 0.5) / (n_Recoll + 1),
    z_HR          = stats::qnorm(HR_corrected),
    z_FAR         = stats::qnorm(FAR_corrected),
    d_prime       = z_HR - z_FAR,
    criterion_c   = -0.5 * (z_HR + z_FAR),
    Group         = assign_group(SERIAL)
  )

drt_gc_dprime    <- compare_groups(drt_sdt$d_prime, drt_sdt$Group)
drt_gc_criterion <- compare_groups(drt_sdt$criterion_c, drt_sdt$Group)

sdt_cohort <- dplyr::if_else(drt_sdt$Group == "Patient", "Patient", assign_control_phase(drt_sdt$SERIAL, drt_sdt$Group))
sdt_phase_table <- purrr::map_dfr(list(c("d_prime", "d' (sensitivity)"), c("criterion_c", "c (criterion)")), function(m) {
  purrr::map_dfr(list(c("Phase1", "Patient", "Phase-1 controls vs. patients"), c("Phase2", "Patient", "Phase-2 controls vs. patients"),
                      c("Phase1", "Phase2", "Phase-1 vs. phase-2 controls")), function(cmp) {
                        sel <- sdt_cohort %in% cmp[1:2]
                        gc <- compare_groups(drt_sdt[[m[1]]][sel], sdt_cohort[sel], group_levels = cmp[1:2])
                        if (is.null(gc)) return(tibble::tibble())
                        d1 <- gc$desc[gc$desc$Group == cmp[1], ]; d2 <- gc$desc[gc$desc$Group == cmp[2], ]
                        tibble::tibble(measure = m[2], Comparison = cmp[3], N_first = d1$N, M_first = d1$M, SD_first = d1$SD,
                                       N_second = d2$N, M_second = d2$M, SD_second = d2$SD, t = gc$t, df = gc$df, p = gc$p,
                                       g = gc$g, g_ci_lo = gc$g_ci[1], g_ci_hi = gc$g_ci[2])
                      })
})
if (nrow(sdt_phase_table) > 0) {
  sdt_phase_table$p_holm <- stats::p.adjust(sdt_phase_table$p, method = "holm")
  readr::write_csv(sdt_phase_table, file.path(output_dir, "SupplTable15_SignalDetection_ByPhase.csv"))
  cat("[DRT SDT] Sensitivity and criterion by control recruitment phase:\n"); print(sdt_phase_table, width = Inf)
}

drt_sdt_summary <- dplyr::bind_rows(
  with(drt_gc_dprime,    tibble::tibble(measure = "d' (sensitivity)", N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], M_ctrl = desc$M[desc$Group == "Control"], SD_ctrl = desc$SD[desc$Group == "Control"], M_pat = desc$M[desc$Group == "Patient"], SD_pat = desc$SD[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2])),
  with(drt_gc_criterion, tibble::tibble(measure = "c (criterion)",    N_ctrl = desc$N[desc$Group == "Control"], N_pat = desc$N[desc$Group == "Patient"], M_ctrl = desc$M[desc$Group == "Control"], SD_ctrl = desc$SD[desc$Group == "Control"], M_pat = desc$M[desc$Group == "Patient"], SD_pat = desc$SD[desc$Group == "Patient"], t = t, df = df, p = p, g = g, g_ci_lo = g_ci[1], g_ci_hi = g_ci[2]))
)
readr::write_csv(drt_sdt_summary, file.path(output_dir, "SupplTable15_SignalDetection.csv"))

cat("[DRT SDT] Signal detection reanalysis of Fam/Recoll (sensitivity vs. criterion):\n")
print(drt_sdt_summary)

# =============================================================================
# 8. CONVERGENT AND DIVERGENT VALIDITY (Figure 2; Supplementary Table 16)
# =============================================================================
# Figure 2: Pearson correlations between all scores and, in patients, the MoCA,
# for the total sample and each group; p values are Benjamini-Hochberg adjusted
# within the total sample and within each group. The diagonal shows split-half
# reliabilities. Supplementary Table 16 contrasts the residualised interference
# scores with non-residualised full-task scores: because baseline and
# interference performance are highly correlated, residualisation removes most
# of the reliable variance shared with other tasks.

fig_cor_labels <- c(
  theta_SART_EI       = "External interference control",
  theta_CI_resilience = "Internal interference control",
  theta_SART_MS       = "Mental speed",
  theta_GRT           = "Reasoning",
  theta_DRT_Fam       = "Familiarity-based recognition",
  theta_DRT_Recoll    = "Recollection-based discrimination",
  theta_PAL           = "Associative memory",
  MoCA                = "MoCA"
)

fig_cor_data <- master_theta %>%
  dplyr::select(SERIAL, Group, dplyr::all_of(names(fig_cor_labels))) %>%
  dplyr::mutate(Group = factor(dplyr::if_else(Group == "Control", "Controls", "Patients"), levels = c("Controls", "Patients")))
names(fig_cor_data)[match(names(fig_cor_labels), names(fig_cor_data))] <- unname(fig_cor_labels)

reliability_diag <- c(
  "External interference control"      = sart_reliability_sb,
  "Internal interference control"      = ci_reliability_sb,
  "Mental speed"                       = ms_reliability_sb,
  "Reasoning"                          = grt_reliability_sb,
  "Familiarity-based recognition"      = drt_reliability_sb,
  "Recollection-based discrimination"  = drt_Recoll_reliability_sb,
  "Associative memory"                 = pal_reliability_sb
)

fig_cor_vars <- unname(fig_cor_labels)
fig_cor_p_adj <- lapply(c(All = "All", Controls = "Controls", Patients = "Patients"), function(sub) {
  d <- if (sub == "All") fig_cor_data else fig_cor_data[fig_cor_data$Group == sub, ]
  pr <- utils::combn(fig_cor_vars, 2, simplify = FALSE)
  pv <- vapply(pr, function(v) {
    ok <- stats::complete.cases(d[[v[1]]], d[[v[2]]])
    if (sum(ok) < 3) NA_real_ else stats::cor.test(d[[v[1]]], d[[v[2]]])$p.value
  }, numeric(1))
  pa <- stats::p.adjust(pv, method = "BH")
  m <- matrix(NA_real_, length(fig_cor_vars), length(fig_cor_vars), dimnames = list(fig_cor_vars, fig_cor_vars))
  for (k in seq_along(pr)) { m[pr[[k]][1], pr[[k]][2]] <- pa[k]; m[pr[[k]][2], pr[[k]][1]] <- pa[k] }
  m
})

fig_cor_upper <- function(data, mapping, ...) {
  x_col <- rlang::as_name(mapping$x); y_col <- rlang::as_name(mapping$y)
  x_val <- data[[x_col]]; y_val <- data[[y_col]]; grp <- data$Group
  cor_lab <- function(x, y, sub) {
    if (sum(complete.cases(x, y)) < 3) return("NA")
    r_ <- stats::cor(x, y, use = "complete.obs"); p_ <- fig_cor_p_adj[[sub]][x_col, y_col]
    stars <- dplyr::case_when(is.na(p_) ~ "", p_ < .001 ~ "***", p_ < .01 ~ "**", p_ < .05 ~ "*", TRUE ~ "")
    sprintf("%.3f%s", r_, stars)
  }
  ggplot() +
    annotate("text", x = .5, y = .65, label = paste("All:", cor_lab(x_val, y_val, "All")), size = 4.5, fontface = "bold") +
    annotate("text", x = .5, y = .45, label = paste("Controls:", cor_lab(x_val[grp == "Controls"], y_val[grp == "Controls"], "Controls")), size = 4, color = ctrl_color) +
    annotate("text", x = .5, y = .25, label = paste("Patients:", cor_lab(x_val[grp == "Patients"], y_val[grp == "Patients"], "Patients")), size = 4, color = pat_color) +
    theme_void() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}
fig_cor_lower <- function(data, mapping, ...) {
  ggplot(data, mapping) + geom_point(alpha = .5, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .8) +
    scale_color_manual(values = c(Controls = ctrl_color, Patients = pat_color)) + theme_corr
}
fig_cor_diag <- function(data, mapping, ...) {
  x_col <- rlang::as_name(mapping$x)
  if (x_col == "MoCA") {
    ggplot(data, mapping) + geom_density(alpha = .6, aes(fill = Group)) +
      scale_fill_manual(values = c(Controls = ctrl_color, Patients = pat_color)) + theme_corr
  } else {
    rel_val <- reliability_diag[x_col]
    ggplot() + annotate("text", x = .5, y = .5, label = if (is.na(rel_val)) "" else sprintf("%.2f", rel_val), size = 6, fontface = "bold") +
      theme_void() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  }
}

fig_correlations <- GGally::ggpairs(
  data         = dplyr::select(fig_cor_data, -SERIAL),
  mapping      = aes(color = Group, fill = Group),
  columns      = unname(fig_cor_labels),
  upper        = list(continuous = fig_cor_upper),
  lower        = list(continuous = fig_cor_lower),
  diag         = list(continuous = fig_cor_diag),
  title        = "Person Score Correlations by Group (Final Latent Constructs)",
  columnLabels = stringr::str_wrap(unname(fig_cor_labels), width = 14)
) + theme_corr + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
                       axis.text.y = element_text(color = "black"))

ggsave(file.path(output_dir, "Figure2_CorrelationReliability.png"), fig_correlations, width = 12, height = 12, dpi = 300, bg = "transparent")

theta_only_cols <- setdiff(unname(fig_cor_labels), "MoCA")
cor_mat <- cor(dplyr::select(fig_cor_data, dplyr::all_of(theta_only_cols)), use = "pairwise.complete.obs")
readr::write_csv(
  tibble::as_tibble(cbind(variable = rownames(cor_mat), as.data.frame(cor_mat))),
  file.path(output_dir, "Figure2_Correlations.csv")
)

# --- Supplementary Table 16: full-task vs. residualised interference scores --

ci_fulltask_theta <- tryCatch({
  ci_mod_fulltask <- TAM::tam.mml(resp = ci_mat_trim, irtmodel = "1PL",
                                  control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
  ci_wle_fulltask <- TAM::tam.wle(ci_mod_fulltask, progress = FALSE)
  tibble::tibble(SERIAL = rownames(ci_mat_trim), theta_CI_full = ci_wle_fulltask$theta)
}, error = function(e) { cat(sprintf("[Suppl. Table 16] CI full-task model failed: %s\n", conditionMessage(e))); NULL })

sart_fulltask_theta <- tryCatch({
  full_rt <- dat_rt %>% dplyr::mutate(logRT_w = winsorize_iterative(logRT, k = 3.5)) %>%
    dplyr::group_by(SERIAL) %>% dplyr::summarise(full_mean_logRT = mean(logRT_w, na.rm = TRUE), .groups = "drop")
  full_acc <- dat_acc %>% dplyr::group_by(SERIAL) %>% dplyr::summarise(full_acc = mean(correct_go, na.rm = TRUE), .groups = "drop")
  full_rt %>% dplyr::full_join(full_acc, by = "SERIAL") %>%
    dplyr::mutate(z_speed = -1 * as.numeric(scale(full_mean_logRT)),
                  z_acc   = as.numeric(scale(full_acc)),
                  theta_SART_full = rowMeans(cbind(z_speed, z_acc), na.rm = TRUE)) %>%
    dplyr::filter(is.finite(theta_SART_full)) %>%
    dplyr::select(SERIAL, theta_SART_full)
}, error = function(e) { cat(sprintf("[Suppl. Table 16] SART full-task score failed: %s\n", conditionMessage(e))); NULL })

if (is.null(ci_fulltask_theta) || is.null(sart_fulltask_theta)) {
  cat("[Suppl. Table 16] Not computed -- see error(s) above.\n")
} else {
  fulltask_data <- master_theta %>%
    dplyr::select(SERIAL, theta_GRT, theta_DRT_Fam, theta_DRT_Recoll, theta_PAL, theta_SART_MS,
                  theta_SART_EI, theta_CI_resilience) %>%
    dplyr::left_join(ci_fulltask_theta, by = "SERIAL") %>%
    dplyr::left_join(sart_fulltask_theta, by = "SERIAL")
  
  fulltask_label_map <- c(
    theta_GRT = "Reasoning", theta_DRT_Fam = "Familiarity-based recognition",
    theta_DRT_Recoll = "Recollection-based discrimination", theta_PAL = "Associative memory",
    theta_SART_MS = "Mental speed",
    theta_SART_EI = "External interference (residualized)", theta_SART_full = "External interference (full task)",
    theta_CI_resilience = "Internal interference (residualized)", theta_CI_full = "Internal interference (full task)"
  )
  fulltask_vars <- names(fulltask_label_map)
  fulltask_mat <- as.matrix(dplyr::select(fulltask_data, dplyr::all_of(fulltask_vars)))
  colnames(fulltask_mat) <- unname(fulltask_label_map[fulltask_vars])
  
  fulltask_ct <- psych::corr.test(fulltask_mat, use = "pairwise", method = "pearson", adjust = "BH")
  fulltask_cor_table <- tibble::as_tibble(cbind(Variable = rownames(fulltask_ct$r), as.data.frame(round(fulltask_ct$r, 3))))
  readr::write_csv(fulltask_cor_table, file.path(output_dir, "SupplTable16_FullTaskVsResidualised.csv"))
  cat("[Suppl. Table 16] Full-task vs. residualized theta correlations:\n")
  print(fulltask_cor_table, n = Inf, width = Inf)
  fulltask_desc <- tibble::tibble(Variable = colnames(fulltask_mat), N = colSums(!is.na(fulltask_mat)),
                                  M = round(colMeans(fulltask_mat, na.rm = TRUE), 3), SD = round(apply(fulltask_mat, 2, stats::sd, na.rm = TRUE), 3))
  fulltask_p_bh <- fulltask_ct$p; fulltask_p_bh[lower.tri(fulltask_p_bh)] <- t(fulltask_ct$p)[lower.tri(fulltask_p_bh)]
  cat("[Suppl. Table 16] N, mean and SD per score:\n"); print(fulltask_desc, n = Inf, width = Inf)
  cat("[Suppl. Table 16] BH-adjusted p-values (symmetric):\n"); print(round(fulltask_p_bh, 4))
  readr::write_csv(fulltask_desc, file.path(output_dir, "SupplTable16_Descriptives.csv"))
  
  ability_labels <- c("Reasoning", "Familiarity-based recognition", "Recollection-based discrimination",
                      "Associative memory", "Mental speed")
  r <- fulltask_ct$r
  ei_resid_mean <- mean(abs(r["External interference (residualized)", ability_labels]))
  ei_full_mean  <- mean(abs(r["External interference (full task)", ability_labels]))
  ci_resid_mean <- mean(abs(r["Internal interference (residualized)", ability_labels]))
  ci_full_mean  <- mean(abs(r["Internal interference (full task)", ability_labels]))
  cat(sprintf("[Suppl. Table 16] Mean |r| with the other ability scores -- external interference: residualised %.2f, full task %.2f; internal interference: residualised %.2f, full task %.2f.\n",
              ei_resid_mean, ei_full_mean, ci_resid_mean, ci_full_mean))
}

# =============================================================================
# 9. EXPLORATORY CRITERION VALIDITY (Supplementary Table 17)
# =============================================================================
# Classification of MoCA-defined impairment (< 26/30) by each score in
# patients: AUC with a Mann-Whitney based test, and the cut-off maximising the
# Youden index. Cut-offs are optimised in the same sample, so the estimates are
# optimistic; the MoCA is a screening instrument, not a diagnostic reference.
compute_auc <- function(theta, impaired) {
  d <- tibble::tibble(theta = theta, impaired = impaired) %>% dplyr::filter(is.finite(theta), !is.na(impaired))
  n_impaired <- sum(d$impaired); n_unimpaired <- sum(!d$impaired)
  if (n_impaired < 5 || n_unimpaired < 5) {
    return(list(auc = NA_real_, p = NA_real_, n = nrow(d), n_impaired = n_impaired, n_unimpaired = n_unimpaired, status = "insufficient data (need >=5 per class)"))
  }
  wt <- tryCatch(stats::wilcox.test(theta ~ impaired, data = d), error = function(e) NULL)
  if (is.null(wt)) return(list(auc = NA_real_, p = NA_real_, n = nrow(d), n_impaired = n_impaired, n_unimpaired = n_unimpaired, status = "test failed"))
  
  auc <- unname(wt$statistic) / (n_unimpaired * n_impaired)
  list(auc = auc, p = wt$p.value, n = nrow(d), n_impaired = n_impaired, n_unimpaired = n_unimpaired, status = "ok")
}

youden_optimal <- function(theta, impaired) {
  d <- tibble::tibble(theta = theta, impaired = impaired) %>% dplyr::filter(is.finite(theta), !is.na(impaired))
  if (dplyr::n_distinct(d$impaired) < 2 || nrow(d) < 10) return(list(cutoff = NA_real_, sens = NA_real_, spec = NA_real_))
  best <- list(J = -Inf, cutoff = NA_real_, sens = NA_real_, spec = NA_real_)
  for (cut in sort(unique(d$theta))) {
    pred_impaired <- d$theta < cut
    sens <- mean(pred_impaired[d$impaired])
    spec <- mean(!pred_impaired[!d$impaired])
    J <- sens + spec - 1
    if (is.finite(J) && J > best$J) best <- list(J = J, cutoff = cut, sens = sens, spec = spec)
  }
  best
}

moca_impaired_vec <- master_theta$MoCA / 100 * 30 < moca_impairment_cutoff
cat(sprintf("[MoCA classification] %d people have usable MoCA data (cutoff: MoCA < %d/30); %d classified impaired, %d not. By group: %s\n",
            sum(!is.na(moca_impaired_vec)), moca_impairment_cutoff,
            sum(moca_impaired_vec, na.rm = TRUE), sum(!moca_impaired_vec, na.rm = TRUE),
            paste(names(table(master_theta$Group[!is.na(moca_impaired_vec)])), table(master_theta$Group[!is.na(moca_impaired_vec)]), sep = "=", collapse = "; ")))

moca_class_tasks <- list(
  list(theta = master_theta$theta_GRT,           label = "GRT (reasoning)"),
  list(theta = master_theta$theta_DRT_Fam,       label = "DRT-Fam (familiarity-based recognition)"),
  list(theta = master_theta$theta_DRT_Recoll,    label = "DRT-Recoll (recollection-based discrimination)"),
  list(theta = master_theta$theta_PAL,           label = "PAL (associative memory)"),
  list(theta = master_theta$theta_CI_resilience, label = "CI (internal interference control)"),
  list(theta = master_theta$theta_SART_EI,       label = "SART-ED (external interference control)"),
  list(theta = master_theta$theta_SART_MS,       label = "SART-ED (mental speed)")
)

moca_classification_results <- purrr::map_dfr(moca_class_tasks, function(t) {
  auc_res <- compute_auc(t$theta, moca_impaired_vec)
  yj <- youden_optimal(t$theta, moca_impaired_vec)
  tibble::tibble(
    Task = t$label, N = auc_res$n, N_impaired = auc_res$n_impaired, N_unimpaired = auc_res$n_unimpaired,
    AUC = auc_res$auc, p = auc_res$p,
    youden_cutoff_theta = yj$cutoff, sensitivity = yj$sens, specificity = yj$spec,
    status = auc_res$status
  )
})

moca_classification_results <- moca_classification_results %>%
  dplyr::mutate(note = "Exploratory criterion validity; cutoff optimised in-sample (estimates optimistic)")
readr::write_csv(moca_classification_results, file.path(output_dir, "SupplTable17_MoCAClassification.csv"))
cat("[MoCA classification] Exploratory criterion validity (AUC classifying MoCA-defined impairment) by domain:\n")
print(moca_classification_results)

# =============================================================================
# 10. STRUCTURAL VALIDITY: S-1 BIFACTOR MODEL (Figure 4; Supplementary Table 18)
# =============================================================================
# Indicators: z-standardised split-half parcels of each score.
# Model: S-1 bifactor model (Eid et al., 2017). Reasoning serves as the
# reference domain without a specific factor, so the general factor g is
# anchored in reasoning. Specific factors for external and internal
# interference control, mental speed, familiarity-based recognition and
# associative memory capture variance beyond g; their loadings are constrained
# equal within each pair of parcels. Specific factors are orthogonal to g and to
# each other, except for a covariance between the familiarity and associative
# memory factors (both figural memory tasks). Residual covariances between
# corresponding mental speed and external interference parcels account for
# their derivation from the same SART-ED trials. The residual variance of one
# associative memory parcel is fixed to zero (Heywood case).
# Estimation: robust maximum likelihood (MLR) with full-information maximum
# likelihood for missing data.
# Alternatives (Supplementary Table 18): single general factor; six correlated
# factors (equal loadings per parcel pair, as in the bifactor model); the
# bifactor model extended by recollection-based discrimination; complete-case
# refit. Figure 4 is drawn directly as SVG from the standardised solution.

df_master <- purrr::reduce(
  list(grt_parcels, drt_parcels, pal_parcels, sart_parcels, ci_parcels),
  dplyr::full_join, by = "SERIAL"
)

parcel_vars <- c("SART_A", "SART_B", "CI_A", "CI_B", "MS_A", "MS_B",
                 "GRT_A", "GRT_B", "DRT_Fam_A", "DRT_Fam_B", "PAL_A", "PAL_B")

df_cfa <- df_master
df_cfa[parcel_vars] <- scale(df_cfa[parcel_vars])

model_bifactor <- '
  g =~ SART_A + SART_B +
       CI_A + CI_B +
       MS_A + MS_B +
       GRT_A + GRT_B +
       DRT_Fam_A + DRT_Fam_B +
       PAL_A + PAL_B

  ExternalInterference =~ ei_s*SART_A    + ei_s*SART_B
  InternalInterference =~ ii_c*CI_A      + ii_c*CI_B
  MentalSpeed          =~ mp_m*MS_A      + mp_m*MS_B
  Familiarity          =~ fm_d*DRT_Fam_A + fm_d*DRT_Fam_B
  AssociativeMemory    =~ fs_p*PAL_A     + fs_p*PAL_B

  # Heywood-case resolution: PAL parcel B residual variance fixed to zero
  PAL_B ~~ 0*PAL_B

  # Local dependence between the two memory-specific factors (freed)
  Familiarity ~~ AssociativeMemory

  # All other between-factor covariances fixed to 0 (orthogonal bifactor structure)
  g ~~ 0*ExternalInterference
  g ~~ 0*InternalInterference
  g ~~ 0*MentalSpeed
  g ~~ 0*Familiarity
  g ~~ 0*AssociativeMemory
  ExternalInterference ~~ 0*InternalInterference
  ExternalInterference ~~ 0*MentalSpeed
  ExternalInterference ~~ 0*Familiarity
  ExternalInterference ~~ 0*AssociativeMemory
  InternalInterference ~~ 0*MentalSpeed
  InternalInterference ~~ 0*Familiarity
  InternalInterference ~~ 0*AssociativeMemory
  MentalSpeed ~~ 0*Familiarity
  MentalSpeed ~~ 0*AssociativeMemory

  # Shared method variance between SART- and MS-derived parcels (matched a/b)
  MS_A ~~ SART_A
  MS_B ~~ SART_B
'

fit_bifactor <- lavaan::cfa(model_bifactor, data = df_cfa, estimator = "MLR", missing = "fiml", std.lv = TRUE)

cat("\n[Figure 4] S-1 bifactor model -- fit indices\n")
fit_idx <- c("cfi.robust", "tli.robust", "rmsea.robust", "rmsea.ci.lower.robust",
             "rmsea.ci.upper.robust", "srmr", "chisq.scaled", "df", "pvalue.scaled")
print(round(lavaan::fitMeasures(fit_bifactor, fit_idx), 3))

fit_table <- as.data.frame(as.list(round(lavaan::fitMeasures(fit_bifactor, fit_idx), 3)))
readr::write_csv(fit_table, file.path(output_dir, "Figure4_BifactorFit.csv"))

std_loadings <- lavaan::standardizedSolution(fit_bifactor) %>% dplyr::filter(op == "=~")
readr::write_csv(std_loadings, file.path(output_dir, "Figure4_StandardisedLoadings.csv"))

model_onefactor <- '
  g =~ SART_A + SART_B + CI_A + CI_B + MS_A + MS_B + GRT_A + GRT_B + DRT_Fam_A + DRT_Fam_B + PAL_A + PAL_B
  PAL_B ~~ 0*PAL_B
  MS_A ~~ SART_A
  MS_B ~~ SART_B
'

model_correlated <- '
  Reasoning            =~ r_c*GRT_A     + r_c*GRT_B
  ExternalInterference =~ e_c*SART_A    + e_c*SART_B
  InternalInterference =~ i_c*CI_A      + i_c*CI_B
  MentalSpeed          =~ m_c*MS_A      + m_c*MS_B
  Familiarity          =~ f_c*DRT_Fam_A + f_c*DRT_Fam_B
  AssociativeMemory    =~ s_c*PAL_A     + s_c*PAL_B
  PAL_B ~~ 0*PAL_B
  MS_A ~~ SART_A
  MS_B ~~ SART_B
'
cfa_fit_row <- function(fit, label) {
  conv <- isTRUE(lavaan::lavInspect(fit, "converged"))
  if (!conv) {
    cat(sprintf("[Figure 4] '%s' did not converge -- fit indices not available.\n", label))
    return(tibble::tibble(model = label, N = lavaan::lavInspect(fit, "nobs"), converged = FALSE))
  }
  fm <- lavaan::fitMeasures(fit, c("chisq.scaled", "df", "pvalue.scaled", "cfi.robust", "tli.robust", "rmsea.robust", "srmr", "aic", "bic"))
  tibble::tibble(model = label, N = lavaan::lavInspect(fit, "nobs"), converged = TRUE, !!!as.list(round(fm, 3)))
}
fit_onefactor  <- tryCatch(lavaan::cfa(model_onefactor,  data = df_cfa, estimator = "MLR", missing = "fiml", std.lv = TRUE), error = function(e) NULL)
fit_correlated <- tryCatch(lavaan::cfa(model_correlated, data = df_cfa, estimator = "MLR", missing = "fiml", std.lv = TRUE), error = function(e) NULL)

df_cfa_recoll <- dplyr::left_join(df_cfa, drt_Recoll_parcels, by = "SERIAL")
df_cfa_recoll[c("DRT_Recoll_A", "DRT_Recoll_B")] <- scale(df_cfa_recoll[c("DRT_Recoll_A", "DRT_Recoll_B")])
model_bifactor_recoll <- paste0(model_bifactor, "
  g =~ DRT_Recoll_A + DRT_Recoll_B
  Recollection =~ rd_d*DRT_Recoll_A + rd_d*DRT_Recoll_B
  Recollection ~~ Familiarity
  g ~~ 0*Recollection
  ExternalInterference ~~ 0*Recollection
  InternalInterference ~~ 0*Recollection
  MentalSpeed ~~ 0*Recollection
  AssociativeMemory ~~ 0*Recollection
")
fit_recoll <- function(model) tryCatch(lavaan::cfa(model, data = df_cfa_recoll, estimator = "MLR", missing = "fiml", std.lv = TRUE),
                                       error = function(e) { cat(sprintf("[Figure 4] Sensitivity model failed: %s\n", conditionMessage(e))); NULL })
fit_bifactor_recoll <- fit_recoll(model_bifactor_recoll)
recoll_label <- "S-1 bifactor including recollection-based discrimination"
if (is.null(fit_bifactor_recoll) || !isTRUE(lavaan::lavInspect(fit_bifactor_recoll, "converged"))) {
  cat("[Figure 4] Sensitivity model did not converge -- refitting without the covariance between the two DRT specific factors.\n")
  model_bifactor_recoll <- sub("  Recollection ~~ Familiarity\n", "  Recollection ~~ 0*Familiarity\n", model_bifactor_recoll, fixed = TRUE)
  fit_bifactor_recoll <- fit_recoll(model_bifactor_recoll)
  recoll_label <- "S-1 bifactor including recollection-based discrimination (DRT factors orthogonal)"
}
recoll_converged <- !is.null(fit_bifactor_recoll) && isTRUE(lavaan::lavInspect(fit_bifactor_recoll, "converged"))
cat(sprintf("[Figure 4] Sensitivity model incl. recollection-based discrimination: %s.\n", if (recoll_converged) "converged" else "did NOT converge -- estimates not interpretable"))
if (recoll_converged) {
  ss_recoll <- lavaan::standardizedSolution(fit_bifactor_recoll)
  cat("[Figure 4] Sensitivity model incl. recollection-based discrimination -- standardised loadings of the Recoll parcels and factor covariance:\n")
  print(ss_recoll[(ss_recoll$op == "=~" & ss_recoll$rhs %in% c("DRT_Recoll_A", "DRT_Recoll_B")) |
                    (ss_recoll$op == "~~" & ss_recoll$lhs %in% c("Recollection", "Familiarity") & ss_recoll$rhs %in% c("Recollection", "Familiarity") & ss_recoll$lhs != ss_recoll$rhs),
                  c("lhs", "op", "rhs", "est.std", "se", "pvalue")])
  ld_rep <- lavaan::standardizedSolution(fit_bifactor) %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs, rhs, est_rep = est.std)
  ld_recoll  <- ss_recoll %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs, rhs, est_recoll = est.std)
  ld_cmp <- dplyr::inner_join(ld_rep, ld_recoll, by = c("lhs", "rhs"))
  cat(sprintf("[Figure 4] Largest change in the reported model's standardised loadings after adding recollection-based discrimination: %.3f\n",
              max(abs(ld_cmp$est_rep - ld_cmp$est_recoll), na.rm = TRUE)))
}
cfa_alt_table <- dplyr::bind_rows(
  cfa_fit_row(fit_bifactor, "S-1 bifactor (reported)"),
  if (!is.null(fit_bifactor_recoll)) cfa_fit_row(fit_bifactor_recoll, recoll_label) else NULL,
  if (!is.null(fit_onefactor))  cfa_fit_row(fit_onefactor,  "Single general factor") else NULL,
  if (!is.null(fit_correlated)) cfa_fit_row(fit_correlated, "Six correlated factors") else NULL
)
readr::write_csv(cfa_alt_table, file.path(output_dir, "SupplTable18_AlternativeModels.csv"))
cat("[Figure 4] Alternative structural models:\n"); print(cfa_alt_table, width = Inf)
if (!is.null(fit_onefactor)) {
  cfa_lrt <- tryCatch(lavaan::lavTestLRT(fit_onefactor, fit_bifactor), error = function(e) NULL)
  if (!is.null(cfa_lrt)) { cat("[Figure 4] Scaled chi-square difference test, single factor vs. S-1 bifactor:\n"); print(cfa_lrt) }
}

fit_bifactor_listwise <- tryCatch(lavaan::cfa(model_bifactor, data = df_cfa, estimator = "MLR", missing = "listwise", std.lv = TRUE),
                                  error = function(e) { cat(sprintf("[Figure 4] Listwise refit failed: %s\n", conditionMessage(e))); NULL })
if (!is.null(fit_bifactor_listwise)) {
  fiml_fit_table <- dplyr::bind_rows(cfa_fit_row(fit_bifactor, "FIML (reported)"), cfa_fit_row(fit_bifactor_listwise, "Complete cases (listwise)"))
  fiml_loadings <- std_loadings %>% dplyr::select(lhs, rhs, est_fiml = est.std) %>%
    dplyr::left_join(lavaan::standardizedSolution(fit_bifactor_listwise) %>% dplyr::filter(op == "=~") %>%
                       dplyr::select(lhs, rhs, est_listwise = est.std), by = c("lhs", "rhs")) %>%
    dplyr::mutate(difference = est_listwise - est_fiml)
  readr::write_csv(fiml_fit_table, file.path(output_dir, "SupplTable18_CompleteCases.csv"))
  readr::write_csv(fiml_loadings, file.path(output_dir, "SupplTable18_CompleteCases_Loadings.csv"))
  cat("[Figure 4] FIML vs. complete-case fit:\n"); print(fiml_fit_table, width = Inf)
  cat(sprintf("[Figure 4] Largest absolute loading difference (listwise - FIML): %.3f\n", max(abs(fiml_loadings$difference), na.rm = TRUE)))
}

# --- Figure 4: path diagram ---------------------------------------------------
# Drawn as SVG with fixed coordinates: g on top, parcels in the middle,
# specific factors at the bottom. Labels are standardised estimates.

std_diag <- lavaan::standardizedSolution(fit_bifactor)

get_loading <- function(lhs, rhs) {
  v <- std_diag$est.std[std_diag$lhs == lhs & std_diag$op == "=~" & std_diag$rhs == rhs]
  if (!length(v) || is.na(v[1])) return("")
  sprintf("%.2f", v[1])
}

get_cov <- function(l, r) {
  v <- std_diag$est.std[(std_diag$lhs == l & std_diag$op == "~~" & std_diag$rhs == r) | (std_diag$lhs == r & std_diag$op == "~~" & std_diag$rhs == l)]
  if (!length(v) || is.na(v[1])) return("")
  sprintf("%.2f", v[1])
}

COL_BLUE_FILL   <- "#6695c5";  COL_BLUE_STROKE  <- "#2e5a9c"
COL_RED_FILL    <- "#eca594";  COL_RED_STROKE   <- "#c95f3f"
COL_GREEN_FILL  <- "#d5f5e3";  COL_GREEN_STROKE <- "#3a9e6f"
COL_TEXT        <- "#1a1a2e"

W       <- 1700L;  H       <- 580L
MAN_W   <- 110L;   MAN_H   <- 44L
LAT_RX  <- 70L;    LAT_RY  <- 36L
G_RX    <- 66L;    G_RY    <- 36L
Y_G     <- 68;     Y_MAN   <- 265;   Y_LAT   <- 480
ARROW_GAP   <- 5;  MARKER_SIZE <- 7L

man_ids <- c("SART_A","SART_B","CI_A","CI_B", "MS_A", "MS_B", "GRT_A","GRT_B","DRT_Fam_A","DRT_Fam_B","PAL_A","PAL_B")
man_labels <- c("EIC a","EIC b","IIC a","IIC b", "MS a", "MS b", "R a","R b","Rec a","Rec b","AM a","AM b")

n_man  <- length(man_ids)
x_pad  <- 70
x_gap  <- (W - 2 * x_pad) / (n_man - 1)
man_x  <- x_pad + (0:(n_man - 1)) * x_gap
g_x    <- mean(man_x);  g_y <- Y_G

man_fill <- c(SART_A=COL_RED_FILL, SART_B=COL_RED_FILL, CI_A=COL_RED_FILL, CI_B=COL_RED_FILL, MS_A=COL_BLUE_FILL, MS_B=COL_BLUE_FILL, GRT_A=COL_BLUE_FILL, GRT_B=COL_BLUE_FILL, DRT_Fam_A=COL_GREEN_FILL, DRT_Fam_B=COL_GREEN_FILL, PAL_A=COL_GREEN_FILL, PAL_B=COL_GREEN_FILL)
man_stroke <- c(SART_A=COL_RED_STROKE, SART_B=COL_RED_STROKE, CI_A=COL_RED_STROKE, CI_B=COL_RED_STROKE, MS_A=COL_BLUE_STROKE, MS_B=COL_BLUE_STROKE, GRT_A=COL_BLUE_STROKE, GRT_B=COL_BLUE_STROKE, DRT_Fam_A=COL_GREEN_STROKE, DRT_Fam_B=COL_GREEN_STROKE, PAL_A=COL_GREEN_STROKE, PAL_B=COL_GREEN_STROKE)
man_tcol <- setNames(rep(COL_TEXT, n_man), man_ids)

lat_info <- list(
  EI  = list(label="External\ninterference\ncontrol", lav="ExternalInterference", ids=c("SART_A","SART_B"), fill=COL_RED_FILL, stroke=COL_RED_STROKE, text=COL_TEXT),
  II  = list(label="Internal\ninterference\ncontrol",  lav="InternalInterference", ids=c("CI_A","CI_B"), fill=COL_RED_FILL, stroke=COL_RED_STROKE, text=COL_TEXT),
  MP  = list(label="Mental\nspeed", lav="MentalSpeed", ids=c("MS_A","MS_B"), fill=COL_BLUE_FILL, stroke=COL_BLUE_STROKE, text=COL_TEXT),
  FM  = list(label="Rec",    lav="Familiarity", ids=c("DRT_Fam_A","DRT_Fam_B"), fill=COL_GREEN_FILL, stroke=COL_GREEN_STROKE, text=COL_TEXT),
  FSM = list(label="Associative\nmemory", lav="AssociativeMemory", ids=c("PAL_A","PAL_B"), fill=COL_GREEN_FILL, stroke=COL_GREEN_STROKE, text=COL_TEXT)
)
for (nm in names(lat_info)) {
  idx <- match(lat_info[[nm]]$ids, man_ids)
  lat_info[[nm]]$x <- mean(man_x[idx])
  lat_info[[nm]]$y <- Y_LAT
}

ellipse_pt <- function(cx, cy, rx, ry, px, py, gap = 0) {
  angle <- atan2((py - cy) * rx, (px - cx) * ry)
  c(cx + (rx + gap) * cos(angle), cy + (ry + gap) * sin(angle))
}

svg_arrow <- function(x1, y1, x2, y2, col, width=1.6, marker="arrow", label=NULL, lcol=COL_TEXT, fsize=12, lx=NULL, ly=NULL) {
  s <- sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.1f" marker-end="url(#%s)"/>', x1, y1, x2, y2, col, width, marker)
  if (!is.null(label) && nchar(label) > 0) {
    plx <- if (!is.null(lx)) lx else (x1 + x2) / 2
    ply <- if (!is.null(ly)) ly else (y1 + y2) / 2
    s <- paste0(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" dy="0.35em" font-family="Arial" font-size="%d" fill="%s">%s</text>', plx, ply, fsize, lcol, label))
  }
  s
}

svg_cov_arc <- function(x1, y1, rx1, x2, y2, rx2, col, label=NULL, fsize=13, y_offset=16) {
  sx  <- x1 + rx1 + ARROW_GAP
  ex  <- x2 - rx2 - ARROW_GAP
  cxc <- (sx + ex) / 2
  cyc <- y1 + 55
  s <- sprintf('<path d="M %.1f %.1f Q %.1f %.1f %.1f %.1f" fill="none" stroke="%s" stroke-width="1.9" stroke-dasharray="6,4" marker-start="url(#arrowBiStart)" marker-end="url(#arrowBiEnd)"/>', sx, y1, cxc, cyc, ex, y2, col)
  if (!is.null(label) && nchar(label) > 0) {
    apex_y <- (y1 + cyc) / 2
    s <- paste0(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-family="Arial" font-size="%d" fill="%s">%s</text>', cxc, apex_y + y_offset, fsize, COL_TEXT, label))
  }
  s
}

svg_resid_arc <- function(x1, x2, y, col, label="", fsize=13, y_offset=16) {
  cy <- y + MAN_H/2 + ARROW_GAP
  cxc <- (x1 + x2) / 2
  cyc <- cy + 65 + abs(x1 - x2)*0.08
  s <- sprintf('<path d="M %.1f %.1f Q %.1f %.1f %.1f %.1f" fill="none" stroke="%s" stroke-width="1.5" stroke-dasharray="5,4" marker-start="url(#arrowBiStart)" marker-end="url(#arrowBiEnd)"/>', x1, cy, cxc, cyc, x2, cy, col)
  if (label != "") {
    apex_y <- (cy + cyc) / 2
    s <- paste0(s, sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" font-family="Arial" font-size="%d" fill="%s">%s</text>', cxc, apex_y + y_offset, fsize, COL_TEXT, label))
  }
  s
}

mk <- function(id, col, reverse=FALSE) {
  if (!reverse) sprintf('<marker id="%s" markerWidth="%d" markerHeight="%d" refX="%d" refY="3.5" orient="auto"><path d="M0,0 L0,7 L%d,3.5 z" fill="%s"/></marker>', id, MARKER_SIZE, MARKER_SIZE, MARKER_SIZE, MARKER_SIZE, col)
  else sprintf('<marker id="%s" markerWidth="%d" markerHeight="%d" refX="0" refY="3.5" orient="auto"><path d="M%d,0 L%d,7 L0,3.5 z" fill="%s"/></marker>', id, MARKER_SIZE, MARKER_SIZE, MARKER_SIZE, MARKER_SIZE, col)
}

els <- character(0)
add <- function(x) els <<- c(els, x)

add(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', W, H, W, H))
add("<defs>")
add(mk("arrowBlue",    COL_BLUE_STROKE))
add(mk("arrowRed",     COL_RED_STROKE))
add(mk("arrowGreen",   COL_GREEN_STROKE))
add(mk("arrowBiEnd",   COL_GREEN_STROKE))
add(mk("arrowBiStart", COL_GREEN_STROKE, reverse=TRUE))
add("</defs>")
add(sprintf('<rect width="%d" height="%d" fill="white"/>', W, H))

PERP_BASE_G <- 45;  PERP_MAX_G <- 200
Y_LABEL_G   <- (Y_G + G_RY + Y_MAN - MAN_H / 2) / 2 + 18
x_distances <- abs(man_x - g_x)
x_max       <- max(x_distances)

for (i in seq_along(man_ids)) {
  mid   <- man_ids[i]
  src   <- ellipse_pt(g_x, g_y, G_RX, G_RY, man_x[i], Y_MAN, ARROW_GAP)
  dst_x <- man_x[i];  dst_y <- Y_MAN - MAN_H / 2 - ARROW_GAP
  mx    <- (src[1] + dst_x) / 2
  t      <- x_distances[i] / x_max
  offset <- PERP_BASE_G * (PERP_MAX_G / PERP_BASE_G)^t
  h_sign <- if (dst_x >= g_x) 1 else -1
  add(svg_arrow(src[1], src[2], dst_x, dst_y, col=COL_BLUE_STROKE, width=1.5, marker="arrowBlue", label=get_loading("g", mid), lcol=COL_TEXT, fsize=22, lx=mx + h_sign * offset, ly=Y_LABEL_G))
}

Y_LABEL_SF     <- (Y_MAN + MAN_H / 2 + Y_LAT - LAT_RY) / 2 + 30
PERP_OFFSET_SF <- 30
spec_col <- c(EI=COL_RED_STROKE, II=COL_RED_STROKE, MP=COL_BLUE_STROKE, FM=COL_GREEN_STROKE, FSM=COL_GREEN_STROKE)
spec_mk  <- c(EI="arrowRed", II="arrowRed", MP="arrowBlue", FM="arrowGreen", FSM="arrowGreen")

for (nm in names(lat_info)) {
  info <- lat_info[[nm]]
  lx   <- info$x;  ly <- info$y
  for (mid in info$ids) {
    i      <- match(mid, man_ids)
    src    <- ellipse_pt(lx, ly, LAT_RX, LAT_RY, man_x[i], Y_MAN, ARROW_GAP)
    dst_x  <- man_x[i];  dst_y <- Y_MAN + MAN_H / 2 + ARROW_GAP
    mx     <- (src[1] + dst_x) / 2
    h_sign <- if (dst_x >= lx) 1 else -1
    add(svg_arrow(src[1], src[2], dst_x, dst_y, col=spec_col[nm], width=1.9, marker=spec_mk[nm], label=get_loading(info$lav, mid), lcol=COL_TEXT, fsize=24, lx=mx + h_sign * PERP_OFFSET_SF, ly=Y_LABEL_SF))
  }
}

add(svg_cov_arc(lat_info$FM$x,  lat_info$FM$y,  LAT_RX, lat_info$FSM$x, lat_info$FSM$y, LAT_RX, col=COL_GREEN_STROKE, label=get_cov("Familiarity","AssociativeMemory")))

add(svg_resid_arc(man_x[match("SART_A", man_ids)], man_x[match("MS_A", man_ids)], Y_MAN, col="#7F8C8D", label=get_cov("MS_A", "SART_A")))
add(svg_resid_arc(man_x[match("SART_B", man_ids)], man_x[match("MS_B", man_ids)], Y_MAN, col="#7F8C8D", label=get_cov("MS_B", "SART_B")))

for (i in seq_along(man_ids)) {
  mid <- man_ids[i]
  cx  <- man_x[i];  cy <- Y_MAN
  x0  <- cx - MAN_W / 2;  y0 <- cy - MAN_H / 2
  add(sprintf('<rect x="%.1f" y="%.1f" width="%d" height="%d" rx="6" ry="6" fill="%s" stroke="%s" stroke-width="1.8"/>', x0, y0, MAN_W, MAN_H, man_fill[mid], man_stroke[mid]))
  add(sprintf('<text x="%.1f" y="%.1f" text-anchor="middle" dy="0.35em" font-family="Arial" font-size="15" font-weight="600" fill="%s">%s</text>', cx, cy, man_tcol[mid], man_labels[i]))
}

add(sprintf('<ellipse cx="%.1f" cy="%d" rx="%d" ry="%d" fill="%s" stroke="%s" stroke-width="2.5"/>', g_x, Y_G, G_RX, G_RY, COL_BLUE_FILL, COL_BLUE_STROKE))
add(sprintf('<text x="%.1f" y="%d" text-anchor="middle" dominant-baseline="central" font-family="Arial" font-size="20" font-weight="bold" fill="%s">g</text>', g_x, Y_G, COL_TEXT))

for (nm in names(lat_info)) {
  info  <- lat_info[[nm]]
  cx    <- info$x
  parts <- strsplit(info$label, "\n")[[1]]
  add(sprintf('<ellipse cx="%.1f" cy="%d" rx="%d" ry="%d" fill="%s" stroke="%s" stroke-width="2.0"/>', cx, Y_LAT, LAT_RX, LAT_RY, info$fill, info$stroke))
  line_h <- 15
  add(sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-family="Arial" font-size="12" font-weight="bold" fill="%s">%s</text>', cx, Y_LAT, info$text, paste0(sprintf('<tspan x="%.1f" dy="%.1f">%s</tspan>', cx, c(-(length(parts) - 1) * line_h / 2, rep(line_h, length(parts) - 1)), parts), collapse = "")))
}

add("</svg>")

svg_string <- paste(els, collapse = "\n")
svg_path <- file.path(output_dir, "Figure4_Bifactor_PathDiagram.svg")
writeLines(svg_string, svg_path)
cat(sprintf("[Figure 4] Path diagram SVG written: %s\n", svg_path))

png_path <- file.path(output_dir, "Figure4_Bifactor_PathDiagram.png")
tryCatch({
  rsvg::rsvg_png(charToRaw(svg_string), file = png_path, width = 2800)
  cat(sprintf("[Figure 4] Path diagram PNG written: %s\n", png_path))
}, error = function(e) {
  cat(sprintf("[Figure 4] Could not render PNG from SVG (%s); the SVG file is unaffected.\n", conditionMessage(e)))
})

# =============================================================================
# 11. EXPORT OF REPORTED TABLES
# =============================================================================
# All reported tables are compiled into one HTML file (readable in Word or a
# browser) in the order of the manuscript. Tables not computed in a run are
# listed as skipped.

html_tables <- character(0)
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE)
}
table_to_html <- function(df, caption) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  cells <- vapply(df, function(col) { v <- as.character(col); v[is.na(v)] <- ""; html_escape(v) }, character(nrow(df)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = nrow(df))
  header <- paste0("<tr>", paste0("<th>", html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- paste0("<tr>", apply(cells, 1, function(r) paste0("<td>", r, "</td>", collapse = "")), "</tr>", collapse = "\n")
  paste0("<h2>", html_escape(caption), "</h2>\n<table>\n", header, "\n", body, "\n</table>\n")
}

format_for_report <- function(df, digits = 3) {
  p_cols <- grep("^p$|^p_|pvalue", names(df), ignore.case = TRUE, value = TRUE)
  df %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric) & !dplyr::any_of(p_cols), ~ round(., digits))) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(p_cols), ~ ifelse(. < 0.001, "<.001", sprintf(paste0("%.", digits, "f"), .))))
}

report_table <- function(df_name, caption, digits = 3) {
  if (!exists(df_name, inherits = TRUE)) {
    cat(sprintf("[Reporting] SKIPPED (not computed in this run): %s\n", caption))
    return(invisible(NULL))
  }
  df <- get(df_name, inherits = TRUE)
  
  if (is.data.frame(df) && !tibble::is_tibble(df)) df <- tibble::as_tibble(df)
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    cat(sprintf("[Reporting] SKIPPED (empty in this run): %s\n", caption))
    return(invisible(NULL))
  }
  df_fmt <- format_for_report(df, digits)
  cat("\n============================================================\n")
  cat(caption, "\n")
  cat("============================================================\n")
  print(df_fmt, n = Inf, width = Inf)
  html_tables <<- c(html_tables, table_to_html(df_fmt, caption))
}

report_table("item_allocation_table",       "Table 1: item sets answered by patients (planned-missingness allocation)")
report_table("participant_characteristics", "Supplementary Table 1: participant characteristics")
report_table("stroke_characteristics",      "Supplementary Table 1: clinical characteristics of patients")
report_table("sample_size_table",           "Supplementary Table 2: analysed sample sizes by task and group")
report_table("reliability_table",           "Supplementary Table 2: split-half reliability (Spearman-Brown corrected; bootstrap 95% CI)")
report_table("sem_table",                   "Supplementary Table 2: cross-sectional standard error of measurement")
report_table("completer_results",           "Supplementary Table 3: completers vs. non-completers (analysed patients)")
report_table("item_fit_full",               "Supplementary Table 4: item difficulty and fit before and after trimming")
report_table("model_comparison_table",      "Supplementary Table 5: IRT model comparisons")
report_table("rasch_vs_2pl_table",          "Supplementary Table 5: Rasch vs. 2PL -- practical impact")
report_table("targeting_table",             "Supplementary Table 6: targeting, test information and conditional SE")
report_table("ld_summary_table",            "Supplementary Table 6: local dependence (Yen's Q3)")
report_table("ic_summary",                  "Supplementary Table 6: interchangeability of random disjoint item halves")
report_table("dif_summary",                 "Supplementary Table 7: differential item functioning")
report_table("dif_sensitivity_table",       "Supplementary Table 8: DIF sensitivity re-scoring")
report_table("effect_size_table",           "Supplementary Table 9: known-groups comparisons")
report_table("latent_mean_table",           "Supplementary Table 10: latent mean differences (multiple-group IRT)")
report_table("sart_group_interaction",      "Supplementary Table 11: SART-ED condition x group interaction")
report_table("joint_covariate_table",       "Supplementary Table 12: joint covariate models")
report_table("sart_severity_table",         "Supplementary Table 13: severity gradient within patients")
report_table("sart_setting_table",          "Supplementary Table 13: setting-independent SART-ED measures")
report_table("phase_table",                 "Supplementary Table 14: comparisons by control recruitment phase")
report_table("drt_common_items_table",      "Supplementary Table 14: DRT phase comparisons on common items")
report_table("drt_sdt_summary",             "Supplementary Table 15: DRT signal detection, controls vs. patients")
report_table("sdt_phase_table",             "Supplementary Table 15: DRT signal detection by control recruitment phase")
report_table("fulltask_cor_table",          "Supplementary Table 16: full-task vs. residualised interference scores")
report_table("moca_classification_results", "Supplementary Table 17: exploratory criterion validity (MoCA)")
report_table("cfa_alt_table",               "Supplementary Table 18: alternative structural models")
report_table("fiml_fit_table",              "Supplementary Table 18: complete-case refit of the S-1 bifactor model")
report_table("fit_table",                   "Figure 4: S-1 bifactor model fit")
report_table("std_loadings",                "Figure 4: standardised loadings")
if (length(html_tables) > 0) {
  html_doc <- c(
    "<!DOCTYPE html>", "<html><head><meta charset=\"utf-8\"><title>ORPheoS reporting tables</title>",
    "<style>body{font-family:Arial,sans-serif;font-size:10pt;margin:2em;} h2{font-size:12pt;margin-top:2em;}",
    "table{border-collapse:collapse;margin-bottom:1em;} th,td{border:1px solid #999;padding:3px 6px;text-align:left;}",
    "th{background:#eee;}</style></head><body>",
    "<h1>ORPheoS / EMA4Stroke -- reporting tables</h1>",
    sprintf("<p>Generated %s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M")),
    html_tables, "</body></html>")
  html_out_path <- file.path(output_dir, "Manuscript_Reporting_Tables.html")
  tryCatch({
    con <- file(html_out_path, open = "w", encoding = "UTF-8"); writeLines(html_doc, con); close(con)
    cat(sprintf("\n[Reporting] Tables written: %s (%d tables; open in Word or a browser)\n", html_out_path, length(html_tables)))
  }, error = function(e) cat(sprintf("[Reporting] HTML document failed to save: %s\n", conditionMessage(e))))
} else {
  cat("[Reporting] No tables were available to compile -- check the SKIPPED messages above.\n")
}
