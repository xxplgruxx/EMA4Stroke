# =============================================================================
# ORPheoS / EMA4Stroke — Reproducible Analysis Script
# Grupe, Hildebrandt, Witt, Kastrup, Thiel, Roheger & Hildebrandt
# Beyond the Clinic: Bridging the Gap in Post-Stroke Cognitive Monitoring
# with the Oldenburg Test Battery for Remote Digital Phenotyping of
# Post-Stroke Cognition (ORPheoS)
#
# This script reproduces every analysis, figure, and in-text statistic
# reported in the main text of the manuscript:
#   - Table 1 completion / sample statistics
#   - Item-level IRT calibration and final item pools per task
#     (GRT, DRT, PAL, CI, SART-ED)
#   - Figure 2  Latent ability (theta) by group, all domains
#   - Figure 3  Inter-task correlation matrix + split-half reliability
#   - Figure 4  S-1 bifactor CFA of the ORPheoS taxonomy
#   - Known-groups validity statistics (Hedges' g, t-tests) per domain
#
# An optional appendix section (education & sex sensitivity analyses,
# Supplementary Tables 5-6) can be toggled off below.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. CONFIGURATION
# -----------------------------------------------------------------------------
# Set this to the folder containing Final_Merged_Dataset.feather and
# EMA4Stroke_Baseline_GoNoGo_Raw.csv. All outputs (figures, tables) are
# written to output_dir.

data_dir   <- "\\\\smb.uni-oldenburg.de\\psychologie$\\PMuS\\GraKo\\Grupe\\GRAKO Stroke & Neuromodulation\\Analyse\\Test Sets"
output_dir <- "./output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Toggle the education/sex sensitivity appendix (Supplementary Tables 5-6).
# Set to FALSE to run only the main-text analyses.
RUN_APPENDIX_SENSITIVITY <- TRUE


# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------
packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble", "glue",
  "arrow", "ggplot2", "gghalves", "GGally", "patchwork", "scales",
  "mirt", "TAM", "psych", "GPArotation", "Gifi", "lavaan",
  "lme4", "effsize", "Matrix", "rlang"
)
new_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

set.seed(20260821)

# -----------------------------------------------------------------------------
# 2. SHARED HELPERS
# -----------------------------------------------------------------------------
apa_colors <- c(hc = "#004E9F", pat = "#E06A4E")
hc_color  <- unname(apa_colors["hc"])
pat_color <- unname(apa_colors["pat"])

theme_apa <- theme_minimal(base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.minor = element_blank()
  )

# Figures 2 and 3 come from a separate original script with its own theme
# (fully blanked gridlines, transparent backgrounds) distinct from the
# per-task theme above.
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

fmt_num <- function(x, d = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmt_p   <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "< .001", paste0("= ", sub("^0\\.", ".", sprintf("%.3f", p)))))

assign_group <- function(serial) {
  dplyr::case_when(
    grepl("^[A-Za-z]{2}\\d{2}$", serial) ~ "Patient",
    !grepl("test|dummy", tolower(serial)) ~ "Healthy",
    TRUE ~ NA_character_
  )
}

# 1PL Rasch theta (EAP) for a binary item matrix, returning NA for rows with
# no data. Shared across GRT / DRT / CI / CFA parcel scoring.
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

# Iterative MNSQ item-fit trimming (TAM), threshold 0.5-1.5 per manuscript
# Section 2.5. Drops the single worst-fitting item at each step until all
# remaining items satisfy both Infit and Outfit criteria.
iterative_mnsq_trim <- function(item_mat, mnsq_lower = 0.5, mnsq_upper = 1.5, min_items = 3) {
  current_items <- colnames(item_mat)
  dropped <- character(0)
  repeat {
    if (length(current_items) < min_items) {
      warning("Iterative MNSQ trimming stopped early (< min_items remaining).")
      break
    }
    sub <- item_mat[, current_items, drop = FALSE]
    mod <- TAM::tam.mml(as.matrix(sub), irtmodel = "1PL", control = list(snodes = 1000, qmc = TRUE), verbose = FALSE)
    fit <- as.data.frame(TAM::msq.itemfit(mod)$itemfit)
    infit_col  <- grep("^Infit$",  names(fit), value = TRUE)[1]
    outfit_col <- grep("^Outfit$", names(fit), value = TRUE)[1]
    item_col   <- grep("^item$|^parameter$", names(fit), value = TRUE)[1]
    fit$misfit <- fit[[infit_col]] > mnsq_upper | fit[[infit_col]] < mnsq_lower |
      fit[[outfit_col]] > mnsq_upper | fit[[outfit_col]] < mnsq_lower
    misfits <- fit[fit$misfit, ]
    if (nrow(misfits) == 0) break
    misfits$dev <- pmax(abs(misfits[[infit_col]] - 1), abs(misfits[[outfit_col]] - 1))
    worst <- misfits[[item_col]][which.max(misfits$dev)]
    dropped <- c(dropped, worst)
    current_items <- setdiff(current_items, worst)
  }
  list(keep = current_items, dropped = dropped)
}

# Group comparison: descriptives, Welch t-test, Hedges' g with 95% CI.
compare_groups <- function(theta, group, group_levels = c("Healthy", "Patient")) {
  d <- tibble::tibble(theta = theta, Group = factor(group, levels = group_levels)) %>%
    dplyr::filter(is.finite(theta), !is.na(Group))
  if (dplyr::n_distinct(d$Group) < 2) return(NULL)
  desc <- d %>% dplyr::group_by(Group) %>% dplyr::summarise(N = dplyr::n(), M = mean(theta), SD = sd(theta), .groups = "drop")
  tt <- t.test(theta ~ Group, data = d)
  gg <- effsize::cohen.d(theta ~ Group, data = d, hedges.correction = TRUE)
  list(desc = desc, t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
       g = unname(gg$estimate), g_ci = gg$conf.int)
}

# Spearman-Brown corrected split-half reliability from two parcel scores.
spearman_brown <- function(x, y) {
  r <- suppressWarnings(cor(x, y, use = "complete.obs"))
  (2 * r) / (1 + r)
}

# =============================================================================
# 3. DATA LOADING
# =============================================================================
# Group is retained as NA for test/dummy SERIALs here; whether those rows are
# excluded from a given task's calibration or reporting differs by task
# below and is handled locally in each section.
df <- arrow::read_feather(file.path(data_dir, "Final_Merged_Dataset.feather"))

df <- df %>% dplyr::mutate(Group = assign_group(SERIAL))

# Group label assignment for controls and patients based on ID.
assign_group_simple <- function(serial) dplyr::if_else(grepl("^[A-Za-z]{2}\\d{2}$", serial), "Patient", "Healthy")

get_union_items <- function(d1, d2, cols) {
  i1 <- cols[colMeans(!is.na(d1[, cols, drop = FALSE])) > 0]
  i2 <- cols[colMeans(!is.na(d2[, cols, drop = FALSE])) > 0]
  union(i1, i2)
}

# =============================================================================
# 4. GRT — GEOMETRIC REASONING TASK
# =============================================================================
# Pipeline: variance screen -> tetrachoric PCA (dimensionality check) ->
# manual removal of items with secondary loadings -> iterative MNSQ trimming
# (0.5-1.5) -> final 1PL Rasch model (1PL vs 2PL compared, 1PL retained on
# BIC parsimony grounds) -> EAP theta scores.

grt_cols_raw <- grep("^GRT", names(df), value = TRUE)
df[grt_cols_raw] <- lapply(df[grt_cols_raw], function(v) suppressWarnings(as.numeric(as.character(v))))

is_binary_col <- function(x) { ux <- unique(na.omit(x)); length(ux) > 0 && all(ux %in% c(0, 1)) }
grt_item_cols <- grt_cols_raw[vapply(df[grt_cols_raw], is_binary_col, logical(1))]

grt_df       <- dplyr::filter(df, !is.na(Group))
grt_patients <- dplyr::filter(grt_df, Group == "Patient")
grt_healthy  <- dplyr::filter(grt_df, Group == "Healthy")
grt_common   <- get_union_items(grt_patients, grt_healthy, grt_item_cols)

grt_combined <- dplyr::bind_rows(
  dplyr::select(grt_patients, SERIAL, dplyr::all_of(grt_common)),
  dplyr::select(grt_healthy,  SERIAL, dplyr::all_of(grt_common))
) %>% dplyr::arrange(SERIAL)

grt_item_var <- sapply(grt_combined[, grt_common, drop = FALSE], var, na.rm = TRUE)
grt_keep_var <- names(grt_item_var[grt_item_var > 0])
grt_X <- dplyr::select(grt_combined, SERIAL, dplyr::all_of(grt_keep_var))
grt_Y <- dplyr::select(grt_X, -SERIAL)

# --- Dimensionality: tetrachoric PCA ---
grt_tet <- suppressWarnings(psych::tetrachoric(grt_Y, na.rm = TRUE))$rho
grt_pca <- psych::pca(grt_tet, nfactors = 2, rotate = "none")
grt_eig <- grt_pca$values
grt_var_exp1 <- grt_pca$Vaccounted["Proportion Var", 1]

# --- Manual exclusion of items with secondary loadings / heterogeneous discrimination ---
# (per manuscript Section 3.1: k=4 secondary loadings, k=6 heterogeneous discrimination)
grt_manual_drop <- c("GRT30", "GRT10", "GRT24", "GRT13")
grt_stage2 <- setdiff(colnames(grt_Y), grt_manual_drop)
grt_Y2 <- grt_Y[, grt_stage2, drop = FALSE]

# --- Iterative MNSQ item-fit trimming ---
grt_rows_any <- rowSums(!is.na(grt_Y2)) > 0
grt_Y2_fit   <- grt_Y2[grt_rows_any, , drop = FALSE]
grt_trim <- iterative_mnsq_trim(grt_Y2_fit)
grt_final_items <- grt_trim$keep
grt_Y_final <- grt_Y2_fit[, grt_final_items, drop = FALSE]

cat(sprintf("[GRT] %d items entered screening -> %d after manual exclusion -> %d after MNSQ trim (%d dropped).\n",
            ncol(grt_Y), ncol(grt_Y2), length(grt_final_items), length(grt_trim$dropped)))

# --- Final 1PL vs 2PL comparison (parsimony: 1PL retained per manuscript) ---
grt_mod_1pl_cmp <- mirt::mirt(grt_Y_final, 1, itemtype = "Rasch", verbose = FALSE)
grt_mod_2pl     <- mirt::mirt(grt_Y_final, 1, itemtype = "2PL", verbose = FALSE)
invisible(utils::capture.output(grt_comp <- anova(grt_mod_1pl_cmp, grt_mod_2pl, verbose = FALSE)))
grt_delta_bic <- grt_comp$BIC[2] - grt_comp$BIC[1]

# Separate SE-enabled model for scoring and reliability
grt_mod_1pl <- mirt::mirt(grt_Y_final, 1, itemtype = "Rasch", SE = TRUE, verbose = FALSE)
grt_rel <- mirt::marginal_rxx(grt_mod_1pl)

# --- Scoring (EAP theta, mapped back to full sample) ---
grt_theta_all <- rep(NA_real_, nrow(grt_X))
grt_theta_all[grt_rows_any] <- mirt::fscores(grt_mod_1pl, method = "EAP", full.scores = TRUE)[, 1]

grt_scores <- tibble::tibble(SERIAL = grt_X$SERIAL, theta = grt_theta_all) %>%
  dplyr::left_join(df %>% dplyr::select(SERIAL, Group) %>% dplyr::distinct(), by = "SERIAL") %>%
  dplyr::filter(!is.na(theta))

grt_group_comp <- compare_groups(grt_scores$theta, grt_scores$Group)

cat(sprintf("[GRT] Marginal reliability rxx = %.2f. Group comparison: t(%.1f) = %.2f, p %s, Hedges g = %.2f [%.2f, %.2f].\n",
            grt_rel, grt_group_comp$df, grt_group_comp$t, fmt_p(grt_group_comp$p),
            grt_group_comp$g, grt_group_comp$g_ci[1], grt_group_comp$g_ci[2]))

# --- CFA parcels: odd/even split-half, fitted independently on FINAL item pool ---
grt_nums <- as.integer(gsub("[^0-9]", "", grt_final_items))
grt_A_items <- grt_final_items[grt_nums %% 2 == 1]
grt_B_items <- grt_final_items[grt_nums %% 2 == 0]

grt_parcel_theta_A <- rasch_theta_eap(grt_X[grt_A_items])
grt_parcel_theta_B <- rasch_theta_eap(grt_X[grt_B_items])
grt_parcels <- tibble::tibble(SERIAL = grt_X$SERIAL, GRT_A = grt_parcel_theta_A, GRT_B = grt_parcel_theta_B)
grt_reliability_sb <- spearman_brown(grt_parcels$GRT_A, grt_parcels$GRT_B)

# =============================================================================
# 5. DRT — DELAYED RECOGNITION TASK (two dimensions: D1 familiarity, D2 recollection)
# =============================================================================
# Pipeline: variance screen -> tetrachoric PCA (supports 2D structure) ->
# manual exclusion -> exploratory 2D IRT loading assignment -> iterative
# Infit-MSQ trimming per dimension (0.5-1.5) -> 1PL Rasch per dimension ->
# EAP theta scores for D1 and D2.

drt_cols_raw <- grep("^DRT", names(df), value = TRUE)
df[drt_cols_raw] <- lapply(df[drt_cols_raw], function(v) suppressWarnings(as.numeric(v)))
is_012 <- function(x) { ux <- unique(na.omit(x)); length(ux) > 0 && all(ux %in% c(0, 1, 2)) }
drt_cols <- drt_cols_raw[vapply(df[drt_cols_raw], is_012, logical(1))]

# DRT never separates out test/dummy SERIALs -- they are retained as
# "Healthy" throughout calibration and reporting, matching the original.
drt_patients <- dplyr::filter(df, grepl("^[A-Za-z]{2}\\d{2}$", SERIAL))
drt_healthy  <- dplyr::filter(df, !grepl("^[A-Za-z]{2}\\d{2}$", SERIAL))
drt_common   <- get_union_items(drt_patients, drt_healthy, drt_cols)
drt_combined <- dplyr::bind_rows(
  dplyr::select(drt_patients, SERIAL, dplyr::all_of(drt_common)),
  dplyr::select(drt_healthy,  SERIAL, dplyr::all_of(drt_common))
) %>% dplyr::arrange(SERIAL)

drt_item_var <- sapply(drt_combined[, drt_common, drop = FALSE], var, na.rm = TRUE)
drt_keep_var <- names(drt_item_var[drt_item_var > 0])
drt_X <- dplyr::select(drt_combined, SERIAL, dplyr::all_of(drt_keep_var))

# Response code 2 ("object unknown") recoded to 0 (incorrect) for binary IRT;
# raw 0/1/2 responses are retained separately where needed.
drt_Y_bin <- dplyr::select(drt_X, -SERIAL)
drt_Y_bin[drt_Y_bin == 2] <- 0

# --- Manual exclusion after item biplot inspection ---
drt_manual_drop <- intersect(c("DRT19", "DRT26", "DRT29"), colnames(drt_Y_bin))
drt_items_stage2 <- setdiff(colnames(drt_Y_bin), drt_manual_drop)
drt_Y_bin_rm <- drt_Y_bin[, drt_items_stage2, drop = FALSE]

# --- Dimensionality: tetrachoric PCA (eigenvalues reported in manuscript) ---
drt_tet <- psych::tetrachoric(as.matrix(drt_Y_bin_rm))
drt_R <- drt_tet$rho; drt_R[is.na(drt_R)] <- 0; diag(drt_R) <- 1
drt_R_stable <- as.matrix(Matrix::nearPD(drt_R, corr = TRUE)$mat)
drt_eig <- eigen(drt_R_stable, symmetric = TRUE, only.values = TRUE)$values

# --- Exploratory 2D assignment (|loading| >= .30) ---
drt_fit_idx <- rowSums(!is.na(drt_Y_bin_rm)) > 0
drt_Y_fit   <- drt_Y_bin_rm[drt_fit_idx, , drop = FALSE]
drt_mod_2d  <- mirt::mirt(drt_Y_fit, 2, itemtype = "2PL", SE = FALSE, verbose = FALSE)
drt_rot     <- summary(drt_mod_2d, rotate = "oblimin", verbose = FALSE)
drt_loadmat <- as.matrix(drt_rot$rotF)
rownames(drt_loadmat) <- colnames(drt_Y_fit)

drt_maxdim <- apply(abs(drt_loadmat), 1, which.max)
drt_maxval <- apply(abs(drt_loadmat), 1, max)
drt_D1_items <- rownames(drt_loadmat)[drt_maxdim == 1 & drt_maxval >= 0.30]
drt_D2_items <- rownames(drt_loadmat)[drt_maxdim == 2 & drt_maxval >= 0.30]

# --- Iterative Infit-MSQ item trimming per dimension (0.5-1.5) ---
iterative_infit_trim <- function(item_mat, lower = 0.5, upper = 1.5) {
  current <- colnames(item_mat)
  dropped <- character(0)
  repeat {
    if (length(current) < 3) break
    mod <- TAM::tam.mml(as.matrix(item_mat[, current, drop = FALSE]), verbose = FALSE)
    fit <- as.data.frame(TAM::msq.itemfit(mod)$itemfit)
    infit_col <- grep("Infit$", names(fit), value = TRUE)[1]
    misfit <- fit[fit[[infit_col]] > upper | fit[[infit_col]] < lower, ]
    if (nrow(misfit) == 0) break
    misfit$dev <- abs(misfit[[infit_col]] - 1)
    worst <- misfit$item[which.max(misfit$dev)]
    dropped <- c(dropped, worst)
    current <- setdiff(current, worst)
  }
  list(keep = current, dropped = dropped)
}

drt_Ybin_D1 <- drt_Y_bin_rm[, intersect(drt_D1_items, colnames(drt_Y_bin_rm)), drop = FALSE]
drt_Ybin_D2 <- drt_Y_bin_rm[, intersect(drt_D2_items, colnames(drt_Y_bin_rm)), drop = FALSE]
drt_sel_D1 <- iterative_infit_trim(drt_Ybin_D1)
drt_sel_D2 <- iterative_infit_trim(drt_Ybin_D2)
drt_final_D1 <- drt_sel_D1$keep
drt_final_D2 <- drt_sel_D2$keep

cat(sprintf("[DRT] %d items entered dimensionality screening -> %d after manual exclusion. D1 final: %d items (%d dropped). D2 final: %d items (%d dropped).\n",
            ncol(drt_Y_bin), ncol(drt_Y_bin_rm),
            length(drt_final_D1), length(drt_sel_D1$dropped),
            length(drt_final_D2), length(drt_sel_D2$dropped)))

# --- Refit 1PL Rasch per dimension on final item pools, EAP scoring ---
fit_dim_theta <- function(Ybin_full, keep_items) {
  Ysub <- Ybin_full[, keep_items, drop = FALSE]
  idx <- rowSums(!is.na(Ysub)) > 0
  mod_1pl <- mirt::mirt(Ysub[idx, , drop = FALSE], 1, itemtype = "Rasch", SE = TRUE, verbose = FALSE)
  mod_2pl <- mirt::mirt(Ysub[idx, , drop = FALSE], 1, itemtype = "2PL", SE = FALSE, verbose = FALSE)
  invisible(utils::capture.output(cmp <- anova(mod_1pl, mod_2pl, verbose = FALSE)))
  dbic <- cmp$BIC[1] - cmp$BIC[2]
  rel <- mirt::marginal_rxx(mod_1pl)
  theta <- rep(NA_real_, nrow(Ybin_full))
  theta[idx] <- mirt::fscores(mod_1pl, method = "EAP", full.scores = TRUE)[, 1]
  list(mod = mod_1pl, rel = rel, dbic = dbic, theta = theta)
}

drt_fit1 <- fit_dim_theta(drt_Y_bin_rm, drt_final_D1)
drt_fit2 <- fit_dim_theta(drt_Y_bin_rm, drt_final_D2)

drt_scores <- tibble::tibble(
  SERIAL = drt_X$SERIAL,
  Group  = assign_group_simple(drt_X$SERIAL),
  theta_DRT_D1 = drt_fit1$theta,
  theta_DRT_D2 = drt_fit2$theta
)

drt_gc_D1 <- compare_groups(drt_scores$theta_DRT_D1, drt_scores$Group)
drt_gc_D2 <- compare_groups(drt_scores$theta_DRT_D2, drt_scores$Group)

cat(sprintf("[DRT D1] rxx = %.2f | Hedges g = %.2f [%.2f, %.2f]\n",
            drt_fit1$rel, drt_gc_D1$g, drt_gc_D1$g_ci[1], drt_gc_D1$g_ci[2]))
cat(sprintf("[DRT D2] rxx = %.2f | Hedges g = %.2f [%.2f, %.2f]\n",
            drt_fit2$rel, drt_gc_D2$g, drt_gc_D2$g_ci[1], drt_gc_D2$g_ci[2]))

# --- CFA parcel: D1 (familiarity-based recognition) only, per Figure 4 caption ---
# odd/even split on the final D1 item pool, fitted independently
drt_D1_nums <- as.integer(gsub("[^0-9]", "", drt_final_D1))
drt_A_items <- drt_final_D1[drt_D1_nums %% 2 == 1]
drt_B_items <- drt_final_D1[drt_D1_nums %% 2 == 0]

drt_parcel_theta_A <- rasch_theta_eap(drt_Y_bin_rm[drt_A_items])
drt_parcel_theta_B <- rasch_theta_eap(drt_Y_bin_rm[drt_B_items])
drt_parcels <- tibble::tibble(SERIAL = drt_X$SERIAL, DRT_A = drt_parcel_theta_A, DRT_B = drt_parcel_theta_B)
drt_reliability_sb <- spearman_brown(drt_parcels$DRT_A, drt_parcels$DRT_B)

# =============================================================================
# 6. PAL — PAIRED ASSOCIATION LEARNING TASK
# =============================================================================
# Pipeline: coalesce 4 parallel versions onto 5 hierarchical levels (L1-L5) ->
# hierarchical stopping-rule propagation -> ordinal discretisation ->
# unidimensionality check (PRINCALS) -> Partial Credit Model (PCM) vs Graded
# Response Model (GRM); PCM retained on BIC parsimony grounds -> EAP theta.

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

# Ordinal discretisation (10 bins) for polytomous IRT
pal_Y_ord <- pal_Y_raw %>% dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(round(.x * 10, 0))))
pal_keep_idx <- rowSums(!is.na(pal_Y_ord)) > 0
pal_Y_ord <- pal_Y_ord[pal_keep_idx, ]
pal_serial_kept <- pal_X$SERIAL[pal_keep_idx]
pal_group_kept <- assign_group_simple(pal_serial_kept)

# --- Dimensionality: PRINCALS ---
pal_unique_counts <- apply(pal_Y_ord, 2, function(x) length(unique(na.omit(x))))
pal_princals_cols <- names(pal_Y_ord)[pal_unique_counts > 1]
pal_Y_princals <- pal_Y_ord[, pal_princals_cols, drop = FALSE]
pal_princals_fit <- Gifi::princals(as.data.frame(lapply(pal_Y_princals, as.ordered)), ndim = 2, ordinal = TRUE, ties = "s")
pal_var_exp1 <- pal_princals_fit$evals[1] / ncol(pal_Y_princals) * 100

# --- Model comparison: PCM (equal discrimination) vs GRM (free discrimination) ---
pal_spec_pcm <- mirt::mirt.model(paste0("F = 1-", ncol(pal_Y_ord), "\nCONSTRAIN = (1-", ncol(pal_Y_ord), ", a1)"))
pal_mod_pcm  <- mirt::mirt(pal_Y_ord, model = pal_spec_pcm, itemtype = "gpcm", SE = TRUE, verbose = FALSE)
pal_mod_grm  <- mirt::mirt(pal_Y_ord, 1, itemtype = "graded", SE = FALSE, verbose = FALSE)
invisible(utils::capture.output(pal_anova <- anova(pal_mod_pcm, pal_mod_grm, verbose = FALSE)))
pal_delta_aic <- pal_anova$AIC[2] - pal_anova$AIC[1]
pal_delta_bic <- pal_anova$BIC[2] - pal_anova$BIC[1]
pal_chi_df <- pal_anova$df[2]
pal_chi_sq <- abs(pal_anova$X2[2])
pal_chi_p  <- pchisq(pal_chi_sq, df = pal_chi_df, lower.tail = FALSE)

pal_rel <- mirt::marginal_rxx(pal_mod_pcm)

# --- Scoring ---
pal_fscores <- mirt::fscores(pal_mod_pcm, method = "EAP", full.scores = TRUE, full.scores.SE = TRUE)
pal_scores <- tibble::tibble(SERIAL = pal_serial_kept, theta = pal_fscores[, "F"], Group = pal_group_kept)

pal_group_comp <- compare_groups(pal_scores$theta, pal_scores$Group)

cat(sprintf("[PAL] PCM retained (delta AIC favours GRM = %.1f, chi2(%d) = %.2f, p %s; delta BIC favours PCM = %+.1f). rxx = %.2f. Hedges g = %.2f [%.2f, %.2f]\n",
            pal_delta_aic, pal_chi_df, pal_chi_sq, fmt_p(pal_chi_p), pal_delta_bic,
            pal_rel, pal_group_comp$g, pal_group_comp$g_ci[1], pal_group_comp$g_ci[2]))

# --- CFA parcels: interleaved levels A = {L1,L3,L5}, B = {L2,L4} ---
# Binarised (any partial credit > 0 counts as pass) with the stopping rule
# re-applied within each parcel, per manuscript Section 2.5.
pal_bin <- pal_Y_raw
pal_bin[pal_level_names] <- lapply(pal_bin[pal_level_names], function(x) as.integer(!is.na(x) & x > 0))
pal_bin <- pal_bin[pal_keep_idx, , drop = FALSE]

pal_A_cols <- c("L1", "L3", "L5")
pal_B_cols <- c("L2", "L4")
pal_A_mat <- t(apply(as.matrix(pal_bin[pal_A_cols]), 1, propagate_failure)); colnames(pal_A_mat) <- pal_A_cols
pal_B_mat <- t(apply(as.matrix(pal_bin[pal_B_cols]), 1, propagate_failure)); colnames(pal_B_mat) <- pal_B_cols

# Falls back to a z-scored sum score if IRT calibration is not stable on a
# 2-3 item parcel (documented fallback for degenerate cases).
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
# 7. CI — COGNITIVE INTERFERENCE TASK
# =============================================================================
# Pipeline: joint 2-dimensional MIRT (Q-matrix: baseline AC vs. interference
# AD) -> iterative item trimming (Outfit 0.5-1.5) -> confirmatory 1D-vs-2D
# comparison -> WLE theta per dimension -> CI resilience = AD residualised
# on AC using the joint model's disattenuated latent slope.

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

# Response coding differs by whether "1" or "0" denotes the correct
# recognition response for a given item; resolved empirically per item.
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

# --- Joint 2D MIRT (Q-matrix) calibration ---
ci_Q <- matrix(0, nrow = ncol(ci_mat), ncol = 2, dimnames = list(colnames(ci_mat), c("Dim_AC", "Dim_AD")))
ci_Q[grepl("^AC", rownames(ci_Q)), 1] <- 1
ci_Q[grepl("^AD", rownames(ci_Q)), 2] <- 1

ci_mod_init <- TAM::tam.mml(resp = ci_mat, Q = ci_Q, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
ci_fit_init <- TAM::tam.fit(ci_mod_init, progress = FALSE)
ci_misfits  <- ci_fit_init$itemfit %>% dplyr::filter(Outfit > 1.5 | Outfit < 0.5 | Infit > 1.5)

ci_pair_key <- function(items) sub("^(AC|AD)_", "", items)

if (nrow(ci_misfits) > 0) {
  ci_drop_items <- as.character(ci_misfits$parameter)
  ci_mat_trim <- ci_mat[, !(colnames(ci_mat) %in% ci_drop_items), drop = FALSE]
} else {
  ci_mat_trim <- ci_mat
}

# Item-fit trimming can remove one member of an AC/AD pair without the
# other; retain only pairs where both the baseline (AC) and interference
# (AD) item survived, so the two item sets stay matched for the
# difference-score parcels below.
ci_ac_all <- grep("^AC", colnames(ci_mat_trim), value = TRUE)
ci_ad_all <- grep("^AD", colnames(ci_mat_trim), value = TRUE)
ci_paired_keys <- intersect(ci_pair_key(ci_ac_all), ci_pair_key(ci_ad_all))
ci_keep_items <- intersect(colnames(ci_mat_trim), c(paste0("AC_", ci_paired_keys), paste0("AD_", ci_paired_keys)))

ci_mat_final <- ci_mat_trim[, ci_keep_items, drop = FALSE]
ci_Q_final   <- ci_Q[rownames(ci_Q) %in% ci_keep_items, , drop = FALSE]
ci_mod_final <- TAM::tam.mml(resp = ci_mat_final, Q = ci_Q_final, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)

cat(sprintf("[CI] %d items entered joint calibration -> %d retained after item-fit trimming and AC/AD pairing (%d dropped).\n",
            ncol(ci_mat), ncol(ci_mat_final), ncol(ci_mat) - ncol(ci_mat_final)))

# --- Confirmatory 1D vs 2D comparison ---
ci_mod_1d <- TAM::tam.mml(resp = ci_mat_final, irtmodel = "1PL", control = list(snodes = 1500, qmc = TRUE), verbose = FALSE)
ci_bic_1d <- ci_mod_1d$ic$BIC
ci_bic_2d <- ci_mod_final$ic$BIC

# --- WLE thetas + latent residualisation ---
ci_wle <- TAM::tam.wle(ci_mod_final, progress = FALSE)
ci_rel_ac <- TAM::WLErel(ci_wle$theta.Dim01, ci_wle$error.Dim01)
ci_rel_ad <- TAM::WLErel(ci_wle$theta.Dim02, ci_wle$error.Dim02)

ci_latent_cov <- ci_mod_final$variance
ci_beta <- ci_latent_cov[1, 2] / ci_latent_cov[1, 1]

ci_theta <- tibble::tibble(
  SERIAL = rownames(ci_mat_final),
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

# --- CFA parcels: odd/even split of trial-matched difference scores (AD-AC),
# fitted independently on the final item pool ---
ci_ac_items <- grep("^AC", colnames(ci_mat_final), value = TRUE)
ci_ad_items <- grep("^AD", colnames(ci_mat_final), value = TRUE)
ci_ac_items <- ci_ac_items[order(ci_pair_key(ci_ac_items))]
ci_ad_items <- ci_ad_items[order(ci_pair_key(ci_ad_items))]
stopifnot(identical(ci_pair_key(ci_ac_items), ci_pair_key(ci_ad_items)))

ci_diff <- as.data.frame(ci_mat_final[, ci_ad_items, drop = FALSE] - ci_mat_final[, ci_ac_items, drop = FALSE])
names(ci_diff) <- sub("^AD", "DIFF", ci_ad_items)
ci_diff_cols <- names(ci_diff)
ci_diff_A <- ci_diff_cols[seq(1, length(ci_diff_cols), by = 2)]
ci_diff_B <- ci_diff_cols[seq(2, length(ci_diff_cols), by = 2)]

ci_bin_diff <- function(mat, d_min, d_max, nbins = 10L) {
  mat01 <- (mat - d_min) / (d_max - d_min)
  mat01[mat01 < 0] <- 0; mat01[mat01 > 1] <- 1
  as.data.frame(matrix(as.integer(round(mat01 * nbins, 0)), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
}
ci_diff_range <- range(as.matrix(ci_diff[ci_diff_cols]), na.rm = TRUE)
ci_A_binned <- ci_bin_diff(as.matrix(ci_diff[ci_diff_A]), ci_diff_range[1], ci_diff_range[2])
ci_B_binned <- ci_bin_diff(as.matrix(ci_diff[ci_diff_B]), ci_diff_range[1], ci_diff_range[2])

ci_parcels <- tibble::tibble(
  SERIAL    = rownames(ci_mat_final),
  CI_DIFF_A = rasch_theta_eap(ci_A_binned),
  CI_DIFF_B = rasch_theta_eap(ci_B_binned)
)
ci_reliability_sb <- spearman_brown(ci_parcels$CI_DIFF_A, ci_parcels$CI_DIFF_B)

# =============================================================================
# 8. SART-ED — SUSTAINED ATTENTION TO RESPONSE TASK WITH EXTERNAL DISTRACTION
# =============================================================================
# Trial-level data are stored as an embedded CSV string per participant and
# must be expanded first. Two outcomes are derived on the FULL trial set:
#   - External Interference control (two-step mixed-effects + MRFA, per
#     manuscript Section 2.5)
#   - Mental Speed (baseline Go-trial log-RT, outlier-screened, single factor)
# Odd/even trial parcels (SART_A/B, MS_A/B) are fitted independently for
# split-half reliability and the CFA.

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

# Single-factor score (minres FA) from a set of resilience indicators, with a
# reflection check so higher scores always mean "more resilient".
latent_score_1f <- function(df_ind, method = c("fa", "pca")) {
  method <- match.arg(method)
  X <- df_ind %>% dplyr::select(-SERIAL) %>% as.data.frame()
  ok_var <- sapply(X, function(v) stats::sd(v, na.rm = TRUE)) > 0
  X <- X[, ok_var, drop = FALSE]
  cc <- stats::complete.cases(X)
  if (ncol(X) < 2 || sum(cc) < 30) return(list(ok = FALSE, theta = tibble::tibble(SERIAL = df_ind$SERIAL, theta = NA_real_)))
  Xcc <- X[cc, , drop = FALSE]
  evals <- eigen(cor(Xcc, use = "pairwise.complete.obs"))$values
  var_exp <- evals[1] / sum(evals) * 100
  fit <- if (method == "fa") psych::fa(Xcc, nfactors = 1, fm = "minres", rotate = "none") else psych::principal(Xcc, nfactors = 1, rotate = "none")
  scores <- as.numeric(fit$scores[, 1])
  loadmean <- mean(as.numeric(unclass(fit$loadings)), na.rm = TRUE)
  if (loadmean < 0) scores <- scores * -1
  theta <- rep(NA_real_, nrow(df_ind)); theta[cc] <- scores
  list(ok = TRUE, theta = tibble::tibble(SERIAL = df_ind$SERIAL, theta = theta), evals = evals, var_exp = var_exp)
}

# --- Load and expand trial-level SART data ---
sart_raw <- readr::read_csv(file.path(data_dir, "EMA4Stroke_Baseline_GoNoGo_Raw.csv"), show_col_types = FALSE, progress = FALSE)
sart_raw$SERIAL[sart_raw$SERIAL %in% c("ZZ90", "ZZ91")] <- "ZZ90_91"

parse_trials <- function(serial, csv_string) {
  if (is.na(csv_string) || !nzchar(csv_string)) return(tibble::tibble())
  txt <- gsub("\\\\n", "\n", csv_string)
  out <- suppressWarnings(readr::read_csv(I(txt), show_col_types = FALSE, progress = FALSE,
                                          col_types = readr::cols(.default = readr::col_character())))
  if (!nrow(out)) return(tibble::tibble())
  out <- out %>% dplyr::rename_with(make.names)
  out %>% dplyr::transmute(
    SERIAL  = as.character(serial),
    Trial   = suppressWarnings(as.integer(Trial)),
    GoNoGo  = suppressWarnings(as.integer(Go.NoGo)),
    Face    = suppressWarnings(as.numeric(Face)),
    Correct = suppressWarnings(as.integer(Correct)),
    RT      = suppressWarnings(as.numeric(RT))
  )
}

sart_trials <- purrr::map2_dfr(sart_raw$SERIAL, sart_raw$GN03_08, parse_trials) %>%
  dplyr::filter(!is.na(SERIAL), SERIAL != "") %>%
  dplyr::mutate(
    GoNoGo    = as01(GoNoGo),
    Correct   = as01(Correct),
    condition = dplyr::if_else(Face == 0, 0L, 1L),
    parcel    = dplyr::if_else(Trial %% 2 == 1, "A", "B"),
    Group     = assign_group(SERIAL)
  )

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

# =============================================================================
# 8a. MAIN ANALYSIS — External Interference (full trial set)
# =============================================================================
sart_elig <- make_cov(sart_go_rt, "go_rt") %>%
  dplyr::full_join(make_cov(sart_go_acc, "go_acc"), by = "SERIAL") %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("cond"), ~ dplyr::coalesce(.x, 0L))) %>%
  dplyr::mutate(elig_go_rt  = cond0_n_go_rt  >= sart_MIN_TRIALS_FULL & cond1_n_go_rt  >= sart_MIN_TRIALS_FULL,
                elig_go_acc = cond0_n_go_acc >= sart_MIN_TRIALS_FULL & cond1_n_go_acc >= sart_MIN_TRIALS_FULL)

# Model 1: log-RT on correct Go trials
dat_rt <- sart_go_rt %>%
  dplyr::inner_join(dplyr::select(sart_elig, SERIAL, elig_go_rt), by = "SERIAL") %>%
  dplyr::filter(elig_go_rt) %>%
  dplyr::left_join(dplyr::select(sart_baseline, SERIAL, baseline_logRT_c), by = "SERIAL") %>%
  dplyr::filter(!is.na(baseline_logRT_c))
fit_go_rt <- lme4::lmer(logRT ~ condition * baseline_logRT_c + (condition | SERIAL), data = dat_rt, REML = TRUE,
                        control = lme4::lmerControl(optimizer = "bobyqa"))
eb_rt <- extract_eb_slope(fit_go_rt) %>% dplyr::rename(eb_slope_rt = eb_slope)

# Model 2: Go-trial accuracy (logit link)
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

# Minimum Residual Factor Analysis -> single External Interference theta
ei_fit <- latent_score_1f(dplyr::select(sart_indicators, SERIAL, ind_rt_resilience, ind_acc_resilience), method = "fa")
sart_ei_theta <- ei_fit$theta %>% dplyr::rename(SART_ExternalInterference_Theta = theta)

# =============================================================================
# 8b. MAIN ANALYSIS — Mental Speed (baseline Go-trial log-RT, outlier-screened)
# =============================================================================
sart_ms_ind <- sart_go_rt %>%
  dplyr::filter(condition == 0L) %>%
  dplyr::group_by(SERIAL) %>%
  dplyr::summarise(ms_mean_logRT = mean(logRT, na.rm = TRUE), ms_sd_logRT = stats::sd(logRT, na.rm = TRUE), ms_n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(ms_speed_ind = -1 * as.numeric(scale(ms_mean_logRT)), ms_cons_ind = -1 * as.numeric(scale(ms_sd_logRT)))

# Univariate (MAD-based, threshold = 3) and robust bivariate (MCD Mahalanobis, alpha=.001) outlier screen
out_uni <- function(x, threshold = 3) {
  med <- stats::median(x, na.rm = TRUE)
  abs_dev <- abs(x - med)
  mad_val <- stats::median(abs_dev, na.rm = TRUE)
  abs_dev > (threshold * mad_val)
}
X_ms <- dplyr::filter(sart_ms_ind, stats::complete.cases(ms_mean_logRT, ms_sd_logRT))
rob <- suppressWarnings(MASS::cov.rob(dplyr::select(X_ms, ms_mean_logRT, ms_sd_logRT), method = "mcd"))
m_dist <- stats::mahalanobis(dplyr::select(X_ms, ms_mean_logRT, ms_sd_logRT), center = rob$center, cov = rob$cov)
bivariate_outlier_ids <- X_ms$SERIAL[m_dist > stats::qchisq(0.999, df = 2)]

ms_clean <- sart_ms_ind %>%
  dplyr::filter(!out_uni(ms_mean_logRT), !out_uni(ms_sd_logRT), !(SERIAL %in% bivariate_outlier_ids))

ms_fit <- latent_score_1f(dplyr::select(ms_clean, SERIAL, ms_speed_ind, ms_cons_ind), method = "fa")
sart_ms_theta <- ms_fit$theta %>% dplyr::rename(SART_MentalSpeed_Theta = theta)

sart_scores <- tibble::tibble(SERIAL = unique(sart_trials$SERIAL), Group = assign_group(unique(sart_trials$SERIAL))) %>%
  dplyr::left_join(sart_ei_theta, by = "SERIAL") %>%
  dplyr::left_join(sart_ms_theta, by = "SERIAL")

sart_gc_ei <- compare_groups(sart_scores$SART_ExternalInterference_Theta, sart_scores$Group)
sart_gc_ms <- compare_groups(sart_scores$SART_MentalSpeed_Theta, sart_scores$Group)

cat(sprintf("[SART-ED External Interference] Hedges g = %.2f [%.2f, %.2f]\n", sart_gc_ei$g, sart_gc_ei$g_ci[1], sart_gc_ei$g_ci[2]))
cat(sprintf("[SART-ED Mental Speed] Hedges g = %.2f [%.2f, %.2f]\n", sart_gc_ms$g, sart_gc_ms$g_ci[1], sart_gc_ms$g_ci[2]))

# =============================================================================
# 8c. CFA PARCELS — independent per-parcel EB slope models (odd vs. even trials)
# =============================================================================
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

# Mental Speed parcels (baseline Go-trial log-RT, odd vs. even trials)
compute_ms_half <- function(p) {
  half <- sart_go_rt %>% dplyr::filter(condition == 0L, parcel == p) %>%
    dplyr::group_by(SERIAL) %>%
    dplyr::summarise(h_mean_logRT = mean(logRT, na.rm = TRUE), h_sd_logRT = stats::sd(logRT, na.rm = TRUE), h_n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(h_n >= 3, is.finite(h_sd_logRT)) %>%
    dplyr::mutate(h_speed_ind = -1 * as.numeric(scale(h_mean_logRT)), h_cons_ind = -1 * as.numeric(scale(h_sd_logRT)))
  fit <- latent_score_1f(dplyr::select(half, SERIAL, h_speed_ind, h_cons_ind), method = "fa")
  fit$theta
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
# 9. MoCA (stroke sample only) & MASTER THETA TABLE
# =============================================================================
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

# The master table's Group label is recomputed fresh here (permissive,
# suffix-tolerant, no test/dummy exclusion), independent of whatever group
# definition a given task used internally for its own calibration/reporting.
assign_group_master <- function(serial) dplyr::if_else(grepl("^[A-Z]{2}\\d{2}(_\\d{2})?$", serial), "Patient", "Healthy")

master_theta <- df %>%
  dplyr::select(SERIAL) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Group = assign_group_master(SERIAL)) %>%
  dplyr::left_join(dplyr::select(grt_scores, SERIAL, theta_GRT = theta), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(drt_scores, SERIAL, theta_DRT_D1, theta_DRT_D2), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(pal_scores, SERIAL, theta_PAL = theta), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(ci_theta, SERIAL, theta_CI_resilience = CI_Resilience), by = "SERIAL") %>%
  dplyr::left_join(dplyr::select(sart_scores, SERIAL, theta_SART_EI = SART_ExternalInterference_Theta, theta_SART_MS = SART_MentalSpeed_Theta), by = "SERIAL") %>%
  dplyr::left_join(moca_df, by = "SERIAL")

readr::write_csv(master_theta, file.path(output_dir, "EMA4Stroke_Master_Theta_Scores.csv"))

# --- Diagnostic: effect of test/dummy SERIALs on the GRT-PAL correlation ---
is_test_dummy <- grepl("test|dummy", tolower(master_theta$SERIAL))
cat(sprintf("\n[Diagnostic] %d of %d SERIALs match a test/dummy pattern.\n", sum(is_test_dummy), nrow(master_theta)))
cat(sprintf("[Diagnostic] GRT-PAL r, all included    : %.3f (n = %d)\n",
            cor(master_theta$theta_GRT, master_theta$theta_PAL, use = "complete.obs"),
            sum(complete.cases(master_theta$theta_GRT, master_theta$theta_PAL))))
cat(sprintf("[Diagnostic] GRT-PAL r, test/dummy excl. : %.3f (n = %d)\n",
            cor(master_theta$theta_GRT[!is_test_dummy], master_theta$theta_PAL[!is_test_dummy], use = "complete.obs"),
            sum(complete.cases(master_theta$theta_GRT[!is_test_dummy], master_theta$theta_PAL[!is_test_dummy]))))

# =============================================================================
# 10. FIGURE 2 — LATENT ABILITY (THETA) BY GROUP, ALL DOMAINS
# =============================================================================
domain_labels <- c(
  theta_SART_EI       = "External interference control",
  theta_CI_resilience = "Internal interference control",
  theta_SART_MS       = "Mental speed",
  theta_GRT           = "Reasoning",
  theta_DRT_D1        = "Familiarity-based recognition",
  theta_DRT_D2        = "Recollection-based discrimination",
  theta_PAL           = "Associative memory"
)
domain_order <- unname(domain_labels)

fig2_data <- master_theta %>%
  dplyr::select(SERIAL, Group, dplyr::all_of(names(domain_labels))) %>%
  tidyr::pivot_longer(cols = -c(SERIAL, Group), names_to = "raw_test", values_to = "score") %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::mutate(
    Test  = factor(domain_labels[raw_test], levels = domain_order),
    Group = factor(dplyr::if_else(Group == "Healthy", "Controls", "Patients"), levels = c("Controls", "Patients"))
  )
levels(fig2_data$Test) <- stringr::str_wrap(levels(fig2_data$Test), width = 18)

sig_data <- fig2_data %>%
  dplyr::group_by(Test) %>%
  dplyr::summarise(p_val = tryCatch(t.test(score ~ Group)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  dplyr::mutate(
    asterisks = dplyr::case_when(is.na(p_val) ~ "", p_val < .001 ~ "***", p_val < .01 ~ "**", p_val < .05 ~ "*", TRUE ~ "ns"),
    y_pos  = max(fig2_data$score, na.rm = TRUE) + 0.10 * diff(range(fig2_data$score, na.rm = TRUE)),
    y_text = max(fig2_data$score, na.rm = TRUE) + 0.16 * diff(range(fig2_data$score, na.rm = TRUE))
  )

fig2 <- ggplot(fig2_data, aes(x = Group, y = score, fill = Group)) +
  gghalves::geom_half_violin(side = "r", alpha = 0.6, color = NA, adjust = 1.2) +
  gghalves::geom_half_point(aes(color = Group), side = "l", range_scale = 0.4, alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8, position = position_nudge(x = 0.03)) +
  scale_fill_manual(values = c(Controls = hc_color, Patients = pat_color)) +
  scale_color_manual(values = c(Controls = hc_color, Patients = pat_color)) +
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
  geom_segment(data = sig_data, aes(x = 1, xend = 2, y = y_pos, yend = y_pos), inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = sig_data, aes(x = 1.5, y = y_text, label = asterisks), inherit.aes = FALSE, size = 5, vjust = 0.5, fontface = "bold")

ggsave(file.path(output_dir, "Figure2_Group_Comparison_Raincloud.png"), fig2, width = 12, height = 14, dpi = 300, bg = "transparent")

# =============================================================================
# 11. FIGURE 3 — INTER-TASK CORRELATION MATRIX & SPLIT-HALF RELIABILITY
# =============================================================================
fig3_labels <- c(
  theta_SART_EI       = "External interference control",
  theta_CI_resilience = "Internal interference control",
  theta_SART_MS       = "Mental speed",
  theta_GRT           = "Reasoning",
  theta_DRT_D1        = "Familiarity-based recognition",
  theta_DRT_D2        = "Recollection-based discrimination",
  theta_PAL           = "Associative memory",
  MoCA                = "MoCA"
)

fig3_data <- master_theta %>%
  dplyr::select(SERIAL, Group, dplyr::all_of(names(fig3_labels))) %>%
  dplyr::mutate(Group = factor(dplyr::if_else(Group == "Healthy", "Controls", "Patients"), levels = c("Controls", "Patients")))
names(fig3_data)[match(names(fig3_labels), names(fig3_data))] <- unname(fig3_labels)

# Split-half reliability per domain, computed from the independently-fitted
# odd/even parcel thetas (Sections 4-8 above).
reliability_diag <- c(
  "External interference control"      = sart_reliability_sb,
  "Internal interference control"      = ci_reliability_sb,
  "Mental speed"                       = ms_reliability_sb,
  "Reasoning"                          = grt_reliability_sb,
  "Familiarity-based recognition"      = drt_reliability_sb,
  "Recollection-based discrimination"  = drt_reliability_sb,
  "Associative memory"                 = pal_reliability_sb
)

fig3_upper <- function(data, mapping, ...) {
  x_col <- rlang::as_name(mapping$x); y_col <- rlang::as_name(mapping$y)
  x_val <- data[[x_col]]; y_val <- data[[y_col]]; grp <- data$Group
  cor_lab <- function(x, y) {
    if (sum(complete.cases(x, y)) < 3) return("")
    ct <- cor.test(x, y); stars <- dplyr::case_when(ct$p.value < .001 ~ "***", ct$p.value < .01 ~ "**", ct$p.value < .05 ~ "*", TRUE ~ "")
    sprintf("%.3f%s", ct$estimate, stars)
  }
  ggplot() +
    annotate("text", x = .5, y = .65, label = paste("All:", cor_lab(x_val, y_val)), size = 4.5, fontface = "bold") +
    annotate("text", x = .5, y = .45, label = paste("Controls:", cor_lab(x_val[grp == "Controls"], y_val[grp == "Controls"])), size = 4, color = hc_color) +
    annotate("text", x = .5, y = .25, label = paste("Patients:", cor_lab(x_val[grp == "Patients"], y_val[grp == "Patients"])), size = 4, color = pat_color) +
    theme_void() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}
fig3_lower <- function(data, mapping, ...) {
  ggplot(data, mapping) + geom_point(alpha = .5, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .8) +
    scale_color_manual(values = c(Controls = hc_color, Patients = pat_color)) + theme_corr
}
fig3_diag <- function(data, mapping, ...) {
  x_col <- rlang::as_name(mapping$x)
  if (x_col == "MoCA") {
    ggplot(data, mapping) + geom_density(alpha = .6, aes(fill = Group)) +
      scale_fill_manual(values = c(Controls = hc_color, Patients = pat_color)) + theme_corr
  } else {
    rel_val <- reliability_diag[x_col]
    ggplot() + annotate("text", x = .5, y = .5, label = if (is.na(rel_val)) "" else sprintf("%.2f", rel_val), size = 6, fontface = "bold") +
      theme_void() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  }
}

fig3 <- GGally::ggpairs(
  data         = dplyr::select(fig3_data, -SERIAL),
  mapping      = aes(color = Group, fill = Group),
  columns      = unname(fig3_labels),
  upper        = list(continuous = fig3_upper),
  lower        = list(continuous = fig3_lower),
  diag         = list(continuous = fig3_diag),
  title        = "Person Score Correlations by Group (Final Latent Constructs)",
  columnLabels = stringr::str_wrap(unname(fig3_labels), width = 14)
) + theme_corr + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
                       axis.text.y = element_text(color = "black"))

ggsave(file.path(output_dir, "Figure3_Correlation_Reliability_Matrix.png"), fig3, width = 12, height = 12, dpi = 300, bg = "transparent")

theta_only_cols <- setdiff(unname(fig3_labels), "MoCA")
cor_mat <- cor(dplyr::select(fig3_data, dplyr::all_of(theta_only_cols)), use = "pairwise.complete.obs")
readr::write_csv(
  tibble::as_tibble(cbind(variable = rownames(cor_mat), as.data.frame(cor_mat))),
  file.path(output_dir, "Table_ThetaCorrelations.csv")
)

# =============================================================================
# 12. FIGURE 4 — S-1 BIFACTOR CFA (MODEL 3: g + 5 SPECIFIC FACTORS)
# =============================================================================
# Manifest indicators are the parcel-level theta scores computed above
# (Sections 4-8), each fitted independently on that task's final calibrated
# item pool. The DRT indicator uses D1 (familiarity-based recognition)
# items only.

df_master <- purrr::reduce(
  list(grt_parcels, drt_parcels, pal_parcels, sart_parcels, ci_parcels),
  dplyr::full_join, by = "SERIAL"
)

parcel_vars <- c("SART_A", "SART_B", "CI_DIFF_A", "CI_DIFF_B", "MS_A", "MS_B",
                 "GRT_A", "GRT_B", "DRT_A", "DRT_B", "PAL_A", "PAL_B")

df_cfa <- df_master
df_cfa[parcel_vars] <- scale(df_cfa[parcel_vars])

model_bifactor <- '
  g =~ SART_A + SART_B +
       CI_DIFF_A + CI_DIFF_B +
       MS_A + MS_B +
       GRT_A + GRT_B +
       DRT_A + DRT_B +
       PAL_A + PAL_B

  ExternalInterference =~ ei_s*SART_A    + ei_s*SART_B
  InternalInterference =~ ii_c*CI_DIFF_A + ii_c*CI_DIFF_B
  MentalProcessing     =~ mp_m*MS_A      + mp_m*MS_B
  FiguralMemory        =~ fm_d*DRT_A     + fm_d*DRT_B
  FiguralSpatialMemory =~ fs_p*PAL_A     + fs_p*PAL_B

  # Heywood-case resolution: PAL parcel B residual variance fixed to zero
  PAL_B ~~ 0*PAL_B

  # Local dependence between the two memory-specific factors (freed)
  FiguralMemory ~~ FiguralSpatialMemory

  # All other between-factor covariances fixed to 0 (orthogonal bifactor structure)
  g ~~ 0*ExternalInterference
  g ~~ 0*InternalInterference
  g ~~ 0*MentalProcessing
  g ~~ 0*FiguralMemory
  g ~~ 0*FiguralSpatialMemory
  ExternalInterference ~~ 0*InternalInterference
  ExternalInterference ~~ 0*MentalProcessing
  ExternalInterference ~~ 0*FiguralMemory
  ExternalInterference ~~ 0*FiguralSpatialMemory
  InternalInterference ~~ 0*MentalProcessing
  InternalInterference ~~ 0*FiguralMemory
  InternalInterference ~~ 0*FiguralSpatialMemory
  MentalProcessing ~~ 0*FiguralMemory
  MentalProcessing ~~ 0*FiguralSpatialMemory

  # Shared method variance between SART- and MS-derived parcels (matched a/b)
  MS_A ~~ SART_A
  MS_B ~~ SART_B
'

fit_bifactor <- lavaan::cfa(model_bifactor, data = df_cfa, estimator = "MLR", missing = "fiml", std.lv = TRUE)

cat("\n[Figure 4] S-1 Bifactor Model — Fit Indices\n")
fit_idx <- c("cfi.robust", "tli.robust", "rmsea.robust", "rmsea.ci.lower.robust",
             "rmsea.ci.upper.robust", "srmr", "chisq.scaled", "df", "pvalue.scaled")
print(round(lavaan::fitMeasures(fit_bifactor, fit_idx), 3))

fit_table <- as.data.frame(as.list(round(lavaan::fitMeasures(fit_bifactor, fit_idx), 3)))
readr::write_csv(fit_table, file.path(output_dir, "Table_Figure4_BifactorFit.csv"))

std_loadings <- lavaan::standardizedSolution(fit_bifactor) %>% dplyr::filter(op == "=~")
readr::write_csv(std_loadings, file.path(output_dir, "Table_Figure4_StandardizedLoadings.csv"))

# Standardized loadings and fit indices underlying the Figure 4 path
# diagram are exported above.

# =============================================================================
# 13. APPENDIX — EDUCATION & SEX SENSITIVITY ANALYSES
#     (Supplementary Tables 5-6; toggle with RUN_APPENDIX_SENSITIVITY)
# =============================================================================
if (RUN_APPENDIX_SENSITIVITY) {
  
  # Demographics for the earlier healthy-control cohort are stored
  # separately from the main dataset and are merged in here.
  demographics_hc_old_file <- file.path(data_dir, "demographics_mCognito.ods")
  
  harmonise_demo <- function(var_main, recode_fun, old_file = demographics_hc_old_file) {
    main <- df %>% dplyr::select(SERIAL, dplyr::all_of(var_main)) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(var_main), as.character))
    if (file.exists(old_file) && requireNamespace("readODS", quietly = TRUE)) {
      old <- readODS::read_ods(old_file) %>%
        dplyr::select(SERIAL, dplyr::all_of(var_main)) %>%
        dplyr::mutate(SERIAL = as.character(SERIAL), dplyr::across(dplyr::all_of(var_main), as.character)) %>%
        dplyr::rename(.old = dplyr::all_of(var_main))
      main <- main %>% dplyr::left_join(old, by = "SERIAL") %>%
        dplyr::mutate(!!var_main := dplyr::coalesce(.data[[var_main]], .data[[".old"]])) %>%
        dplyr::select(-.old)
    }
    main %>% dplyr::mutate(label = recode_fun(.data[[var_main]]), Group = assign_group_simple(SERIAL))
  }
  
  # --- Education: binary Abitur contrast ---
  df_edu <- harmonise_demo("school", function(school) dplyr::case_when(
    school %in% c("1", "No degree")                              ~ "No degree",
    school %in% c("2", "Secondary (Hauptschule)")                 ~ "Hauptschule",
    school %in% c("3", "Intermediate (Realschule)")                ~ "Realschule",
    school %in% c("4", "University entrance (Abitur)")             ~ "Abitur",
    TRUE ~ NA_character_
  )) %>% dplyr::mutate(D_Abitur = dplyr::case_when(label == "Abitur" ~ 1L, label %in% c("No degree", "Hauptschule", "Realschule") ~ 0L, TRUE ~ NA_integer_))
  
  # --- Sex: binary Male contrast ---
  df_sex <- harmonise_demo("sex", function(sex) dplyr::case_when(
    tolower(sex) %in% c("m", "male", "1", "mannlich", "m\u00e4nnlich") ~ "Male",
    tolower(sex) %in% c("f", "w", "female", "2", "weiblich")          ~ "Female",
    TRUE ~ NA_character_
  )) %>% dplyr::mutate(D_Male = dplyr::case_when(label == "Male" ~ 1L, label == "Female" ~ 0L, TRUE ~ NA_integer_))
  
  # Both groups are entered jointly into the regression; residuals are
  # grand-mean re-centred to preserve the original theta scale.
  residualize_and_compare <- function(theta, group, dummy, task_label, score_label) {
    d <- tibble::tibble(theta = theta, Group = group, dummy = dummy, .row = seq_along(theta)) %>% dplyr::filter(!is.na(theta))
    fit_idx <- !is.na(d$dummy)
    grand_mean <- mean(d$theta[fit_idx], na.rm = TRUE)
    lm_fit <- lm(theta ~ dummy, data = d[fit_idx, ])
    d$theta_adj <- NA_real_
    d$theta_adj[fit_idx] <- residuals(lm_fit) + grand_mean
    
    gc_orig <- compare_groups(d$theta, d$Group)
    gc_adj  <- compare_groups(d$theta_adj, d$Group)
    if (is.null(gc_orig) || is.null(gc_adj)) return(tibble::tibble())
    
    tibble::tibble(
      Task = task_label, Score = score_label,
      N_hc = gc_orig$desc$N[gc_orig$desc$Group == "Healthy"], N_pat = gc_orig$desc$N[gc_orig$desc$Group == "Patient"],
      g_unadjusted = gc_orig$g, g_unadjusted_p = gc_orig$p,
      g_adjusted   = gc_adj$g,  g_adjusted_p   = gc_adj$p,
      delta_g      = gc_adj$g - gc_orig$g
    )
  }
  
  appendix_tasks <- list(
    list(theta = master_theta$theta_SART_EI,       label = "SART-ED (external interference)"),
    list(theta = master_theta$theta_CI_resilience,  label = "CI (internal interference)"),
    list(theta = master_theta$theta_SART_MS,        label = "SART-ED (mental speed)"),
    list(theta = master_theta$theta_GRT,            label = "GRT (reasoning)"),
    list(theta = master_theta$theta_DRT_D1,         label = "DRT-D1 (familiarity-based recognition)"),
    list(theta = master_theta$theta_DRT_D2,         label = "DRT-D2 (recollection-based discrimination)"),
    list(theta = master_theta$theta_PAL,            label = "PAL (associative memory)")
  )
  
  edu_dummy <- df_edu$D_Abitur[match(master_theta$SERIAL, df_edu$SERIAL)]
  sex_dummy <- df_sex$D_Male[match(master_theta$SERIAL, df_sex$SERIAL)]
  
  suppl_table5_education <- purrr::map_dfr(appendix_tasks, ~ residualize_and_compare(
    .x$theta, master_theta$Group, edu_dummy, .x$label, "theta"))
  suppl_table6_sex <- purrr::map_dfr(appendix_tasks, ~ residualize_and_compare(
    .x$theta, master_theta$Group, sex_dummy, .x$label, "theta"))
  
  readr::write_csv(suppl_table5_education, file.path(output_dir, "SupplTable5_Education_Sensitivity.csv"))
  readr::write_csv(suppl_table6_sex,       file.path(output_dir, "SupplTable6_Sex_Sensitivity.csv"))
  
  cat("\n[Appendix] Education and sex sensitivity tables written to output_dir.\n")
}

# =============================================================================
# END OF SCRIPT
# Outputs written to output_dir:
#   EMA4Stroke_Master_Theta_Scores.csv       - per-participant theta, all domains
#   Figure2_Group_Comparison_Raincloud.png   - Figure 2
#   Figure3_Correlation_Reliability_Matrix.png / Table_ThetaCorrelations.csv - Figure 3
#   Table_Figure4_BifactorFit.csv / Table_Figure4_StandardizedLoadings.csv  - Figure 4
#   SupplTable5_Education_Sensitivity.csv, SupplTable6_Sex_Sensitivity.csv  - appendix
# =============================================================================