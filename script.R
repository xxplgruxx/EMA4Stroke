# ─────────────────────────────────────────────────────────────────────────────
# IRT Analysis for GRT (General Recognition Test) – EMA4Stroke
# Clean, robust template with Rasch-based item selection using TAM::tam.fit
# Supports planned missingness (matrix sampling) without imputation.
# Outputs: theta CSV, kept/dropped item lists, item-fit table, APA-style figures,
# group comparison, and README logging.
# ─────────────────────────────────────────────────────────────────────────────

# --- Step 0: Library Section (install + load; UTF-8 throughout) ---
packages <- c(
  "readr","dplyr","tidyr","ggplot2","mirt","psych","GPArotation",
  "effsize","arrow","ggrepel","patchwork","stringr","purrr","tibble","Gifi",
  "TAM",      # for Rasch item fit under missingness
  "ggdist","gghalves"  # for half-violin rainclouds
)
new_pkgs <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

set.seed(20250821)  # reproducibility for EM starts/optimizers

# --- Step 0a: APA plotting theme & project palette (Arial + specified colors) ---
apa_colors <- c(
  rgb(0,78,159, maxColorValue = 255),   # Pantone 2945
  rgb(0,62,107, maxColorValue = 255),   # Pantone 7693
  rgb(0,120,120, maxColorValue = 255),  # Pantone 3295
  rgb(213,59,10, maxColorValue = 255),  # Pantone 1665
  rgb(0,171,217, maxColorValue = 255),  # Pantone 2995
  rgb(0,168,121, maxColorValue = 255),  # Pantone 339
  rgb(238,113,0, maxColorValue = 255),  # Pantone 152
  rgb(91,197,242, maxColorValue = 255), # Pantone 298
  rgb(148,193,28, maxColorValue = 255), # Pantone 376
  rgb(243,145,0, maxColorValue = 255),  # Pantone 144
  rgb(161,217,248, maxColorValue = 255),# Pantone 291
  rgb(199,211,0, maxColorValue = 255),  # Pantone 382
  rgb(253,195,0, maxColorValue = 255)   # Pantone 7408
)

theme_apa <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.minor = element_blank()
  )

# --- Step 1: Config & Data ---
data_folder <- "\\\\smb.uni-oldenburg.de\\psychologie$\\PMuS\\GraKo\\Grupe\\GRAKO Stroke & Neuromodulation\\Analyse\\Test Sets\\"
input_file  <- file.path(data_folder, "Final_Merged_Dataset.feather")   # Feather (UTF-8)
output_tag  <- "GRT"
test_prefix <- "^GRT"   # adjust when reusing template for other tests

# Load
df <- arrow::read_feather(input_file)

# Identify groups from SERIAL
patients <- dplyr::filter(df, grepl("^[A-Za-z]{2}\\d{2}$", SERIAL))
healthy  <- dplyr::filter(df, !grepl("^[A-Za-z]{2}\\d{2}$", SERIAL))

# Keep only columns that look like items (start with prefix) and are binary 0/1 (or NA)
all_cols <- names(df)
item_cols_raw <- grep(test_prefix, all_cols, value = TRUE)
df[item_cols_raw] <- lapply(df[item_cols_raw], function(v) suppressWarnings(as.numeric(v)))
is_binary_col <- function(x) {
  ux <- unique(na.omit(x))
  all(ux %in% c(0,1)) && length(ux) > 0
}
item_cols <- item_cols_raw[vapply(df[item_cols_raw], is_binary_col, logical(1))]

# Items observed in at least one respondent per group (planned missingness aware)
get_union_items <- function(d1, d2, cols) {
  i1 <- cols[colMeans(!is.na(d1[, cols, drop = FALSE])) > 0]
  i2 <- cols[colMeans(!is.na(d2[, cols, drop = FALSE])) > 0]
  union(i1, i2)
}
common_items <- get_union_items(patients, healthy, item_cols)

combined <- bind_rows(
  select(patients, SERIAL, all_of(common_items)),
  select(healthy,  SERIAL, all_of(common_items))
) %>% arrange(SERIAL)

# Remove zero-variance items
item_vars <- sapply(combined[, common_items, drop = FALSE], var, na.rm = TRUE)
keep_var  <- names(item_vars[item_vars > 0])
X <- dplyr::select(combined, SERIAL, all_of(keep_var))
Y <- dplyr::select(X, -SERIAL)

# ─────────────────────────────────────────────────────────────────────────────
# NEW: Stage 1 — Response-rate filter (≥ 45 non-missing respondents)
# ─────────────────────────────────────────────────────────────────────────────
min_respondents <- 45L
item_nobs <- colSums(!is.na(Y))
stage1_keep <- names(item_nobs[item_nobs >= min_respondents])
stage1_drop <- setdiff(colnames(Y), stage1_keep)

cat("\n[Stage 1] Response-rate filter (>= ", min_respondents, "):\n", sep = "")
cat("Kept (", length(stage1_keep), "): ", paste(stage1_keep, collapse = ", "), "\n", sep = "")
if (length(stage1_drop)) {
  cat("Dropped (", length(stage1_drop), "): ", paste(stage1_drop, collapse = ", "), "\n", sep = "")
}

# Restrict Y to Stage 1 kept items for downstream analyses
Y1 <- Y[, stage1_keep, drop = FALSE]

# --- Tetrachoric-based dimensionality screen (optional; unchanged logic) ---
if (ncol(Y1) >= 3) {
  message("Running tetrachoric correlations (Stage 1 set)…")
  tetra <- psych::tetrachoric(as.matrix(Y1))
  Rtet  <- tetra$rho
  Rtet[!is.finite(Rtet)] <- 0
  Rtet[is.na(Rtet)] <- 0
  suppressWarnings({
    fa.parallel(Rtet, fa = "pc", n.iter = 100,
                main = paste0(output_tag, ": Parallel Analysis (tetrachoric PC; Stage 1 set)"))
  })
} else {
  message("< 3 items after Stage 1; skipping dimensionality screen.")
}

# --- Princals biplot (items only, with labels) on Stage 1 set ---
Y1_fac <- as.data.frame(lapply(Y1, function(x) factor(x, levels = c(0,1))))
ok_cols <- vapply(Y1_fac, function(f) length(unique(stats::na.omit(f))) >= 2, logical(1))
Y1_fac <- Y1_fac[ok_cols]

if (ncol(Y1_fac) >= 3) {
  pc_out <- Gifi::princals(Y1_fac, ndim = 2, ordinal = TRUE)
  os_name <- names(pc_out)[grepl("object", names(pc_out), ignore.case = TRUE)][1]
  if (is.na(os_name)) stop("princals output lacked object scores.")
  S <- as.data.frame(pc_out[[os_name]])
  colnames(S) <- c("Dim1","Dim2")
  
  # Item vectors = correlations with princals dimensions
  Y1_num <- as.data.frame(lapply(Y1_fac, function(f) as.numeric(as.character(f))))
  item_vec <- purrr::map_dfr(names(Y1_num), function(j) {
    x <- Y1_num[[j]]
    tibble::tibble(
      item = j,
      r1 = suppressWarnings(cor(x, S$Dim1, use = "pairwise.complete.obs")),
      r2 = suppressWarnings(cor(x, S$Dim2, use = "pairwise.complete.obs"))
    )
  })
  item_vec[is.na(item_vec)] <- 0
  arrows_df <- item_vec %>% mutate(x0 = 0, y0 = 0, x1 = r1, y1 = r2)
  
  p_princals <- ggplot(arrows_df) +
    geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1),
                 arrow = arrow(length = unit(0.15, "cm")),
                 linewidth = 0.5, color = "black") +
    ggrepel::geom_text_repel(
      aes(x = x1, y = y1, label = item),
      size = 3, seed = 20250821,
      min.segment.length = 0, segment.size = 0.2
    ) +
    labs(title = paste0(output_tag, ": Princals item biplot (Stage 1 set)"),
         subtitle = "Item vectors (correlations with princals dimensions)",
         x = "Dimension 1", y = "Dimension 2") +
    theme_apa
  print(p_princals)
  ggsave(file.path(data_folder, paste0(output_tag, "_Princals_ItemsOnly_Stage1.png")),
         p_princals, width = 7, height = 7, dpi = 300)
} else {
  message("Princals skipped: fewer than 3 usable items after Stage 1.")
}

# ─────────────────────────────────────────────────────────────────────────────
# NEW: Stage 2 — Manual exclusions after looking at the biplot
# Fill this vector with item names you wish to drop based on the biplot.
# Example: manual_drop <- c("GRT03","GRT12")
# Leave as character(0) if none.
# ─────────────────────────────────────────────────────────────────────────────
manual_drop <- c("GRT30", "GRT10", "GRT24", "GRT13")  # <-- EDIT HERE after inspecting the saved biplot
manual_drop <- intersect(manual_drop, colnames(Y1))  # safety

stage2_keep <- setdiff(colnames(Y1), manual_drop)
stage2_drop <- manual_drop

cat("\n[Stage 2] Manual exclusions after biplot:\n")
cat("Kept (", length(stage2_keep), "): ", paste(stage2_keep, collapse = ", "), "\n", sep = "")
if (length(stage2_drop)) {
  cat("Dropped (", length(stage2_drop), "): ", paste(stage2_drop, collapse = ", "), "\n", sep = "")
}

# Prepare matrix for fitting based on Stage 2 kept items
Y2 <- Y1[, stage2_keep, drop = FALSE]

# Remove persons with all-NA before model fitting
Y2_idx <- rowSums(!is.na(Y2)) > 0
if (!any(Y2_idx)) stop("All rows are all-NA after Stage 2 item selection.")
Y2_fit <- Y2[Y2_idx, , drop = FALSE]

# --- Optional model comparison (Rasch vs 2PL) on Stage 2 set ---
mod_rasch_s2 <- mirt::mirt(Y2_fit, 1, itemtype = "Rasch", SE = FALSE, verbose = FALSE)
mod_2pl_s2   <- mirt::mirt(Y2_fit, 1, itemtype = "2PL",   SE = FALSE, verbose = FALSE)
mod_2pl_s2   <- tryCatch(mirt::optim.mirt(mod_2pl_s2, method = "BFGS"), error = function(e) mod_2pl_s2)

compare_models <- function(m1, m2, label1 = "Rasch", label2 = "2PL") {
  safe_num <- function(x) if (inherits(x, "try-error") || is.null(x)) NA_real_ else as.numeric(x)[1]
  out <- tibble::tibble(
    Model  = c(label1, label2),
    AIC    = c(safe_num(AIC(m1)), safe_num(AIC(m2))),
    BIC    = c(safe_num(BIC(m1)), safe_num(BIC(m2))),
    LogLik = c(safe_num(logLik(m1)), safe_num(logLik(m2)))
  )
  print(out)
  cat("\nLikelihood Ratio Test (", label2, " vs ", label1, "):\n", sep = "")
  print(tryCatch(anova(m1, m2), error = function(e) e$message))
  invisible(out)
}
cat("\nModel comparison on Stage 2 set:\n")
invisible(compare_models(mod_rasch_s2, mod_2pl_s2, "Rasch (S2)", "2PL (S2)"))

# ─────────────────────────────────────────────────────────────────────────────
# NEW: Stage 3 — Rasch item fit (TAM) with Holm correction on Stage 2 kept items
# ─────────────────────────────────────────────────────────────────────────────
tam_mod <- TAM::tam.mml(resp = as.matrix(Y2_fit), irtmodel = "1PL")
tam_fit <- TAM::tam.fit(tam_mod)
fit_tab <- as.data.frame(tam_fit$itemfit)

# Normalize to a unified p column "p.S_X2"
cn <- names(fit_tab)
if (all(c("Infit_p","Outfit_p") %in% cn)) {
  fit_tab$p.S_X2 <- pmin(fit_tab$Infit_p, fit_tab$Outfit_p, na.rm = TRUE)
} else {
  p_candidates <- intersect(cn, c("p.S_X2","p.X2","p_value","p","pX2","pv","p.value"))
  if (length(p_candidates) >= 1) {
    names(fit_tab)[match(p_candidates[1], names(fit_tab))] <- "p.S_X2"
  } else {
    x2_candidates <- intersect(cn, c("X2","X2.item","stat","Statistic","chi2","Chi2","AX2"))
    df_candidates <- intersect(cn, c("df","df.item","DF","df.X2","AX2.df"))
    if (length(x2_candidates) >= 1 && length(df_candidates) >= 1) {
      fit_tab$"p.S_X2" <- 1 - pchisq(as.numeric(fit_tab[[x2_candidates[1]]]),
                                     df = as.numeric(fit_tab[[df_candidates[1]]]))
    } else stop("Could not locate p or (X2, df) in TAM item-fit output.")
  }
}
if (!"item" %in% names(fit_tab)) {
  item_candidates <- intersect(names(fit_tab), c("item","Item","parameter","par","name"))
  fit_tab$item <- if (length(item_candidates)) fit_tab[[item_candidates[1]]] else rownames(fit_tab)
}

alpha <- 0.05
fit_tab$padj_holm <- p.adjust(fit_tab$"p.S_X2", method = "holm")

stage3_keep <- fit_tab %>% dplyr::filter(padj_holm >= alpha) %>% dplyr::pull(item)
stage3_keep <- intersect(stage3_keep, colnames(Y2_fit))
stage3_drop <- setdiff(colnames(Y2_fit), stage3_keep)

cat("\n[Stage 3] Rasch item fit (Holm-corrected, α=.05) on Stage 2 kept items:\n")
print(fit_tab[order(fit_tab$padj_holm), c("item","p.S_X2","padj_holm")])
cat("Kept (", length(stage3_keep), "): ", paste(stage3_keep, collapse = ", "), "\n", sep = "")
if (length(stage3_drop)) {
  cat("Dropped (", length(stage3_drop), "): ", paste(stage3_drop, collapse = ", "), "\n", sep = "")
}

# Combined decisions + reason tag (for plotting/report)
decisions <- bind_rows(
  tibble(item = stage1_drop, reason = "Dropped: <45 respondents"),
  tibble(item = stage2_drop, reason = "Dropped: manual (biplot)"),
  tibble(item = stage3_drop, reason = "Dropped: Rasch misfit (Holm)"),
  tibble(item = stage3_keep, reason = "Kept: all stages passed")
) %>%
  distinct(item, .keep_all = TRUE)

# Save decision log
readr::write_csv(decisions, file.path(data_folder, paste0(output_tag, "_Item_Decisions.csv")))

# --- Step 6: Refit Rasch on final retained items (Stage 3 keep) ---
Y_final <- Y2_fit[, stage3_keep, drop = FALSE]
if (ncol(Y_final) < 3) stop("Fewer than 3 items retained after Stage 3; aborting theta estimation.")
mod_rasch_final <- mirt::mirt(Y_final, 1, itemtype = "Rasch", SE = FALSE, verbose = FALSE)
mod_2pl_final   <- mirt::mirt(Y_final, 1, itemtype = "2PL",   SE = FALSE, verbose = FALSE)
mod_2pl_final   <- tryCatch(mirt::optim.mirt(mod_2pl_final, method = "BFGS"), error = function(e) mod_2pl_final)

cat("\nModel comparison on FINAL item set:\n")
invisible(compare_models(mod_rasch_final, mod_2pl_final, "Rasch (final)", "2PL (final)"))

# --- Step 7: ICCs for retained items (APA style) ---
plot_iccs <- function(model, title = "Item Characteristic Curves") {
  pars <- coef(model, IRTpars = TRUE, simplify = TRUE)$items
  theta_seq <- seq(-4, 4, length.out = 400)
  icc_df <- purrr::map_dfr(rownames(pars), function(it) {
    a <- pars[it, if ("a1" %in% colnames(pars)) "a1" else "a"]
    b <- pars[it, "b"]
    tibble::tibble(Item = it, theta = theta_seq, P = plogis(a * (theta_seq - b)))
  })
  ggplot(icc_df, aes(theta, P, group = Item)) +
    geom_line(linewidth = 0.8, alpha = 0.9, color = apa_colors[1]) +
    labs(title = title, x = expression(theta), y = "P(correct)") +
    theme_apa
}
p_icc <- plot_iccs(mod_rasch_final, title = paste0(output_tag, ": ICCs (Final Retained Items)"))
print(p_icc)

# --- Step 8: Ability estimates (EAP) for all rows (respect planned missingness) ---
theta_all <- rep(NA_real_, nrow(Y))
# Map: only rows used in Stage 2 fitting are eligible (Y2_idx)
# and among those, only rows with at least one final item observed are scored
has_final_obs <- rowSums(!is.na(Y2[ , stage3_keep, drop = FALSE])) > 0
scored <- mirt::fscores(mod_rasch_final, method = "EAP")[,1]
theta_all[Y2_idx & has_final_obs] <- scored

scores <- tibble::tibble(
  SERIAL = X$SERIAL,
  !!paste0("theta_", tolower(output_tag)) := theta_all
)
readr::write_csv(scores, file.path(data_folder, paste0("EMA4Stroke_Theta_", output_tag, ".csv")))

# --- Step 9: Theta distributions & group comparison (APA style) ---
scores2 <- scores %>%
  mutate(Group = ifelse(grepl("^[A-Za-z]{2}\\d{2}$", SERIAL), "Patient", "Healthy")) %>%
  filter(!is.na(.data[[paste0("theta_", tolower(output_tag))]]))

# Density
p_theta <- ggplot(scores2, aes(x = .data[[paste0("theta_", tolower(output_tag))]], fill = Group, color = Group)) +
  geom_density(alpha = 0.35, adjust = 1.0) +
  scale_fill_manual(values = c("Healthy" = apa_colors[1], "Patient" = apa_colors[4])) +
  scale_color_manual(values = c("Healthy" = apa_colors[1], "Patient" = apa_colors[4])) +
  labs(title = paste0(output_tag, ": Ability Distributions by Group"),
       x = expression(theta), y = "Density") +
  theme_apa
print(p_theta)

# Welch t-test + effect size
if (nrow(scores2) > 1) {
  stats <- scores2 %>%
    group_by(Group) %>%
    summarize(M = mean(.data[[paste0("theta_", tolower(output_tag))]], na.rm = TRUE),
              SD = sd(.data[[paste0("theta_", tolower(output_tag))]],  na.rm = TRUE),
              N = dplyr::n(), .groups = "drop")
  print(stats)
  t_out <- tryCatch(
    t.test(stats::reformulate("Group", response = paste0("theta_", tolower(output_tag))), data = scores2),
    error = function(e) NULL
  )
  if (!is.null(t_out)) print(t_out)
  
  cat("\n== Effect size (Cohen's d, Hedges correction) ==\n")
  d_out <- effsize::cohen.d(
    stats::reformulate("Group", response = paste0("theta_", tolower(output_tag))),
    data = scores2,
    hedges.correction = TRUE, conf.level = 0.95
  )
  print(d_out)
}

# --- Optional: "Manhattan"-style item-fit view reflecting final decision reasons ---
if ("p.S_X2" %in% names(fit_tab)) {
  df_manh_raw <- fit_tab %>%
    mutate(item = as.character(item)) %>%
    left_join(decisions, by = "item") %>%
    mutate(reason = ifelse(is.na(reason), "Kept/Dropped earlier", reason)) %>%
    arrange(p.S_X2) %>%
    mutate(rank = row_number())
  
  m <- nrow(df_manh_raw)
  alpha <- 0.05
  holm_step <- tibble::tibble(rank = 1:m, crit = alpha / (m - rank + 1))
  
  p_manh_step <- ggplot(df_manh_raw, aes(x = rank, y = p.S_X2, color = reason)) +
    geom_step(data = holm_step, aes(x = rank, y = crit), inherit.aes = FALSE,
              linetype = "dashed", linewidth = 0.6) +
    geom_point(size = 2.0) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(df_manh_raw, grepl("^Dropped", reason)),
      aes(label = item), size = 3, seed = 20250821,
      min.segment.length = 0, segment.size = 0.2
    ) +
    scale_color_manual(values = c(
      "Kept: all stages passed"           = apa_colors[1],
      "Dropped: <45 respondents"          = apa_colors[11],
      "Dropped: manual (biplot)"          = apa_colors[7],
      "Dropped: Rasch misfit (Holm)"      = apa_colors[4],
      "Kept/Dropped earlier"              = "gray50"
    ), name = "Decision reason") +
    labs(
      title = paste0(output_tag, ": Item-fit raw p with Holm step boundary (Stage 3)"),
      subtitle = "Dashed step = Holm critical value α/(m−i+1); annotated = dropped",
      x = "Item rank (by raw p)", y = "Raw item-fit p-value (TAM)"
    ) +
    theme_apa + theme(legend.position = "bottom")
  
  print(p_manh_step)
  ggsave(file.path(data_folder, paste0(output_tag, "_ItemFit_Manhattan_RawP_HolmBoundary_FINAL.png")),
         p_manh_step, width = 10, height = 6, dpi = 300)
}

# ─────────────────────────────────────────────────────────────────────────────
# Minimal sample characteristics (optional, unchanged)
# ─────────────────────────────────────────────────────────────────────────────
df$group <- ifelse(grepl("^ZZ", df$SERIAL) | grepl("^[A-Za-z]{2}\\d{2}$", df$SERIAL),
                   "Patient", "Healthy")
cat("\n== Counts by group ==\n"); print(table(df$group)); cat("\n")

# Helper: columns with any entries for healthy (quick check)
cols_healthy_nonmissing <- colnames(df)[colSums(!is.na(df[df$group == "Healthy", ])) > 0]
# print(cols_healthy_nonmissing)

# ─────────────────────────────────────────────────────────────────────────────
# Notes:
# • Stage 1 threshold is set by `min_respondents` (default 45).
# • Stage 2 manual list: edit `manual_drop` after you review the saved biplot.
# • Stage 3 Holm adjustment applied only to Stage 2 kept items, as requested.
# • Decision log written to: *_Item_Decisions.csv
# • Theta saved to: EMA4Stroke_Theta_GRT.csv
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# APPEND-ONLY PLOTTING CODE (Manhattan: Holm-corrected; Princals: manual drops;
# Raincloud: group comparison). Assumes `theme_apa`, `apa_colors` already exist.
# ─────────────────────────────────────────────────────────────────────────────

# === 1) Manhattan plot using Holm-corrected item-fit (fixed α = .05 line) ===
# Requires: fit_tab with columns: item, padj_holm (produced in Stage 3)
stopifnot(all(c("item","padj_holm") %in% names(fit_tab)))
df_manh_holm <- fit_tab |>
  dplyr::mutate(
    item   = as.character(item),
    status = ifelse(padj_holm < 0.05, "Dropped (Holm < .05)", "Kept (Holm ≥ .05)")
  ) |>
  dplyr::arrange(padj_holm) |>
  dplyr::mutate(rank = dplyr::row_number())

p_manh_holm <- ggplot(df_manh_holm, aes(x = rank, y = padj_holm, color = status)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.7) +
  geom_point(size = 2.1) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(df_manh_holm, padj_holm < 0.05),
    aes(label = item), size = 3, seed = 20250821,
    min.segment.length = 0, segment.size = 0.2, show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Kept (Holm ≥ .05)"     = apa_colors[1],
    "Dropped (Holm < .05)"  = apa_colors[4]
  ), name = "Decision") +
  labs(
    title = paste0(output_tag, ": Rasch item fit (Holm-corrected p)"),
    subtitle = "Dashed horizontal line = α = .05",
    x = "Item rank (by Holm-adjusted p)", y = "Holm-adjusted p-value"
  ) +
  theme_apa + theme(legend.position = "bottom")

print(p_manh_holm)
ggsave(file.path(data_folder, paste0(output_tag, "_ItemFit_Manhattan_HolmFixed05.png")),
       p_manh_holm, width = 9, height = 5.5, dpi = 300)

# === 2) Princals biplot with manually dropped items colored in red ===========
# Requires: Stage 1 dataset `Y1_fac` and object scores `S` (from earlier),
# and `manual_drop` character vector. If not present, rebuild arrows_df.

if (!exists("manual_drop")) manual_drop <- character(0)

# Rebuild arrows_df defensively if needed
if (!exists("arrows_df") || !is.data.frame(arrows_df)) {
  if (!exists("Y1_fac") || !exists("S")) {
    # Attempt to rebuild from Stage 1 inputs
    # (Assumes Y1 is available; else stop with an informative message)
    if (!exists("Y1")) stop("Princals objects not found. Ensure Stage 1 (Y1/Y1_fac/S) ran before this block.")
    Y1_fac <- as.data.frame(lapply(Y1, function(x) factor(x, levels = c(0,1))))
    ok_cols <- vapply(Y1_fac, function(f) length(unique(stats::na.omit(f))) >= 2, logical(1))
    Y1_fac <- Y1_fac[ok_cols]
    stopifnot(ncol(Y1_fac) >= 3)
    pc_out <- Gifi::princals(Y1_fac, ndim = 2, ordinal = TRUE)
    os_name <- names(pc_out)[grepl("object", names(pc_out), ignore.case = TRUE)][1]
    stopifnot(!is.na(os_name))
    S <- as.data.frame(pc_out[[os_name]])
    colnames(S) <- c("Dim1","Dim2")
  }
  Y1_num <- as.data.frame(lapply(Y1_fac, function(f) as.numeric(as.character(f))))
  item_vec <- purrr::map_dfr(names(Y1_num), function(j) {
    x <- Y1_num[[j]]
    tibble::tibble(
      item = j,
      r1 = suppressWarnings(cor(x, S$Dim1, use = "pairwise.complete.obs")),
      r2 = suppressWarnings(cor(x, S$Dim2, use = "pairwise.complete.obs"))
    )
  })
  item_vec[is.na(item_vec)] <- 0
  arrows_df <- item_vec %>% dplyr::mutate(x0 = 0, y0 = 0, x1 = r1, y1 = r2)
}

arrows_df$flag <- ifelse(arrows_df$item %in% manual_drop, "Dropped (manual)", "Kept")

p_princals_flag <- ggplot(arrows_df) +
  geom_segment(
    aes(x = x0, y = y0, xend = x1, yend = y1, color = flag),
    arrow = arrow(length = unit(0.15, "cm")),
    linewidth = 0.6
  ) +
  ggrepel::geom_text_repel(
    aes(x = x1, y = y1, label = item, color = flag),
    size = 3, seed = 20250821,
    min.segment.length = 0, segment.size = 0.2
  ) +
  scale_color_manual(values = c(
    "Kept"             = apa_colors[1],  # blue
    "Dropped (manual)" = apa_colors[4]   # red/orange
  ), name = "Stage 2 flag") +
  labs(
    title = paste0(output_tag, ": Princals item biplot (Stage 1 set)"),
    subtitle = "Manually dropped items highlighted in red",
    x = "Dimension 1", y = "Dimension 2"
  ) +
  theme_apa + theme(legend.position = "bottom")

print(p_princals_flag)
ggsave(file.path(data_folder, paste0(output_tag, "_Princals_ItemsOnly_Stage1_ManualFlag.png")),
       p_princals_flag, width = 7, height = 7, dpi = 300)

# === 3) Raincloud plot for group comparison (Healthy vs Patient) =============
# Requires: `scores2` with columns Group ∈ {Healthy, Patient} and theta column.
theta_col <- paste0("theta_", tolower(output_tag))
stopifnot(all(c("Group", theta_col) %in% names(scores2)))

scores2$Group <- factor(scores2$Group, levels = c("Healthy","Patient"))
set.seed(20250821)

p_rain <- ggplot(scores2, aes(x = Group, y = .data[[theta_col]], fill = Group, color = Group)) +
  gghalves::geom_half_violin(
    side = "l", alpha = 0.6, color = "black", linewidth = 0.4, width = 0.9, trim = FALSE
  ) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.4) +
  gghalves::geom_half_point(
    side = "r", range_scale = 0.4, alpha = 0.55, size = 1.3,
    position = position_jitter(width = 0.07, height = 0)
  ) +
  scale_fill_manual(values = c("Healthy" = apa_colors[1], "Patient" = apa_colors[4])) +
  scale_color_manual(values = c("Healthy" = apa_colors[1], "Patient" = apa_colors[4])) +
  labs(title = paste0(output_tag, ": Raincloud Theta Plot by Group (Half-violin)"),
       x = NULL, y = expression(theta)) +
  theme_apa +
  theme(legend.position = "none")

print(p_rain)
ggsave(file.path(data_folder, paste0(output_tag, "_Raincloud_Theta_ByGroup.png")),
       p_rain, width = 7.5, height = 5.5, dpi = 300)
