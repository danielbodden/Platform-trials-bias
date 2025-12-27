#!/usr/bin/env Rscript
# ============================================================
#  Nonconcurrent controls – plotting (PDF only)
#
#  Fixes included:
#   1) Y-scale:
#      - Type I error is ALWAYS 0..0.2 (all columns)
#      - Power is ALWAYS 0..1.0 (all columns)
#      Achieved via facet_grid(scales="free_y") + geom_blank limits.
#
#   2) Chronological-bias reference detection:
#      - NEVER uses eta=0 (since chronobias often has eta=0 everywhere)
#      - Uses explicit "Reference/no bias" strings OR beta_time==0
#
#   3) Patchwork stacking:
#      - Robust: plot_alloc_2x3_save() returns list(plot=..., file=...)
#      - wrap_plots() always receives ggplot objects
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) suppressPackageStartupMessages(library(patchwork))

# ---- Style sheet ----
source("PT_bias/plot_style_bimj.R")

# Ensure the reference entry is labeled consistently (even if the style sheet is older).
if (exists("BIMJ_RAND_LABELS")) {
  BIMJ_RAND_LABELS$PBR0 <- expression("Reference (no bias)")
}

in_dir  <- "PT_bias/results_nonconcurrent"
out_dir <- "plots_nonconc_both"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

BGRID    <- as.integer(seq(1, 49, by = 12))
X_LAB_TB <- expression(paste("Arm B start ", t[B]))

safe_save_pdf <- function(plot, filename, width, height) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  force(plot)
  
  dev_open <- FALSE
  tryCatch({
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(filename, width = width, height = height, onefile = TRUE)
    } else {
      grDevices::pdf(filename, width = width, height = height, onefile = TRUE)
    }
    dev_open <- TRUE
    print(plot)
  }, finally = {
    if (dev_open && grDevices::dev.cur() != 1) grDevices::dev.off()
  })
  
  message("Written: ", normalizePath(filename, mustWork = FALSE))
  invisible(filename)
}

parse_int <- function(x) suppressWarnings(as.integer(as.character(x)))

get_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

ensure_series_col <- function(df) {
  if (!("series" %in% names(df)) && ("ser" %in% names(df))) df$series <- df$ser
  df
}

measure_factor <- function(x) {
  factor(as.character(x), levels = c("Type I error", "Power"))
}

series_factor_parsed <- function(x) {
  x <- as.character(x)
  x <- dplyr::case_when(
    x %in% c("t1e_A", "A", "armA", "ArmA") ~ "A",
    x %in% c("t1e_B", "B", "armB", "ArmB") ~ "B",
    x %in% c("fwer", "FWER", "P_reject")   ~ "FWER",
    TRUE ~ x
  )
  bimj_series_factor(x)
}

x_scale_tb <- function() {
  scale_x_continuous(
    breaks = BGRID,
    labels = as.character(BGRID),
    limits = c(1, 49),
    expand = expansion(add = c(2, 2))
  )
}

legend_layout_theme <- function() {
  theme(
    legend.box       = "vertical",
    legend.direction = "horizontal",
    legend.box.just  = "left",
    legend.position  = "bottom"
  )
}

# ------------------------------------------------------------
# Reference extraction from allocation-bias files:
# config == "PBR (eta=0)"
# ------------------------------------------------------------
extract_ref_no_bias_alloc <- function(in_dir) {
  
  t1e_path <- file.path(in_dir, "allocbias_nonconc_t1e_preferB_bothModels.csv")
  pow_path <- file.path(in_dir, "allocbias_nonconc_power_preferB_bothModels.csv")
  
  if (!file.exists(t1e_path) || !file.exists(pow_path)) {
    stop("Missing allocation-bias files for reference extraction: ", t1e_path, " / ", pow_path)
  }
  
  t1e_wide <- read.csv(t1e_path, check.names = FALSE)
  cfg_col  <- get_col(t1e_wide, c("config", "cfg"))
  if (is.na(cfg_col)) stop("T1E file has no config column (expected 'config').")
  
  ref_t1e <- t1e_wide %>%
    filter(.data[[cfg_col]] == "PBR (eta=0)") %>%
    mutate(b_start = parse_int(b_start)) %>%
    pivot_longer(cols = c(t1e_A, t1e_B, fwer), names_to = "tmp_series", values_to = "rate") %>%
    mutate(
      series  = recode(tmp_series, "t1e_A" = "A", "t1e_B" = "B", "fwer" = "FWER"),
      measure = "Type I error"
    ) %>%
    select(b_start, series, analysis_model, measure, rate)
  
  pow_long <- read.csv(pow_path, check.names = FALSE) %>% ensure_series_col()
  cfg_col2 <- get_col(pow_long, c("config", "cfg"))
  if (is.na(cfg_col2)) stop("Power file has no config column (expected 'config').")
  
  ref_pow <- pow_long %>%
    filter(.data[[cfg_col2]] == "PBR (eta=0)") %>%
    mutate(
      b_start = parse_int(b_start),
      measure = "Power"
    ) %>%
    select(b_start, series, analysis_model, measure, rate)
  
  ref <- bind_rows(ref_t1e, ref_pow) %>%
    filter(!is.na(b_start), !is.na(rate))
  
  if (nrow(ref) == 0) {
    stop("Reference extraction returned 0 rows for config == 'PBR (eta=0)'.")
  }
  
  ref %>% mutate(panel = "Reference (no bias)")
}

ref_no_bias <- extract_ref_no_bias_alloc(in_dir)

# ------------------------------------------------------------
# Common factor enrichment
# ------------------------------------------------------------
add_common_factors <- function(dat) {
  if (!("panel" %in% names(dat))) dat$panel <- NA_character_
  
  dat <- dat %>%
    mutate(
      b_start = parse_int(b_start),
      panel   = as.character(panel)
    )
  
  rk <- suppressWarnings(as.character(bimj_rand_key(dat$panel)))
  
  is_ref <- grepl("reference|no[[:space:]]*bias", dat$panel, ignore.case = TRUE) |
    dat$panel %in% c("Reference (no bias)", "PBR0")
  
  rk[is.na(rk) & is_ref] <- "PBR0"
  rk[dat$panel %in% c("Reference (no bias)", "PBR0")] <- "PBR0"
  
  if (exists("BIMJ_RAND_LEVELS")) {
    rand_key <- factor(rk, levels = BIMJ_RAND_LEVELS)
  } else {
    rand_key <- factor(rk)
  }
  
  dat %>%
    mutate(
      rand_key           = rand_key,
      analysis_model_fct = bimj_model_factor(analysis_model),
      measure_fct        = measure_factor(measure),
      series_fct         = series_factor_parsed(series)
    ) %>%
    filter(!is.na(b_start), !is.na(rate), !is.na(rand_key), !is.na(analysis_model_fct))
}

# ------------------------------------------------------------
# Core 2x3 plot (T1E row + Power row) with reference overlay
# ------------------------------------------------------------
make_2x3_plot <- function(main_dat, ref_dat, title) {
  
  plot_dat <- bind_rows(main_dat, ref_dat) %>% distinct()
  if (nrow(plot_dat) == 0) stop("make_2x3_plot(): plot_dat has 0 rows.")
  
  series_vals <- if (is.factor(plot_dat$series_fct)) levels(plot_dat$series_fct) else unique(plot_dat$series_fct)
  meas_lvls   <- if (is.factor(plot_dat$measure_fct)) levels(plot_dat$measure_fct) else unique(plot_dat$measure_fct)
  
  # Force y-range for each row:
  # - Type I error: 0..0.2
  # - Power:       0..1.0
  blank_t1e <- tidyr::expand_grid(
    b_start     = 1L,
    series_fct  = series_vals,
    measure_fct = factor("Type I error", levels = meas_lvls),
    rate        = c(0, 0.2)
  )
  
  blank_pow <- tidyr::expand_grid(
    b_start     = 1L,
    series_fct  = series_vals,
    measure_fct = factor("Power", levels = meas_lvls),
    rate        = c(0, 1.0)
  )
  
  if (is.factor(plot_dat$series_fct)) {
    blank_t1e$series_fct <- factor(blank_t1e$series_fct, levels = levels(plot_dat$series_fct))
    blank_pow$series_fct <- factor(blank_pow$series_fct, levels = levels(plot_dat$series_fct))
  }
  
  alpha_df <- data.frame(
    measure_fct = factor("Type I error", levels = meas_lvls),
    yint = 0.025
  )
  
  ggplot(
    plot_dat,
    aes(
      x = b_start, y = rate,
      color = rand_key, shape = rand_key,
      linetype = analysis_model_fct,
      group = interaction(rand_key, analysis_model_fct)
    )
  ) +
    geom_blank(data = blank_t1e, aes(x = b_start, y = rate), inherit.aes = FALSE) +
    geom_blank(data = blank_pow, aes(x = b_start, y = rate), inherit.aes = FALSE) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_point(size = 1.9, na.rm = TRUE) +
    geom_hline(
      data = alpha_df,
      aes(yintercept = yint),
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.8
    ) +
    # IMPORTANT: free_y so each row can have its own fixed range via geom_blank
    facet_grid(measure_fct ~ series_fct, scales = "free_y",
               labeller = labeller(series_fct = label_parsed)) +
    labs(title = title, x = X_LAB_TB, y = NULL) +
    x_scale_tb() +
    coord_cartesian(xlim = c(1, 49), clip = "off") +
    bimj_scale_color_rand(include_ref = TRUE) +
    bimj_scale_shape_rand(include_ref = TRUE) +
    bimj_scale_linetype_model() +
    guides(
      color    = guide_legend(order = 1, nrow = 2, byrow = TRUE),
      shape    = guide_legend(order = 1, nrow = 2, byrow = TRUE),
      linetype = guide_legend(order = 2, nrow = 1, byrow = TRUE)
    ) +
    bimj_theme(base_size = 13) +
    legend_layout_theme()
}

# ============================================================
# Allocation bias plots
# ============================================================

alloc_t1e_path <- file.path(in_dir, "allocbias_nonconc_t1e_preferB_bothModels.csv")
alloc_pow_path <- file.path(in_dir, "allocbias_nonconc_power_preferB_bothModels.csv")

alloc_t1e_wide <- read.csv(alloc_t1e_path, check.names = FALSE)
alloc_pow_long <- read.csv(alloc_pow_path, check.names = FALSE) %>% ensure_series_col()

alloc_t1e_long <- alloc_t1e_wide %>%
  mutate(b_start = parse_int(b_start)) %>%
  pivot_longer(cols = c(t1e_A, t1e_B, fwer), names_to = "tmp_series", values_to = "rate") %>%
  mutate(
    series  = recode(tmp_series, "t1e_A" = "A", "t1e_B" = "B", "fwer" = "FWER"),
    measure = "Type I error"
  ) %>%
  select(b_start, panel, analysis_model, config, series, rate, measure)

alloc_pow_long <- alloc_pow_long %>%
  mutate(
    b_start = parse_int(b_start),
    measure = "Power"
  ) %>%
  select(b_start, panel, analysis_model, config, series, rate, measure)

build_alloc_main <- function(main_config) {
  bind_rows(
    alloc_t1e_long %>% filter(config == main_config),
    alloc_pow_long %>% filter(config == main_config)
  )
}

# Returns list(plot=ggplot, file=path) so stacking can never receive a filename
plot_alloc_2x3_save <- function(cfg, title, outfile, width = 11, height = 8) {
  main_dat <- build_alloc_main(cfg) %>% add_common_factors() %>% distinct()
  ref_dat  <- ref_no_bias %>% add_common_factors() %>% distinct()
  
  p <- make_2x3_plot(main_dat, ref_dat, title)
  
  safe_save_pdf(p, outfile, width = width, height = height)
  
  list(plot = p, file = outfile)
}

TITLE_PREFER_B  <- "Prefer B"
TITLE_PREFER_AB <- "Prefer A and B"

p_favorB_1 <- plot_alloc_2x3_save(
  cfg     = "favor_B (eta>0)",
  title   = TITLE_PREFER_B,
  outfile = file.path(out_dir, "nonconc_allocbias_preferB_2x3.pdf")
)

p_favorB_2 <- plot_alloc_2x3_save(
  cfg     = "favor_B_2step (eta>0)",
  title   = TITLE_PREFER_B,
  outfile = file.path(out_dir, "nonconc_allocbias_preferB_2step_2x3.pdf")
)

p_allExp_1 <- plot_alloc_2x3_save(
  cfg     = "favor_all_exp (eta>0)",
  title   = TITLE_PREFER_AB,
  outfile = file.path(out_dir, "nonconc_allocbias_favorAllExp_2x3.pdf")
)

p_allExp_2 <- plot_alloc_2x3_save(
  cfg     = "favor_all_exp_2step (eta>0)",
  title   = TITLE_PREFER_AB,
  outfile = file.path(out_dir, "nonconc_allocbias_favorAllExp_2step_2x3.pdf")
)

if (has_patchwork) {
  suppressPackageStartupMessages(library(patchwork))
  
  p_stack_1 <- patchwork::wrap_plots(
    list(p_favorB_1$plot, p_allExp_1$plot),
    ncol = 1,
    guides = "collect"
  ) & legend_layout_theme()
  
  safe_save_pdf(
    p_stack_1,
    file.path(out_dir, "nonconc_allocbias_1step_favorB_then_favorAllExp.pdf"),
    width = 11, height = 16
  )
  
  p_stack_2 <- patchwork::wrap_plots(
    list(p_favorB_2$plot, p_allExp_2$plot),
    ncol = 1,
    guides = "collect"
  ) & legend_layout_theme()
  
  safe_save_pdf(
    p_stack_2,
    file.path(out_dir, "nonconc_allocbias_2step_favorB_then_favorAllExp.pdf"),
    width = 11, height = 16
  )
}

# ============================================================
# Chronological bias plot (2x3)
#
# Reference detection here must NOT use eta=0.
# ============================================================

# Locate chronobias CSVs (robust: allow nested folders)
find_first_file <- function(root, filename_regex) {
  hits <- list.files(root, pattern = filename_regex, recursive = TRUE, full.names = TRUE)
  if (length(hits) > 0) return(hits[[1]])
  NA_character_
}

chrono_t1e_path <- find_first_file(in_dir, "chronobias_nonconc_t1e_bothModels\\.csv$")
chrono_pow_path <- find_first_file(in_dir, "chronobias_nonconc_power_bothModels\\.csv$")

split_main_ref_chrono <- function(df) {
  cfg_col  <- get_col(df, c("config", "cfg", "scenario"))
  beta_col <- get_col(df, c("beta_time", "beta", "betaTime"))
  
  is_ref <- rep(FALSE, nrow(df))
  
  # Chronobias reference rows: explicit markers, or beta==0
  if (!is.na(cfg_col)) {
    cfg <- as.character(df[[cfg_col]])
    
    # Explicit markers
    is_ref <- is_ref | grepl("no[[:space:]]*bias|reference", cfg, ignore.case = TRUE)
    
    # Explicit beta==0 in text (exact)
    is_ref <- is_ref | grepl("beta(?:_time)?[[:space:]]*=[[:space:]]*0(?:\\.0+)?\\b", cfg, ignore.case = TRUE)
  }
  
  # Numeric beta column if present
  if (!is.na(beta_col)) {
    bt <- suppressWarnings(as.numeric(df[[beta_col]]))
    is_ref <- is_ref | (!is.na(bt) & bt == 0)
  }
  
  # If panel already encodes reference, keep it
  if ("panel" %in% names(df)) {
    pnl <- as.character(df$panel)
    is_ref <- is_ref | pnl %in% c("Reference (no bias)", "PBR0")
    is_ref <- is_ref | grepl("reference|no[[:space:]]*bias", pnl, ignore.case = TRUE)
  }
  
  # Safety: if everything classified as reference but there are multiple configs,
  # fall back to splitting only by explicit "no bias/reference" markers.
  if (all(is_ref) && !is.na(cfg_col)) {
    cfg <- as.character(df[[cfg_col]])
    is_ref2 <- grepl("no[[:space:]]*bias|reference", cfg, ignore.case = TRUE)
    if (any(!is_ref2) && any(is_ref2)) is_ref <- is_ref2
  }
  
  list(
    main = df[!is_ref, , drop = FALSE],
    ref  = df[ is_ref, , drop = FALSE]
  )
}

if (!is.na(chrono_t1e_path) && !is.na(chrono_pow_path) &&
    file.exists(chrono_t1e_path) && file.exists(chrono_pow_path)) {
  
  chrono_t1e_raw <- read.csv(chrono_t1e_path, check.names = FALSE) %>% ensure_series_col()
  chrono_pow_raw <- read.csv(chrono_pow_path, check.names = FALSE) %>% ensure_series_col()
  
  chrono_t1e_raw <- chrono_t1e_raw %>%
    mutate(b_start = parse_int(b_start), measure = "Type I error")
  
  chrono_pow_raw <- chrono_pow_raw %>%
    mutate(b_start = parse_int(b_start), measure = "Power")
  
  chrono_all <- bind_rows(chrono_t1e_raw, chrono_pow_raw)
  sp <- split_main_ref_chrono(chrono_all)
  
  chrono_main <- sp$main %>%
    select(b_start, panel, analysis_model, series, rate, measure)
  
  chrono_ref_from_file <- sp$ref %>%
    filter(
      panel == "Reference (no bias)" |
        grepl("Permuted block randomization", panel, fixed = TRUE) |
        grepl("^PBR", panel)
    ) %>%
    select(b_start, series, analysis_model, measure, rate) %>%
    mutate(panel = "Reference (no bias)")
  
  ref_for_chrono <- if (nrow(chrono_ref_from_file) > 0) chrono_ref_from_file else ref_no_bias
  
  chrono_main2 <- chrono_main %>% add_common_factors() %>% distinct()
  
  if (nrow(chrono_main2) == 0) {
    stop(
      "Chronobias main data is empty after reference-splitting.\n",
      "This indicates your chronobias config/panel strings are being misclassified as reference.\n",
      "Check that only the true no-bias rows have beta_time==0 or explicit 'no bias/reference' markers."
    )
  }
  
  p_chrono <- make_2x3_plot(
    main_dat = chrono_main2,
    ref_dat  = ref_for_chrono %>% add_common_factors() %>% distinct(),
    title    = ""
  )
  
  safe_save_pdf(
    p_chrono,
    file.path(out_dir, "nonconc_chronobias_2x3.pdf"),
    width = 11, height = 8
  )
  
} else {
  message("Chronobias CSVs not found -> skipped nonconc_chronobias_2x3.pdf")
  message("Expected: ", chrono_t1e_path)
  message("Expected: ", chrono_pow_path)
}

message("Done. Plots written to: ", normalizePath(out_dir, mustWork = FALSE))
