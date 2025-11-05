#!/usr/bin/env Rscript
# ============================================================
#  Per-arm violin plots of bias metrics by randomization procedure
#  Arms: A, B vs shared control D (plots show A & B separately)
#
#  Metrics (per arm):
#    - MSE (outcome − alpha)
#    - Mean allocation bias
#    - Mean chronological bias
#
#  Scenarios (3):
#    1) Concurrent (t-test)
#    2) Nonconcurrent (t-test)
#    3) Nonconcurrent (ANOVA with period factor)
#
#  Procedures (3):
#    - Complete randomization
#    - Block randomization (1*arms)   [block_factor = 1]
#    - Block randomization (4*arms)   [block_factor = 4]
#
#  Fixed settings:
#    alpha (one-sided, greater) = 0.025
#    max_group_size = 24 (per arm), expected_total = 96
#    B starts after 16 patients
#    allocation bias η = 0.08, chronological bias β = 0.16 (both ON)
#
#  ENV:
#    OUT_DIR (default "results"), N_RUNS (default 200), SEED_BASE (default 20251031)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------- knobs ----------
n_runs    <- as.integer(Sys.getenv("N_RUNS",    "2000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", "20251031"))

# ---------- fixed trial settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96
B_start_fixed   <- 16L

exp_arms <- c("A","B")
mu_null  <- c(A=0, B=0, D=0)

# Bias magnitudes (both ON)
alloc_bias_val  <- 0.08
chrono_beta_val <- 0.16

# Output dir
out_dir <- Sys.getenv("OUT_DIR", unset = "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# UTF-8 safe devices (headless nodes)
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
use_cairo <- TRUE
png_device <- function(...) if (use_cairo) grDevices::png(type = "cairo", ...) else grDevices::png(...)
pdf_device <- function(...) if (use_cairo) grDevices::cairo_pdf(...) else grDevices::pdf(...)

# ---------- simulator ----------
# Assumes simulation_functions.R is in the working dir (or set your path)
source("simulation_functions.R")

# ---------- one trial -> per-arm bias metrics --------------
run_one_trial_bias <- function(rand_mode = c("complete","block"),
                               block_factor = 1L,
                               concurrent_only = TRUE,
                               analysis_model = c("ttest","anova_period"),
                               rep_seed = NULL) {
  rand_mode      <- match.arg(rand_mode)
  analysis_model <- match.arg(analysis_model)
  if (!is.null(rep_seed)) set.seed(rep_seed)
  
  res <- platform_trials_simulation(
    max_group_size  = max_group_size,
    mu              = mu_null,
    alpha           = alpha_one_sided,          # also MSE ref in patched simulator
    arm_start       = c(A = 0L, B = B_start_fixed),
    concurrent_only = concurrent_only,
    expected_total  = expected_total,
    beta_time       = chrono_beta_val,          # ON
    rand_mode       = rand_mode,
    block_factor    = block_factor,
    alloc_bias      = alloc_bias_val,           # ON
    exp_arms        = exp_arms,
    test_side       = test_side,
    alternative     = alternative,
    bias_policy     = "favor_all_exp",
    analysis_model  = analysis_model
  )
  
  bm <- res$bias_metrics
  if (is.null(bm) || length(bm) == 0) {
    return(data.frame(arm=character(0), metric=character(0), value=numeric(0)))
  }
  data.frame(
    arm    = rownames(bm),
    metric = colnames(bm),
    value  = as.numeric(bm),
    stringsAsFactors = FALSE
  )
}

# ---------- scenarios & procedures ----------
scenarios <- list(
  list(key = "Concurrent (t-test)",         concurrent_only = TRUE,  analysis_model = "ttest"),
  list(key = "Nonconcurrent (t-test)",      concurrent_only = FALSE, analysis_model = "ttest"),
  list(key = "Nonconcurrent (ANOVA)",       concurrent_only = FALSE, analysis_model = "anova_period")
)

procedures <- list(
  list(key = "Complete randomization",           rand_mode="complete", block_factor=1L),
  list(key = "Block randomization (1*arms)",     rand_mode="block",    block_factor=1L),
  list(key = "Block randomization (4*arms)",     rand_mode="block",    block_factor=4L)
)

# ---------- run all combos; collect per-run, per-arm metrics ----------
collect_runs <- function() {
  rows <- list(); idx <- 0L
  for (sc_i in seq_along(scenarios)) {
    sc <- scenarios[[sc_i]]
    for (pr_i in seq_along(procedures)) {
      pr <- procedures[[pr_i]]
      seeds <- seed_base + seq_len(n_runs) + 10000L*sc_i + 1000L*pr_i
      for (r in seq_len(n_runs)) {
        df <- run_one_trial_bias(
          rand_mode       = pr$rand_mode,
          block_factor    = pr$block_factor,
          concurrent_only = sc$concurrent_only,
          analysis_model  = sc$analysis_model,
          rep_seed        = seeds[r]
        )
        if (!nrow(df)) next
        df$scenario  <- sc$key
        df$procedure <- pr$key
        df$run       <- r
        idx <- idx + 1L
        rows[[idx]] <- df
      }
    }
  }
  dat <- do.call(rbind, rows)
  
  # ---- Robust metric detection & renaming (per arm) -------------
  all_metrics <- unique(dat$metric)
  
  # MSE: prefer exact alpha; fall back to any mse_outcome_vs_*
  mse_expected <- sprintf("mse_outcome_vs_%.3f", alpha_one_sided)  # "mse_outcome_vs_0.025"
  mse_any <- grep("^mse_outcome_vs_", all_metrics, value = TRUE)
  if (!length(mse_any)) warning("No MSE column found in bias_metrics. Present: ", paste(all_metrics, collapse=", "))
  mse_col <- if (mse_expected %in% all_metrics) mse_expected else if (length(mse_any)) mse_any[1] else NA_character_
  
  # Allocation bias: look for exact + regex backup
  alloc_col <- if ("mean_allocation_bias" %in% all_metrics) {
    "mean_allocation_bias"
  } else {
    hit <- grep("alloc", all_metrics, ignore.case = TRUE, value = TRUE)
    if (length(hit)) hit[1] else NA_character_
  }
  
  # Chronological bias: exact + regex backup (fix for your missing column)
  chrono_col <- if ("mean_chronological_bias" %in% all_metrics) {
    "mean_chronological_bias"
  } else {
    hit <- grep("chron", all_metrics, ignore.case = TRUE, value = TRUE)
    if (length(hit)) hit[1] else NA_character_
  }
  
  # Keep + rename nicely
  keep <- c(mse_col, alloc_col, chrono_col)
  keep <- keep[!is.na(keep)]
  nice_names <- c(
    setNames("MSE (outcome \u2212 0.025)", mse_col),
    setNames("Mean allocation bias",       alloc_col),
    setNames("Mean chronological bias",    chrono_col)
  )
  
  dat <- subset(dat, arm %in% c("A","B") & metric %in% keep)
  dat$metric <- unname(nice_names[dat$metric])
  
  # Diagnostics
  message("Detected metrics -> ",
          paste(sprintf("%s = '%s'",
                        c("MSE","Allocation","Chronological"),
                        c(mse_col, alloc_col, chrono_col)), collapse = " | "))
  
  # factor orders
  dat$scenario  <- factor(dat$scenario,
                          levels = c("Concurrent (t-test)",
                                     "Nonconcurrent (t-test)",
                                     "Nonconcurrent (ANOVA)"))
  dat$procedure <- factor(dat$procedure,
                          levels = c("Complete randomization",
                                     "Block randomization (1*arms)",
                                     "Block randomization (4*arms)"))
  dat$metric    <- factor(dat$metric,
                          levels = c("MSE (outcome \u2212 0.025)",
                                     "Mean allocation bias",
                                     "Mean chronological bias"))
  dat$arm <- factor(dat$arm, levels = c("A","B"))
  dat
}

# ---------- plotting (VIOLIN) ----------
make_violins <- function(dat, file_stub = "violins_bias_metrics_by_proc_per_arm") {
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  write.csv(dat, csv_path, row.names = FALSE)
  
  p <- ggplot(dat, aes(x = procedure, y = value, fill = arm)) +
    geom_violin(position = position_dodge(width = 0.85), trim = FALSE, alpha = 0.7) +
    # median point per group
    stat_summary(fun = median, geom = "point",
                 position = position_dodge(width = 0.85), size = 1.3, color = "black") +
    facet_grid(rows = vars(scenario), cols = vars(metric), scales = "free_y") +
    labs(
      title = "Per-arm bias metrics by randomization procedure (violin plots)",
      subtitle = "A & B shown separately; B starts at 16; both biases (η = 0.08, β = 0.16); MSE ref = α",
      x = "Randomization procedure",
      y = "Metric value",
      fill = "Arm"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 15, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  png_path <- file.path(out_dir, paste0(file_stub, ".png"))
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  
  png_device(filename = png_path, width = 14, height = 7.8, units = "in", res = 150)
  print(p); dev.off()
  pdf_device(file = pdf_path, width = 14, height = 7.8)
  print(p); dev.off()
  
  message("Wrote:")
  message("  Data: ", csv_path)
  message("  PNG : ", png_path)
  message("  PDF : ", pdf_path)
}

# ---------- run ----------
set.seed(seed_base)
dat_all <- collect_runs()
make_violins(dat_all)
