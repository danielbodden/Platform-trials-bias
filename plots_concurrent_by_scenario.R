#!/usr/bin/env Rscript
# ============================================================
#  Rejection probability vs. timing of Arm B (A,B vs D)
#  Concurrent case: 4 separate figures (No bias; Alloc; Chrono; Both) – t-test
#  Nonconcurrent case: 4 separate figures, each overlaying t-test vs ANOVA
#
#  Fixed settings:
#    alpha (one-sided, greater) = 0.025
#    max_group_size = 24 (per arm), expected_total = 96
#    allocation bias (eta) = 0.08, chronological bias (beta) = 0.16
#
#  Cluster-friendly:
#    - N_SIM_PER_POINT, SEED_BASE, SLURM_CPUS_PER_TASK, OUT_DIR
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------- cluster knobs ----------
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "200"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE",       "20251031"))

# ---------- fixed trial settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

exp_arms <- c("A","B")
mu_null  <- c(A=0, B=0, D=0)

# Bias magnitudes
alloc_bias_val  <- 0.08   # eta
chrono_beta_val <- 0.16   # beta

# Sweep grid for timing of Arm B opening (patients)
# (latest start point is 48 patients)
b_grid <- seq(0, 48, by = 12)

# I/O
out_dir <- Sys.getenv("OUT_DIR", unset = "PT_bias/results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# --- robust UTF-8 plotting on headless nodes (Cairo devices) ---
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
use_cairo <- TRUE
pdf_device <- function(...) if (use_cairo) grDevices::cairo_pdf(...) else grDevices::pdf(...)

# --- simulator ---
source("PT_bias/simulation_functions.R")

# --- core sweep for one panel ---
.run_sweep <- function(
    rand_mode       = c("complete","block","bigstick"),
    block_factor    = 1L,
    alloc_on        = FALSE,
    chrono_on       = FALSE,
    concurrent_only = TRUE,
    analysis_model  = c("ttest","lm_time"),
    seed_bump       = 0L,
    two_step        = FALSE
) {
  rand_mode      <- match.arg(rand_mode)
  analysis_model <- match.arg(analysis_model)
  
  rows <- lapply(seq_along(b_grid), function(i) {
    bs <- b_grid[i]
    set.seed(seed_base + seed_bump + i +
               1000L * block_factor +
               ifelse(rand_mode == "complete", 9999L, 0L) +
               ifelse(analysis_model == "lm_time", 7777L, 0L))
    
    out <- calc_rejection_summary(
      n_sim          = n_sim_per_point,
      n_cores        = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      max_group_size = max_group_size,
      mu             = mu_null,
      alpha          = alpha_one_sided,
      arm_start      = c(A = 0, B = bs),
      concurrent_only= concurrent_only,
      expected_total = expected_total,
      beta_time      = if (chrono_on) chrono_beta_val else 0,
      rand_mode      = rand_mode,
      block_factor   = block_factor,
      alloc_bias     = if (alloc_on) alloc_bias_val else 0,
      exp_arms       = exp_arms,
      test_side      = test_side,
      alternative    = alternative,
      bias_policy    = "favor_B",
      analysis_model = analysis_model,
      two_step       = two_step
    )
    
    pr   <- out$per_comparison_rejection_rate
    want <- intersect(c("A_vs_D","B_vs_D"), names(pr))
    
    df_ab <- data.frame(
      B_start = bs,
      series  = sub("_vs_.*$", "", want),  # "A","B"
      rate    = as.numeric(pr[want]),
      stringsAsFactors = FALSE
    )
    df_fwer <- data.frame(
      B_start = bs,
      series  = "FWER",
      rate    = out$fwer,
      stringsAsFactors = FALSE
    )
    rbind(df_ab, df_fwer)
  })
  
  do.call(rbind, rows)
}

# --- assemble 4 panels for one scenario + one model (1-step) ---
build_panel_data_single_model <- function(alloc_on, chrono_on, concurrent_only, scenario_label, analysis_model = "ttest") {
  d_b1   <- .run_sweep("block",    1L, alloc_on, chrono_on, concurrent_only, analysis_model,
                       seed_bump = 10L, two_step = FALSE);  d_b1$panel   <- "Block randomization (1*arms, 1-step)"
                       d_b4   <- .run_sweep("block",    8L, alloc_on, chrono_on, concurrent_only, analysis_model,
                                            seed_bump = 40L, two_step = FALSE);  d_b4$panel   <- "Block randomization (8*arms, 1-step)"
                       d_bs   <- .run_sweep("bigstick", 1L, alloc_on, chrono_on, concurrent_only, analysis_model,
                                            seed_bump = 20L, two_step = FALSE);  d_bs$panel   <- "Big-stick design (a = 2, 1-step)"
                       d_comp <- .run_sweep("complete", 1L, alloc_on, chrono_on, concurrent_only, analysis_model,
                                            seed_bump = 0L,  two_step = FALSE);  d_comp$panel <- "Complete randomization (1-step)"
                       
                       dat <- rbind(d_b1, d_b4, d_bs, d_comp)
                       dat$series   <- factor(dat$series, levels = c("A","B","FWER"))
                       dat$scenario <- scenario_label
                       dat$model    <- analysis_model
                       dat
}

# --- assemble 4 panels for one scenario + one model (2-step) ---
build_panel_data_single_model_2step <- function(alloc_on, chrono_on, concurrent_only, scenario_label, analysis_model = "ttest") {
  d_b1   <- .run_sweep("block",    1L, alloc_on, chrono_on, concurrent_only, analysis_model,
                       seed_bump = 110L, two_step = TRUE);  d_b1$panel   <- "Block randomization (1*arms, 2-step)"
                       d_b4   <- .run_sweep("block",    8L, alloc_on, chrono_on, concurrent_only, analysis_model,
                                            seed_bump = 140L, two_step = TRUE);  d_b4$panel   <- "Block randomization (8*arms, 2-step)"
                       d_bs   <- .run_sweep("bigstick", 1L, alloc_on, chrono_on, concurrent_only, analysis_model,
                                            seed_bump = 120L, two_step = TRUE);  d_bs$panel   <- "Big-stick design (a = 2, 2-step)"
                       d_comp <- .run_sweep("complete", 1L, alloc_on, chrono_on, concurrent_only, analysis_model,
                                            seed_bump = 100L, two_step = TRUE);  d_comp$panel <- "Complete randomization (2-step)"
                       
                       dat <- rbind(d_b1, d_b4, d_bs, d_comp)
                       dat$series   <- factor(dat$series, levels = c("A","B","FWER"))
                       dat$scenario <- scenario_label
                       dat$model    <- analysis_model
                       dat
}

# --- assemble nonconcurrent overlay (t-test vs ANOVA, 1-step) ---
build_panel_data_overlay_models <- function(alloc_on, chrono_on, scenario_label) {
  d_t  <- build_panel_data_single_model(alloc_on, chrono_on, concurrent_only = FALSE,
                                        scenario_label = scenario_label, analysis_model = "ttest")
  d_an <- build_panel_data_single_model(alloc_on, chrono_on, concurrent_only = FALSE,
                                        scenario_label = scenario_label, analysis_model = "lm_time")
  dat <- rbind(d_t, d_an)
  dat$model <- factor(dat$model, levels = c("ttest","lm_time"), labels = c("t-test","ANOVA"))
  dat
}

# --- assemble nonconcurrent overlay (t-test vs ANOVA, 2-step) ---
build_panel_data_overlay_models_2step <- function(alloc_on, chrono_on, scenario_label) {
  d_t  <- build_panel_data_single_model_2step(alloc_on, chrono_on, concurrent_only = FALSE,
                                              scenario_label = scenario_label, analysis_model = "ttest")
  d_an <- build_panel_data_single_model_2step(alloc_on, chrono_on, concurrent_only = FALSE,
                                              scenario_label = scenario_label, analysis_model = "lm_time")
  dat <- rbind(d_t, d_an)
  dat$model <- factor(dat$model, levels = c("ttest","lm_time"), labels = c("t-test","ANOVA"))
  dat
}

# --- common plot pieces ---
add_reference_line <- function() {
  geom_hline(yintercept = 0.025, color = "red", linetype = "dashed", linewidth = 0.8)
}
limits_cartesian <- function() {
  coord_cartesian(ylim = c(0, 0.10), xlim = c(0, 48))
}

# --- plotting helpers (CSV + PDF only) ---
plot_and_save_single_model <- function(dat, title_top, subtitle, file_stub) {
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  
  write.csv(dat, csv_path, row.names = FALSE)
  
  p <- ggplot(dat, aes(x = B_start, y = rate, color = series)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    add_reference_line() +
    facet_wrap(~ panel, nrow = 1) +
    scale_color_manual(values = c("A"="#1b9e77","B"="#d95f02","FWER"="#7570b3"),
                       name = NULL, breaks = c("A","B","FWER"),
                       labels = c("Arm A","Arm B","FWER")) +
    limits_cartesian() +
    labs(
      title = title_top,
      subtitle = subtitle,
      x = "Timing of addition of Arm B (patients)",
      y = "Rejection probability"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  pdf_device(file = pdf_path, width = 12, height = 3.8)
  print(p); dev.off()
  
  message("Wrote:")
  message("  Data: ", csv_path)
  message("  PDF : ", pdf_path)
}

plot_and_save_overlay_models <- function(dat, title_top, subtitle, file_stub) {
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  
  write.csv(dat, csv_path, row.names = FALSE)
  
  p <- ggplot(dat, aes(x = B_start, y = rate, color = series,
                       linetype = model, shape = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    add_reference_line() +
    facet_wrap(~ panel, nrow = 1) +
    scale_color_manual(values = c("A"="#1b9e77","B"="#d95f02","FWER"="#7570b3"),
                       name = NULL, breaks = c("A","B","FWER"),
                       labels = c("Arm A","Arm B","FWER")) +
    scale_linetype_manual(values = c("t-test" = "solid", "ANOVA" = "dashed"), name = "Analysis") +
    scale_shape_manual(values = c("t-test" = 16, "ANOVA" = 17), name = "Analysis") +
    limits_cartesian() +
    labs(
      title = title_top,
      subtitle = subtitle,
      x = "Timing of addition of Arm B (patients)",
      y = "Rejection probability"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  pdf_device(file = pdf_path, width = 12, height = 3.8)
  print(p); dev.off()
  
  message("Wrote (overlay models):")
  message("  Data: ", csv_path)
  message("  PDF : ", pdf_path)
}

# ============================================================
# Concurrent (t-test), 4 scenarios — ASCII subtitles, 1-step
# ============================================================
dat_conc_nb   <- build_panel_data_single_model(FALSE, FALSE, TRUE, "No bias",                                   "ttest")
plot_and_save_single_model(dat_conc_nb,   "Concurrent (1-step): Rejection probability vs timing of Arm B", "No bias",                                 "conc_reject_vs_Bstart_no_bias")

dat_conc_ab   <- build_panel_data_single_model(TRUE,  FALSE, TRUE, "Allocation bias only (eta = 0.08)",        "ttest")
plot_and_save_single_model(dat_conc_ab,   "Concurrent (1-step): Rejection probability vs timing of Arm B", "Allocation bias only (eta = 0.08)",        "conc_reject_vs_Bstart_alloc_bias")

dat_conc_cb   <- build_panel_data_single_model(FALSE, TRUE,  TRUE, "Chronological bias only (beta = 0.16)",     "ttest")
plot_and_save_single_model(dat_conc_cb,   "Concurrent (1-step): Rejection probability vs timing of Arm B", "Chronological bias only (beta = 0.16)",     "conc_reject_vs_Bstart_chrono_bias")

dat_conc_both <- build_panel_data_single_model(TRUE,  TRUE,  TRUE, "Allocation + chronological (eta = 0.08, beta = 0.16)", "ttest")
plot_and_save_single_model(dat_conc_both, "Concurrent (1-step): Rejection probability vs timing of Arm B", "Allocation + chronological (eta = 0.08, beta = 0.16)", "conc_reject_vs_Bstart_alloc_plus_chrono")

# ============================================================
# Concurrent (t-test), 4 scenarios — 2-step randomization
# ============================================================
dat_conc_nb_2   <- build_panel_data_single_model_2step(FALSE, FALSE, TRUE, "No bias",                                   "ttest")
plot_and_save_single_model(dat_conc_nb_2,   "Concurrent (2-step): Rejection probability vs timing of Arm B", "No bias",                                 "conc_reject_vs_Bstart_no_bias_2step")

dat_conc_ab_2   <- build_panel_data_single_model_2step(TRUE,  FALSE, TRUE, "Allocation bias only (eta = 0.08)",        "ttest")
plot_and_save_single_model(dat_conc_ab_2,   "Concurrent (2-step): Rejection probability vs timing of Arm B", "Allocation bias only (eta = 0.08)",        "conc_reject_vs_Bstart_alloc_bias_2step")

dat_conc_cb_2   <- build_panel_data_single_model_2step(FALSE, TRUE,  TRUE, "Chronological bias only (beta = 0.16)",     "ttest")
plot_and_save_single_model(dat_conc_cb_2,   "Concurrent (2-step): Rejection probability vs timing of Arm B", "Chronological bias only (beta = 0.16)",     "conc_reject_vs_Bstart_chrono_bias_2step")

dat_conc_both_2 <- build_panel_data_single_model_2step(TRUE,  TRUE,  TRUE, "Allocation + chronological (eta = 0.08, beta = 0.16)", "ttest")
plot_and_save_single_model(dat_conc_both_2, "Concurrent (2-step): Rejection probability vs timing of Arm B", "Allocation + chronological (eta = 0.08, beta = 0.16)", "conc_reject_vs_Bstart_alloc_plus_chrono_2step")

# ============================================================
# Nonconcurrent (overlay t-test vs ANOVA), 4 scenarios — 1-step
# ============================================================
dat_non_nb   <- build_panel_data_overlay_models(FALSE, FALSE, "No bias")
plot_and_save_overlay_models(
  dat_non_nb,
  title_top = "Nonconcurrent (1-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — No bias",
  file_stub = "nonconc_reject_vs_Bstart_no_bias_overlay_models"
)

dat_non_ab   <- build_panel_data_overlay_models(TRUE,  FALSE, "Allocation bias only (eta = 0.08)")
plot_and_save_overlay_models(
  dat_non_ab,
  title_top = "Nonconcurrent (1-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — Allocation bias only (eta = 0.08)",
  file_stub = "nonconc_reject_vs_Bstart_alloc_bias_overlay_models"
)

dat_non_cb   <- build_panel_data_overlay_models(FALSE, TRUE,  "Chronological bias only (beta = 0.16)")
plot_and_save_overlay_models(
  dat_non_cb,
  title_top = "Nonconcurrent (1-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — Chronological bias only (beta = 0.16)",
  file_stub = "nonconc_reject_vs_Bstart_chrono_bias_overlay_models"
)

dat_non_both <- build_panel_data_overlay_models(TRUE,  TRUE,  "Allocation + chronological (eta = 0.08, beta = 0.16)")
plot_and_save_overlay_models(
  dat_non_both,
  title_top = "Nonconcurrent (1-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — Allocation + chronological (eta = 0.08, beta = 0.16)",
  file_stub = "nonconc_reject_vs_Bstart_alloc_plus_chrono_overlay_models"
)

# ============================================================
# Nonconcurrent (overlay t-test vs ANOVA), 4 scenarios — 2-step
# ============================================================
dat_non_nb_2   <- build_panel_data_overlay_models_2step(FALSE, FALSE, "No bias")
plot_and_save_overlay_models(
  dat_non_nb_2,
  title_top = "Nonconcurrent (2-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — No bias",
  file_stub = "nonconc_reject_vs_Bstart_no_bias_overlay_models_2step"
)

dat_non_ab_2   <- build_panel_data_overlay_models_2step(TRUE,  FALSE, "Allocation bias only (eta = 0.08)")
plot_and_save_overlay_models(
  dat_non_ab_2,
  title_top = "Nonconcurrent (2-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — Allocation bias only (eta = 0.08)",
  file_stub = "nonconc_reject_vs_Bstart_alloc_bias_overlay_models_2step"
)

dat_non_cb_2   <- build_panel_data_overlay_models_2step(FALSE, TRUE,  "Chronological bias only (beta = 0.16)")
plot_and_save_overlay_models(
  dat_non_cb_2,
  title_top = "Nonconcurrent (2-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — Chronological bias only (beta = 0.16)",
  file_stub = "nonconc_reject_vs_Bstart_chrono_bias_overlay_models_2step"
)

dat_non_both_2 <- build_panel_data_overlay_models_2step(TRUE,  TRUE,  "Allocation + chronological (eta = 0.08, beta = 0.16)")
plot_and_save_overlay_models(
  dat_non_both_2,
  title_top = "Nonconcurrent (2-step): Rejection probability vs timing of Arm B",
  subtitle  = "Overlay: t-test vs ANOVA — Allocation + chronological (eta = 0.08, beta = 0.16)",
  file_stub = "nonconc_reject_vs_Bstart_alloc_plus_chrono_overlay_models_2step"
)

message("All figures written to: ", normalizePath(out_dir, mustWork = FALSE))
