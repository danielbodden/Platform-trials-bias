#!/usr/bin/env Rscript
# ============================================================
#  Type I error vs allocation bias (eta) for different
#  randomization schemes and allocation bias policies
#
#  - Arms: A, B vs D (control)
#  - Start times: A at 0, B at 10
#  - Concurrent controls only, t-test analysis
#  - 1-step: PBR (bf=1,8), big-stick (a=2), CR;
#            allocation policies 1–3 (in one figure)
#  - 2-step: same randomization schemes;
#            allocation policies 1–5 (in one figure)
#
#  Figures:
#    - facets: 3 rows (A, B, FWER) × 4 cols (rand. schemes)
#    - lines: allocation bias policies
#
#  Extended:
#    - Also generate nonconcurrent figures (concurrent_only = FALSE)
#      overlaying t-test vs ANOVA (lm_time), analogous to the
#      "plots_concurrent_by_scenario" script.
#
#  Cluster-friendly:
#    - N_SIM_PER_POINT, SEED_BASE, SLURM_CPUS_PER_TASK, OUT_DIR
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------- cluster knobs ----------
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "20000"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE",       "20251120"))

# ---------- fixed trial settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

exp_arms <- c("A","B")
mu_null  <- c(A=0, B=0, D=0)

# Arm start times: A at 0, B at 10
arm_start <- c(A = 0, B = 10)

# Sweep grid for allocation bias (eta)
eta_grid <- seq(0, 0.20, by = 0.04)

# I/O
out_dir <- Sys.getenv("OUT_DIR", unset = "PT_bias/results_t1e_vs_alloc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# --- robust UTF-8 plotting on headless nodes (Cairo devices) ---
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
use_cairo <- TRUE
pdf_device <- function(...) if (use_cairo) grDevices::cairo_pdf(...) else grDevices::pdf(...)

# --- simulator ---
source("PT_bias/simulation_functions.R")

# ---------- helper: human-readable policy labels ----------
# Short, ordinal labels 1–5 used in the legend
.policy_ordinal_label <- function(bias_policy) {
  switch(
    bias_policy,
    "favor_B"             = "Policy 1: Prefer B",
    "favor_all_exp"       = "Policy 2: Prefer all experimental",
    "average"             = "Policy 3: Average",
    "favor_B_2step"       = "Policy 4: Prefer B (2-step cohort)",
    "favor_all_exp_2step" = "Policy 5: Prefer all experimental (2-step cohort)",
    bias_policy
  )
}

# (kept for potential reuse; not essential now)
.policy_stub <- function(bias_policy, step_mode) {
  paste0(
    if (step_mode == "1-step") "1step_" else "2step_",
    bias_policy
  )
}

# ---------- core sweep: one rand_mode / block_factor / two_step / policy ----------
.run_sweep_alloc <- function(
    rand_mode      = c("complete","block","bigstick"),
    block_factor   = 1L,
    two_step       = FALSE,
    bias_policy    = "favor_B",
    seed_bump      = 0L,
    concurrent_only= TRUE,
    analysis_model = c("ttest","lm_time")
) {
  rand_mode      <- match.arg(rand_mode)
  analysis_model <- match.arg(analysis_model)
  
  rows <- lapply(seq_along(eta_grid), function(i) {
    eta <- eta_grid[i]
    
    set.seed(seed_base + seed_bump + i +
               1000L * block_factor +
               ifelse(rand_mode == "complete", 9999L, 0L) +
               ifelse(two_step, 5000L, 0L) +
               ifelse(analysis_model == "lm_time", 7777L, 0L))
    
    out <- calc_rejection_summary(
      n_sim          = n_sim_per_point,
      n_cores        = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      max_group_size = max_group_size,
      mu             = mu_null,
      alpha          = alpha_one_sided,
      arm_start      = arm_start,
      concurrent_only= concurrent_only,
      expected_total = expected_total,
      beta_time      = 0,
      rand_mode      = rand_mode,
      block_factor   = block_factor,
      alloc_bias     = eta,
      exp_arms       = exp_arms,
      test_side      = test_side,
      alternative    = alternative,
      bias_policy    = bias_policy,
      analysis_model = analysis_model,
      two_step       = two_step,
      bigstick_a     = 2L
    )
    
    pr   <- out$per_comparison_rejection_rate
    # should contain A_vs_D, B_vs_D
    want <- intersect(c("A_vs_D","B_vs_D"), names(pr))
    
    df_ab <- data.frame(
      alloc_bias = eta,
      series     = sub("_vs_.*$", "", want),  # "A","B"
      rate       = as.numeric(pr[want]),
      stringsAsFactors = FALSE
    )
    df_fwer <- data.frame(
      alloc_bias = eta,
      series     = "FWER",
      rate       = out$fwer,
      stringsAsFactors = FALSE
    )
    rbind(df_ab, df_fwer)
  })
  
  do.call(rbind, rows)
}

# ---------- assemble panels for 1-step (policies 1–3) ----------
# now generic: can be used for concurrent/nonconcurrent and t-test/ANOVA
build_panel_data_1step <- function(bias_policy,
                                   concurrent_only = TRUE,
                                   analysis_model  = "ttest") {
  # 1-step: only policies 1–3 are allowed:
  #   "favor_B", "favor_all_exp", "average"
  step_mode <- "1-step"
  
  d_b1   <- .run_sweep_alloc("block",    1L, two_step = FALSE, bias_policy = bias_policy,
                             seed_bump =  10L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_b1$panel <- "PBR (block factor = 1, 1-step)"
  
  d_b8   <- .run_sweep_alloc("block",    8L, two_step = FALSE, bias_policy = bias_policy,
                             seed_bump =  80L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_b8$panel <- "PBR (block factor = 8, 1-step)"
  
  d_bs   <- .run_sweep_alloc("bigstick", 1L, two_step = FALSE, bias_policy = bias_policy,
                             seed_bump =  20L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_bs$panel <- "Big-stick (a = 2, 1-step)"
  
  d_comp <- .run_sweep_alloc("complete", 1L, two_step = FALSE, bias_policy = bias_policy,
                             seed_bump =   0L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_comp$panel <- "Complete randomization (1-step)"
  
  dat <- rbind(d_b1, d_b8, d_bs, d_comp)
  dat$series    <- factor(dat$series, levels = c("A","B","FWER"))
  dat$step_mode <- step_mode
  dat$policy    <- bias_policy
  dat$model     <- analysis_model
  dat
}

# ---------- assemble panels for 2-step (policies 1–5) ----------
# now generic: can be used for concurrent/nonconcurrent and t-test/ANOVA
build_panel_data_2step <- function(bias_policy,
                                   concurrent_only = TRUE,
                                   analysis_model  = "ttest") {
  # 2-step: policies 1–5:
  #   "favor_B", "favor_all_exp", "average",
  #   "favor_B_2step", "favor_all_exp_2step"
  step_mode <- "2-step"
  
  d_b1   <- .run_sweep_alloc("block",    1L, two_step = TRUE, bias_policy = bias_policy,
                             seed_bump = 110L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_b1$panel <- "PBR (block factor = 1, 2-step)"
  
  d_b8   <- .run_sweep_alloc("block",    8L, two_step = TRUE, bias_policy = bias_policy,
                             seed_bump = 180L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_b8$panel <- "PBR (block factor = 8, 2-step)"
  
  d_bs   <- .run_sweep_alloc("bigstick", 1L, two_step = TRUE, bias_policy = bias_policy,
                             seed_bump = 120L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_bs$panel <- "Big-stick (a = 2, 2-step)"
  
  d_comp <- .run_sweep_alloc("complete", 1L, two_step = TRUE, bias_policy = bias_policy,
                             seed_bump = 100L,
                             concurrent_only = concurrent_only,
                             analysis_model  = analysis_model); 
  d_comp$panel <- "Complete randomization (2-step)"
  
  dat <- rbind(d_b1, d_b8, d_bs, d_comp)
  dat$series    <- factor(dat$series, levels = c("A","B","FWER"))
  dat$step_mode <- step_mode
  dat$policy    <- bias_policy
  dat$model     <- analysis_model
  dat
}

# ---------- common plot pieces ----------
add_reference_line <- function() {
  geom_hline(yintercept = alpha_one_sided, color = "red", linetype = "dashed", linewidth = 0.8)
}

limits_cartesian <- function() {
  coord_cartesian(ylim = c(0, 0.10), xlim = c(min(eta_grid), max(eta_grid)))
}

# ---------- plotting: one figure per step_mode with all policies (concurrent, t-test) ----------
plot_and_save_step_mode <- function(dat, step_mode, policies, file_stub) {
  # Map policies to ordinal labels
  dat$policy_label <- factor(
    vapply(dat$policy, .policy_ordinal_label, character(1)),
    levels = vapply(policies, .policy_ordinal_label, character(1))
  )
  
  # Save raw data
  csv_path <- file.path(out_dir, paste0("t1e_vs_alloc_", file_stub, ".csv"))
  pdf_path <- file.path(out_dir, paste0("t1e_vs_alloc_", file_stub, ".pdf"))
  write.csv(dat, csv_path, row.names = FALSE)
  
  p <- ggplot(dat, aes(x = alloc_bias, y = rate, color = policy_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    add_reference_line() +
    facet_grid(series ~ panel) +  # 3 rows (A,B,FWER) × 4 cols (rand schemes)
    limits_cartesian() +
    labs(
      title     = sprintf("Type I error vs allocation bias (η) – %s", step_mode),
      subtitle  = "Rows: A, B, FWER; Columns: randomization procedures; Lines: allocation policies",
      x         = expression(paste("Allocation bias parameter ", eta)),
      y         = "Type I error",
      color     = "Allocation bias policy"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box      = "horizontal",
      strip.text      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold")
    )
  
  pdf_device(file = pdf_path, width = 11.5, height = 7.0)
  print(p); dev.off()
  
  message("Wrote:")
  message("  Data: ", csv_path)
  message("  PDF : ", pdf_path)
}

# ---------- plotting: nonconcurrent overlay t-test vs ANOVA ----------
plot_and_save_step_mode_overlay <- function(dat, step_mode, policies, file_stub) {
  # Policy labels
  dat$policy_label <- factor(
    vapply(dat$policy, .policy_ordinal_label, character(1)),
    levels = vapply(policies, .policy_ordinal_label, character(1))
  )
  # Model factor
  dat$model <- factor(dat$model, levels = c("ttest","lm_time"),
                      labels = c("t-test","ANOVA"))
  
  csv_path <- file.path(out_dir, paste0("t1e_vs_alloc_", file_stub, ".csv"))
  pdf_path <- file.path(out_dir, paste0("t1e_vs_alloc_", file_stub, ".pdf"))
  write.csv(dat, csv_path, row.names = FALSE)
  
  p <- ggplot(dat, aes(x = alloc_bias, y = rate,
                       color = policy_label,
                       linetype = model,
                       shape = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    add_reference_line() +
    facet_grid(series ~ panel) +
    limits_cartesian() +
    labs(
      title     = sprintf("Type I error vs allocation bias (η) – %s (nonconcurrent)", step_mode),
      subtitle  = "Rows: A, B, FWER; Columns: randomization procedures; Lines: policies, line type: analysis model",
      x         = expression(paste("Allocation bias parameter ", eta)),
      y         = "Type I error",
      color     = "Allocation bias policy",
      linetype  = "Analysis model",
      shape     = "Analysis model"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box      = "horizontal",
      strip.text      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold")
    )
  
  pdf_device(file = pdf_path, width = 11.5, height = 7.0)
  print(p); dev.off()
  
  message("Wrote (nonconcurrent overlay):")
  message("  Data: ", csv_path)
  message("  PDF : ", pdf_path)
}

# ============================================================
# Run sweeps and generate figures
# ============================================================

# 1-step policies: 1–3
policies_1step <- c("favor_B", "favor_all_exp", "average")

message("Running 1-step (policies 1–3, concurrent, t-test) ...")
dat_1_list <- lapply(policies_1step, function(bp) {
  message("  1-step, policy = ", bp, " ...")
  build_panel_data_1step(bp, concurrent_only = TRUE, analysis_model = "ttest")
})
dat_1_all <- do.call(rbind, dat_1_list)

# Ensure consistent ordering
dat_1_all$panel  <- factor(dat_1_all$panel,
                           levels = c("PBR (block factor = 1, 1-step)",
                                      "PBR (block factor = 8, 1-step)",
                                      "Big-stick (a = 2, 1-step)",
                                      "Complete randomization (1-step)"))
dat_1_all$series <- factor(dat_1_all$series, levels = c("A","B","FWER"))

plot_and_save_step_mode(
  dat      = dat_1_all,
  step_mode= "1-step",
  policies = policies_1step,
  file_stub= "1step_all_policies"
)

# 2-step policies: 1–5
policies_2step <- c("favor_B", "favor_all_exp", "average", "favor_B_2step", "favor_all_exp_2step")

message("Running 2-step (policies 1–5, concurrent, t-test) ...")
dat_2_list <- lapply(policies_2step, function(bp) {
  message("  2-step, policy = ", bp, " ...")
  build_panel_data_2step(bp, concurrent_only = TRUE, analysis_model = "ttest")
})
dat_2_all <- do.call(rbind, dat_2_list)

# Ensure consistent ordering
dat_2_all$panel  <- factor(dat_2_all$panel,
                           levels = c("PBR (block factor = 1, 2-step)",
                                      "PBR (block factor = 8, 2-step)",
                                      "Big-stick (a = 2, 2-step)",
                                      "Complete randomization (2-step)"))
dat_2_all$series <- factor(dat_2_all$series, levels = c("A","B","FWER"))

plot_and_save_step_mode(
  dat      = dat_2_all,
  step_mode= "2-step",
  policies = policies_2step,
  file_stub= "2step_all_policies"
)

# ============================================================
# Nonconcurrent (concurrent_only = FALSE) – overlay t-test vs ANOVA
# ============================================================

# 1-step, nonconcurrent
message("Running 1-step (policies 1–3, NONCONCURRENT, t-test vs ANOVA) ...")
dat_1_nc_list <- lapply(policies_1step, function(bp) {
  message("  1-step NONCONC, policy = ", bp, " ...")
  d_t  <- build_panel_data_1step(bp, concurrent_only = FALSE, analysis_model = "ttest")
  d_an <- build_panel_data_1step(bp, concurrent_only = FALSE, analysis_model = "lm_time")
  rbind(d_t, d_an)
})
dat_1_nc_all <- do.call(rbind, dat_1_nc_list)

dat_1_nc_all$panel  <- factor(dat_1_nc_all$panel,
                              levels = c("PBR (block factor = 1, 1-step)",
                                         "PBR (block factor = 8, 1-step)",
                                         "Big-stick (a = 2, 1-step)",
                                         "Complete randomization (1-step)"))
dat_1_nc_all$series <- factor(dat_1_nc_all$series, levels = c("A","B","FWER"))

plot_and_save_step_mode_overlay(
  dat      = dat_1_nc_all,
  step_mode= "1-step",
  policies = policies_1step,
  file_stub= "1step_all_policies_nonconc_overlay_models"
)

# 2-step, nonconcurrent
message("Running 2-step (policies 1–5, NONCONCURRENT, t-test vs ANOVA) ...")
dat_2_nc_list <- lapply(policies_2step, function(bp) {
  message("  2-step NONCONC, policy = ", bp, " ...")
  d_t  <- build_panel_data_2step(bp, concurrent_only = FALSE, analysis_model = "ttest")
  d_an <- build_panel_data_2step(bp, concurrent_only = FALSE, analysis_model = "lm_time")
  rbind(d_t, d_an)
})
dat_2_nc_all <- do.call(rbind, dat_2_nc_list)

dat_2_nc_all$panel  <- factor(dat_2_nc_all$panel,
                              levels = c("PBR (block factor = 1, 2-step)",
                                         "PBR (block factor = 8, 2-step)",
                                         "Big-stick (a = 2, 2-step)",
                                         "Complete randomization (2-step)"))
dat_2_nc_all$series <- factor(dat_2_nc_all$series, levels = c("A","B","FWER"))

plot_and_save_step_mode_overlay(
  dat      = dat_2_nc_all,
  step_mode= "2-step",
  policies = policies_2step,
  file_stub= "2step_all_policies_nonconc_overlay_models"
)

message("All Type I error vs allocation bias figures written to: ",
        normalizePath(out_dir, mustWork = FALSE))
