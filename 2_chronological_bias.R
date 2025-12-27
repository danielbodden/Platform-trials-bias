#!/usr/bin/env Rscript
# ============================================================
#  Chronological bias (linear time trend), concurrent controls
#  Three B-start scenarios:
#    1) Arm B starts at 0
#    2) Arm B starts at 25 patients
#    3) Arm B starts at 48 patients
#
#  NO allocation bias (alloc_bias = 0).
#
#  For each scenario, we compute:
#    - Type I error (null: mu_A = mu_B = mu_D = 0)
#    - Power (alternative: mu_A = mu_B = effect_size, mu_D = 0)
#
#  Randomization procedures:
#    - Permuted block randomization (factor = 1)
#    - Big stick design (a = 3)
#    - Complete randomization
#
#  Output:
#    PT_bias/results_chronobias_threeScenarios/
#      chronobias_concurrent_threeScenarios_t1e_power.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ---------- cluster knobs ----------
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "2000"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE",       "20251205"))
n_cores_default <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# ---------- trial settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

# effect size as in your previous chronobias script
effect_size     <- 0.653

exp_arms <- c("A", "B")

# B-start scenarios (patients)
bstart_scenarios <- c(1L, 25L, 49L)
scenario_labels  <- c(
  "Scenario 1: Arm B starts at 1 patients",
  "Scenario 2: Arm B starts at 25 patients",
  "Scenario 3: Arm B starts at 49 patients"
)
names(scenario_labels) <- as.character(bstart_scenarios)

# Chronological bias strength grid (linear trend)
chrono_strength_grid <- seq(0, 6, by = 1)

# I/O
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = "PT_bias/results_chronobias_threeScenarios"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# Use same simulator as your existing code
source("PT_bias/simulation_functions.R")

# ---------- helper: panel labels ----------
panel_label <- function(rand_mode, block_factor = 1L, bigstick_a = 3L) {
  if (rand_mode == "block") {
    sprintf("Permuted block randomization (factor = %d)", block_factor)
  } else if (rand_mode == "bigstick") {
    sprintf("Big stick design (a = %d)", bigstick_a)
  } else if (rand_mode == "complete") {
    "Complete randomization"
  } else {
    rand_mode
  }
}

# ---------- sweep over chrono_strength for one scenario + rand_mode ----------

run_sweep_chrono <- function(mu,
                             arm_start,
                             rand_mode      = c("block", "bigstick", "complete"),
                             block_factor   = 1L,
                             bigstick_a     = 3L,
                             two_step       = TRUE,
                             seed_bump      = 0L) {
  rand_mode <- match.arg(rand_mode)
  
  rows <- lapply(seq_along(chrono_strength_grid), function(i) {
    strength <- chrono_strength_grid[i]
    
    beta_time       <- strength
    chronobias_type <- "linear"
    chronobias_incr <- 0
    
    set.seed(seed_base + seed_bump + i +
               1000L * block_factor +
               ifelse(rand_mode == "complete", 9999L, 0L) +
               ifelse(two_step, 5000L, 0L))
    
    out <- calc_rejection_summary(
      n_sim           = n_sim_per_point,
      n_cores         = n_cores_default,
      max_group_size  = max_group_size,
      mu              = mu,
      alpha           = alpha_one_sided,
      arm_start       = arm_start,
      concurrent_only = TRUE,
      expected_total  = expected_total,
      beta_time       = beta_time,
      chronobias_type = chronobias_type,
      chronobias_incr = chronobias_incr,
      rand_mode       = rand_mode,
      block_factor    = block_factor,
      alloc_bias      = 0,          # NO allocation bias
      exp_arms        = exp_arms,
      test_side       = test_side,
      alternative     = alternative,
      bias_policy     = "favor_B",  # irrelevant since alloc_bias=0
      analysis_model  = "ttest",
      two_step        = two_step,
      bigstick_a      = bigstick_a
    )
    
    pr   <- out$per_comparison_rejection_rate
    want <- intersect(c("A_vs_D", "B_vs_D"), names(pr))
    
    df_ab <- data.frame(
      chrono_strength = strength,
      series          = sub("_vs_.*$", "", want),  # "A","B"
      rate            = as.numeric(pr[want]),
      stringsAsFactors = FALSE
    )
    
    df_fwer <- data.frame(
      chrono_strength = strength,
      series          = "FWER",
      rate            = out$fwer,
      stringsAsFactors = FALSE
    )
    
    rbind(df_ab, df_fwer)
  })
  
  dat <- do.call(rbind, rows)
  dat$series <- factor(dat$series, levels = c("A", "B", "FWER"))
  dat
}

# ---------- build data for one scenario & measure ----------

build_data_for_scenario <- function(mu, measure_label, b_start) {
  arm_start <- c(A = 0L, B = b_start)
  scen_lab  <- scenario_labels[as.character(b_start)]
  
  # Permuted block randomization
  d_pbr <- run_sweep_chrono(
    mu          = mu,
    arm_start   = arm_start,
    rand_mode   = "block",
    block_factor= 1L,
    bigstick_a  = 3L,
    two_step    = TRUE,
    seed_bump   = 10L
  )
  d_pbr$panel <- panel_label("block", block_factor = 1L, bigstick_a = 3L)
  
  # Big stick
  d_bs <- run_sweep_chrono(
    mu          = mu,
    arm_start   = arm_start,
    rand_mode   = "bigstick",
    block_factor= 1L,
    bigstick_a  = 3L,
    two_step    = TRUE,
    seed_bump   = 20L
  )
  d_bs$panel <- panel_label("bigstick", block_factor = 1L, bigstick_a = 3L)
  
  # Complete randomization
  d_cr <- run_sweep_chrono(
    mu          = mu,
    arm_start   = arm_start,
    rand_mode   = "complete",
    block_factor= 1L,
    bigstick_a  = 3L,
    two_step    = TRUE,
    seed_bump   = 30L
  )
  d_cr$panel <- panel_label("complete", block_factor = 1L, bigstick_a = 3L)
  
  dat <- bind_rows(d_pbr, d_bs, d_cr)
  dat$measure        <- measure_label
  dat$scenario_b_start <- b_start
  dat$scenario_label   <- scen_lab
  
  dat
}

# ============================================================
#  Main: build Type I error and Power for all scenarios
# ============================================================

mu_null <- c(A = 0, B = 0, D = 0)
mu_alt  <- c(A = effect_size, B = effect_size, D = 0)

rows_t1e   <- list()
rows_power <- list()

message("Running chronological bias sweeps for three B-start scenarios ...")

for (b_start in bstart_scenarios) {
  message("  Scenario with B start = ", b_start, " ...")
  
  rows_t1e[[length(rows_t1e) + 1L]] <-
    build_data_for_scenario(mu = mu_null, measure_label = "Type I error", b_start = b_start)
  
  rows_power[[length(rows_power) + 1L]] <-
    build_data_for_scenario(mu = mu_alt, measure_label = "Power",        b_start = b_start)
}

dat_t1e   <- bind_rows(rows_t1e)
dat_power <- bind_rows(rows_power)

dat_all <- bind_rows(dat_t1e, dat_power)

outfile <- file.path(out_dir, "chronobias_concurrent_threeScenarios_t1e_power.csv")
write.csv(dat_all, outfile, row.names = FALSE)

message("Done. Data written to: ", normalizePath(outfile, mustWork = FALSE))
