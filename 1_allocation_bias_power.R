#!/usr/bin/env Rscript
# ============================================================
#  Allocation bias (concurrent controls) – POWER
#  Same sweep/grid as 1_allocation_bias_... (Type I error),
#  but under alternative with effect size 0.653 for both arms.
#
#  Output (CSV):
#    power_vs_alloc_concurrent_threeScenarios_steps.csv
#
#  Notes for plotting:
#    - y-axis should be 0..1
#    - y-label should be "Power"
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

source("PT_bias/simulation_functions.R")

# ---------- global settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

exp_arms <- c("A","B")

# Alternative (power): effect size 0.653 for both experimental arms, control D = 0
mu_alt <- c(A = 0.653, B = 0.653, D = 0)

# Sweep grid for allocation bias strength (lambda/eta in your code)
eta_grid <- seq(0, 0.20, by = 0.04)

# B-start scenarios (as used in your concurrent plots)
bstart_scenarios <- c(1L, 25L, 49L)
scenario_labels  <- c(
  "Arm B start at 1",
  "Arm B start at 25",
  "Arm B start at 49"
)
names(scenario_labels) <- as.character(bstart_scenarios)

# Biasing policies (base)
bias_policies_base <- c("favor_all_exp","favor_B")

# Randomization procedures
rand_modes <- c("block","bigstick","complete")

rand_label <- function(rm, block_factor = 1L, bigstick_a = 3L) {
  if (rm == "block") {
    sprintf("Permuted block randomization (factor = %d)", block_factor)
  } else if (rm == "bigstick") {
    sprintf("Big stick design (a = %d)", bigstick_a)
  } else if (rm == "complete") {
    "Complete randomization"
  } else {
    rm
  }
}

# Step schemes (same as 1_)
scheme_df <- data.frame(
  scheme_id    = c("step1", "step2_total", "step2_cohort"),
  scheme_label = c(
    "1-step randomization",
    "2-step randomization (total counts only)",
    "2-step randomization (total counts by cohort)"
  ),
  two_step     = c(FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Helper: map base bias policy + scheme to actual bias_policy string
get_bias_policy_full <- function(bias_base, scheme_id) {
  if (scheme_id == "step2_cohort") {
    # 2-step: total counts by cohort → *_2step
    if (bias_base == "favor_all_exp") {
      return("favor_all_exp_2step")
    } else if (bias_base == "favor_B") {
      return("favor_B_2step")
    } else {
      stop("Unknown bias_base for step2_cohort: ", bias_base)
    }
  } else {
    # 1-step and 2-step (total counts only) use the same policy name
    return(bias_base)
  }
}

# Output dir (same default as 1_)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = "PT_bias/results_allocbias_concurrent_threeScenarios"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# Reproducibility controls
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "1000"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE", "20251204"))
n_cores_default <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# ============================================================
#  Main sweep (POWER)
# ============================================================

rows_all <- list()

for (s_idx in seq_along(bstart_scenarios)) {
  b_start  <- bstart_scenarios[s_idx]
  scen_lab <- scenario_labels[as.character(b_start)]
  
  arm_start <- c(A = 0L, B = b_start)
  
  message("=== ", scen_lab, " (B start = ", b_start, ") ===")
  
  for (sch_idx in seq_len(nrow(scheme_df))) {
    scheme_id    <- scheme_df$scheme_id[sch_idx]
    scheme_label <- scheme_df$scheme_label[sch_idx]
    two_step     <- scheme_df$two_step[sch_idx]
    
    message("  Scheme: ", scheme_label, " [", scheme_id, "]")
    
    for (bias_base in bias_policies_base) {
      bias_full <- get_bias_policy_full(bias_base, scheme_id)
      
      message("    Bias base: ", bias_base, " -> ", bias_full)
      
      for (rm in rand_modes) {
        rlab <- rand_label(rm, block_factor = 1L, bigstick_a = 3L)
        
        for (eta in eta_grid) {
          # deterministic seeding across grid points
          set.seed(seed_base +
                     100000L * s_idx +
                     10000L  * sch_idx +
                     1000L   * match(bias_base, bias_policies_base) +
                     100L    * match(rm, rand_modes) +
                     as.integer(round(eta * 1000)))
          
          out <- calc_rejection_summary(
            n_sim           = n_sim_per_point,
            n_cores         = n_cores_default,
            max_group_size  = max_group_size,
            mu              = mu_alt,
            alpha           = alpha_one_sided,
            arm_start       = arm_start,
            concurrent_only = TRUE,
            expected_total  = expected_total,
            
            # no time trend for allocation-bias plots
            beta_time       = 0,
            chronobias_type = "linear",
            chronobias_incr = 0,
            
            rand_mode       = rm,
            block_factor    = 1L,
            bigstick_a      = 3L,
            
            alloc_bias      = eta,
            exp_arms        = exp_arms,
            
            test_side       = test_side,
            alternative     = alternative,
            
            bias_policy     = bias_full,
            analysis_model  = "ttest",
            two_step        = two_step
          )
          
          pr   <- out$per_comparison_rejection_rate
          fwer <- out$fwer
          
          rate_A <- unname(pr["A_vs_D"])
          rate_B <- unname(pr["B_vs_D"])
          
          # Arm A
          rows_all[[length(rows_all) + 1L]] <- data.frame(
            scenario_b_start  = b_start,
            scenario_label    = scen_lab,
            scheme_id         = scheme_id,
            scheme_label      = scheme_label,
            rand_mode         = rm,
            rand_label        = rlab,
            bias_policy_base  = bias_base,
            bias_policy_full  = bias_full,
            alloc_bias        = eta,
            series            = "A",
            rate              = rate_A,
            stringsAsFactors  = FALSE
          )
          
          # Arm B
          rows_all[[length(rows_all) + 1L]] <- data.frame(
            scenario_b_start  = b_start,
            scenario_label    = scen_lab,
            scheme_id         = scheme_id,
            scheme_label      = scheme_label,
            rand_mode         = rm,
            rand_label        = rlab,
            bias_policy_base  = bias_base,
            bias_policy_full  = bias_full,
            alloc_bias        = eta,
            series            = "B",
            rate              = rate_B,
            stringsAsFactors  = FALSE
          )
          
          # Disjunctive power (kept as "FWER" series for compatibility with your facet logic)
          rows_all[[length(rows_all) + 1L]] <- data.frame(
            scenario_b_start  = b_start,
            scenario_label    = scen_lab,
            scheme_id         = scheme_id,
            scheme_label      = scheme_label,
            rand_mode         = rm,
            rand_label        = rlab,
            bias_policy_base  = bias_base,
            bias_policy_full  = bias_full,
            alloc_bias        = eta,
            series            = "FWER",
            rate              = fwer,
            stringsAsFactors  = FALSE
          )
        }
      }
    }
  }
}

dat_all <- dplyr::bind_rows(rows_all)

outfile <- file.path(out_dir, "power_vs_alloc_concurrent_threeScenarios_steps.csv")
write.csv(dat_all, outfile, row.names = FALSE)

message("Done. Data written to: ", normalizePath(outfile, mustWork = FALSE))
