#!/usr/bin/env Rscript
# ============================================================
#  Nonconcurrent controls: data generation (chronological + allocation bias)
#  Outputs:
#    Chronological bias:
#      chronobias_nonconc_t1e_bothModels.csv
#      chronobias_nonconc_power_bothModels.csv
#
#    Allocation bias (multiple policies + PBR eta=0 reference):
#      allocbias_nonconc_t1e_preferB_bothModels.csv
#      allocbias_nonconc_power_preferB_bothModels.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

source("PT_bias/simulation_functions.R")

## ---- Global settings ----
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "1000"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE",       "20251204"))
n_cores_default <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

exp_arms <- c("A","B")

# For power
mu_null <- c(A = 0,     B = 0,     D = 0)
mu_alt  <- c(A = 0.653, B = 0.653, D = 0)

# ------------------------------------------------------------
# B-start grid for PLOTS/CSV: 1–49 (by 12)
# Internal simulation uses 0-based indexing -> subtract 1
# ------------------------------------------------------------
bstart_grid <- as.integer(seq(1, 49, by = 12))   # 1,13,25,37,49

# Chronological bias settings
beta_time_chrono <- as.numeric(Sys.getenv("BETA_TIME_CHRONO", "0.826"))
chronobias_type  <- "linear"

# Allocation bias settings
eta_alloc <- as.numeric(Sys.getenv("ETA_ALLOC", "0.0826"))

# Allocation-bias configs (policy + 1-step/2-step)
alloc_configs <- list(
  list(label = "favor_B (eta>0)",             bias_policy = "favor_B",             two_step = FALSE),
  list(label = "favor_B_2step (eta>0)",       bias_policy = "favor_B_2step",       two_step = TRUE),
  list(label = "favor_all_exp (eta>0)",       bias_policy = "favor_all_exp",       two_step = FALSE),
  list(label = "favor_all_exp_2step (eta>0)", bias_policy = "favor_all_exp_2step", two_step = TRUE)
)

# Output dir
out_dir <- Sys.getenv("OUT_DIR", unset = "PT_bias/results_nonconcurrent")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

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

# ============================================================
#  Chronological bias: Type I error & Power (nonconcurrent)
# ============================================================

run_chronobias_nonconc <- function(mu_vec, measure_label) {
  chrono_cfg <- sprintf("Chronological bias (beta=%.3f)", beta_time_chrono)
  rand_modes        <- c("block","bigstick","complete")
  ana_models        <- c("ttest","lm_time")
  block_factor_main <- 1L
  bigstick_a_main   <- 3L
  
  rows <- list()
  
  for (rm in rand_modes) {
    lab_panel <- panel_label(rm, block_factor_main, bigstick_a_main)
    
    for (am in ana_models) {
      for (i in seq_along(bstart_grid)) {
        
        bstart_disp <- bstart_grid[i]
        bstart_int  <- as.integer(bstart_disp - 1L)  # IMPORTANT: internal 0-based
        
        set.seed(seed_base +
                   1000L * i +
                   match(rm, rand_modes) * 100L +
                   match(am, ana_models))
        
        out <- calc_rejection_summary(
          n_sim           = n_sim_per_point,
          n_cores         = n_cores_default,
          max_group_size  = max_group_size,
          mu              = mu_vec,
          alpha           = alpha_one_sided,
          arm_start       = c(A = 0L, B = bstart_int),
          concurrent_only = FALSE,
          expected_total  = expected_total,
          beta_time       = beta_time_chrono,
          chronobias_type = chronobias_type,
          chronobias_incr = 0,
          rand_mode       = rm,
          block_factor    = block_factor_main,
          alloc_bias      = 0,
          exp_arms        = exp_arms,
          test_side       = test_side,
          alternative     = alternative,
          bias_policy     = "average",
          analysis_model  = am,
          two_step        = FALSE,
          bigstick_a      = bigstick_a_main
        )
        
        pr   <- out$per_comparison_rejection_rate
        fwer <- out$fwer
        
        rows[[length(rows) + 1L]] <- data.frame(
          b_start        = bstart_disp,
          series         = "A",
          rate           = unname(pr["A_vs_D"]),
          panel          = lab_panel,
          measure        = measure_label,
          analysis_model = am,
          config         = chrono_cfg,
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- data.frame(
          b_start        = bstart_disp,
          series         = "B",
          rate           = unname(pr["B_vs_D"]),
          panel          = lab_panel,
          measure        = measure_label,
          analysis_model = am,
          config         = chrono_cfg,
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- data.frame(
          b_start        = bstart_disp,
          series         = "FWER",
          rate           = fwer,
          panel          = lab_panel,
          measure        = measure_label,
          analysis_model = am,
          config         = chrono_cfg,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  
  # ------------------------------------------------------------
  # Reference (no bias): block randomization, beta_time = 0, alloc_bias = 0
  # Added so each chronological-bias CSV also contains a no-bias reference series.
  # ------------------------------------------------------------
  rm_ref      <- "block"
  lab_ref     <- panel_label(rm_ref, block_factor_main, bigstick_a_main)
  ref_cfg     <- "PBR (no bias)"
  
  for (am in ana_models) {
    for (i in seq_along(bstart_grid)) {
      
      bstart_disp <- bstart_grid[i]
      bstart_int  <- as.integer(bstart_disp - 1L)
      
      set.seed(seed_base +
                 20000L * i +
                 900L +
                 match(am, ana_models))
      
      out_ref <- calc_rejection_summary(
        n_sim           = n_sim_per_point,
        n_cores         = n_cores_default,
        max_group_size  = max_group_size,
        mu              = mu_vec,
        alpha           = alpha_one_sided,
        arm_start       = c(A = 0L, B = bstart_int),
        concurrent_only = FALSE,
        expected_total  = expected_total,
        beta_time       = 0,
        chronobias_type = chronobias_type,
        chronobias_incr = 0,
        rand_mode       = rm_ref,
        block_factor    = block_factor_main,
        alloc_bias      = 0,
        exp_arms        = exp_arms,
        test_side       = test_side,
        alternative     = alternative,
        bias_policy     = "average",
        analysis_model  = am,
        two_step        = FALSE,
        bigstick_a      = bigstick_a_main
      )
      
      pr_ref   <- out_ref$per_comparison_rejection_rate
      fwer_ref <- out_ref$fwer
      
      rows[[length(rows) + 1L]] <- data.frame(
        b_start        = bstart_disp,
        series         = "A",
        rate           = unname(pr_ref["A_vs_D"]),
        panel          = lab_ref,
        measure        = measure_label,
        analysis_model = am,
        config         = ref_cfg,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        b_start        = bstart_disp,
        series         = "B",
        rate           = unname(pr_ref["B_vs_D"]),
        panel          = lab_ref,
        measure        = measure_label,
        analysis_model = am,
        config         = ref_cfg,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        b_start        = bstart_disp,
        series         = "FWER",
        rate           = fwer_ref,
        panel          = lab_ref,
        measure        = measure_label,
        analysis_model = am,
        config         = ref_cfg,
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(rows)
}

message("Running nonconcurrent chronobias (Type I error + Power).")

dat_chrono_t1e   <- run_chronobias_nonconc(mu_null, "Type I error")
dat_chrono_power <- run_chronobias_nonconc(mu_alt,  "Power")

write.csv(dat_chrono_t1e,
          file = file.path(out_dir, "chronobias_nonconc_t1e_bothModels.csv"),
          row.names = FALSE)
write.csv(dat_chrono_power,
          file = file.path(out_dir, "chronobias_nonconc_power_bothModels.csv"),
          row.names = FALSE)

# ============================================================
#  Allocation bias: Type I error + allocation index (B vs D)
#  Includes: favor_B, favor_B_2step, favor_all_exp, favor_all_exp_2step, plus PBR eta=0 reference
# ============================================================

run_allocbias_nonconc_t1e <- function() {
  rand_modes        <- c("block","bigstick","complete")
  ana_models        <- c("ttest","lm_time")
  block_factor_main <- 1L
  bigstick_a_main   <- 3L
  
  rows <- list()
  
  for (rm in rand_modes) {
    lab_panel <- panel_label(rm, block_factor_main, bigstick_a_main)
    
    for (am in ana_models) {
      for (i in seq_along(bstart_grid)) {
        
        bstart_disp <- bstart_grid[i]
        bstart_int  <- as.integer(bstart_disp - 1L)  # IMPORTANT: internal 0-based
        
        # Biased configs
        for (cfg_i in seq_along(alloc_configs)) {
          cfg <- alloc_configs[[cfg_i]]
          
          set.seed(seed_base +
                     5000L * i +
                     cfg_i * 100000L +
                     match(rm, rand_modes) * 100L +
                     match(am, ana_models))
          
          out <- calc_rejection_summary(
            n_sim           = n_sim_per_point,
            n_cores         = n_cores_default,
            max_group_size  = max_group_size,
            mu              = mu_null,
            alpha           = alpha_one_sided,
            arm_start       = c(A = 0L, B = bstart_int),
            concurrent_only = FALSE,
            expected_total  = expected_total,
            beta_time       = 0,
            chronobias_type = "linear",
            chronobias_incr = 0,
            rand_mode       = rm,
            block_factor    = block_factor_main,
            alloc_bias      = eta_alloc,
            exp_arms        = exp_arms,
            test_side       = test_side,
            alternative     = alternative,
            bias_policy     = cfg$bias_policy,
            analysis_model  = am,
            two_step        = cfg$two_step,
            bigstick_a      = bigstick_a_main
          )
          
          pr   <- out$per_comparison_rejection_rate
          fwer <- out$fwer
          
          # allocation index for B vs D, if available
          alloc_B <- NA_real_
          ss <- out$sizes_summary
          if (!is.null(ss) && is.data.frame(ss) && "comp" %in% names(ss)) {
            row_B <- ss[ss$comp == "B_vs_D", , drop = FALSE]
            if (nrow(row_B) >= 1L) {
              arm_mean  <- row_B$arm_n_mean[1L]
              ctrl_mean <- row_B$ctrl_n_mean[1L]
              tot_mean  <- arm_mean + ctrl_mean
              if (!is.na(tot_mean) && tot_mean > 0) {
                alloc_B <- (arm_mean - ctrl_mean) / tot_mean
              }
            }
          }
          
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            analysis_model = am,
            config         = cfg$label,
            t1e_A          = unname(pr["A_vs_D"]),
            t1e_B          = unname(pr["B_vs_D"]),
            fwer           = fwer,
            alloc_bias_B   = alloc_B,
            stringsAsFactors = FALSE
          )
        }
        
        # Reference: PBR without allocation bias (eta=0), only for permuted block
        if (rm == "block") {
          set.seed(seed_base + 7000L * i + match(am, ana_models))
          
          out_ref <- calc_rejection_summary(
            n_sim           = n_sim_per_point,
            n_cores         = n_cores_default,
            max_group_size  = max_group_size,
            mu              = mu_null,
            alpha           = alpha_one_sided,
            arm_start       = c(A = 0L, B = bstart_int),
            concurrent_only = FALSE,
            expected_total  = expected_total,
            beta_time       = 0,
            chronobias_type = "linear",
            chronobias_incr = 0,
            rand_mode       = "block",
            block_factor    = block_factor_main,
            alloc_bias      = 0,
            exp_arms        = exp_arms,
            test_side       = test_side,
            alternative     = alternative,
            bias_policy     = "average",
            analysis_model  = am,
            two_step        = FALSE,
            bigstick_a      = bigstick_a_main
          )
          
          pr_ref   <- out_ref$per_comparison_rejection_rate
          fwer_ref <- out_ref$fwer
          
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            analysis_model = am,
            config         = "PBR (eta=0)",
            t1e_A          = unname(pr_ref["A_vs_D"]),
            t1e_B          = unname(pr_ref["B_vs_D"]),
            fwer           = fwer_ref,
            alloc_bias_B   = NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  dplyr::bind_rows(rows)
}

message("Running nonconcurrent allocation bias (Type I error; multiple policies).")
dat_alloc_t1e <- run_allocbias_nonconc_t1e()

write.csv(dat_alloc_t1e,
          file = file.path(out_dir, "allocbias_nonconc_t1e_preferB_bothModels.csv"),
          row.names = FALSE)

# ============================================================
#  Allocation bias: Power
#  Includes: favor_B, favor_B_2step, favor_all_exp, favor_all_exp_2step, plus PBR eta=0 reference
# ============================================================

run_allocbias_nonconc_power <- function() {
  rand_modes        <- c("block","bigstick","complete")
  ana_models        <- c("ttest","lm_time")
  block_factor_main <- 1L
  bigstick_a_main   <- 3L
  
  rows <- list()
  
  for (rm in rand_modes) {
    lab_panel <- panel_label(rm, block_factor_main, bigstick_a_main)
    
    for (am in ana_models) {
      for (i in seq_along(bstart_grid)) {
        
        bstart_disp <- bstart_grid[i]
        bstart_int  <- as.integer(bstart_disp - 1L)  # IMPORTANT: internal 0-based
        
        # Biased configs
        for (cfg_i in seq_along(alloc_configs)) {
          cfg <- alloc_configs[[cfg_i]]
          
          set.seed(seed_base +
                     9000L * i +
                     cfg_i * 100000L +
                     match(rm, rand_modes) * 100L +
                     match(am, ana_models))
          
          out_biased <- calc_rejection_summary(
            n_sim           = n_sim_per_point,
            n_cores         = n_cores_default,
            max_group_size  = max_group_size,
            mu              = mu_alt,
            alpha           = alpha_one_sided,
            arm_start       = c(A = 0L, B = bstart_int),
            concurrent_only = FALSE,
            expected_total  = expected_total,
            beta_time       = 0,
            chronobias_type = "linear",
            chronobias_incr = 0,
            rand_mode       = rm,
            block_factor    = block_factor_main,
            alloc_bias      = eta_alloc,
            exp_arms        = exp_arms,
            test_side       = test_side,
            alternative     = alternative,
            bias_policy     = cfg$bias_policy,
            analysis_model  = am,
            two_step        = cfg$two_step,
            bigstick_a      = bigstick_a_main
          )
          
          pr_biased   <- out_biased$per_comparison_rejection_rate
          fwer_biased <- out_biased$fwer
          
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            config         = cfg$label,
            series         = "A",
            analysis_model = am,
            rate           = unname(pr_biased["A_vs_D"]),
            stringsAsFactors = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            config         = cfg$label,
            series         = "B",
            analysis_model = am,
            rate           = unname(pr_biased["B_vs_D"]),
            stringsAsFactors = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            config         = cfg$label,
            series         = "FWER",
            analysis_model = am,
            rate           = fwer_biased,
            stringsAsFactors = FALSE
          )
        }
        
        # Reference only for permuted block without bias
        if (rm == "block") {
          set.seed(seed_base + 12000L * i + match(am, ana_models))
          
          out_ref <- calc_rejection_summary(
            n_sim           = n_sim_per_point,
            n_cores         = n_cores_default,
            max_group_size  = max_group_size,
            mu              = mu_alt,
            alpha           = alpha_one_sided,
            arm_start       = c(A = 0L, B = bstart_int),
            concurrent_only = FALSE,
            expected_total  = expected_total,
            beta_time       = 0,
            chronobias_type = "linear",
            chronobias_incr = 0,
            rand_mode       = "block",
            block_factor    = block_factor_main,
            alloc_bias      = 0,
            exp_arms        = exp_arms,
            test_side       = test_side,
            alternative     = alternative,
            bias_policy     = "average",
            analysis_model  = am,
            two_step        = FALSE,
            bigstick_a      = bigstick_a_main
          )
          
          pr_ref   <- out_ref$per_comparison_rejection_rate
          fwer_ref <- out_ref$fwer
          
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            config         = "PBR (eta=0)",
            series         = "A",
            analysis_model = am,
            rate           = unname(pr_ref["A_vs_D"]),
            stringsAsFactors = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            config         = "PBR (eta=0)",
            series         = "B",
            analysis_model = am,
            rate           = unname(pr_ref["B_vs_D"]),
            stringsAsFactors = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            b_start        = bstart_disp,
            panel          = lab_panel,
            config         = "PBR (eta=0)",
            series         = "FWER",
            analysis_model = am,
            rate           = fwer_ref,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  dplyr::bind_rows(rows)
}

message("Running nonconcurrent allocation bias (Power; multiple policies).")
dat_alloc_power <- run_allocbias_nonconc_power()

write.csv(dat_alloc_power,
          file = file.path(out_dir, "allocbias_nonconc_power_preferB_bothModels.csv"),
          row.names = FALSE)

message("All nonconcurrent data written to: ",
        normalizePath(out_dir, mustWork = FALSE))
