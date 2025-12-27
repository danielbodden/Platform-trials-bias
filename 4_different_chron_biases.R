#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)  # not strictly needed here, but harmless
})

# ---------- cluster knobs ----------
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "1000"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE",       "20251126"))

# ---------- fixed trial settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

exp_arms <- c("A","B")
mu_null  <- c(A = 0, B = 0, D = 0)

# >>> requested: standard chronological bias magnitude
chrono_beta_val <- 0.826

# Sweep grid for timing of Arm B opening (patients)
b_grid <- seq(0, 48, by = 12)

# Randomization parameters (aligned with other scripts)
block_factor_main <- 1L
bigstick_a_main   <- 3L

rand_label <- function(rm,
                       block_factor = block_factor_main,
                       bigstick_a   = bigstick_a_main) {
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

# I/O
out_dir <- Sys.getenv("OUT_DIR", unset = "PT_bias/results_t1e_vs_Bstart_chrono_grid")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# --- simulator ---
source("PT_bias/simulation_functions.R")

# ============================================================
# Chronological bias TYPES and mapping to parameters
# ============================================================
chrono_types <- c("none", "linear", "stepwise", "inv_u", "seasonal")

chrono_label <- function(key, beta = chrono_beta_val) {
  switch(key,
         none     = "No chronological bias",
         linear   = sprintf("Linear trend (beta = %.3f)", beta),
         stepwise = sprintf("Stepwise trend (beta = %.3f)", beta),
         inv_u    = sprintf("Inverted-U trend (beta = %.3f)", beta),
         seasonal = sprintf("Seasonal trend (beta = %.3f)", beta)
  )
}

get_chrono_params_fixed <- function(key) {
  beta <- chrono_beta_val
  
  if (key == "none") {
    return(list(
      chronobias_type = "linear",
      beta_time       = 0,
      chronobias_incr = 0
    ))
  }
  if (key == "linear") {
    return(list(
      chronobias_type = "linear",
      beta_time       = beta,
      chronobias_incr = 0
    ))
  }
  if (key == "stepwise") {
    step_pattern <- c(0, beta/3, 2*beta/3, beta)
    return(list(
      chronobias_type = "stepwise",
      beta_time       = 0,
      chronobias_incr = rep(step_pattern, length.out = 200L)
    ))
  }
  if (key == "inv_u") {
    return(list(
      chronobias_type = "inv_u",
      beta_time       = beta,
      chronobias_incr = 0
    ))
  }
  if (key == "seasonal") {
    return(list(
      chronobias_type = "seasonal",
      beta_time       = beta,
      chronobias_incr = 0
    ))
  }
  stop("Unknown chrono bias type: ", key)
}

# ============================================================
# Core sweep over b_grid for one rand_mode × chrono_type × model
# ============================================================
.run_sweep_Bstart_type_model <- function(
    rand_mode       = c("complete","block","bigstick"),
    block_factor    = 1L,
    rand_proc_label = "Complete randomization",
    chrono_key      = "linear",
    analysis_model  = c("ttest","lm_time"),
    concurrent_only = FALSE,
    two_step        = FALSE,
    seed_bump       = 0L
) {
  rand_mode      <- match.arg(rand_mode)
  analysis_model <- match.arg(analysis_model, c("ttest","lm_time"))
  
  rows <- lapply(seq_along(b_grid), function(i) {
    bs <- b_grid[i]
    params <- get_chrono_params_fixed(chrono_key)
    
    set.seed(seed_base + seed_bump + i +
               1000L * block_factor +
               ifelse(rand_mode == "complete", 9999L, 0L) +
               ifelse(two_step, 5000L, 0L) +
               ifelse(analysis_model == "lm_time", 7777L, 0L))
    
    out <- calc_rejection_summary(
      n_sim           = n_sim_per_point,
      n_cores         = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      max_group_size  = max_group_size,
      mu              = mu_null,
      alpha           = alpha_one_sided,
      arm_start       = c(A = 0, B = bs),
      concurrent_only = concurrent_only,
      expected_total  = expected_total,
      beta_time       = params$beta_time,
      chronobias_type = params$chronobias_type,
      chronobias_incr = params$chronobias_incr,
      rand_mode       = rand_mode,
      block_factor    = block_factor,
      alloc_bias      = 0,
      exp_arms        = exp_arms,
      test_side       = test_side,
      alternative     = alternative,
      bias_policy     = "favor_B",
      analysis_model  = analysis_model,
      two_step        = two_step,
      bigstick_a      = bigstick_a_main
    )
    
    pr <- out$per_comparison_rejection_rate
    want <- intersect(c("A_vs_D","B_vs_D"), names(pr))
    
    df_ab <- data.frame(
      B_start        = bs,
      series         = sub("_vs_.*$", "", want),  # "A","B"
      t1e            = as.numeric(pr[want]),
      chrono_key     = chrono_key,
      rand_proc      = rand_proc_label,
      analysis_model = analysis_model,
      stringsAsFactors = FALSE
    )
    
    df_fwer <- data.frame(
      B_start        = bs,
      series         = "FWER",
      t1e            = out$fwer,
      chrono_key     = chrono_key,
      rand_proc      = rand_proc_label,
      analysis_model = analysis_model,
      stringsAsFactors = FALSE
    )
    
    rbind(df_ab, df_fwer)
  })
  
  do.call(rbind, rows)
}

build_panel_data_Bstart_series_by_rand_model <- function(concurrent_only = FALSE) {
  dat_list <- list()
  
  for (ck in chrono_types) {
    
    run_for_model <- function(analysis_model, model_bump) {
      d_block <- .run_sweep_Bstart_type_model(
        rand_mode       = "block",
        block_factor    = block_factor_main,
        rand_proc_label = rand_label("block"),
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 10L + model_bump
      )
      
      d_bs <- .run_sweep_Bstart_type_model(
        rand_mode       = "bigstick",
        block_factor    = block_factor_main,
        rand_proc_label = rand_label("bigstick"),
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 20L + model_bump
      )
      
      d_comp <- .run_sweep_Bstart_type_model(
        rand_mode       = "complete",
        block_factor    = block_factor_main,
        rand_proc_label = rand_label("complete"),
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 0L + model_bump
      )
      
      rbind(d_block, d_bs, d_comp)
    }
    
    dat_ttest <- run_for_model("ttest",   model_bump = 0L)
    dat_lm    <- run_for_model("lm_time", model_bump = 1000L)
    
    dat_list[[ck]] <- rbind(dat_ttest, dat_lm)
  }
  
  dat <- do.call(rbind, dat_list)
  
  dat$series <- factor(dat$series, levels = c("A","B","FWER"))
  dat$chronobias <- factor(
    vapply(dat$chrono_key, chrono_label, character(1)),
    levels = vapply(chrono_types, chrono_label, character(1))
  )
  dat$rand_proc <- factor(
    dat$rand_proc,
    levels = c(rand_label("block"), rand_label("bigstick"), rand_label("complete"))
  )
  dat$analysis_model <- factor(
    dat$analysis_model,
    levels = c("ttest","lm_time"),
    labels = c("t-test", "ANOVA (lm_time)")
  )
  
  dat
}

# ============================================================
# Run + save
# ============================================================
message("Generating data (chrono beta = ", chrono_beta_val, ") ...")

dat_Bstart <- build_panel_data_Bstart_series_by_rand_model(concurrent_only = FALSE)

file_stub <- Sys.getenv("FILE_STUB", "t1e_vs_Bstart_chronobias_beta0826_by_rand_and_model")
csv_path  <- file.path(out_dir, paste0(file_stub, ".csv"))
rds_path  <- file.path(out_dir, paste0(file_stub, ".rds"))

write.csv(dat_Bstart, csv_path, row.names = FALSE)
saveRDS(dat_Bstart, rds_path)

message("Wrote:")
message("  CSV: ", csv_path)
message("  RDS: ", rds_path)
