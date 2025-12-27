#!/usr/bin/env Rscript
# ============================================================
#  Power vs timing of Arm B (A,B vs D) with chronological bias
#
#  Layout:
#    Rows    = Series (Arm A, Arm B, FWER)
#    Columns = Chronological bias TYPE
#    Color   = Randomization procedure
#    Linetype = Analysis model (t-test vs ANOVA [lm_time])
#
#  Scenario:
#    - Power: mu = (0.6, 0.6, 0) for (A, B, D)
#    - Allocation bias: 0
#    - Chronological bias: fixed strength beta = 0.16
#    - alpha (one-sided, greater) = 0.025
#    - max_group_size  = 24 (per arm)
#    - expected_total  = 96
#    - arm_start: A = 0, B varies over b_grid
#    - concurrent_only = FALSE (nonconcurrent analyses)
#    - randomization procedures:
#        * Block (1*arms, 1-step)
#        * Block (8*arms, 1-step)
#        * Big-stick (a = 2, 1-step)
#        * Complete randomization (1-step)
#
#  Time trend parameterization (fixed size = 0.16):
#    Let beta := 0.16.
#    - None:
#        chronobias_type = "linear"
#        beta_time       = 0
#        chronobias_incr = 0
#    - Linear:
#        chronobias_type = "linear"
#        beta_time       = beta
#        chronobias_incr = 0
#    - Stepwise:
#        chronobias_type = "stepwise"
#        beta_time       = 0
#        chronobias_incr = repeating pattern (0, beta/3, 2*beta/3, beta)
#    - Inverted-U:
#        chronobias_type = "inv_u"
#        beta_time       = beta
#        chronobias_incr = 0
#    - Seasonal:
#        chronobias_type = "seasonal"
#        beta_time       = beta
#        chronobias_incr = 0
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------- cluster knobs ----------
n_sim_per_point <- as.integer(Sys.getenv("N_SIM_PER_POINT", "20000"))
seed_base       <- as.integer(Sys.getenv("SEED_BASE",       "20251126"))

# ---------- fixed trial settings ----------
alpha_one_sided <- 0.025
test_side       <- "one.sided"
alternative     <- "greater"

max_group_size  <- 24
expected_total  <- 96

exp_arms <- c("A","B")
mu_alt   <- c(A = 0.6, B = 0.6, D = 0)

# Chronological bias magnitude (same as chrono_beta_val in existing code)
chrono_beta_val <- 0.16

# Sweep grid for timing of Arm B opening (patients)
# (same as in plots_concurrent_by_scenario.R)
b_grid <- seq(0, 48, by = 12)

# I/O
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = "PT_bias/results_power_vs_Bstart_chrono_grid"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Writing outputs to: ", normalizePath(out_dir, mustWork = FALSE))

# --- robust UTF-8 plotting on headless nodes (Cairo devices) ---
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
use_cairo <- TRUE
pdf_device <- function(file, ...) {
  if (use_cairo) grDevices::cairo_pdf(file = file, ...) else grDevices::pdf(file = file, ...)
}

# --- simulator ---
source("PT_bias/simulation_functions.R")

# ============================================================
# Chronological bias TYPES and mapping to parameters
# ============================================================

chrono_types <- c("none", "linear", "stepwise", "inv_u", "seasonal")

chrono_label <- function(key) {
  switch(key,
         none     = "No chronological bias",
         linear   = "Linear trend (beta = 0.16)",
         stepwise = "Stepwise trend (beta = 0.16)",
         inv_u    = "Inverted-U trend (beta = 0.16)",
         seasonal = "Seasonal trend (beta = 0.16)")
}

get_chrono_params_fixed <- function(key) {
  # Returns list(chronobias_type, beta_time, chronobias_incr)
  # for given bias type, with fixed size = chrono_beta_val.
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
# analysis_model ∈ {"ttest","lm_time"}
# ============================================================

.run_sweep_Bstart_type_model <- function(
    rand_mode       = c("complete","block","bigstick"),
    block_factor    = 1L,
    rand_label      = "Complete randomization",
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
      mu              = mu_alt,
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
      bigstick_a      = 2L
    )
    
    pr   <- out$per_comparison_rejection_rate
    want <- intersect(c("A_vs_D","B_vs_D"), names(pr))
    
    # A and B power
    df_ab <- data.frame(
      B_start        = bs,
      series         = sub("_vs_.*$", "", want),  # "A","B"
      power          = as.numeric(pr[want]),
      chrono_key     = chrono_key,
      rand_proc      = rand_label,
      analysis_model = analysis_model,
      stringsAsFactors = FALSE
    )
    # "FWER" under the alternative: prob(at least one rejection)
    df_fwer <- data.frame(
      B_start        = bs,
      series         = "FWER",
      power          = out$fwer,
      chrono_key     = chrono_key,
      rand_proc      = rand_label,
      analysis_model = analysis_model,
      stringsAsFactors = FALSE
    )
    
    rbind(df_ab, df_fwer)
  })
  
  do.call(rbind, rows)
}

# ============================================================
# Assemble data across chrono TYPES × rand modes × models
#   Rows   = series (A, B, FWER)
#   Columns= chrono types
#   Color  = randomization procedure
#   Linetype = analysis model
# ============================================================

build_panel_data_Bstart_series_by_rand_model <- function(concurrent_only = FALSE) {
  dat_list <- list()
  
  for (ck in chrono_types) {
    
    run_for_model <- function(analysis_model, model_bump) {
      # Block (1*arms)
      d_b1 <- .run_sweep_Bstart_type_model(
        rand_mode       = "block",
        block_factor    = 1L,
        rand_label      = "Block randomization (1*arms, 1-step)",
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 10L + model_bump
      )
      
      # Block (8*arms)
      d_b8 <- .run_sweep_Bstart_type_model(
        rand_mode       = "block",
        block_factor    = 8L,
        rand_label      = "Block randomization (8*arms, 1-step)",
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 40L + model_bump
      )
      
      # Big-stick
      d_bs <- .run_sweep_Bstart_type_model(
        rand_mode       = "bigstick",
        block_factor    = 1L,
        rand_label      = "Big-stick design (a = 2, 1-step)",
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 20L + model_bump
      )
      
      # Complete randomization
      d_comp <- .run_sweep_Bstart_type_model(
        rand_mode       = "complete",
        block_factor    = 1L,
        rand_label      = "Complete randomization (1-step)",
        chrono_key      = ck,
        analysis_model  = analysis_model,
        concurrent_only = concurrent_only,
        two_step        = FALSE,
        seed_bump       = 0L + model_bump
      )
      
      rbind(d_b1, d_b8, d_bs, d_comp)
    }
    
    dat_ttest <- run_for_model("ttest",   model_bump = 0L)
    dat_lm    <- run_for_model("lm_time", model_bump = 1000L)
    
    dat_ck <- rbind(dat_ttest, dat_lm)
    dat_list[[ck]] <- dat_ck
  }
  
  dat <- do.call(rbind, dat_list)
  
  dat$series <- factor(dat$series, levels = c("A","B","FWER"))
  dat$chronobias <- factor(
    vapply(dat$chrono_key, chrono_label, character(1)),
    levels = vapply(chrono_types, chrono_label, character(1))
  )
  dat$rand_proc <- factor(
    dat$rand_proc,
    levels = c(
      "Block randomization (1*arms, 1-step)",
      "Block randomization (8*arms, 1-step)",
      "Big-stick design (a = 2, 1-step)",
      "Complete randomization (1-step)"
    )
  )
  dat$analysis_model <- factor(
    dat$analysis_model,
    levels = c("ttest","lm_time"),
    labels = c("t-test", "ANOVA (lm_time)")
  )
  
  dat
}

# ============================================================
# Plotting
#   facet_grid(series ~ chronobias)
#   color    = randomization procedure
#   linetype = analysis model
# ============================================================

plot_and_save_Bstart_series_by_rand_model <- function(dat, file_stub) {
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  
  write.csv(dat, csv_path, row.names = FALSE)
  
  p <- ggplot(dat, aes(x = B_start, y = power,
                       color = rand_proc,
                       linetype = analysis_model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 48)) +
    facet_grid(series ~ chronobias) +
    labs(
      title    = "Power vs timing of Arm B (chronological bias beta = 0.16)",
      subtitle = paste(
        "Effect size: A = B = 0.6 vs D = 0.",
        "Rows: Arm A, Arm B, FWER; Columns: time-trend types.",
        "Curves: randomization procedures; Linetype: t-test vs ANOVA (lm_time)."
      ),
      x        = "Timing of addition of Arm B (patients)",
      y        = "Power",
      color    = "Randomization procedure",
      linetype = "Analysis model",
      caption  = paste(
        "Chronological bias parameterization with size beta = 0.16:",
        "Linear / Inverted-U / Seasonal use beta_time = 0.16;",
        "Stepwise uses repeating increments (0, 0.16/3, 2×0.16/3, 0.16);",
        "None uses beta_time = 0 and no increments."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box      = "vertical",
      strip.text      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold")
    )
  
  pdf_device(file = pdf_path, width = 12, height = 8.5)
  print(p)
  dev.off()
  
  message("Wrote:")
  message("  Data: ", csv_path)
  message("  PDF : ", pdf_path)
}

# ============================================================
# Run
# ============================================================

message("Running power vs B-start (chrono beta = 0.16, series rows, time trend columns, rand & model as curves) ...")

dat_Bstart <- build_panel_data_Bstart_series_by_rand_model(concurrent_only = FALSE)

plot_and_save_Bstart_series_by_rand_model(
  dat       = dat_Bstart,
  file_stub = "power_vs_Bstart_chronobias_beta016_by_rand_and_model_es06"
)

message("Done. Written to: ",
        normalizePath(out_dir, mustWork = FALSE))
