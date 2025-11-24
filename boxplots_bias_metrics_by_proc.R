#!/usr/bin/env Rscript
# ============================================================
# FIGURE: Mean Bias Differences (Arm − control, by control set)
# - X: Randomization procedure
# - Rows: Arm A, Arm B
# - Columns: Allocation bias, Chronological bias  × Control set
#   Control sets:
#     • Concurrent only
#     • Nonconcurrent
# - Thin dashed line at y=0; colorful violins + boxplots
# - Value per run: mean(arm) − mean(control)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(parallel)
  library(grid)
})

# ---------- knobs ----------
n_runs    <- as.integer(Sys.getenv("N_RUNS",    "3000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", "20251107"))
out_dir   <- Sys.getenv("OUT_DIR", unset = "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

# Trial constants
max_group_size_exp <- 24L
B_start_fixed      <- 16L
alloc_bias_val     <- 0.08
chrono_beta_val    <- 0.16
expected_total     <- 96L  # <- ensure this exists in workers

# ---------- simulator ----------
if (file.exists("simulation_functions.R")) source("simulation_functions.R")
if (!exists("platform_trials_simulation"))
  stop("platform_trials_simulation() not found. Please provide it via simulation_functions.R")

# ---------- randomization procedures ----------
# Order: Block(1*arms), Block(8*arms), Big-stick(a=2), Complete
procedures <- list(
  list(key = "Block randomization (1*arms)", rand_mode = "block",    block_factor = 1L),
  list(key = "Block randomization (8*arms)", rand_mode = "block",    block_factor = 8L),
  list(key = "Big-stick design (a = 2)",     rand_mode = "bigstick", block_factor = 1L),
  list(key = "Complete randomization",       rand_mode = "complete", block_factor = 1L)
)

# ---------- one run -> mean bias differences (Arm − control) ----------
# Produces rows for both control sets:
#   1) Concurrent only:        D in [t_open(arm), t_close(arm)]
#   2) Concurrent + prior:     D with t <= t_close(arm)
# Value in each row: mean(arm) − mean(control)
one_run_extract <- function(rand_mode, block_factor, rep_seed,
                            B_start_fixed, max_group_size_exp,
                            alloc_bias_val, chrono_beta_val,
                            two_step) {
  
  set.seed(rep_seed)
  res <- platform_trials_simulation(
    exp_arms        = c("A","B"),
    arm_start       = c(A = 0L, B = B_start_fixed),
    max_group_size  = max_group_size_exp,
    expected_total  = expected_total,
    alpha           = 0.025,
    test_side       = "one.sided",
    alternative     = "greater",
    rand_mode       = rand_mode,
    block_factor    = block_factor,
    alloc_bias      = alloc_bias_val,
    beta_time       = chrono_beta_val,
    concurrent_only = TRUE,        # analysis setting; we compute control sets below
    analysis_model  = "ttest",
    bias_policy     = "favor_all_exp",
    return_detail   = TRUE,
    two_step        = two_step
  )
  
  tr <- res$trace_df
  if (is.null(tr) || !nrow(tr)) return(data.frame())
  
  # time trend values per allocation time
  et <- if (!is.null(res$expected_total)) res$expected_total else expected_total
  s_t <- pmin(tr$t, et) / et
  chrono_val <- chrono_beta_val * s_t
  
  make_rows_for_arm <- function(arm_name) {
    t0 <- res$window_open[[arm_name]]
    t1 <- res$window_close[[arm_name]]
    if (is.na(t0) || is.na(t1) || t1 < t0) return(NULL)
    
    # indices for allocations
    in_win   <- tr$t >= t0 & tr$t <= t1
    arm_idx  <- which(in_win & tr$assigned == arm_name)
    
    # control sets
    D_conc_idx   <- which(in_win & tr$assigned == "D")          # concurrent only
    D_conc_prior <- which(tr$assigned == "D" & tr$t <= t1)      # nonconcurrent (≤ close)
    
    # ---- Allocation bias (mean difference) ----
    mean_alloc_arm   <- if (length(arm_idx))      mean(tr$bias_val[arm_idx])      else 0
    mean_alloc_D_c   <- if (length(D_conc_idx))   mean(tr$bias_val[D_conc_idx])   else 0
    mean_alloc_D_cp  <- if (length(D_conc_prior)) mean(tr$bias_val[D_conc_prior]) else 0
    
    # ---- Chronological bias (mean difference) ----
    mean_chron_arm   <- if (length(arm_idx))      mean(chrono_val[arm_idx])      else 0
    mean_chron_D_c   <- if (length(D_conc_idx))   mean(chrono_val[D_conc_idx])   else 0
    mean_chron_D_cp  <- if (length(D_conc_prior)) mean(chrono_val[D_conc_prior]) else 0
    
    rbind(
      data.frame(
        arm         = arm_name,
        metric      = "Allocation bias",
        control_set = "Concurrent only",
        value       = mean_alloc_arm - mean_alloc_D_c,
        stringsAsFactors = FALSE
      ),
      data.frame(
        arm         = arm_name,
        metric      = "Allocation bias",
        control_set = "Nonconcurrent",
        value       = mean_alloc_arm - mean_alloc_D_cp,
        stringsAsFactors = FALSE
      ),
      data.frame(
        arm         = arm_name,
        metric      = "Chronological bias",
        control_set = "Concurrent only",
        value       = mean_chron_arm - mean_chron_D_c,
        stringsAsFactors = FALSE
      ),
      data.frame(
        arm         = arm_name,
        metric      = "Chronological bias",
        control_set = "Nonconcurrent",
        value       = mean_chron_arm - mean_chron_D_cp,
        stringsAsFactors = FALSE
      )
    )
  }
  
  out <- do.call(rbind, Filter(NROW, list(make_rows_for_arm("A"), make_rows_for_arm("B"))))
  if (!NROW(out)) return(data.frame())
  out
}

# ---------- collect all procedures ----------
collect_all <- function(two_step = FALSE) {
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  cores <- max(1L, cores)
  
  tasks <- list()
  for (p_i in seq_along(procedures)) {
    pr    <- procedures[[p_i]]
    seeds <- seed_base + seq_len(n_runs) + 1000L * p_i
    for (r in seq_len(n_runs)) {
      tasks[[length(tasks) + 1L]] <- list(
        rand_mode          = pr$rand_mode,
        block_factor       = pr$block_factor,
        rep_seed           = seeds[r],
        procedure          = pr$key,
        B_start_fixed      = B_start_fixed,
        max_group_size_exp = max_group_size_exp,
        alloc_bias_val     = alloc_bias_val,
        chrono_beta_val    = chrono_beta_val
      )
    }
  }
  
  run_task <- function(task) {
    df <- one_run_extract(
      rand_mode          = task$rand_mode,
      block_factor       = task$block_factor,
      rep_seed           = task$rep_seed,
      B_start_fixed      = task$B_start_fixed,
      max_group_size_exp = task$max_group_size_exp,
      alloc_bias_val     = task$alloc_bias_val,
      chrono_beta_val    = task$chrono_beta_val,
      two_step           = two_step
    )
    if (nrow(df)) df$procedure <- task$procedure
    df
  }
  
  parts <- if (.Platform$OS.type == "windows") {
    cl <- makeCluster(cores, type = "PSOCK")
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                 OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1")
      if (file.exists("simulation_functions.R")) source("simulation_functions.R")
      NULL
    })
    # export all knobs/functions needed by workers
    clusterExport(
      cl,
      varlist = c("one_run_extract","procedures","n_runs","seed_base",
                  "B_start_fixed","max_group_size_exp","alloc_bias_val",
                  "chrono_beta_val","expected_total","two_step"),
      envir = environment()
    )
    parLapply(cl, tasks, run_task)
  } else {
    if (cores > 1L) mclapply(tasks, run_task, mc.cores = cores) else lapply(tasks, run_task)
  }
  
  bias_rows <- do.call(rbind, parts)
  if (!nrow(bias_rows)) stop("No mean-bias rows produced — check simulation output.")
  
  bias_rows$procedure <- factor(
    bias_rows$procedure,
    levels = c("Block randomization (1*arms)",
               "Block randomization (8*arms)",
               "Big-stick design (a = 2)",
               "Complete randomization")
  )
  bias_rows$arm <- factor(bias_rows$arm,
                          levels = c("A","B"),
                          labels = c("Arm A","Arm B"))
  bias_rows$metric <- factor(bias_rows$metric,
                             levels = c("Allocation bias","Chronological bias"))
  bias_rows$control_set <- factor(
    bias_rows$control_set,
    levels = c("Concurrent only","Nonconcurrent")
  )
  bias_rows
}

# ---------- plot (CSV + PDF only) ----------
plot_and_save <- function(bias_rows, file_stub, title_suffix) {
  
  proc_cols <- c(
    "Block randomization (1*arms)" = "#d95f02",
    "Block randomization (8*arms)" = "#e7298a",
    "Big-stick design (a = 2)"     = "#7570b3",
    "Complete randomization"       = "#1b9e77"
  )
  
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  write.csv(bias_rows, csv_path, row.names = FALSE)
  
  # facet columns by metric × control set
  bias_rows$facet_col <- interaction(bias_rows$metric, bias_rows$control_set,
                                     drop = TRUE, sep = " — ")
  
  p <- ggplot(bias_rows, aes(x = procedure, y = value, fill = procedure)) +
    geom_violin(trim = FALSE, alpha = 0.45, color = NA) +
    geom_boxplot(width = 0.20, outlier.shape = 16, outlier.size = 0.8, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "black") +
    facet_grid(
      rows = vars(arm),
      cols = vars(facet_col),
      scales = "fixed",
      switch = "y"
    ) +
    scale_fill_manual(values = proc_cols, name = "Randomization") +
    labs(
      title = paste0(
        "Mean bias difference per run (Arm − control, by control set) — ",
        title_suffix
      ),
      subtitle = sprintf(
        "Mean(arm) − Mean(control); B starts at %d; η = %.2f; β = %.2f; max per arm = %d",
        B_start_fixed, alloc_bias_val, chrono_beta_val, max_group_size_exp
      ),
      x = "Randomization procedure",
      y = "Mean difference"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 15, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.spacing.y = unit(16, "pt")
    )
  
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  grDevices::cairo_pdf(file = pdf_path, width = 14, height = 10)
  print(p); dev.off()
  
  message("Wrote PDF + CSV: ", pdf_path, " (data: ", csv_path, ")")
}

# ---------- run ----------
set.seed(seed_base)

# 1-step randomization
bias_rows_1 <- collect_all(two_step = FALSE)
plot_and_save(
  bias_rows_1,
  file_stub     = "mean_bias_diff_AB_with_nonconc_1step",
  title_suffix  = "1-step randomization"
)

# 2-step randomization
bias_rows_2 <- collect_all(two_step = TRUE)
plot_and_save(
  bias_rows_2,
  file_stub     = "mean_bias_diff_AB_with_nonconc_2step",
  title_suffix  = "2-step randomization"
)
