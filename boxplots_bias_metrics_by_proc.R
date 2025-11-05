#!/usr/bin/env Rscript
# ============================================================
# ONE FIGURE: Bias vs Randomization Procedure
# - X: Arm A, Arm B, Control (D)
# - Columns: Mean allocation bias, Mean chronological bias
# - Rows: Randomization procedures (stacked)
# - Thin dashed line at y=0; colorful violins + boxplots
#
# Env knobs:
#   N_RUNS    (default 300)
#   SEED_BASE (default 20251031)
#   OUT_DIR   (default "results")
#   SLURM_CPUS_PER_TASK respected on Linux
#
# Requires: simulation_functions.R in the same directory
# Outputs: results/bias_by_proc_SINGLE.(csv|png|pdf)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(parallel)
  library(grid)   # unit()
})

# ---------- knobs ----------
n_runs    <- as.integer(Sys.getenv("N_RUNS",    "40000"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", "20251031"))
out_dir   <- Sys.getenv("OUT_DIR", unset = "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Ensure UTF-8 plotting on headless nodes
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

# Fixed trial settings (for subtitle; control varies naturally)
max_group_size_exp <- 24L   # max patients PER EXPERIMENTAL ARM (A,B)
B_start_fixed      <- 16L
alloc_bias_val     <- 0.08
chrono_beta_val    <- 0.16

# ---------- simulator ----------
source("simulation_functions.R")

# ---------- procedures ----------
procedures <- list(
  list(key = "Complete randomization",           rand_mode="complete", block_factor=1L),
  list(key = "Block randomization (1*arms)",     rand_mode="block",    block_factor=1L),
  list(key = "Block randomization (2*arms)",     rand_mode="block",    block_factor=2L),
  list(key = "Block randomization (8*arms)",     rand_mode="block",    block_factor=8L)
)

# ---------- one run -> bias rows ----------
one_run_extract <- function(rand_mode, block_factor, rep_seed,
                            B_start_fixed, max_group_size_exp,
                            alloc_bias_val, chrono_beta_val) {
  set.seed(rep_seed)
  res <- platform_trials_simulation(
    exp_arms        = c("A","B"),
    arm_start       = c(A=0L, B=B_start_fixed),
    max_group_size  = max_group_size_exp,
    expected_total  = 96,
    alpha           = 0.025,
    test_side       = "one.sided",
    alternative     = "greater",
    rand_mode       = rand_mode,
    block_factor    = block_factor,
    alloc_bias      = alloc_bias_val,
    beta_time       = chrono_beta_val,
    concurrent_only = TRUE,            # fixed; not shown in figure
    analysis_model  = "ttest",
    bias_policy     = "favor_all_exp"
  )
  
  bm <- res$bias_metrics
  if (is.null(bm) || !length(bm)) {
    data.frame(arm=character(0), metric=character(0), value=numeric(0))
  } else {
    data.frame(
      arm    = rep(rownames(bm), times = ncol(bm)),
      metric = rep(colnames(bm), each  = nrow(bm)),
      value  = as.vector(bm),
      stringsAsFactors = FALSE
    )
  }
}

# ---------- collect all procedures (parallel-safe) ----------
collect_all <- function() {
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  cores <- max(1L, cores)
  
  tasks <- list()
  for (p_i in seq_along(procedures)) {
    pr <- procedures[[p_i]]
    seeds <- seed_base + seq_len(n_runs) + 1000L * p_i
    for (r in seq_len(n_runs)) {
      tasks[[length(tasks)+1L]] <- list(
        rand_mode        = pr$rand_mode,
        block_factor     = pr$block_factor,
        rep_seed         = seeds[r],
        procedure        = pr$key,
        # pass constants explicitly (avoid missing-object issues)
        B_start_fixed    = B_start_fixed,
        max_group_size_exp = max_group_size_exp,
        alloc_bias_val   = alloc_bias_val,
        chrono_beta_val  = chrono_beta_val
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
      chrono_beta_val    = task$chrono_beta_val
    )
    if (nrow(df)) df$procedure <- task$procedure
    df
  }
  
  parts <- if (.Platform$OS.type == "windows") {
    cl <- makeCluster(cores, type = "PSOCK")
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    # Prevent BLAS oversubscription
    clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1")
      source("simulation_functions.R"); NULL
    })
    # Export the runner and any referenced objects
    clusterExport(
      cl,
      varlist = c("one_run_extract"),
      envir = environment()
    )
    parLapply(cl, tasks, run_task)
  } else {
    # On Unix, forked workers inherit the environment; still safe to use mclapply
    if (cores > 1L) mclapply(tasks, run_task, mc.cores = cores) else lapply(tasks, run_task)
  }
  
  bias_rows <- do.call(rbind, parts)
  if (!nrow(bias_rows)) stop("No bias rows produced — check simulation_functions.R metrics.")
  
  # keep only the two metrics we want
  keep_map <- c(
    "mean_allocation_bias"    = "Mean allocation bias",
    "mean_chronological_bias" = "Mean chronological bias"
  )
  bias_rows <- subset(bias_rows, metric %in% names(keep_map))
  bias_rows$metric <- keep_map[bias_rows$metric]
  
  # Keep only A, B, D; label nicely
  bias_rows <- subset(bias_rows, arm %in% c("A","B","D"))
  bias_rows$arm <- factor(bias_rows$arm,
                          levels = c("A","B","D"),
                          labels = c("Arm A","Arm B","Control (D)"))
  
  # order procedures
  bias_rows$procedure <- factor(
    bias_rows$procedure,
    levels = c("Complete randomization",
               "Block randomization (1*arms)",
               "Block randomization (2*arms)",
               "Block randomization (8*arms)")
  )
  
  bias_rows
}

# ---------- plotting ----------
plot_and_save <- function(bias_rows, file_stub = "bias_by_proc_SINGLE") {
  
  # color palette per procedure (colorful)
  proc_cols <- c(
    "Complete randomization"           = "#1b9e77",
    "Block randomization (1*arms)"     = "#d95f02",
    "Block randomization (2*arms)"     = "#7570b3",
    "Block randomization (8*arms)"     = "#e7298a"
  )
  
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  write.csv(bias_rows, csv_path, row.names = FALSE)
  
  p <- ggplot(bias_rows, aes(x = arm, y = value, fill = procedure)) +
    geom_violin(trim = FALSE, alpha = 0.45, color = NA) +
    geom_boxplot(width = 0.20, outlier.shape = 16, outlier.size = 0.8, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "black") +
    facet_grid(
      rows = vars(procedure),      # stacked procedures
      cols = vars(metric),
      scales = "fixed",
      switch = "y"
    ) +
    scale_fill_manual(values = proc_cols, name = "Randomization") +
    labs(
      title = "Bias vs Randomization Procedure",
      subtitle = sprintf(
        "B starts at %d; η = %.2f, β = %.2f; max patients per analysis = 2*%d",
        B_start_fixed, alloc_bias_val, chrono_beta_val, max_group_size_exp
      ),
      x = "Arm",
      y = "Bias value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 15, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.spacing.y = unit(16, "pt")  # extra whitespace between stacked rows
    )
  
  png_path <- file.path(out_dir, paste0(file_stub, ".png"))
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  
  grDevices::png(filename = png_path, type = "cairo", width = 14, height = 10, units = "in", res = 150)
  print(p); dev.off()
  grDevices::cairo_pdf(file = pdf_path, width = 14, height = 10)
  print(p); dev.off()
  
  message("Wrote: ", png_path, " and ", pdf_path, " (data: ", csv_path, ")")
}

# ---------- run ----------
set.seed(seed_base)
bias_rows <- collect_all()
plot_and_save(bias_rows)
