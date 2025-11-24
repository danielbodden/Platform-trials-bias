# ============================================================
#  Validation for platform trials simulation
#  (1) Correct generation of randomization sequences 
#  (2) Type I error check at alpha = 0.05
#  (3) Correct allocation biasing policy
#  (4) Power check for block and complete randomization
# ============================================================

# --- 0) Setup ------------------------------------------------
set.seed(123)
source("simulation_functions.R")
library(ggplot2)

# Small helpers
code_to_arm <- function(k) c("A","B","C","D")[k]
arm_to_code <- function(ch) match(ch, c("A","B","C","D"))
`%||%` <- function(x, y) if (is.null(x)) y else x


# ============================================================
# (1) RANDOMIZATION SEQUENCE AUDIT 
#     - Prints core trial settings first
#     - Asserts no allocations to closed arms
#     - Prints planned, observed, and used windows + realized sample sizes
# ============================================================
# ============================================================
# (1) RANDOMIZATION SEQUENCE AUDIT 
#     - Prints core trial settings first
#     - Asserts no allocations to closed arms
#     - Prints planned, observed, and used windows + realized sample sizes
# ============================================================
audit_randomization <- function(res, arm_start, N = 100) {
  if (is.null(res$trace_df)) stop("trace_df is missing. Call platform_trials_simulation(return_detail=TRUE).")
  tr <- res$trace_df
  
  # ---- derive arm sets ----
  arms_all <- names(res$final_counts) %||% c("A","B","C","D")
  exp_arms <- names(res$window_open)  %||% setdiff(arms_all, "D")
  n_exp    <- length(exp_arms)
  ctrl_arm <- setdiff(arms_all, exp_arms); if (length(ctrl_arm) != 1L) ctrl_arm <- "D"
  
  # ---- core settings (echo) ----
  concurrent_only_used <- res$concurrent_only %||% TRUE
  rand_mode_used       <- res$rand_mode       %||% "?"
  alloc_bias_used      <- res$alloc_bias      %||% 0
  block_factor_used    <- res$block_factor    %||% NA
  max_group_size_used  <- res$max_group_size  %||% NA
  expected_total_used  <- res$expected_total  %||% NA
  beta_time_used       <- res$beta_time       %||% 0
  
  cat("\n=== TRIAL SETTINGS SUMMARY ===\n")
  cat(sprintf(" Randomization mode : %s\n", rand_mode_used))
  cat(sprintf(" Block factor       : %s\n", block_factor_used))
  cat(sprintf(" Max group size     : %s per arm\n", max_group_size_used))
  cat(sprintf(" Concurrent-only    : %s\n", concurrent_only_used))
  cat(sprintf(" Allocation bias η  : %s\n", alloc_bias_used))
  cat(sprintf(" Time trend β_time  : %s\n", beta_time_used))
  cat(sprintf(" Expected total N   : %s\n\n", expected_total_used))
  
  # ---- windows (observed) ----
  obs_open  <- (res$window_open  %||% rep(NA_integer_, n_exp))
  obs_close <- (res$window_close %||% rep(NA_integer_, n_exp))
  
  # ---- "used" windows (based on actual allocations) ----
  arm_open_used  <- obs_open
  arm_close_used <- obs_close
  for (k in seq_len(n_exp)) {
    arm_lab <- exp_arms[k]
    tt <- tr$t[tr$assigned == arm_lab]
    if (length(tt)) {
      if (is.na(arm_open_used[k]))  arm_open_used[k]  <- min(tt)
      if (is.na(arm_close_used[k])) arm_close_used[k] <- max(tt)
    }
  }
  
  # ---- sequence (first N) ----
  n <- min(N, nrow(tr))
  
  # open set in trace_df is a comma-separated string of arm labels, e.g. "A,B,D"
  parse_open <- function(s) {
    if (is.na(s) || nchar(s) == 0) return(integer(0))
    labs <- strsplit(s, ",", fixed = TRUE)[[1]]
    match(labs, arms_all)
  }
  open_list <- lapply(tr$open[seq_len(n)], parse_open)
  
  assign_seq <- match(tr$assigned[seq_len(n)], arms_all)
  times      <- tr$t[seq_len(n)]
  
  df <- data.frame(
    i        = seq_len(n),
    t        = times,
    arm      = arms_all[assign_seq],
    open_set = I(open_list),
    stringsAsFactors = FALSE
  )
  
  # ---- validate: in open set & within (used) windows ----
  is_valid_time <- logical(n)
  is_in_openset <- logical(n)
  for (i in seq_len(n)) {
    tt <- df$t[i]
    a  <- df$arm[i]
    os <- df$open_set[[i]] %||% integer(0)
    a_code <- match(a, arms_all)
    is_in_openset[i] <- a_code %in% os
    
    if (a %in% exp_arms) {
      k <- match(a, exp_arms)
      t0 <- arm_open_used[k]; t1 <- arm_close_used[k]
      is_valid_time[i] <- !is.na(t0) && tt >= t0 && (is.na(t1) || tt <= t1)
    } else {
      # control valid when any experimental window is active at time t
      any_open <- any(!is.na(arm_open_used) & !is.na(arm_close_used) & tt >= arm_open_used & tt <= arm_close_used)
      still_running <- any(!is.na(arm_open_used) & is.na(arm_close_used) & tt >= arm_open_used)
      is_valid_time[i] <- (any_open || still_running)
    }
  }
  df$valid_in_openset <- is_in_openset
  df$valid_time       <- is_valid_time
  df$valid            <- is_in_openset & is_valid_time
  
  # ---- window summary table ----
  comp_names <- paste0(exp_arms, "_vs_", ctrl_arm)
  rs <- res$realized_sizes
  cols <- if (!is.null(colnames(rs))) match(comp_names, colnames(rs)) else seq_len(n_exp)
  
  win_tbl <- data.frame(
    Arm               = exp_arms,
    Planned_Start     = as.integer(arm_start[exp_arms]),
    Observed_Open     = as.integer(obs_open),
    Close_Time        = as.integer(obs_close),
    Window_Start_Used = as.integer(arm_open_used),
    Window_End_Used   = as.integer(arm_close_used),
    Arm_n_realized    = as.integer(rs["arm",  cols, drop=TRUE]),
    Ctrl_n_realized   = as.integer(rs["ctrl", cols, drop=TRUE]),
    stringsAsFactors = FALSE
  )
  win_tbl$Window_Length_Used <- with(
    win_tbl,
    ifelse(is.na(Window_Start_Used) | is.na(Window_End_Used),
           NA_integer_,
           Window_End_Used - Window_Start_Used + 1L)
  )
  
  cat("--- ARM WINDOW SUMMARY (planned vs. observed vs. used) ---\n")
  print(win_tbl, row.names = FALSE)
  cat("\n--- RANDOMIZATION SEQUENCE (first", n, "allocations) ---\n")
  print(df, row.names = FALSE)
  
  if (!all(df$valid_in_openset)) warning("Found allocations to arms NOT in the open set at that time.")
  if (!all(df$valid_time))       warning("Found allocations outside the analysis windows (with fallback).")
  if (all(df$valid))             cat("All first", n, "allocations are valid. ✅\n")
  
  invisible(list(
    settings = list(
      rand_mode       = rand_mode_used,
      block_factor    = block_factor_used,
      max_group_size  = max_group_size_used,
      concurrent_only = concurrent_only_used,
      alloc_bias      = alloc_bias_used,
      beta_time       = beta_time_used,
      expected_total  = expected_total_used
    ),
    windows = win_tbl,
    seq = df
  ))
}

# ============================================================
# (1.1) PLOTTING FUNCTION FOR TRIAL TIMELINE
#       Shows the allocations to each arm, open and closing windows of arms.
# ============================================================
plot_trial_timeline <- function(res, arm_start, title = "Platform trial timeline",
                                show_allocations = TRUE) {
  if (is.null(res$trace_df)) stop("trace_df is missing. Call platform_trials_simulation(return_detail=TRUE).")
  tr <- res$trace_df
  
  arms_all <- names(res$final_counts) %||% c("A","B","C","D")
  exp_arms <- names(res$window_open)  %||% setdiff(arms_all, "D")
  n_exp    <- length(exp_arms)
  ctrl_arm <- setdiff(arms_all, exp_arms); if (length(ctrl_arm) != 1L) ctrl_arm <- "D"
  
  obs_open  <- (res$window_open  %||% rep(NA_integer_, n_exp))
  obs_close <- (res$window_close %||% rep(NA_integer_, n_exp))
  
  # "used" windows based on actual allocations
  arm_open_used  <- obs_open
  arm_close_used <- obs_close
  for (k in seq_len(n_exp)) {
    arm_lab <- exp_arms[k]
    tt <- tr$t[tr$assigned == arm_lab]
    if (length(tt)) {
      if (is.na(arm_open_used[k]))  arm_open_used[k]  <- min(tt)
      if (is.na(arm_close_used[k])) arm_close_used[k] <- max(tt)
    }
  }
  
  arm_names  <- exp_arms
  arm_factor <- factor(arm_names, levels = rev(arm_names))
  
  df_windows <- data.frame(
    arm   = arm_factor,
    start = as.numeric(arm_open_used),
    end   = as.numeric(arm_close_used)
  )
  df_windows$start_plot <- ifelse(is.na(df_windows$start), df_windows$end, df_windows$start)
  df_windows$end_plot   <- ifelse(is.na(df_windows$end),   df_windows$start, df_windows$end)
  
  df_marks <- do.call(
    rbind,
    c(
      lapply(seq_len(n_exp), function(i)
        data.frame(arm = factor(arm_names[i], levels = rev(arm_names)),
                   t = obs_open[i],  what = "Arm opens",  stringsAsFactors = FALSE)),
      lapply(seq_len(n_exp), function(i)
        data.frame(arm = factor(arm_names[i], levels = rev(arm_names)),
                   t = obs_close[i], what = "Arm closes", stringsAsFactors = FALSE))
    )
  )
  df_marks <- df_marks[!is.na(df_marks$t), , drop = FALSE]
  
  df_planned <- data.frame(
    arm = arm_factor,
    t   = as.numeric(arm_start[exp_arms])
  )
  df_planned$y <- as.numeric(df_planned$arm)
  
  p <- ggplot2::ggplot()
  
  # used windows
  df_seg <- subset(df_windows, !is.na(start_plot) & !is.na(end_plot) & end_plot >= start_plot)
  if (nrow(df_seg) > 0) {
    p <- p + ggplot2::geom_segment(
      data = df_seg,
      ggplot2::aes(x = start_plot, xend = end_plot, y = arm, yend = arm, color = arm),
      linewidth = 8, lineend = "round", alpha = 0.35
    )
  }
  
  # open/close markers
  if (nrow(df_marks) > 0) {
    p <- p + ggplot2::geom_point(
      data = df_marks,
      ggplot2::aes(x = t, y = arm, shape = what),
      size = 3, stroke = 1.1, color = "black", fill = "white"
    )
  }
  
  # planned starts
  p <- p + ggplot2::geom_segment(
    data = df_planned,
    ggplot2::aes(x = t, xend = t, y = y - 0.35, yend = y + 0.35),
    inherit.aes = FALSE,
    linewidth = 0.6, color = "black"
  ) +
    ggplot2::geom_text(
      data = df_planned,
      ggplot2::aes(x = t, y = arm, label = "P"),
      vjust = -1.2, size = 3
    )
  
  # allocations
  if (isTRUE(show_allocations)) {
    df_alloc <- data.frame(
      t   = tr$t,
      arm = tr$assigned
    )
    control_label <- paste0(ctrl_arm, " (control)")
    y_levels <- c(rev(exp_arms), control_label)
    
    p <- p +
      ggplot2::geom_point(
        data = subset(df_alloc, arm %in% exp_arms),
        ggplot2::aes(x = t, y = factor(arm, levels = rev(arm_names)), color = arm),
        size = 1.8, alpha = 0.85,
        position = ggplot2::position_jitter(height = 0.06, width = 0)
      ) +
      ggplot2::geom_point(
        data = subset(df_alloc, arm == ctrl_arm),
        ggplot2::aes(x = t, y = factor(control_label, levels = y_levels)),
        size = 1.6, alpha = 0.6, color = "gray40", shape = 16
      )
  }
  
  # colors as before
  base_cols <- c("A" = "#1b9e77", "B" = "#d95f02", "C" = "#7570b3")
  use_cols  <- base_cols[names(base_cols) %in% exp_arms]
  
  p +
    ggplot2::scale_color_manual(values = use_cols) +
    ggplot2::scale_shape_manual(values = c("Arm opens" = 21, "Arm closes" = 24)) +
    ggplot2::labs(title = title, x = "Patient randomized", y = NULL, color = "Arm", shape = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# ============================================================
# (2) TYPE I ERROR GRID (no bias), with both tests:
# - analysis_model ∈ {"ttest", "anova_period"}
# - ANOVA allowed only if concurrent_only == FALSE
# ============================================================
run_type1_grid <- function(
    nsim = 2000,
    max_group_size = 50,
    block_factor = 1,
    alpha = 0.05,
    test_side = "two.sided"
) {
  base_grid <- expand.grid(
    concurrent_only = c(TRUE, FALSE),
    rand_mode = c("complete", "block"),
    startB = c(0L, 50L, 100L, 150L),
    stringsAsFactors = FALSE
  )
  
  # For each base row, decide which analysis models to run
  # - If concurrent_only == TRUE  -> only "ttest"
  # - If concurrent_only == FALSE -> both "ttest" and "anova_period"
  expanded <- do.call(
    rbind,
    lapply(seq_len(nrow(base_grid)), function(i) {
      cfg <- base_grid[i, ]
      models <- if (isTRUE(cfg$concurrent_only)) "ttest" else c("ttest", "anova_period")
      data.frame(cfg[rep(1, length(models)), , drop = FALSE],
                 analysis_model = models,
                 row.names = NULL,
                 stringsAsFactors = FALSE)
    })
  )
  
  out <- data.frame(
    concurrent_only = expanded$concurrent_only,
    rand_mode       = expanded$rand_mode,
    startB          = expanded$startB,
    analysis_model  = expanded$analysis_model,
    A = NA_real_, B = NA_real_, C = NA_real_, any = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(expanded))) {
    cfg <- expanded[i, ]
    arm_start <- c(A = 1L, B = cfg$startB, C = 1L)
    
    rej_mat <- matrix(FALSE, nrow = nsim, ncol = 3)
    
    for (s in seq_len(nsim)) {
      res <- platform_trials_simulation(
        max_group_size  = max_group_size,
        mu              = c(A = 0, B = 0, C = 0, D = 0),  # Type I (null true)
        alpha           = alpha,
        arm_start       = arm_start,
        concurrent_only = cfg$concurrent_only,
        expected_total  = 200,
        beta_time       = 0,                              # no time trend
        rand_mode       = cfg$rand_mode,                  # CR or PBR
        block_factor    = block_factor,
        alloc_bias      = 0,                              # NO allocation bias
        test_side       = test_side,
        analysis_model  = cfg$analysis_model
      )
      rej_mat[s, ] <- unname(res$reject)
    }
    
    rates <- colMeans(rej_mat)
    out$A[i]   <- rates[1]
    out$B[i]   <- rates[2]
    out$C[i]   <- rates[3]
    out$any[i] <- mean(rowSums(rej_mat) > 0)
    
    cat(sprintf(
      "Done %d/%d: conc=%s, mode=%s, startB=%d, model=%s  |  A=%.3f B=%.3f C=%.3f FWER=%.3f\n",
      i, nrow(expanded), cfg$concurrent_only, cfg$rand_mode, cfg$startB, cfg$analysis_model,
      out$A[i], out$B[i], out$C[i], out$any[i]
    ))
  }
  
  out
}


# ============================================================
# --- Run the validations ------------------------------------
# ============================================================

# 1) Randomization + analysis audit 
ex1 <- platform_trials_simulation(return_detail = TRUE,
  max_group_size  = 30,
  mu              = c(A=0, B=0, C=0, D=0),
  alpha           = 0.05,
  arm_start       = c(A=1, B=10, C=1),
  concurrent_only = TRUE,
  expected_total  = 200,
  beta_time       = 0,
  rand_mode       = "block",  
  block_factor    = 2,
  alloc_bias       = 0
)

ex2 <- platform_trials_simulation(return_detail = TRUE,
  max_group_size  = 30,
  mu              = c(A=0, B=0, C=0, D=0),
  alpha           = 0.05,
  arm_start       = c(A=1, B=10, C=1),
  concurrent_only = TRUE,
  expected_total  = 200,
  beta_time       = 0,
  rand_mode       = "complete",  
  block_factor    = 1,
  alloc_bias       = 0
)

ex3 <- platform_trials_simulation(return_detail = TRUE,
  max_group_size  = 30,
  mu              = c(A=0, B=0, C=0, D=0),
  alpha           = 0.05,
  arm_start       = c(A=1, B=10, C=1),
  concurrent_only = FALSE,
  expected_total  = 200,
  beta_time       = 0,
  rand_mode       = "complete",   
  block_factor    = 1,
  alloc_bias       = 0
)


ex4 <- platform_trials_simulation(return_detail = TRUE,
  max_group_size  = 30,
  mu              = c(A=0, B=0, C=0, D=0),
  alpha           = 0.05,
  arm_start       = c(A=1, B=10, C=1),
  concurrent_only = FALSE,
  expected_total  = 200,
  beta_time       = 0,
  rand_mode       = "block",  
  block_factor    = 1,
  alloc_bias       = 0
)


ex5 <- platform_trials_simulation(return_detail = TRUE,
  max_group_size  = 30,
  mu              = c(A=0, B=0, D=0),
  alpha           = 0.05,
  arm_start       = c(A=1, B=10),
  concurrent_only = FALSE,
  expected_total  = 200,
  beta_time       = 0,
  rand_mode       = "block",  
  block_factor    = 1,
  alloc_bias       = 0,
  exp_arms        = c("A","B") 
)

# (1) Sequence audit
audit_randomization(ex1, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex2, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex3, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex4, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex5, arm_start = c(A=1, B=10), N = 150)

# Plotting timeline of trials
p1 <- plot_trial_timeline(
  res = ex1,
  arm_start = c(A = 1, B = 10, C = 1),
  title = "Trial timeline – PBR, concurrent-only, B opens at t=10",
  show_allocations = TRUE
)
print(p1)

p2 <- plot_trial_timeline(
  res = ex2,
  arm_start = c(A = 1, B = 10, C = 1),
  title = "Trial timeline – CR, concurrent-only, B opens at t=10",
  show_allocations = TRUE
)
print(p2)

p3 <- plot_trial_timeline(
  res = ex3,
  arm_start = c(A = 1, B = 10, C = 1),
  title = "Trial timeline – CR, nonconcurrent, B opens at t=10",
  show_allocations = TRUE
)
print(p3)

p4 <- plot_trial_timeline(
  res = ex4,
  arm_start = c(A = 1, B = 10, C = 1),
  title = "Trial timeline – PBR, nonconcurrent, B opens at t=10",
  show_allocations = TRUE
)
print(p4)

p5 <- plot_trial_timeline(
  res = ex5,
  arm_start = c(A = 1, B = 10),
  title = "Trial timeline – PBR, nonconcurrent, B opens at t=10",
  show_allocations = TRUE
)
print(p5)

# 2) Type I error grid (no bias, one-sided)
type1_results <- run_type1_grid(
  nsim = 2000,         
  max_group_size = 30, 
  block_factor = 1,
  test_side="one.sided",
  alpha=0.025
)
cat("\n=== TYPE I ERROR SUMMARY (target ~ 0.025) ===\n")
print(type1_results, row.names = FALSE)


# 3) Responder type tracking for allocation biasing policy validation
# Note: bias_policy="favor_B" has not been validated, as it is not being used.
res <- platform_trials_simulation(return_detail = TRUE, max_group_size = 24, alloc_bias = 0.2,
                                  beta_time = 0.1, expected_total = 150,
                                  rand_mode = "block", block_factor = 1,
                                  analysis_model = "ttest", exp_arms = c("A","B"), arm_start = c(A=0,B=16), bias_policy="favor_all_exp")
print(res)

res <- platform_trials_simulation(return_detail = TRUE, max_group_size = 24, alloc_bias = 0.2, two_step=TRUE,
                                  beta_time = 0.1, expected_total = 150,
                                  rand_mode = "block", block_factor = 1,
                                  analysis_model = "ttest", exp_arms = c("A","B"), arm_start = c(A=1,B=16), bias_policy="favor_all_exp")
print(res)


# ============================
# (4) Power validation
# For a parallel group design the power is 80% for an effect size of 0.56
# Therefore the power should be ~80%
# ============================
calc_rejection_summary(n_sim = 5000,
                       max_group_size = 50,
                       mu = c(A=0.56,B=0.56,C=0.56,D=0),
                       alpha = 0.05,
                       arm_start = c(A=1,B=100,C=200),
                       concurrent_only = TRUE,
                       expected_total = 300,
                       beta_time = 0,
                       rand_mode = "block",
                       block_factor = 1,
                       alloc_bias = 0)

calc_rejection_summary(n_sim = 5000,
                       max_group_size = 50,
                       mu = c(A=0.56,B=0.56,C=0.56,D=0),
                       alpha = 0.05,
                       arm_start = c(A=1,B=100,C=200),
                       concurrent_only = TRUE,
                       expected_total = 300,
                       beta_time = 0,
                       rand_mode = "complete",
                       block_factor = 1,
                       alloc_bias = 0)

calc_rejection_summary(n_sim = 5000,
                       max_group_size = 24,
                       mu = c(A=0.826,B=0.826,C=0.826,D=0),
                       arm_start = c(A=1,B=100,C=200),
                       concurrent_only = TRUE,
                       expected_total = 300,
                       beta_time = 0,
                       rand_mode = "complete",
                       block_factor = 1,
                       alloc_bias = 0,
                       alpha=0.025,
                       test_side="one.sided")

calc_rejection_summary(n_sim = 5000,
                       max_group_size = 24,
                       mu = c(A=0.826,B=0.826,C=0.826,D=0),
                       arm_start = c(A=1,B=100,C=200),
                       concurrent_only = TRUE,
                       expected_total = 300,
                       beta_time = 0,
                       rand_mode = "block",
                       block_factor = 1,
                       alloc_bias = 0,
                       alpha=0.025,
                       test_side="one.sided")


###############

# ============================================================
# BIGSTICK AUDIT (one-step)
# Mirrors next_assignment() logic for rand_mode = "bigstick"
# Uses per-period local_counts (*_local) in trace_df.
# ============================================================
audit_bigstick_onestep <- function(res, bigstick_a, N = 200) {
  if (is.null(res$trace_df)) stop("trace_df is missing. Call platform_trials_simulation(return_detail=TRUE).")
  tr <- res$trace_df
  
  rand_mode_used <- res$rand_mode %||% "?"
  if (!identical(rand_mode_used, "bigstick")) {
    warning(sprintf("audit_bigstick_onestep(): rand_mode is '%s', not 'bigstick'.", rand_mode_used))
  }
  
  # arm names from *_local columns (these mirror local_counts in platform_trials_simulation)
  local_cols <- grep("_local$", names(tr), value = TRUE)
  if (!length(local_cols)) stop("No *_local columns found in trace_df; cannot audit bigstick.")
  arms_all <- sub("_local$", "", local_cols)
  
  # helper: parse open set labels -> arm labels
  parse_open <- function(s) {
    if (is.na(s) || nchar(s) == 0) return(character(0))
    strsplit(s, ",", fixed = TRUE)[[1]]
  }
  
  n <- min(N, nrow(tr))
  imbalance  <- rep(NA_integer_, n)
  max_open   <- rep(NA_integer_, n)
  min_open   <- rep(NA_integer_, n)
  rule_type  <- rep(NA_character_, n)  # "free" or "restricted"
  ok         <- rep(NA, n)
  allowed    <- vector("list", n)
  
  for (i in seq_len(n)) {
    open_arms <- parse_open(tr$open[i])
    open_arms <- open_arms[open_arms %in% arms_all]
    if (!length(open_arms)) next
    
    # local_counts BEFORE assignment i for the open set
    lc_row <- setNames(
      as.numeric(tr[i, paste0(open_arms, "_local"), drop = TRUE]),
      open_arms
    )
    
    max_c <- max(lc_row)
    min_c <- min(lc_row)
    imb   <- max_c - min_c
    
    imbalance[i] <- imb
    max_open[i]  <- max_c
    min_open[i]  <- min_c
    
    # bigstick rule as in next_assignment()
    if (imb > bigstick_a) {
      # allocate randomly among arms with smallest per-period counts
      rule_type[i] <- "restricted"
      allowed_arms <- names(lc_row)[lc_row == min_c]
    } else {
      # equal randomization among all open arms
      rule_type[i] <- "free"
      allowed_arms <- names(lc_row)
    }
    allowed[[i]] <- allowed_arms
    
    assigned_arm <- tr$assigned[i]
    ok[i] <- assigned_arm %in% allowed_arms
  }
  
  audit_df <- data.frame(
    i         = seq_len(n),
    t         = tr$t[seq_len(n)],
    assigned  = tr$assigned[seq_len(n)],
    imbalance = imbalance,
    max_open  = max_open,
    min_open  = min_open,
    rule_type = rule_type,
    ok        = ok,
    stringsAsFactors = FALSE
  )
  audit_df$allowed_arms <- vapply(
    allowed,
    function(x) if (is.null(x)) "" else paste(x, collapse = ","),
    character(1)
  )
  
  cat("\n=== BIGSTICK ONE-STEP AUDIT (first", n, "allocations) ===\n")
  cat(" rand_mode  =", rand_mode_used, "\n")
  cat(" bigstick_a =", bigstick_a, "\n")
  
  valid_rows <- !is.na(audit_df$ok)
  if (any(valid_rows)) {
    viol_rate <- mean(!audit_df$ok[valid_rows])
    cat(sprintf(" Checked %d allocations with a defined rule.\n", sum(valid_rows)))
    cat(sprintf(" Violations of bigstick rule: %d (%.2f%%)\n",
                sum(!audit_df$ok[valid_rows]),
                100 * viol_rate))
  } else {
    cat(" No allocations with a defined rule.\n")
  }
  
  if (any(valid_rows) && any(!audit_df$ok[valid_rows])) {
    cat("\n--- EXAMPLES OF VIOLATIONS (up to 20 rows) ---\n")
    print(utils::head(audit_df[valid_rows & !audit_df$ok, ], 20), row.names = FALSE)
  } else if (any(valid_rows)) {
    cat(" No bigstick rule violations detected. ✅\n")
  }
  
  invisible(audit_df)
}



ex_bs <- platform_trials_simulation(
  return_detail  = TRUE,
  max_group_size = 30,
  mu             = c(A=0,B=0,D=0),
  alpha          = 0.05,
  arm_start      = c(A=1,B=10,D=1),
  concurrent_only= TRUE,
  expected_total = 200,
  beta_time      = 0,
  rand_mode      = "bigstick",
  block_factor   = 1,
  alloc_bias     = 0,
  bigstick_a     = 2L,
  exp_arms       = c("A","B")         # <- important
)

audit_randomization(ex_bs, arm_start = c(A=1,B=10,D=1), N = 20)
audit_bigstick_onestep(ex_bs, bigstick_a = 2L, N = 150)


set.seed(1)
ex_bs <- platform_trials_simulation(
  return_detail  = TRUE,
  max_group_size = 30,
  mu             = c(A=0,B=0,D=0),
  alpha          = 0.05,
  arm_start      = c(A=1,B=10,D=1),
  concurrent_only= TRUE,
  expected_total = 200,
  beta_time      = 0,
  rand_mode      = "bigstick",   # or "block", "complete"
  block_factor   = 1,
  alloc_bias     = 0,
  bigstick_a     = 2L,
  exp_arms       = c("A","B"),
  two_step       = TRUE          # to force the two-step branch
)

