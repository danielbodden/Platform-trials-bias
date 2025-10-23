# ============================================================
#  Validation for platform trials simulation
#  (1) Correct generation of randomization sequences 
#  (2) Type I error check at alpha = 0.05
#  (3) Correct allocation biasing policy
# ============================================================

# --- 0) Setup ------------------------------------------------
set.seed(123)
source("simulation_functions.R")
source("main.R")


# Small helpers
code_to_arm <- function(k) c("A","B","C","D")[k]
arm_to_code <- function(ch) match(ch, c("A","B","C","D"))
`%||%` <- function(x, y) if (is.null(x)) y else x


# ============================================================
# 0.1) Logged copy of the simulator (no logic changes)
#      Adds logging arrays so we can validate sequences & windows.
#      This function has not been updated with newest changes.
# ============================================================
platform_trials_simulation_logged <- function(
    max_group_size  = 50,
    mu              = c(A=0, B=0, C=0, D=0),
    alpha           = 0.05,
    arm_start       = c(A=1, B=100, C=1),
    concurrent_only = TRUE,
    expected_total  = 200,
    beta_time       = 0,
    rand_mode       = c("block", "complete"),
    block_factor    = 1,
    alloc_bias       = 0
) {
  rand_mode <- match.arg(rand_mode)
  
  # Fixed order A,B,C,D
  mu_vec         <- mu[c("A","B","C","D")]
  arm_starts_vec <- c(arm_start["A"], arm_start["B"], arm_start["C"], 0L)
  
  # Counters
  arm_counts   <- integer(4L)
  local_counts <- integer(4L)
  target_exp   <- rep(max_group_size, 3L)
  
  # Storage
  cap      <- max(1024L, as.integer(2.5 * expected_total))
  assign_i <- integer(cap)
  times_i  <- integer(cap)
  bias_i   <- numeric(cap)
  
  # LOG buffers
  open_set_log       <- vector("list", cap)                 # open codes per tick
  finished_mask_log  <- matrix(FALSE, nrow = cap, ncol = 3) # finished per arm per tick
  
  idx <- 0L; t <- 0L
  
  # Randomization
  current_block   <- integer(0)
  block_pos       <- 0L
  block_signature <- -1L
  .sig_bitmask <- function(open_codes) sum(bitwShiftL(1L, open_codes - 1L))
  build_block   <- function(open_codes) sample(rep(open_codes, each = block_factor))
  next_assignment <- function(open_codes) {
    if (rand_mode == "complete") return(sample(open_codes, 1L))
    sig <- .sig_bitmask(open_codes)
    need_new <- (length(current_block) == 0L) ||
      (block_pos >= length(current_block)) ||
      (sig != block_signature)
    if (need_new) {
      current_block   <<- build_block(open_codes)
      block_pos       <<- 0L
      block_signature <<- sig
      local_counts[open_codes] <<- 0L
    }
    block_pos <<- block_pos + 1L
    current_block[[block_pos]]
  }
  
  t_open            <- rep(NA_integer_, 3L)  # A,B,C
  close_time        <- rep(NA_integer_, 3L)
  ctrl_n            <- integer(3L)
  global_ctrl_total <- 0L
  
  # Finished flags
  finished_mask <- rep(FALSE, 3L)
  
  repeat {
    if (all(finished_mask)) break
    t <- t + 1L
    
    # Which arms are open by time t?
    open_codes <- which(arm_starts_vec <= t)
    
    # Drop finished experimental arms
    if (any(finished_mask)) {
      fin_idx   <- which(finished_mask)
      open_codes <- open_codes[!(open_codes %in% fin_idx)]
    }
    
    # Ensure control present whenever any experimental is open
    if (any(open_codes %in% 1:3) && !(4L %in% open_codes)) open_codes <- c(open_codes, 4L)
    
    newly_eligible <- intersect(which(is.na(t_open)), intersect(1:3, open_codes))
    if (length(newly_eligible)) {
      t_open[newly_eligible] <- t
      if (!concurrent_only) {
        # non-concurrent: past controls count from opening
        ctrl_n[newly_eligible] <- global_ctrl_total
      }
    }
    
    # If no experimental arm is open yet, just log & advance
    if (!any(open_codes %in% 1:3)) {
      if (idx + 1L <= cap) {
        open_set_log[[idx + 1L]]   <- open_codes
        finished_mask_log[idx + 1L, ] <- finished_mask
      }
      next
    }
    
    # Grow storage if needed
    if (idx >= cap) {
      grow     <- max(512L, as.integer(0.5 * cap))
      assign_i <- c(assign_i, integer(grow))
      times_i  <- c(times_i,  integer(grow))
      bias_i   <- c(bias_i,   numeric(grow))
      open_set_log       <- c(open_set_log,    vector("list", grow))
      finished_mask_log  <- rbind(finished_mask_log, matrix(FALSE, nrow = grow, ncol = 3))
      cap      <- cap + grow
    }
    
    # Allocation bias only
    counts_open <- local_counts[open_codes]
    m           <- min(counts_open)
    eq_min      <- open_codes[counts_open == m]
    if (length(eq_min) == 1L) {
      expected_code <- eq_min
      bias_next     <- if (expected_code == 4L) -alloc_bias else +alloc_bias
    } else {
      expected_code <- sample(eq_min, 1L)
      bias_next     <- 0
    }
    
    # Actual assignment
    arm_code <- next_assignment(open_codes)
    
    # Record assignment & logs
    idx                      <- idx + 1L
    assign_i[idx]            <- arm_code
    times_i [idx]            <- t
    bias_i  [idx]            <- bias_next
    arm_counts[arm_code]     <- arm_counts[arm_code]   + 1L
    local_counts[arm_code]   <- local_counts[arm_code] + 1L
    open_set_log[[idx]]      <- open_codes
    finished_mask_log[idx, ] <- finished_mask
    
    if (arm_code %in% 1:3) {
      k <- arm_code
      if (!finished_mask[k] && arm_counts[k] >= max_group_size && ctrl_n[k] >= max_group_size) {
        close_time[k]    <- t
        finished_mask[k] <- TRUE
      }
    } else {
      # control
      global_ctrl_total <- global_ctrl_total + 1L
      
      # increment ctrl tallies for any open, unfinished arm
      open_mask <- !is.na(t_open) & !finished_mask
      if (any(open_mask)) ctrl_n[open_mask] <- ctrl_n[open_mask] + 1L
      
      # check for newly finished arms
      if (any(!finished_mask)) {
        newly_done <- which(!finished_mask &
                              (arm_counts[1:3] >= max_group_size) &
                              (ctrl_n >= max_group_size))
        if (length(newly_done)) {
          close_time[newly_done]    <- t
          finished_mask[newly_done] <- TRUE
        }
      }
    }
  }
  
  # Trim
  if (idx > 0L) {
    assign_i <- assign_i[seq_len(idx)]
    times_i  <- times_i [seq_len(idx)]
    bias_i   <- bias_i  [seq_len(idx)]
    finished_mask_log <- finished_mask_log[seq_len(idx), , drop = FALSE]
    open_set_log      <- open_set_log[seq_len(idx)]
  }
  
  # Outcomes
  s_t     <- if (idx > 0L) pmin(1, times_i / expected_total) else numeric(0)
  mu_pat  <- if (idx > 0L) (mu_vec[assign_i] + beta_time * s_t + bias_i) else numeric(0)
  outcomes <- if (idx > 0L) rnorm(idx, mean = mu_pat, sd = 1) else numeric(0)
  
  reject   <- logical(3)
  realized <- matrix(0L, 2, 3, dimnames = list(c("arm","ctrl"), c("A_vs_D","B_vs_D","C_vs_D")))
  for (k in 1:3) {
      arm_open <- t_open[k]
      t_end    <- close_time[k]
    x_idx <- which(assign_i == k & times_i >= arm_open & times_i <= t_end)
    nx    <- length(x_idx)
    if (concurrent_only) {
      y_idx <- which(assign_i == 4L & times_i >= arm_open & times_i <= t_end)
    } else {
      y_idx <- which(assign_i == 4L & times_i <= t_end)
    }
    ny <- length(y_idx)
    realized["arm",  k] <- nx
    realized["ctrl", k] <- ny
    if (nx < 2 || ny < 2) { reject[k] <- FALSE; next }
    res <- two_sample_t_pooled(outcomes[x_idx], outcomes[y_idx])
    reject[k] <- abs(res$t) > qt(1 - alpha/2, df = res$df)
  }
  names(reject) <- c("A_vs_D","B_vs_D","C_vs_D")
  
  return(list(
    reject            = reject,
    num_rej           = sum(reject),
    realized_sizes    = realized,
    total_randomized  = idx,
    final_counts      = setNames(arm_counts, c("A","B","C","D")),
    window_open       = setNames(t_open,     c("A","B","C")),
    window_close      = setNames(close_time, c("A","B","C")),
    assign_seq        = assign_i,
    assign_times      = times_i,
    open_set_log      = open_set_log,
    finished_mask_log = finished_mask_log,
    # --- NEW: echo back settings so the audit can print them ---
    concurrent_only   = concurrent_only,
    rand_mode         = rand_mode,
    block_factor      = block_factor,
    max_group_size    = max_group_size,
    expected_total    = expected_total,
    beta_time         = beta_time,
    alloc_bias         = alloc_bias
  ))
}
# ============================================================
# (1) RANDOMIZATION SEQUENCE AUDIT 
#     - Prints core trial settings first
#     - Asserts no allocations to closed arms
#     - Prints planned, observed, and used windows + realized sample sizes
# ============================================================
audit_randomization <- function(res, arm_start, N = 100) {
  n <- min(N, length(res$assign_seq))
  
  # Extract raw logs
  assign <- res$assign_seq
  times  <- res$assign_times
  concurrent_only_used <- res$concurrent_only %||% attr(res, "concurrent_only") %||% TRUE
  rand_mode_used       <- res$rand_mode       %||% attr(res, "rand_mode")       %||% "?"
  alloc_bias_used       <- res$alloc_bias       %||% attr(res, "alloc_bias")       %||% 0
  block_factor_used    <- res$block_factor    %||% attr(res, "block_factor")    %||% NA
  max_group_size_used  <- res$max_group_size  %||% attr(res, "max_group_size")  %||% NA
  expected_total_used  <- res$expected_total  %||% attr(res, "expected_total")  %||% NA
  beta_time_used       <- res$beta_time       %||% attr(res, "beta_time")       %||% 0
  
  # ============================================================
  # 1. Print trial settings summary
  # ============================================================
  cat("\n=== TRIAL SETTINGS SUMMARY ===\n")
  cat(sprintf(" Randomization mode : %s\n", rand_mode_used))
  cat(sprintf(" Block factor       : %s\n", block_factor_used))
  cat(sprintf(" Max group size     : %s per arm\n", max_group_size_used))
  cat(sprintf(" Concurrent-only    : %s\n", concurrent_only_used))
  cat(sprintf(" Allocation bias η  : %s\n", alloc_bias_used))
  cat(sprintf(" Time trend β_time  : %s\n", beta_time_used))
  cat(sprintf(" Expected total N   : %s\n", expected_total_used))
  cat("\n")
  
  # ============================================================
  # 2. Basic sequence info
  # ============================================================
  df <- data.frame(
    i        = seq_len(n),
    t        = times[seq_len(n)],
    arm      = c("A","B","C","D")[assign[seq_len(n)]],
    open_set = I(res$open_set_log[seq_len(n)]),
    stringsAsFactors = FALSE
  )
  
  obs_open  <- res$window_open
  obs_close <- res$window_close
  
  # ============================================================
  # 3. Build "used" window times
  # ============================================================
  arm_labels <- c("A","B","C")
  arm_open  <- obs_open
  arm_close <- obs_close
  for (k in 1:3) {
    x_times <- times[assign == k]
    if (length(x_times)) {
      if (is.na(arm_open[k]))  arm_open[k]  <- min(x_times)
      if (is.na(arm_close[k])) arm_close[k] <- max(x_times)
    }
  }
  
  # ============================================================
  # 4. Validate allocations
  # ============================================================
  is_valid_time <- logical(n)
  is_in_openset <- logical(n)
  for (i in seq_len(n)) {
    tt <- df$t[i]
    a  <- df$arm[i]
    os <- df$open_set[[i]] %||% integer(0)
    a_code <- match(a, c("A","B","C","D"))
    is_in_openset[i] <- a_code %in% os
    if (a %in% c("A","B","C")) {
      k <- match(a, arm_labels)
      t0 <- arm_open[k]; t1 <- arm_close[k]
      is_valid_time[i] <- !is.na(t0) && tt >= t0 && (is.na(t1) || tt <= t1)
    } else {
      any_open <- any(!is.na(arm_open) & !is.na(arm_close) & tt >= arm_open & tt <= arm_close)
      still_running <- any(!is.na(arm_open) & is.na(arm_close) & tt >= arm_open)
      is_valid_time[i] <- (any_open || still_running)
    }
  }
  df$valid_in_openset <- is_in_openset
  df$valid_time       <- is_valid_time
  df$valid            <- is_in_openset & is_valid_time
  
  # ============================================================
  # 5. Window summary table (planned / observed / used)
  # ============================================================
  win_tbl <- data.frame(
    Arm               = arm_labels,
    Planned_Start     = as.integer(c(arm_start["A"], arm_start["B"], arm_start["C"])),
    Observed_Open     = as.integer(obs_open),
    Close_Time        = as.integer(obs_close),
    Window_Start_Used = as.integer(arm_open),
    Window_End_Used   = as.integer(arm_close),
    Arm_n_realized    = as.integer(res$realized_sizes["arm", 1:3]),
    Ctrl_n_realized   = as.integer(res$realized_sizes["ctrl", 1:3]),
    stringsAsFactors = FALSE
  )
  win_tbl$Window_Length_Used <- with(win_tbl,
                                     ifelse(is.na(Window_Start_Used) | is.na(Window_End_Used),
                                            NA_integer_,
                                            Window_End_Used - Window_Start_Used + 1L)
  )
  
  cat("--- ARM WINDOW SUMMARY (planned vs. observed vs. used) ---\n")
  print(win_tbl, row.names = FALSE)
  
  cat("\n--- RANDOMIZATION SEQUENCE (first", n, "allocations) ---\n")
  print(df, row.names = FALSE)
  
  # ============================================================
  # 6. Validation messages
  # ============================================================
  if (!all(df$valid_in_openset)) {
    warning("Found allocations to arms that were NOT in the open set at that time.")
  }
  if (!all(df$valid_time)) {
    warning("Found allocations outside the analysis windows (with fallback).")
  }
  if (all(df$valid)) {
    cat("All first", n, "allocations are within the open set AND within the analysis windows. ✅\n")
  }
  
  invisible(list(settings = list(
    rand_mode       = rand_mode_used,
    block_factor    = block_factor_used,
    max_group_size  = max_group_size_used,
    concurrent_only = concurrent_only_used,
    alloc_bias       = alloc_bias_used,
    beta_time       = beta_time_used,
    expected_total  = expected_total_used),
    windows = win_tbl,
    seq = df))
}

# ============================================================
# (1.1) PLOTTING FUNCTION FOR TRIAL TIMELINE
#       Shows the allocations to each arm, open and closing windows of arms.
# ============================================================

plot_trial_timeline <- function(res, arm_start, title = "Platform trial timeline",
                                show_allocations = TRUE) {
  assign <- res$assign_seq
  times  <- res$assign_times
  if (is.null(assign) || length(assign) == 0L) {
    stop("Result has no assignment log (assign_seq). Did the trial randomize anyone?")
  }
  
  obs_open  <- (res$window_open  %||% rep(NA_integer_, 3L))
  obs_close <- (res$window_close %||% rep(NA_integer_, 3L))
  
  arm_open_used  <- obs_open
  arm_close_used <- obs_close
  for (k in 1:3) {
    x_times <- times[assign == k]
    if (length(x_times)) {
      if (is.na(arm_open_used[k]))  arm_open_used[k]  <- min(x_times)
      if (is.na(arm_close_used[k])) arm_close_used[k] <- max(x_times)
    }
  }
  
  arm_names <- c("A","B","C")
  arm_factor <- factor(arm_names, levels = rev(arm_names))
  
  df_windows <- data.frame(
    arm   = arm_factor,
    start = as.numeric(arm_open_used),
    end   = as.numeric(arm_close_used)
  )
  df_windows$start_plot <- ifelse(is.na(df_windows$start), df_windows$end, df_windows$start)
  df_windows$end_plot   <- ifelse(is.na(df_windows$end),   df_windows$start, df_windows$end)
  
  df_marks <- rbind(
    data.frame(arm = factor("A", levels = rev(arm_names)), t = obs_open[1],  what = "Arm opens"),
    data.frame(arm = factor("A", levels = rev(arm_names)), t = obs_close[1], what = "Arm closes"),
    data.frame(arm = factor("B", levels = rev(arm_names)), t = obs_open[2],  what = "Arm opens"),
    data.frame(arm = factor("B", levels = rev(arm_names)), t = obs_close[2], what = "Arm closes"),
    data.frame(arm = factor("C", levels = rev(arm_names)), t = obs_open[3],  what = "Arm opens"),
    data.frame(arm = factor("C", levels = rev(arm_names)), t = obs_close[3], what = "Arm closes")
  )
  df_marks <- df_marks[!is.na(df_marks$t), , drop = FALSE]
  
  # Planned starts + precompute numeric y for a small vertical tick
  df_planned <- data.frame(
    arm = arm_factor,
    t   = as.numeric(c(arm_start["A"], arm_start["B"], arm_start["C"]))
  )
  df_planned$y <- as.numeric(df_planned$arm)
  
  p <- ggplot()
  
  # Used window (thick bar)
  df_seg <- subset(df_windows, !is.na(start_plot) & !is.na(end_plot) & end_plot >= start_plot)
  if (nrow(df_seg) > 0) {
    p <- p + geom_segment(
      data = df_seg,
      aes(x = start_plot, xend = end_plot, y = arm, yend = arm, color = arm),
      linewidth = 8, lineend = "round", alpha = 0.35
    )
  }
  
  # Arm opens/close markers
  if (nrow(df_marks) > 0) {
    p <- p + geom_point(
      data = df_marks,
      aes(x = t, y = arm, shape = what),
      size = 3, stroke = 1.1, color = "black", fill = "white"
    )
  }
  
  # Planned start as a small vertical tick
  p <- p + geom_segment(
    data = df_planned,
    aes(x = t, xend = t, y = y - 0.35, yend = y + 0.35),
    inherit.aes = FALSE,
    linewidth = 0.6, color = "black"
  ) +
    geom_text(
      data = df_planned,
      aes(x = t, y = arm, label = "P"),
      vjust = -1.2, size = 3
    )
  
  if (isTRUE(show_allocations)) {
    df_alloc <- data.frame(
      t   = times,
      arm = c("A","B","C","D")[as.integer(assign)]
    )
    df_alloc$arm_f <- factor(ifelse(df_alloc$arm == "D", "D (control)", df_alloc$arm),
                             levels = c("C","B","A","D (control)"))
    
    p <- p +
      geom_point(
        data = subset(df_alloc, arm != "D"),
        aes(x = t, y = factor(arm, levels = rev(arm_names)), color = arm),
        size = 1.8, alpha = 0.85, position = position_jitter(height = 0.06, width = 0)
      ) +
      geom_point(
        data = subset(df_alloc, arm == "D"),
        aes(x = t, y = factor("D (control)", levels = c("C","B","A","D (control)"))),
        size = 1.6, alpha = 0.6, color = "gray40", shape = 16
      ) 
  }
  
  p +
    scale_color_manual(values = c("A" = "#1b9e77", "B" = "#d95f02", "C" = "#7570b3")) +
    scale_shape_manual(values = c("Arm opens" = 21, "Arm closes" = 24)) +
    labs(title = title, x = "Patient randomized", y = NULL, color = "Arm", shape = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
}


# ============================================================
# (2) TYPE I ERROR GRID (no bias): expect ≈ 0.05
#     concurrent_only ∈ {TRUE,FALSE}
#     rand_mode ∈ {"complete"(CR),"block"(PBR)}
#     arm_start B ∈ {0, 50, 100, 150}
# ============================================================
run_type1_grid <- function(
    nsim = 2000,  
    max_group_size = 50,
    block_factor = 1
) {
  grid <- expand.grid(
    concurrent_only = c(TRUE, FALSE),
    rand_mode = c("complete","block"),
    startB = c(0L, 50L, 100L, 150L),
    stringsAsFactors = FALSE
  )
  out <- data.frame(grid, A=NA_real_, B=NA_real_, C=NA_real_, any=NA_real_)
  for (i in seq_len(nrow(grid))) {
    cfg <- grid[i, ]
    arm_start <- c(A=1L, B=cfg$startB, C=1L)
    rej_mat <- matrix(FALSE, nrow = nsim, ncol = 3)
    for (s in 1:nsim) {
      res <- platform_trials_simulation(
        max_group_size  = max_group_size,
        mu              = c(A=0, B=0, C=0, D=0),  # null true
        alpha           = 0.05,
        arm_start       = arm_start,
        concurrent_only = cfg$concurrent_only,
        expected_total  = 200,
        beta_time       = 0,                      # no time trend
        rand_mode       = cfg$rand_mode,          # CR or PBR
        block_factor    = block_factor,
        alloc_bias       = 0                       # NO bias
      )
      rej_mat[s, ] <- unname(res$reject)
    }
    rates <- colMeans(rej_mat)
    out$A[i]   <- rates[1]
    out$B[i]   <- rates[2]
    out$C[i]   <- rates[3]
    out$any[i] <- mean(rowSums(rej_mat) > 0)
    cat(sprintf("Done %d/%d: conc=%s, mode=%s, startB=%d  |  A=%.3f B=%.3f C=%.3f FWER=%.3f\n",
                i, nrow(grid), cfg$concurrent_only, cfg$rand_mode, cfg$startB,
                out$A[i], out$B[i], out$C[i], out$any[i]))
  }
  out
}

# ============================
# (3) Responder type tracker
#    logs per-patient expected & actual allocation
# ============================
trace_randomization_sequence <- function(
    max_group_size  = 20,
    mu              = c(A=0, B=0, C=0, D=0),
    alpha           = 0.05,
    arm_start       = c(A=0, B=20, C=0),
    concurrent_only = TRUE,
    rand_mode       = c("block","complete"),
    block_factor    = 1,
    expected_total  = 200,
    beta_time       = 0,
    alloc_bias      = 0.25,
    alloc_bias       = NULL,
    max_steps       = 1000,   
    seed            = 1
) {
  set.seed(seed)
  rand_mode <- match.arg(rand_mode)
  if (!is.null(alloc_bias)) alloc_bias <- alloc_bias
  
  # ---- fixed mapping and starts
  mu_vec         <- mu[c("A","B","C","D")]
  arm_starts_vec <- c(arm_start["A"], arm_start["B"], arm_start["C"], 0L)
  
  # ---- counters (global and period-local)
  arm_counts   <- integer(4L)  # total A..D
  local_counts <- integer(4L)  # period-local A..D
  target_exp   <- rep(max_group_size, 3L)
  
  # ---- window bookkeeping
  t_open            <- rep(NA_integer_, 3L)  # A,B,C
  close_time        <- rep(NA_integer_, 3L)
  ctrl_n            <- integer(3L)
  global_ctrl_total <- 0L
  finished_mask     <- rep(FALSE, 3L)
  
  # ---- randomization engine (identical signature logic)
  current_block   <- integer(0)
  block_pos       <- 0L
  block_signature <- -1L
  .sig_bitmask <- function(open_codes) sum(bitwShiftL(1L, open_codes - 1L))
  build_block   <- function(open_codes) sample(rep(open_codes, each = block_factor))
  next_assignment <- function(open_codes) {
    if (rand_mode == "complete") return(sample(open_codes, 1L))
    sig <- .sig_bitmask(open_codes)
    need_new <- (length(current_block) == 0L) ||
      (block_pos >= length(current_block)) ||
      (sig != block_signature)
    if (need_new) {
      current_block   <<- build_block(open_codes)
      block_pos       <<- 0L
      block_signature <<- sig
    }
    block_pos <<- block_pos + 1L
    current_block[[block_pos]]
  }
  
  # ---- period tracking
  period_idx     <- 1L
  period_t_start <- 1L
  .start_new_period <- function(t_now) {
    period_idx     <<- period_idx + 1L
    period_t_start <<- t_now
    local_counts[] <<- 0L
  }
  
  # ---- bias policy
  decide_bias <- function(open_codes, local_counts, alloc_bias) {
    # A=1, B=2, C=3, D=4
    if (length(open_codes) == 0L) return(list(expected=NA_integer_, bias_cat="neu", bias_val=0))
    counts_open <- local_counts[open_codes]
    min_count   <- min(counts_open)
    eq_min      <- open_codes[counts_open == min_count]
    
    if (!(2L %in% open_codes)) {
      # B not open -> researcher does not try to bias
      return(list(expected=NA_integer_, bias_cat="neu", bias_val=0))
    } else if ((4L %in% eq_min) && (2L %in% eq_min)) {
      # control tied with B -> neutral
      return(list(expected=2L, bias_cat="neu", bias_val=0))
    } else if (2L %in% eq_min) {
      # B favored when B least among open codes (ties with A/C allowed)
      return(list(expected=2L, bias_cat="pos", bias_val=+alloc_bias))
    } else if (4L %in% eq_min) {
      # control least -> negative responder
      return(list(expected=4L, bias_cat="neg", bias_val=-alloc_bias))
    } else {
      # all other cases -> neutral
      return(list(expected=eq_min[1L], bias_cat="neu", bias_val=0))
    }
  }
  
  # ---- trace storage
  trace <- data.frame(
    t = integer(0),
    period = integer(0),
    open = character(0),
    A_local = integer(0), B_local = integer(0), C_local = integer(0), D_local = integer(0),
    expected = character(0),
    bias_cat = character(0),
    bias_val = numeric(0),
    assigned = character(0),
    stringsAsFactors = FALSE
  )
  
  # ---- main loop
  t <- 0L
  steps <- 0L
  while (steps < max_steps && !all(finished_mask)) {
    steps <- steps + 1L
    t <- t + 1L
    
    # which arms are eligible (open by time and not finished)?
    open_codes <- which(arm_starts_vec <= t)
    if (any(finished_mask)) {
      fin_idx <- which(finished_mask)
      open_codes <- open_codes[!(open_codes %in% fin_idx)]
    }
    if (any(open_codes %in% 1:3) && !(4L %in% open_codes)) open_codes <- c(open_codes, 4L)
    
    # new experimental opens -> start new period
    newly_eligible <- intersect(which(is.na(t_open)), intersect(1:3, open_codes))
    if (length(newly_eligible)) {
      .start_new_period(t)
      t_open[newly_eligible] <- t
      if (!concurrent_only) ctrl_n[newly_eligible] <- global_ctrl_total
    }
    
    # if no experimental arm open, continue time
    if (!any(open_codes %in% 1:3)) next
    
    # take a snapshot of local counts BEFORE the assignment
    snap_local <- local_counts
    
    # decide bias according to policy
    pol <- decide_bias(open_codes, local_counts, alloc_bias)
    
    # actual assignment via randomization
    arm_code <- next_assignment(open_codes)
    
    # update counts
    arm_counts[arm_code]   <- arm_counts[arm_code]   + 1L
    local_counts[arm_code] <- local_counts[arm_code] + 1L
    
    # control accounting for windows
    if (arm_code %in% 1:3) {
      k <- arm_code
      # check close if both arm and matching controls hit target
      if (!finished_mask[k] && arm_counts[k] >= max_group_size && ctrl_n[k] >= max_group_size) {
        close_time[k]    <- t
        finished_mask[k] <- TRUE
        .start_new_period(t + 1L)  # next time step starts new period
      }
    } else {
      # assigned control
      global_ctrl_total <- global_ctrl_total + 1L
      open_mask <- !is.na(t_open) & !finished_mask
      if (any(open_mask)) ctrl_n[open_mask] <- ctrl_n[open_mask] + 1L
      
      # when controls push some arm(s) over the line
      if (any(!finished_mask)) {
        newly_done <- which(!finished_mask &
                              (arm_counts[1:3] >= max_group_size) &
                              (ctrl_n >= max_group_size))
        if (length(newly_done)) {
          close_time[newly_done]    <- t
          finished_mask[newly_done] <- TRUE
          .start_new_period(t + 1L)
        }
      }
    }
    
    # write a trace row
    open_lbl <- paste(c("A","B","C","D")[open_codes], collapse=",")
    expected_lbl <- if (is.na(pol$expected)) "none" else c("A","B","C","D")[pol$expected]
    assigned_lbl <- c("A","B","C","D")[arm_code]
    trace <- rbind(trace, data.frame(
      t = t,
      period = period_idx,
      open = open_lbl,
      A_local = snap_local[1], B_local = snap_local[2],
      C_local = snap_local[3], D_local = snap_local[4],
      expected = expected_lbl,
      bias_cat = pol$bias_cat,
      bias_val = pol$bias_val,
      assigned = assigned_lbl,
      stringsAsFactors = FALSE
    ))
  }
  
  # return the full trace and a small summary
  list(
    trace = trace,
    final_counts = setNames(arm_counts, c("A","B","C","D")),
    closed_at = setNames(close_time, c("A","B","C")),
    periods_recorded = unique(trace$period)
  )
}

# ============================
# (4) Power validation
#    to be finished...
# ============================
summarize_runs(n_sim = 500,
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

summarize_runs(n_sim = 500,
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

# ============================================================
# --- Run the validations ------------------------------------
# ============================================================

# 1) Randomization + analysis audit 
ex1 <- platform_trials_simulation_logged(
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

ex2 <- platform_trials_simulation_logged(
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

ex3 <- platform_trials_simulation_logged(
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


ex4 <- platform_trials_simulation_logged(
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

# (1) Sequence audit
audit_randomization(ex1, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex2, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex3, arm_start = c(A=1, B=10, C=1), N = 150)
audit_randomization(ex4, arm_start = c(A=1, B=10, C=1), N = 150)

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



# 2) Type I error grid (no bias)
type1_results <- run_type1_grid(
  nsim = 60000,         
  max_group_size = 30, 
  block_factor = 1
)

cat("\n=== TYPE I ERROR SUMMARY (target ~ 0.05) ===\n")
print(type1_results, row.names = FALSE)


# 3) Responder type tracking for allocation biasing policy validation
tr <- trace_randomization_sequence(
  max_group_size = 5,
  arm_start = c(A=0, B=8, C=0),
  rand_mode = "block",
  block_factor = 2,
  alloc_bias = 0.3,
  seed = 42
)
print(tr)
