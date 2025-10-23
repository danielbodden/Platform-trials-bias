# ============================================================
#  Platform trials simulation with allocation bias & time trend
#  for different randomization schemes: Complete randomization & permuted blocks.
#  Randomization process restarts whenever a new arm opens or closes.
#  Treatment arm closes when both it and its matching controls >= max_group_size.
# ============================================================


# Manual implementation of standard deviation for speed
fast_sd <- function(x, mean_x) sqrt(sum((x - mean_x)^2) / (length(x) - 1))

# Two-sample t-test with pooled variance
# returns list with t statistic and degrees of freedom
two_sample_t_pooled <- function(x, y) {
  nx <- length(x); ny <- length(y)
  mx <- mean(x);  my <- mean(y)
  sx <- fast_sd(x, mx); sy <- fast_sd(y, my)
  sp2 <- ((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2)
  se  <- sqrt(sp2 * (1/nx + 1/ny))
  list(t = (mx - my) / se, df = nx + ny - 2)
}

# Main simulation function
platform_trials_simulation <- function(
    max_group_size  = 50,                      # per experimental arm target n
    mu              = c(A=0, B=0, C=0, D=0),   # true means by arm
    alpha           = 0.05,                    # two-sided alpha
    arm_start       = c(A=0, B=100, C=0),      # Patient timing at which experimental arm opens (D always 0)
    concurrent_only = TRUE,                    # controls must be within [t_open, t_end] if TRUE
    rand_mode       = c("block", "complete"),    # randomization mode
    block_factor    = 1,                       # repeats per code in a block (only for "block")
    expected_total  = 200,                     # used for time trend scaling
    beta_time       = 0,                       # linear time trend coefficient
    alloc_bias      = 0                       # allocation bias
) {
  rand_mode <- match.arg(rand_mode)

  # --- fixed order A,B,C,D (codes 1..4)
  mu_vec         <- mu[c("A","B","C","D")]
  arm_starts_vec <- c(arm_start["A"], arm_start["B"], arm_start["C"], 0L)
  
  # --- counters
  arm_counts   <- integer(4L)  # A,B,C,D
  local_counts <- integer(4L)
  target_exp   <- rep(max_group_size, 3L)
  
  # --- storage (grown if needed)
  cap          <- max(1024L, as.integer(2.5 * expected_total))
  assign_i     <- integer(cap)
  times_i      <- integer(cap)
  alloc_bias_i <- numeric(cap)   
  
  idx <- 0L; t <- 0L
  
  
  # --- Randomization
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
  
  # --- Arm opening and closing times
  t_open            <- rep(NA_integer_, 3L)  # A,B,C
  close_time        <- rep(NA_integer_, 3L)
  ctrl_n            <- integer(3L)
  global_ctrl_total <- 0L
  
  # --- finished flags (speed)
  finished_mask <- rep(FALSE, 3L)
  
  # ============================
  # NEW: period tracking (each time an arm opens or closes we start a new period)
  #      We keep a single data frame with counts per arm (A,B,C,D) and bias category (pos/neg/neu).
  # ============================
  period_idx <- 1L
  period_t_start <- 1L
  current_counts <- matrix(0L, nrow = 4L, ncol = 3L,
                           dimnames = list(c("A","B","C","D"), c("pos","neg","neu")))
  period_counts_df <- data.frame(
    period = integer(0), t_start = integer(0), t_end = integer(0),
    arm = factor(character(0), levels = c("A","B","C","D")),
    pos = integer(0), neg = integer(0), neu = integer(0),
    stringsAsFactors = FALSE
  )
  # helper: flush current period counts into the data frame and start a new period
  .flush_and_start_new_period <- function(t_now) {
    # append only if there has been any allocation in the current period
    if (sum(current_counts) > 0L) {
      # wide rows per arm
      add <- data.frame(
        period = period_idx,
        t_start = period_t_start,
        t_end = t_now - 1L,
        arm = factor(rownames(current_counts), levels = c("A","B","C","D")),
        pos = current_counts[, "pos"],
        neg = current_counts[, "neg"],
        neu = current_counts[, "neu"],
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      period_counts_df <<- rbind(period_counts_df, add)
    }
    # reset for next period
    period_idx <<- period_idx + 1L
    period_t_start <<- t_now
    current_counts[,] <<- 0L
  }
  
  # ============================
  # main loop
  # ============================
  repeat {
    if (all(finished_mask)) break
    t <- t + 1L
    
    # which arms have started by now?
    open_codes <- which(arm_starts_vec <= t)
    
    # drop finished experimental arms
    if (any(finished_mask)) {
      fin_idx <- which(finished_mask)
      open_codes <- open_codes[!(open_codes %in% fin_idx)]
    }
    
    # ensure control available whenever any experimental is open
    if (any(open_codes %in% 1:3) && !(4L %in% open_codes)) open_codes <- c(open_codes, 4L)
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Nset t_open at the moment an experimental arm first becomes ELIGIBLE
    newly_eligible <- intersect(which(is.na(t_open)), intersect(1:3, open_codes))
    if (length(newly_eligible)) {
      # opening an arm starts a new period
      .flush_and_start_new_period(t)
      local_counts[] <- 0L   # <-- add this line
      
      t_open[newly_eligible] <- t
      # for non-concurrent rule: all past controls count at opening
      if (!concurrent_only) ctrl_n[newly_eligible] <- global_ctrl_total
    }
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # if no experimental arm is open yet, just advance time
    if (!any(open_codes %in% 1:3)) next
    
    # grow storage if needed
    if (idx >= cap) {
      grow     <- max(512L, as.integer(0.5 * cap))
      assign_i <- c(assign_i, integer(grow))
      times_i  <- c(times_i,  integer(grow))
      alloc_bias_i <- c(alloc_bias_i, numeric(grow))
      cap      <- cap + grow
    }
    # expected-allocation (for allocation biasing policy)
    counts_open <- local_counts[open_codes]
    min_count   <- min(counts_open)
    eq_min      <- open_codes[counts_open == min_count]
    
    # Allocation bias policy:
    # - If Arm B (code 2) is NOT open -> neutral only (no attempt to bias)
    # - If control (code 4) AND arm B are both tied in this period -> neutral
    # - If arm B alone is among the least -> positive responder (+alloc_bias)
    # - If control alone is among the least -> negative responder (-alloc_bias)
    # - Otherwise -> neutral
    if (!(2L %in% open_codes)) {
      alloc_bias_next <- 0;            bias_cat <- "neu"
    } else if ((4L %in% eq_min) && (2L %in% eq_min)) {
      alloc_bias_next <- 0;            bias_cat <- "neu"
    } else if (2L %in% eq_min) {
      alloc_bias_next <- +alloc_bias;  bias_cat <- "pos"
    } else if (4L %in% eq_min) {
      alloc_bias_next <- -alloc_bias;  bias_cat <- "neg"
    } else {
      alloc_bias_next <- 0;            bias_cat <- "neu"
    }
    
    # actual assignment
    arm_code <- next_assignment(open_codes)
    
    # record
    idx                    <- idx + 1L
    assign_i[idx]          <- arm_code
    times_i [idx]          <- t
    alloc_bias_i[idx]      <- alloc_bias_next
    arm_counts[arm_code]   <- arm_counts[arm_code]   + 1L
    local_counts[arm_code] <- local_counts[arm_code] + 1L
    
    # period-category bookkeeping by assigned arm
    current_counts[c("A","B","C","D")[arm_code], bias_cat] <-
      current_counts[c("A","B","C","D")[arm_code], bias_cat] + 1L
    
    # fast incremental window logic
    if (arm_code %in% 1:3) {
      k <- arm_code
      if (!finished_mask[k] && arm_counts[k] >= max_group_size && ctrl_n[k] >= max_group_size) {
        close_time[k]    <- t
        finished_mask[k] <- TRUE
        # NEW: closing an arm starts a new period
        .flush_and_start_new_period(t)
        local_counts[] <- 0L   # <-- add this line
      }
    } else {
      # control
      global_ctrl_total <- global_ctrl_total + 1L
      
      # increment matching controls for any arm whose window is open and not finished
      open_mask <- !is.na(t_open) & !finished_mask
      if (any(open_mask)) ctrl_n[open_mask] <- ctrl_n[open_mask] + 1L
      
      # any arms now meet BOTH quotas?
      if (any(!finished_mask)) {
        newly_done <- which(!finished_mask &
                              (arm_counts[1:3] >= max_group_size) &
                              (ctrl_n >= max_group_size))
        if (length(newly_done)) {
          close_time[newly_done]    <- t
          finished_mask[newly_done] <- TRUE
          # one or more closings at this time -> single new period
          .flush_and_start_new_period(t)
        }
      }
    }
  }
  
  # edge case: nobody randomized
  if (idx == 0L) {
    # flush empty and return
    return(list(
      reject            = setNames(rep(FALSE, 3), c("A_vs_D","B_vs_D","C_vs_D")),
      num_rej           = 0L,
      realized_sizes    = rbind(arm = c(0,0,0), ctrl = c(0,0,0)),
      total_randomized  = 0L,
      final_counts      = setNames(arm_counts, c("A","B","C","D")),
      window_open       = setNames(t_open,  c("A","B","C")),
      window_close      = setNames(close_time, c("A","B","C")),
      period_counts_df  = period_counts_df
    ))
  }
  
  # trim
  assign_i     <- assign_i[seq_len(idx)]
  times_i      <- times_i [seq_len(idx)]
  alloc_bias_i <- alloc_bias_i[seq_len(idx)]
  
  # outcomes
  s_t     <- pmin(times_i, expected_total) / expected_total
  mu_pat  <- mu_vec[assign_i] + beta_time * s_t + alloc_bias_i
  outcomes <- rnorm(idx, mean = mu_pat, sd = 1)
  
  # testing on windows [t_open[k], close_time[k]]
  reject   <- logical(3)
  realized <- matrix(0L, 2, 3, dimnames = list(c("arm","ctrl"), c("A_vs_D","B_vs_D","C_vs_D")))
  for (k in 1:3) {
    arm_open <- t_open[k]
    t_end    <- close_time[k]
    
    # arm observations within window
    x_idx <- which(assign_i == k & times_i >= arm_open & times_i <= t_end)
    nx    <- length(x_idx)
    
    # matching controls per rule
    if (concurrent_only) {
      y_idx <- which(assign_i == 4L & times_i >= arm_open & times_i <= t_end)
    } else {
      y_idx <- which(assign_i == 4L & times_i <= t_end)
    }
    ny <- length(y_idx)
    
    realized["arm",  k] <- nx
    realized["ctrl", k] <- ny
    
    # t-test needs at least two allocation for each group
    if (nx < 2 || ny < 2) { reject[k] <- FALSE; next }
    
    res <- two_sample_t_pooled(outcomes[x_idx], outcomes[y_idx])
    reject[k] <- abs(res$t) > qt(1 - alpha/2, df = res$df)
  }
  names(reject) <- c("A_vs_D","B_vs_D","C_vs_D")
  
  .flush_and_start_new_period(t + 1L)  # end the last period at t
  
  list(
    reject           = reject,
    num_rej          = sum(reject),
    realized_sizes   = realized,
    total_randomized = idx,
    final_counts     = setNames(arm_counts, c("A","B","C","D")),
    window_open      = setNames(t_open,     c("A","B","C")),
    window_close     = setNames(close_time, c("A","B","C")),
    period_counts_df = period_counts_df     
  )
}

platform_trials_simulation(rand_mode = "block")  # test run
platform_trials_simulation(rand_mode = "complete")  # test run

