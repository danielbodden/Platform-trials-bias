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

# ---------- Bias policy helper (with explicit 'no arm' on control–exp ties) ----------
.decide_bias <- function(open_codes, local_counts, alloc_bias,
                         policy = c("favor_B","favor_all_exp"), n_exp) {
  policy <- match.arg(policy)
  if (length(open_codes) == 0L) {
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  ctrl_code   <- n_exp + 1L
  counts_open <- local_counts[open_codes]
  min_count   <- min(counts_open)
  eq_min      <- open_codes[counts_open == min_count]
  
  if (policy == "favor_B") {
    B_code <- if (n_exp >= 2L) 2L else NA_integer_
    
    # ---- your rule: if control and B are tied for least -> EXPECT NO ARM ----
    if (!is.na(B_code) && (ctrl_code %in% eq_min) && (B_code %in% eq_min)) {
      return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
    }
    
    # B uniquely among least -> positive
    if (!is.na(B_code) && (B_code %in% eq_min) && !(ctrl_code %in% eq_min)) {
      return(list(expected = B_code, bias_cat = "pos", bias_val = +alloc_bias))
    }
    
    # Control uniquely among least -> negative
    if (ctrl_code %in% eq_min && (is.na(B_code) || !(B_code %in% eq_min))) {
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    }
    
    # All other situations -> neutral, pick among least OPEN arms
    return(list(expected = sample(eq_min, 1L), bias_cat = "neu", bias_val = 0))
  }
  
  # policy == "favor_all_exp"
  exp_in_min  <- any(eq_min %in% seq_len(n_exp))
  ctrl_in_min <- (ctrl_code %in% eq_min)
  
  # ---- your rule: control & any experimental tied among least -> EXPECT NO ARM ----
  if (ctrl_in_min && exp_in_min) {
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  # Control least alone -> negative
  if (ctrl_in_min) {
    return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
  }
  
  # One or more experimental least (no control) -> positive (pick among exp ties)
  if (exp_in_min) {
    exp_eq <- intersect(eq_min, seq_len(n_exp))
    return(list(expected = sample(exp_eq, 1L), bias_cat = "pos", bias_val = +alloc_bias))
  }
  
  # Fallback neutral (shouldn’t really happen)
  return(list(expected = sample(eq_min, 1L), bias_cat = "neu", bias_val = 0))
}

# ============================================================
# Bias policy depending on cohort (for two-step randomization)
# ============================================================
.decide_bias_by_cohort <- function(cohort_code, local_counts, n_exp, alloc_bias) {
  # cohort_code = 1L (A) or 2L (B)
  # control is always n_exp + 1
  ctrl_code <- n_exp + 1L
  
  if (!(cohort_code %in% c(1L, 2L))) {
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  arm_count  <- local_counts[cohort_code]
  ctrl_count <- local_counts[ctrl_code]
  
  if (arm_count < ctrl_count) {
    # arm has fewer -> positive (favor arm)
    return(list(expected = cohort_code, bias_cat = "pos", bias_val = +alloc_bias))
  } else if (arm_count > ctrl_count) {
    # arm has more -> negative (favor control)
    return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
  } else {
    # equal -> neutral
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
}


.fit_lm_time <- function(y, trt, per, alpha = 0.05,
                         test_side = c("two.sided","one.sided"),
                         alternative = c("greater","less")) {
  test_side   <- match.arg(test_side)
  alternative <- match.arg(alternative)
  
  # Ensure contrasts encode "first period as baseline"
  trt <- factor(trt, levels = c("ctrl","arm"))
  per <- droplevels(factor(per))
  if (nlevels(per) < 2L) {
    fit <- lm(y ~ trt)
  } else {
    # default treatment contrasts: per’s first level is baseline -> ∑_{s=2}^{S} ν_s I(s_j=s)
    fit <- lm(y ~ trt + per)
  }
  sm <- summary(fit)$coefficients
  coef_name <- if ("trtarm" %in% rownames(sm)) "trtarm" else grep("^trt", rownames(sm), value = TRUE)[1]
  if (!length(coef_name)) return(list(reject = FALSE))
  
  beta_hat <- sm[coef_name, "Estimate"]
  p_two    <- sm[coef_name, "Pr(>|t|)"]
  
  if (test_side == "two.sided") {
    list(reject = is.finite(p_two) && (p_two < alpha))
  } else if (alternative == "greater") {
    list(reject = is.finite(p_two) && beta_hat > 0 && (p_two/2 < alpha))
  } else {
    list(reject = is.finite(p_two) && beta_hat < 0 && (p_two/2 < alpha))
  }
}



# Main simulation function
platform_trials_simulation <- function(
    max_group_size  = 50,                      # per experimental arm target n
    mu              = c(A=0, B=0, C=0, D=0),   # true means by arm
    alpha           = 0.05,                    # two-sided alpha
    arm_start       = c(A=0, B=100, C=0),      # Patient timing at which experimental arm opens (D always 0)
    concurrent_only = TRUE,                    # controls must be within [t_open, t_end] if TRUE
    rand_mode       = c("block", "complete"),  # randomization mode
    block_factor    = 1,                       # repeats per code in a block (only for "block")
    expected_total  = 200,                     # used for time trend scaling
    beta_time       = 0,                       # linear time trend coefficient
    alloc_bias      = 0,                       # allocation bias
    exp_arms        = c("A","B","C"),           # <<< ADD: flexible experimental arms
    test_side       = c("two.sided","one.sided"),
    alternative     = c("greater","less"),
    bias_policy     = c("favor_B","favor_all_exp"),
    analysis_model  = c("ttest","lm_time"),   # <<< ADD
    return_detail   = FALSE,         # <<< NEW
    return_log      = FALSE,          # <<< NEW (alias; either triggers the same trace)
    trial_sample_size = 2L * max_group_size,  # NEW: per-arm total (arm + concurrent ctrl)
    two_step        = FALSE              # <<< NEW: enable two-step randomization
    
) {
  if (analysis_model == "lm_time" && isTRUE(concurrent_only)) {
    stop("lm_time is only allowed when concurrent_only = FALSE")
  }
  rand_mode   <- match.arg(rand_mode)
  test_side   <- match.arg(test_side)
  alternative <- match.arg(alternative)
  bias_policy <- match.arg(bias_policy)
  analysis_model <- match.arg(analysis_model)     # <<< ADD  
  # --- fixed order A,B,C,D (codes 1..4)
  # <<< EDIT:
  all_arms <- c(exp_arms, "D")
  mu_vec   <- mu[all_arms]
  n_exp    <- length(exp_arms)
  arm_starts_vec <- c(arm_start[exp_arms], 0L)
  
  
  # --- counters
  arm_counts   <- integer(length(all_arms))  # <<< EDIT
  local_counts <- integer(length(all_arms))  # <<< EDIT
  target_exp   <- rep(max_group_size, n_exp) # <<< EDIT
  
  # --- storage (grown if needed)
  cap          <- max(1024L, as.integer(2.5 * expected_total))
  assign_i     <- integer(cap)
  times_i      <- integer(cap)
  alloc_bias_i <- numeric(cap)
  period_i     <- integer(cap)                    # <<< ADD
  
  
  # --- trace storage (only if return_detail/log = TRUE) ---
  want_trace <- isTRUE(return_detail) || isTRUE(return_log)
  if (want_trace) {
    expected_code_i  <- integer(cap)
    bias_cat_i       <- character(cap)
    open_set_log     <- vector("list", cap)
    local_counts_log <- matrix(NA_integer_, nrow = cap, ncol = length(all_arms),
                               dimnames = list(NULL, all_arms))
  } else {
    expected_code_i  <- bias_cat_i <- NULL
    open_set_log     <- NULL
    local_counts_log <- NULL
  }
  
  
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
  
  # --- Two-step randomization state (independent of the one-step block state) ---
  # Step 1 (cohort: A vs B)
  current_block_cohort   <- integer(0)
  block_pos_cohort       <- 0L
  block_signature_cohort <- -1L
  
  # Step 2 (within cohort: {A,D} or {B,D}) – we keep separate blocks for A-branch and B-branch
  current_block_within   <- list(
    A = integer(0),
    B = integer(0)
  )
  block_pos_within       <- list(A = 0L, B = 0L)
  block_signature_within <- list(A = -1L, B = -1L)
  
  .build_block_simple <- function(codes) sample(rep(codes, each = block_factor))
  
  .next_step1_cohort <- function(open_codes) {
    # Only valid if both A and B are open; otherwise fall back to one-step later
    # A code is 1L, B code is 2L under your fixed order (A,B,C,D)
    cohort_open <- intersect(open_codes, c(1L, 2L))
    if (length(cohort_open) != 2L) return(NA_integer_)
    if (rand_mode == "complete") {
      return(sample(cohort_open, 1L))
    }
    # block mode: signature by which cohorts are open (A,B)
    sig <- .sig_bitmask(cohort_open)
    need_new <- (length(current_block_cohort) == 0L) ||
      (block_pos_cohort >= length(current_block_cohort)) ||
      (sig != block_signature_cohort)
    if (need_new) {
      current_block_cohort   <<- .build_block_simple(cohort_open)
      block_pos_cohort       <<- 0L
      block_signature_cohort <<- sig
    }
    block_pos_cohort <<- block_pos_cohort + 1L
    current_block_cohort[[block_pos_cohort]]
  }
  
  # ============================================================
  # Step-2 allocator: permuted-block randomization (length 3)
  # Pattern: 2 × experimental + 1 × control  →  permute
  # Ensures global ≈ 1 : 1 : 1 allocation when both cohorts active.
  # ============================================================
  .next_step2_within <- function(branch, ctrl_code) {
    if (!(branch %in% c(1L, 2L))) return(branch)  # defensive
    key <- if (branch == 1L) "A" else "B"
    
    if (rand_mode == "complete") {
      # simple fallback: choose arm vs control with 2:1 probability
      probs <- c(2/3, 1/3)
      return(sample(c(branch, ctrl_code), 1L, prob = probs))
    }
    
    # --- block-mode (permuted blocks of length 3) ---
    sig <- .sig_bitmask(c(branch, ctrl_code))
    need_new <- (length(current_block_within[[key]]) == 0L) ||
      (block_pos_within[[key]] >= length(current_block_within[[key]])) ||
      (sig != block_signature_within[[key]])
    
    if (need_new) {
      # build new permuted block: two arms + one control
      block <- c(branch, branch, ctrl_code)
      current_block_within[[key]]   <<- sample(block, length(block))
      block_pos_within[[key]]       <<- 0L
      block_signature_within[[key]] <<- sig
    }
    
    block_pos_within[[key]] <<- block_pos_within[[key]] + 1L
    current_block_within[[key]][[block_pos_within[[key]]]]
  }
  
  
  
  .two_step_next_assignment <- function(open_codes, n_exp) {
    # Only implemented for exactly two experimental arms A and B being open
    exp_open <- intersect(open_codes, seq_len(n_exp))
    if (length(exp_open) != 2L || !all(sort(exp_open) == c(1L,2L))) {
      return(NA_integer_)  # signal to fall back to one-step
    }
    ctrl_code <- n_exp + 1L
    # Step 1: choose cohort (A vs B)
    cohort_code <- .next_step1_cohort(open_codes)
    if (is.na(cohort_code)) return(NA_integer_)  # fallback
    # Step 2: within chosen cohort, choose EXP vs D
    .next_step2_within(cohort_code, ctrl_code)
  }
  
  
  # --- Arm opening and closing times
  t_open     <- rep(NA_integer_, n_exp)        # <<< EDIT
  close_time <- rep(NA_integer_, n_exp)        # <<< EDIT
  ctrl_n     <- integer(n_exp)                 # <<< EDIT
  global_ctrl_total <- 0L
  
  # --- finished flags (speed)
  finished_mask <- rep(FALSE, n_exp)           # <<< EDIT
  
  # ============================
  # NEW: period tracking ...
  # ============================
  period_idx <- 1L
  period_t_start <- 1L
  current_counts <- matrix(0L, nrow = length(all_arms), ncol = 3L,
                           dimnames = list(all_arms, c("pos","neg","neu")))  # <<< EDIT
  period_counts_df <- data.frame(
    period = integer(0), t_start = integer(0), t_end = integer(0),
    arm = factor(character(0), levels = all_arms),  # <<< EDIT
    pos = integer(0), neg = integer(0), neu = integer(0),
    stringsAsFactors = FALSE
  )
  .flush_and_start_new_period <- function(t_now) {
    if (sum(current_counts) > 0L) {
      add <- data.frame(
        period = period_idx,
        t_start = period_t_start,
        t_end = t_now - 1L,
        arm = factor(rownames(current_counts), levels = all_arms),  # <<< EDIT
        pos = current_counts[, "pos"],
        neg = current_counts[, "neg"],
        neu = current_counts[, "neu"],
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      period_counts_df <<- rbind(period_counts_df, add)
    }
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
    
    open_codes <- which(arm_starts_vec <= t)
    if (any(finished_mask)) {
      fin_idx <- which(finished_mask)
      open_codes <- open_codes[!(open_codes %in% fin_idx)]
    }
    if (any(open_codes %in% seq_len(n_exp)) && !((n_exp+1L) %in% open_codes)) open_codes <- c(open_codes, n_exp+1L)  # <<< EDIT control always last
    
    newly_eligible <- intersect(which(is.na(t_open)), intersect(seq_len(n_exp), open_codes))
    if (length(newly_eligible)) {
      .flush_and_start_new_period(t)
      local_counts[] <- 0L
      t_open[newly_eligible] <- t
      # For stopping, we want concurrent control counts only:
      # -> do NOT seed ctrl_n with past controls, always start at 0
      # ctrl_n[newly_eligible] stays 0
    }
    
    
    if (!any(open_codes %in% seq_len(n_exp))) next
    
    if (idx >= cap) {
      grow     <- max(512L, as.integer(0.5 * cap))
      assign_i <- c(assign_i, integer(grow))
      times_i  <- c(times_i,  integer(grow))
      alloc_bias_i <- c(alloc_bias_i, numeric(grow))
      period_i     <- c(period_i,     integer(grow))
      
      # ---- 2.6: grow trace buffers when tracing is enabled ----
      if (want_trace) {
        expected_code_i <- c(expected_code_i, integer(grow))
        bias_cat_i      <- c(bias_cat_i,      character(grow))
        open_set_log    <- c(open_set_log,    vector("list", grow))
        add <- matrix(
          NA_integer_, nrow = grow, ncol = ncol(local_counts_log),
          dimnames = list(NULL, colnames(local_counts_log))
        )
        local_counts_log <- rbind(local_counts_log, add)
      }
      
      cap <- cap + grow
    }
    
    
    # Allocation bias policy (current period only) via helper
    pol <- .decide_bias(open_codes, local_counts, alloc_bias, bias_policy, n_exp)
    alloc_bias_next <- pol$bias_val
    bias_cat        <- pol$bias_cat
    expected_code   <- pol$expected
    
    if (isTRUE(two_step)) {
      # If two-step mode: determine cohort first for biasing policy
      exp_open <- intersect(open_codes, seq_len(n_exp))
      if (length(exp_open) == 2L && all(c(1L,2L) %in% exp_open)) {
        ctrl_code <- n_exp + 1L
        cohort_code <- .next_step1_cohort(open_codes)  # pre-sample cohort
        # compute bias based on that cohort
        pol <- .decide_bias_by_cohort(cohort_code, local_counts, n_exp, alloc_bias)
        alloc_bias_next <- pol$bias_val
        bias_cat        <- pol$bias_cat
        expected_code   <- pol$expected
      } else {
        # fallback to global bias rule if two-step inactive or not both arms open
        pol <- .decide_bias(open_codes, local_counts, alloc_bias, bias_policy, n_exp)
        alloc_bias_next <- pol$bias_val
        bias_cat        <- pol$bias_cat
        expected_code   <- pol$expected
      }
    } else {
      # regular one-step mode
      pol <- .decide_bias(open_codes, local_counts, alloc_bias, bias_policy, n_exp)
      alloc_bias_next <- pol$bias_val
      bias_cat        <- pol$bias_cat
      expected_code   <- pol$expected
    }
    
    arm_code <- if (isTRUE(two_step)) {
      ac2 <- .two_step_next_assignment(open_codes, n_exp)
      if (is.na(ac2)) next_assignment(open_codes) else ac2
    } else {
      next_assignment(open_codes)
    }
    
    
    if (want_trace) {
      
      # snapshot BEFORE making the assignment
      local_counts_log[idx + 1L, ] <- local_counts
    }
    
    
    
    idx                    <- idx + 1L
    assign_i[idx]          <- arm_code
    times_i [idx]          <- t
    alloc_bias_i[idx]      <- alloc_bias_next
    period_i[idx]          <- period_idx          # <<< ADD
    arm_counts[arm_code]   <- arm_counts[arm_code]   + 1L
    local_counts[arm_code] <- local_counts[arm_code] + 1L
    
    if (want_trace) {
      expected_code_i[idx] <- if (is.na(expected_code)) NA_integer_ else as.integer(expected_code)
      bias_cat_i[idx]      <- bias_cat
      open_set_log[[idx]]  <- open_codes
    }
    
    
    current_counts[all_arms[arm_code], bias_cat] <- current_counts[all_arms[arm_code], bias_cat] + 1L  # <<< EDIT
    
    if (arm_code %in% seq_len(n_exp)) {
      k <- arm_code
      if (!finished_mask[k] && (arm_counts[k] + ctrl_n[k]) >= trial_sample_size) {
        close_time[k]    <- t
        finished_mask[k] <- TRUE
        .flush_and_start_new_period(t)
        local_counts[] <- 0L
      }
      
    } else {
      global_ctrl_total <- global_ctrl_total + 1L
      open_mask <- !is.na(t_open) & !finished_mask
      if (any(open_mask)) ctrl_n[open_mask] <- ctrl_n[open_mask] + 1L
      if (any(!finished_mask)) {
        newly_done <- which(!finished_mask &
                              ((arm_counts[seq_len(n_exp)] + ctrl_n) >= trial_sample_size))
        
        if (length(newly_done)) {
          close_time[newly_done]    <- t
          finished_mask[newly_done] <- TRUE
          .flush_and_start_new_period(t)
        }
      }
    }
  }
  
  # edge case: nobody randomized
  if (idx == 0L) {
    # empty responder counts (2 rows per exp arm: Arm, D_for_Arm)
    rc <- matrix(0L, nrow = 2L * n_exp, ncol = 3L,
                 dimnames = list(
                   c(exp_arms, paste0("D_for_", exp_arms)),
                   c("pos","neg","neu")
                 ))
    
    res <- list(
      bias_metrics      = matrix(numeric(0), nrow = 0, ncol = 0),
      reject            = setNames(rep(FALSE, n_exp), paste0(exp_arms, "_vs_D")),
      num_rej           = 0L,
      realized_sizes    = rbind(arm = rep(0, n_exp), ctrl = rep(0, n_exp)),
      total_randomized  = 0L,
      final_counts      = setNames(arm_counts, all_arms),
      window_open       = setNames(t_open,  exp_arms),
      window_close      = setNames(close_time, exp_arms),
      period_counts_df  = period_counts_df,
      responder_counts  = rc
    )
    # --- Echo settings ---
    res$rand_mode       <- rand_mode
    res$block_factor    <- block_factor
    res$max_group_size  <- max_group_size
    res$concurrent_only <- concurrent_only
    res$alloc_bias      <- alloc_bias
    res$beta_time       <- beta_time
    res$expected_total  <- expected_total
    res$alpha           <- alpha
    res$test_side       <- test_side
    res$alternative     <- alternative
    res$bias_policy     <- bias_policy
    res$analysis_model  <- analysis_model
    res$exp_arms        <- exp_arms
    
    return(res)
  }
  
  
  
  # trim + outcomes
  assign_i     <- assign_i[seq_len(idx)]
  
  times_i      <- times_i [seq_len(idx)]
  if (want_trace && idx > 0L) {
    expected_code_i  <- expected_code_i[seq_len(idx)]
    bias_cat_i       <- bias_cat_i[seq_len(idx)]
    open_set_log     <- open_set_log[seq_len(idx)]
    local_counts_log <- local_counts_log[seq_len(idx), , drop = FALSE]
  }
  alloc_bias_i <- alloc_bias_i[seq_len(idx)]
  period_i     <- period_i[seq_len(idx)]          # <<< ADD
  # -- outcomes needed for bias_metrics (fixes 'outcomes not found') --
  s_t     <- pmin(times_i, expected_total) / expected_total
  mu_pat  <- mu_vec[assign_i] + beta_time * s_t + alloc_bias_i
  outcomes <- rnorm(idx, mean = mu_pat, sd = 1)
  
  
  # ---------------------------
  # Responder counts per arm and per concurrent-control window
  # ---------------------------
  # Derive responder category from alloc_bias_i sign (works without return_detail)
  responder_cat <- ifelse(alloc_bias_i > 0, "pos",
                          ifelse(alloc_bias_i < 0, "neg", "neu"))
  # Data for counting
  .df_all <- data.frame(
    t         = times_i,
    code      = assign_i,
    arm_label = all_arms[assign_i],
    responder = factor(responder_cat, levels = c("pos","neg","neu")),
    stringsAsFactors = FALSE
  )
  
  # Prepare result matrix:
  responder_counts <- matrix(
    0L, nrow = 2L * n_exp, ncol = 3L,
    dimnames = list(
      c(exp_arms, paste0("D_for_", exp_arms)),
      c("pos","neg","neu")
    )
  )
  
  # Fill counts for each experimental arm and its concurrent control
  for (k in seq_len(n_exp)) {
    arm_name  <- exp_arms[k]
    dfor_name <- paste0("D_for_", arm_name)
    
    # Arm window [open, close]
    t0 <- t_open[k]
    t1 <- close_time[k]
    
    # Arm A/B/... counts (they’re only assigned while open, but guard with window anyway)
    idx_arm <- which(.df_all$code == k & !is.na(t0) & !is.na(t1) &
                       .df_all$t >= t0 & .df_all$t <= t1)
    if (length(idx_arm)) {
      tab_arm <- table(.df_all$responder[idx_arm])
      responder_counts[arm_name, names(tab_arm)] <- as.integer(tab_arm)
    }
    
    # Concurrent controls for this arm: D during that arm’s window
    idx_ctrl <- which(.df_all$arm_label == "D" & !is.na(t0) & !is.na(t1) &
                        .df_all$t >= t0 & .df_all$t <= t1)
    if (length(idx_ctrl)) {
      tab_ctrl <- table(.df_all$responder[idx_ctrl])
      responder_counts[dfor_name, names(tab_ctrl)] <- as.integer(tab_ctrl)
    }
  }
  
  
  
  
  # --- per-arm bias/performance metrics --------------------
  # Mean squared error (outcome - alpha), mean allocation bias, mean chronological bias
  if (idx > 0L) {
    arm_names <- names(mu_vec)               # z.B. c("A","B","D")
    n_arms    <- length(arm_names)
    alloc_vec <- if (exists("alloc_bias_i")) alloc_bias_i else bias_i
    
    mse_col <- sprintf("mse_outcome_vs_%.3f", alpha)
    
    # Basismatrix: alle "globalen" Arm-Zeilen inkl. D (wie zuvor)
    arm_metrics <- matrix(NA_real_, nrow = n_arms, ncol = 3,
                          dimnames = list(arm_names,
                                          c(mse_col,
                                            "mean_allocation_bias",
                                            "mean_chronological_bias")))
    for (code in seq_len(n_arms)) {
      idxs <- which(assign_i == code)
      if (length(idxs) == 0L) next
      arm_metrics[code, mse_col]                    <- mean((outcomes[idxs] - alpha)^2)
      arm_metrics[code, "mean_allocation_bias"]     <- mean(alloc_vec[idxs])
      arm_metrics[code, "mean_chronological_bias"]  <- mean(beta_time * pmin(times_i[idxs], expected_total) / expected_total)
    }
    
    # >>> NEU: pro Experimentalarm k eine concurrent-Control-Zeile D_for_<Arm> hinzufügen
    # Zeilencontainer erweitern
    add_rows <- paste0("D_for_", exp_arms)
    add_mat  <- matrix(NA_real_, nrow = length(add_rows), ncol = 3,
                       dimnames = list(add_rows, colnames(arm_metrics)))
    
    for (k in seq_len(n_exp)) {
      arm_open <- t_open[k]
      t_end    <- close_time[k]
      if (is.na(arm_open) || is.na(t_end) || t_end < arm_open) next
      
      # Zeitfenster des Arms
      in_win <- times_i >= arm_open & times_i <= t_end
      
      # Concurrent-Controls für diesen Arm: Code = n_exp + 1 (D)
      y_idx <- which(assign_i == (n_exp + 1L) & in_win)
      
      if (length(y_idx) >= 1L) {
        add_mat[k, mse_col]                   <- mean((outcomes[y_idx] - alpha)^2)
        add_mat[k, "mean_allocation_bias"]    <- mean(alloc_vec[y_idx])
        add_mat[k, "mean_chronological_bias"] <- mean(beta_time * pmin(times_i[y_idx], expected_total) / expected_total)
      }
    }
    
    # zusammenführen: A, B, D  +  D_for_A, D_for_B, ...
    arm_metrics <- rbind(arm_metrics, add_mat)
    
  } else {
    mse_col <- sprintf("mse_outcome_vs_%.3f", alpha)
    arm_metrics <- matrix(numeric(0), nrow = 0, ncol = 3,
                          dimnames = list(NULL,
                                          c(mse_col,
                                            "mean_allocation_bias",
                                            "mean_chronological_bias")))
  }
  
  
  # testing
  reject   <- logical(n_exp)
  realized <- matrix(0L, 2, n_exp, dimnames = list(c("arm","ctrl"), paste0(exp_arms, "_vs_D")))  # <<< EDIT
  for (k in seq_len(n_exp)) {
    arm_open <- t_open[k]; t_end <- close_time[k]
    x_idx <- which(assign_i == k & times_i >= arm_open & times_i <= t_end)
    nx    <- length(x_idx)
    if (concurrent_only) {
      y_idx <- which(assign_i == (n_exp+1L) & times_i >= arm_open & times_i <= t_end)
    } else {
      y_idx <- which(assign_i == (n_exp+1L) & times_i <= t_end)
    }
    ny <- length(y_idx)
    realized["arm",  k] <- nx
    realized["ctrl", k] <- ny
    if (nx < 2 || ny < 2) { reject[k] <- FALSE; next }
    res <- two_sample_t_pooled(outcomes[x_idx], outcomes[y_idx])
    res <- two_sample_t_pooled(outcomes[x_idx], outcomes[y_idx])
    if (analysis_model == "ttest") {
      if (test_side == "two.sided") {
        reject[k] <- abs(res$t) > qt(1 - alpha/2, df = res$df)
      } else {
        if (alternative == "greater") {
          reject[k] <- res$t > qt(1 - alpha, df = res$df)
        } else {
          reject[k] <- res$t < qt(alpha, df = res$df)
        }
      }
    } else {  # analysis_model == "lm_time"
      # observations and covariates for the arm’s analysis window
      y   <- c(outcomes[x_idx], outcomes[y_idx])
      trt <- c(rep("arm", length(x_idx)), rep("ctrl", length(y_idx)))
      per <- c(period_i[x_idx],            period_i[y_idx])   # period fixed effects
      
      ans <- .fit_lm_time(
        y, trt, per, alpha = alpha,
        test_side = test_side, alternative = alternative
      )
      reject[k] <- ans$reject
    }
    
  }
  names(reject) <- paste0(exp_arms, "_vs_D")  # <<< EDIT
  
  .flush_and_start_new_period(t + 1L)
  
  
  trace_df <- NULL
  if (want_trace && idx > 0L) {
    # human-readable labels
    code_to_lab <- function(code) ifelse(is.na(code), "none", all_arms[code])
    # open set as labels
    open_labs <- vapply(open_set_log, function(v) paste(all_arms[v], collapse=","), character(1))
    # local counts (prefix with arm name)
    lc_df <- as.data.frame(local_counts_log, stringsAsFactors = FALSE)
    names(lc_df) <- paste0(names(lc_df), "_local")
    
    trace_df <- data.frame(
      i        = seq_len(idx),
      t        = times_i,
      period   = period_i,
      open     = open_labs,
      expected = code_to_lab(expected_code_i),
      responder= bias_cat_i,            # "pos"/"neg"/"neu"
      bias_val = alloc_bias_i,          # numeric added to outcome mean
      assigned = all_arms[assign_i],
      lc_df,
      stringsAsFactors = FALSE
    )
  }
  
  res <- list(
    bias_metrics      = arm_metrics,
    reject            = reject,
    num_rej           = sum(reject),
    realized_sizes    = realized,
    total_randomized  = idx,
    final_counts      = setNames(arm_counts, all_arms),
    window_open       = setNames(t_open,     exp_arms),
    window_close      = setNames(close_time, exp_arms),
    period_counts_df  = period_counts_df,
    responder_counts  = responder_counts,   # <<< NEW
    trace_df          = trace_df
  )
  
  
  # --- Echo settings so validators can print them ---
  res$rand_mode       <- rand_mode
  res$block_factor    <- block_factor
  res$max_group_size  <- max_group_size
  res$concurrent_only <- concurrent_only
  res$alloc_bias      <- alloc_bias
  res$beta_time       <- beta_time
  res$expected_total  <- expected_total
  res$alpha           <- alpha
  res$test_side       <- test_side
  res$alternative     <- alternative
  res$bias_policy     <- bias_policy
  res$analysis_model  <- analysis_model
  res$exp_arms        <- exp_arms
  
  return(res)
  
}

# Function to summarize runs of platform trials locally
summarize_runs <- function(n_sim = 500,
                           max_group_size = 50,
                           mu = c(A=0,B=0,C=0,D=0),
                           alpha = 0.05,
                           arm_start = c(A=0,B=100,C=0),
                           concurrent_only = TRUE,
                           expected_total = 200,
                           beta_time = 0,
                           rand_mode = "block",
                           block_factor = 1,
                           alloc_bias = 0,
                           two_step = FALSE) {         # <<< NEW
  res <- replicate(n_sim,
                   platform_trials_simulation(max_group_size=max_group_size,
                                              mu=mu,
                                              alpha=alpha,
                                              arm_start=arm_start,
                                              concurrent_only=concurrent_only,
                                              expected_total=expected_total,
                                              beta_time=beta_time,
                                              rand_mode=rand_mode,
                                              block_factor=block_factor,
                                              alloc_bias=alloc_bias,
                                              two_step=two_step),
                   simplify = FALSE)
  rej <- do.call(rbind, lapply(res, function(x) as.integer(x$reject)))
  colnames(rej) <- c("A_vs_D","B_vs_D","C_vs_D")
  list(
    per_comparison_rejection_rate = colMeans(rej, na.rm = TRUE),
    mean_final_ctrl_size = mean(sapply(res, function(x) x$final_counts["D"])),
    mean_arm_sizes       = colMeans(do.call(rbind, lapply(res, function(x) as.vector(x$realized_sizes["arm", ]))))
  )
}

# Summarizes run of platform_trials_simulations() by returning mean values of the key metrics
# using parallel processing
summarize_runs_par <- function(n_sim = 500,
                               max_group_size = 50,
                               mu = c(A=0,B=0,C=0,D=0),
                               alpha = 0.05,
                               arm_start = c(A=0,B=100,C=0),
                               concurrent_only = TRUE,
                               expected_total = 200,
                               beta_time = 0,
                               rand_mode = "block",
                               block_factor = 1,
                               alloc_bias = 0,
                               n_cores = NULL,
                               show_progress = TRUE,
                               two_step=FALSE) {
  if (is.null(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
    if (is.na(n_cores) || n_cores < 1) n_cores <- max(1L, parallel::detectCores(logical = FALSE))
  }
  n_cores <- max(1L, n_cores)
  
  # split n_sim into chunks
  chunk_sizes <- rep(n_sim %/% n_cores, n_cores)
  remainder <- n_sim %% n_cores
  if (remainder > 0) chunk_sizes[seq_len(remainder)] <- chunk_sizes[seq_len(remainder)] + 1L
  chunk_sizes <- chunk_sizes[chunk_sizes > 0]
  if (length(chunk_sizes) == 0) stop("n_sim must be >= 1")
  
  common_args <- list(
    max_group_size=max_group_size,
    mu=mu,
    alpha=alpha,
    arm_start=arm_start,
    concurrent_only=concurrent_only,
    expected_total=expected_total,
    beta_time=beta_time,
    rand_mode=rand_mode,
    block_factor=block_factor,
    alloc_bias=alloc_bias,
    two_step =two_step
  )
  
  cl <- parallel::makeCluster(length(chunk_sizes), type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  # avoid oversubscription inside workers
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1")
    NULL
  })
  parallel::clusterExport(cl, varlist = c("common_args"), envir = environment())
  parallel::clusterEvalQ(cl, source("PT_bias/simulation_functions.R"))
  
  seed_master <- as.integer(Sys.time()) %% .Machine$integer.max
  parallel::clusterSetRNGStream(cl, seed_master)
  
  worker_fun <- function(chunk_size, common_args) {
    sum_reject <- c(A_vs_D=0, B_vs_D=0, C_vs_D=0)
    sum_ctrl_D <- 0
    sum_arm_sz <- c(A=0, B=0, C=0)
    N <- 0L
    for (i in seq_len(chunk_size)) {
      x <- do.call(platform_trials_simulation, common_args)
      
      # defensive checks
      if (is.null(x$reject) || length(x$reject) < 3) next
      if (is.null(x$final_counts) || !"D" %in% names(x$final_counts)) next
      if (is.null(x$realized_sizes) || !all(c("A_vs_D","B_vs_D","C_vs_D") %in% colnames(x$realized_sizes))) next
      
      sum_reject <- sum_reject + as.integer(x$reject)
      sum_ctrl_D <- sum_ctrl_D + as.numeric(x$final_counts["D"])
      v <- as.numeric(x$realized_sizes["arm", c("A_vs_D","B_vs_D","C_vs_D")])
      names(v) <- c("A","B","C")
      sum_arm_sz <- sum_arm_sz + v
      N <- N + 1L
    }
    list(sum_reject=sum_reject, sum_ctrl_D=sum_ctrl_D, sum_arm_sz=sum_arm_sz, N=N)
  }
  
  if (show_progress) {
    message(sprintf("Running %d sims on %d core(s): chunks = %s",
                    n_sim, length(chunk_sizes), paste(chunk_sizes, collapse = "+")))
  }
  
  parts <- parallel::parLapply(cl, chunk_sizes, worker_fun, common_args = common_args)
  
  total_N    <- sum(vapply(parts, `[[`, integer(1), "N"))
  sum_reject <- Reduce(`+`, lapply(parts, `[[`, "sum_reject"))
  sum_ctrl_D <- sum(vapply(parts, `[[`, numeric(1), "sum_ctrl_D"))
  sum_arm_sz <- Reduce(`+`, lapply(parts, `[[`, "sum_arm_sz"))
  
  per_comp <- sum_reject / total_N
  names(per_comp) <- c("A_vs_D","B_vs_D","C_vs_D")
  
  list(
    per_comparison_rejection_rate = per_comp,
    mean_final_ctrl_size          = sum_ctrl_D / total_N,
    mean_arm_sizes                = sum_arm_sz / total_N
  )
}

#platform_trials_simulation(exp_arms = c("A","B"), rand_mode = "block")
#platform_trials_simulation(exp_arms = c("A","B"), rand_mode = "complete")


################


# ============================================================
#  Unified summary for both serial and parallel usage
#  - Use n_cores = 1 for serial (validation.r)
#  - Use n_cores > 1 for parallel (main.r / SLURM)
#  Returns the same structure expected by your scripts.
# ============================================================
calc_rejection_summary <- function(
    n_sim = 500,
    n_cores = 1L,
    max_group_size  = 50,
    mu              = c(A=0, B=0, C=0, D=0),
    alpha           = 0.05,
    arm_start       = c(A=0, B=100, C=0),
    concurrent_only = TRUE,
    expected_total  = 200,
    beta_time       = 0,
    rand_mode       = c("block","complete"),
    block_factor    = 1,
    alloc_bias      = 0,
    exp_arms        = NULL,   # optional; inferred if NULL
    seed            = NULL,
    verbose_every   = 0L,      # prints every k iterations in serial; ignored in parallel
    test_side       = c("two.sided","one.sided"),
    alternative     = c("greater","less"),
    bias_policy     = c("favor_B","favor_all_exp"),
    analysis_model  = c("ttest","lm_time"),  # rename from "anova_period"
    two_step = FALSE
) {
  rand_mode   <- match.arg(rand_mode)
  test_side   <- match.arg(test_side)
  alternative <- match.arg(alternative)
  bias_policy <- match.arg(bias_policy)
  analysis_model <- match.arg(analysis_model)
  if (!is.null(seed)) set.seed(seed)
  
  analysis_model <- match.arg(analysis_model)

  # --- helpers (local)
  infer_exp_arms <- function(mu, arm_start) {
    exp <- intersect(names(mu), c("A","B","C"))
    exp <- exp[exp != "D"]
    if (!length(exp)) exp <- c("A","B","C")
    if (!is.null(arm_start)) {
      exp <- intersect(exp, names(arm_start))
      if (!length(exp)) exp <- c("A","B","C")
    }
    exp
  }
  do_one <- function(args) {
    res <- do.call(platform_trials_simulation, args)
    # return both reject vector and realized sizes matrix
    list(reject = res$reject, sizes = res$realized_sizes)
  }
  
  if (is.null(exp_arms)) exp_arms <- infer_exp_arms(mu, arm_start)
  n_exp <- length(exp_arms)
  comp_names <- paste0(exp_arms, "_vs_D")
  
  base_args <- list(
    max_group_size  = max_group_size,
    mu              = mu,
    alpha           = alpha,
    arm_start       = arm_start,
    concurrent_only = concurrent_only,
    expected_total  = expected_total,
    beta_time       = beta_time,
    rand_mode       = rand_mode,
    block_factor    = block_factor,
    alloc_bias      = alloc_bias,
    exp_arms        = exp_arms,
    test_side       = test_side,
    alternative     = alternative,
    bias_policy     = bias_policy,
    analysis_model  = analysis_model,
    two_step = two_step
  )
  
  # storage for results
  rej_mat   <- matrix(FALSE, nrow = n_sim, ncol = n_exp)
  colnames(rej_mat) <- comp_names
  arm_n_mat  <- matrix(NA_integer_, nrow = n_sim, ncol = n_exp)
  ctrl_n_mat <- matrix(NA_integer_, nrow = n_sim, ncol = n_exp)
  colnames(arm_n_mat)  <- comp_names
  colnames(ctrl_n_mat) <- comp_names
  
  # --- serial runner
  run_serial <- function() {
    for (i in seq_len(n_sim)) {
      if (verbose_every > 0L && (i %% verbose_every == 0L)) {
        message(sprintf("[calc_rejection_summary] %d/%d ...", i, n_sim))
      }
      ans <- do_one(base_args)
      
      # rejects
      r <- ans$reject
      if (!is.null(names(r))) {
        rej_mat[i, colnames(rej_mat)] <<- unname(r[colnames(rej_mat)])
      } else {
        rej_mat[i, ] <<- unname(r[seq_len(n_exp)])
      }
      
      # realized sizes (rows: "arm","ctrl"; columns per comparison)
      siz <- ans$sizes
      if (!is.null(colnames(siz))) {
        idx <- match(comp_names, colnames(siz))
      } else {
        idx <- seq_len(n_exp)
      }
      arm_n_mat [i, ] <<- as.integer(siz["arm",  idx])
      ctrl_n_mat[i, ] <<- as.integer(siz["ctrl", idx])
    }
  }
  
  # --- parallel runner
  run_parallel <- function() {
    seeds <- sample.int(.Machine$integer.max, n_sim)
    .one_idx <- function(i) { set.seed(seeds[i]); do_one(base_args) }
    
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl,
                              varlist = c("platform_trials_simulation", "base_args", "seeds"),
                              envir = environment())
      res_list <- parallel::parLapply(cl, seq_len(n_sim), .one_idx)
    } else {
      res_list <- parallel::mclapply(seq_len(n_sim), .one_idx, mc.cores = n_cores)
    }
    
    for (i in seq_len(n_sim)) {
      ans <- res_list[[i]]
      
      r <- ans$reject
      if (!is.null(names(r))) {
        rej_mat[i, colnames(rej_mat)] <<- unname(r[colnames(rej_mat)])
      } else {
        rej_mat[i, ] <<- unname(r[seq_len(n_exp)])
      }
      
      siz <- ans$sizes
      if (!is.null(colnames(siz))) {
        idx <- match(comp_names, colnames(siz))
      } else {
        idx <- seq_len(n_exp)
      }
      arm_n_mat [i, ] <<- as.integer(siz["arm",  idx])
      ctrl_n_mat[i, ] <<- as.integer(siz["ctrl", idx])
    }
  }
  
  if (!is.null(n_cores) && n_cores >= 2L) run_parallel() else run_serial()
  
  # compute summaries
  per_cmp <- colMeans(rej_mat)
  fwer <- mean(rowSums(rej_mat) > 0)
  
  # sizes: means & sds
  sizes_summary <- data.frame(
    comp        = comp_names,
    arm_n_mean  = colMeans(arm_n_mat,  na.rm = TRUE),
    arm_n_sd    = apply(arm_n_mat,  2, function(x) sd(x,  na.rm = TRUE)),
    ctrl_n_mean = colMeans(ctrl_n_mat, na.rm = TRUE),
    ctrl_n_sd   = apply(ctrl_n_mat, 2, function(x) sd(x, na.rm = TRUE)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  list(
    per_comparison_rejection_rate = setNames(per_cmp, comp_names),
    fwer = fwer,
    n_sim = n_sim,
    sizes_summary = sizes_summary,
    settings = list(
      max_group_size=max_group_size, alpha=alpha,
      concurrent_only=concurrent_only, expected_total=expected_total,
      beta_time=beta_time, rand_mode=rand_mode, block_factor=block_factor,
      alloc_bias=alloc_bias, exp_arms=exp_arms, n_cores=n_cores, two_step=two_step
    )
  )
}
