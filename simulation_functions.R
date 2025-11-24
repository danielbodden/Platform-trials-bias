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
                         policy = c("favor_B","favor_all_exp","average"), n_exp) {
  policy <- match.arg(policy)
  if (length(open_codes) == 0L) {
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  ctrl_code <- n_exp + 1L
  if (!(ctrl_code %in% open_codes)) {
    # no control open -> neutral
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  ctrl_count <- local_counts[ctrl_code]
  exp_codes  <- intersect(open_codes, seq_len(n_exp))
  if (!length(exp_codes)) {
    # no experimental arm open -> neutral
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  if (policy == "favor_B") {
    # 1) Prefer B
    # If B open:
    #   Negative, if Control < B  (favor control)
    #   Positive, if Control > B  (favor B)
    #   Neutral,  if Control = B
    # If B not open:
    #   Negative
    B_code <- if (n_exp >= 2L) 2L else NA_integer_
    B_open <- !is.na(B_code) && (B_code %in% exp_codes)
    
    if (!B_open) {
      # B not open -> negative, favor control
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    }
    
    B_count <- local_counts[B_code]
    if (ctrl_count < B_count) {
      # Control < B -> negative, favor control
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    } else if (ctrl_count > B_count) {
      # Control > B -> positive, favor B
      return(list(expected = B_code, bias_cat = "pos", bias_val = +alloc_bias))
    } else {
      # Control = B -> neutral
      return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
    }
  }
  
  if (policy == "favor_all_exp") {
    # 2) Prefer all experimental
    # Negative, when Control ≤ any arm AND an arm exists with Control < arm
    # Positive, when any arm < Control
    # Neutral if all arms and control have equal size.
    # Mixed case (some arms smaller, some larger) -> POSITIVE (favor under-represented exp arms).
    exp_counts <- local_counts[exp_codes]
    
    all_equal   <- all(exp_counts == ctrl_count)
    any_greater <- any(exp_counts > ctrl_count)  # Arm > Control
    any_smaller <- any(exp_counts < ctrl_count)  # Arm < Control
    
    if (all_equal) {
      # all equal -> neutral
      return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
    }
    
    # Negative: at least one experimental arm has more patients than control,
    # and none has fewer -> clearly over-represented exp arms
    if (any_greater && !any_smaller) {
      # favor control
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    }
    
    # Positive: at least one experimental arm has fewer patients than control
    # This includes the mixed case (some >, some <).
    if (any_smaller) {
      # choose among experimental arms with fewer patients than control
      exp_less <- exp_codes[local_counts[exp_codes] < ctrl_count]
      if (!length(exp_less)) exp_less <- exp_codes
      expected <- sample(exp_less, 1L)
      return(list(expected = expected, bias_cat = "pos", bias_val = +alloc_bias))
    }
    
    # Fallback (shouldn't really happen): neutral
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  # policy == "average"
  # 3) Average policy
  # Positive, if control > mean(A,B)
  # Negative, if control < mean(A,B)
  # Neutral,  if control = mean(A,B)
  # (generalized to mean over open experimental arms)
  exp_counts <- local_counts[exp_codes]
  m_exp      <- mean(exp_counts)
  
  if (ctrl_count > m_exp) {
    # Control > mean -> positive, favor experimental (under-represented)
    exp_under <- exp_codes[local_counts[exp_codes] < ctrl_count]
    if (!length(exp_under)) exp_under <- exp_codes
    expected <- sample(exp_under, 1L)
    return(list(expected = expected, bias_cat = "pos", bias_val = +alloc_bias))
  } else if (ctrl_count < m_exp) {
    # Control < mean -> negative, favor control
    return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
  } else {
    # equal -> neutral
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
}

# ============================================================
# Bias policy depending on cohort (for two-step randomization)
# only implemented for a total of 3 arms: A, B, D (D is control)
# Uses cohort-specific control counts: ctrl_by_cohort[A], ctrl_by_cohort[B]
# ============================================================
.decide_bias_by_cohort <- function(cohort_code, open_codes, local_counts,
                                   n_exp, alloc_bias, policy,
                                   ctrl_by_cohort) {
  # cohort_code = 1L (A) or 2L (B)
  # control is always n_exp + 1
  ctrl_code <- n_exp + 1L
  
  if (!(cohort_code %in% c(1L, 2L))) {
    return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
  }
  
  # global arm and control counts (across cohorts)
  arm_count_global  <- local_counts[cohort_code]
  ctrl_count_global <- local_counts[ctrl_code]
  
  # cohort-specific control count: number of control patients that came
  # from this cohort in the two-step allocation
  ctrl_cohort_count <- ctrl_by_cohort[cohort_code]
  
  if (policy == "favor_B_2step") {
    # 4) Prefer B (2-step cohort-based)
    # - When both A and B are open and cohort B is chosen:
    #     Compare 2 * ctrl_B_count (only controls from cohort B)
    #     to B_count (global B count).
    #     If 2 * ctrl_B_count < B_count -> NEGATIVE (favor control)
    #     If 2 * ctrl_B_count > B_count -> POSITIVE (favor B)
    #     otherwise -> NEUTRAL
    # - When cohort A is chosen: always NEGATIVE (favor control),
    #   as long as B is open.
    
    B_code <- if (n_exp >= 2L) 2L else NA_integer_
    B_open <- !is.na(B_code) && (B_code %in% open_codes)
    
    if (is.na(B_code) || !B_open) {
      # B not open -> fall back: negative, favor control
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    }
    
    if (cohort_code == 1L) {
      # Cohort A chosen while B is open:
      # Always negative -> send patient to control, so B is indirectly favored.
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    }
    
    # Cohort B chosen
    B_count      <- local_counts[B_code]
    ctrl_B_count <- ctrl_by_cohort[2L]  # controls from cohort B only
    
    if (2L * ctrl_B_count < B_count) {
      # Too many B vs controls in cohort B -> NEGATIVE (favor control)
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    } else if (2L * ctrl_B_count > B_count) {
      # Too many controls (in B-cohort) vs B -> POSITIVE (favor B)
      return(list(expected = B_code, bias_cat = "pos", bias_val = +alloc_bias))
    } else {
      # balanced -> neutral
      return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
    }
  }
  
  if (policy == "favor_all_exp_2step") {
    # 5) Prefer all experimental (2-step, cohort-based)
    # Cohort X chosen:
    #   Compare 2 * ctrl_X_count (controls from cohort X only)
    #   with Arm_X_count (global):
    #   - If 2 * ctrl_X_count < Arm_X_count -> NEGATIVE (favor control)
    #   - If 2 * ctrl_X_count > Arm_X_count -> POSITIVE (favor that arm)
    #   - If equal -> NEUTRAL
    #
    # For 2-arm situations (only one experimental + control) we fall back
    # to the 1-step "favor_all_exp" policy in the main loop.
    
    ctrl_X_count <- ctrl_cohort_count
    arm_X_count  <- arm_count_global
    
    if (2L * ctrl_X_count < arm_X_count) {
      # Too many experimental patients vs cohort-specific controls -> favor control
      return(list(expected = ctrl_code, bias_cat = "neg", bias_val = -alloc_bias))
    } else if (2L * ctrl_X_count > arm_X_count) {
      # Too many controls in this cohort vs arm -> favor this experimental arm
      return(list(expected = cohort_code, bias_cat = "pos", bias_val = +alloc_bias))
    } else {
      # balanced -> neutral
      return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
    }
  }
  
  # Fallback (should not happen): neutral
  return(list(expected = NA_integer_, bias_cat = "neu", bias_val = 0))
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
    rand_mode       = c("block", "complete", "bigstick"),  # randomization mode
    block_factor    = 1,                       # repeats per code in a block (only for "block")
    expected_total  = 200,                     # used for time trend scaling
    beta_time       = 0,                       # linear time trend coefficient
    chronobias_incr = 0,                       # bias per period
    alloc_bias      = 0,                       # allocation bias 
    chronobias_type = c("linear", "stepwise",
                        "inv_u", "seasonal"), # type of chronological bias
    exp_arms        = c("A","B","C"),           # <<< ADD: flexible experimental arms
    test_side       = c("two.sided","one.sided"),
    alternative     = c("greater","less"),
    bias_policy     = c("favor_B","favor_all_exp","average",
                        "favor_B_2step","favor_all_exp_2step"),
    analysis_model  = c("ttest","lm_time"),
    return_detail   = FALSE,
    return_log      = FALSE,
    trial_sample_size = 2L * max_group_size,   # per-arm total (arm + concurrent ctrl)
    two_step        = FALSE,                   # enable two-step randomization
    bigstick_a      = 2L                       # max tolerated per-period imbalance for 'bigstick'
) {
  
  rand_mode   <- match.arg(rand_mode)
  test_side   <- match.arg(test_side)
  alternative <- match.arg(alternative)
  chronobias_type <- match.arg(chronobias_type)
  bias_policy <- match.arg(bias_policy)
  
  # 2-step-specific policies cannot be used in 1-step randomization
  if (!isTRUE(two_step) && bias_policy %in% c("favor_B_2step","favor_all_exp_2step")) {
    stop("bias_policy 'favor_B_2step' and 'favor_all_exp_2step' can only be used when two_step = TRUE")
  }
  
  # --- fixed order A,B,C,D (codes 1..4)
  # <<< EDIT:
  
  if (analysis_model == "anova_period" && isTRUE(concurrent_only)) {
    stop("anova_period is only allowed when concurrent_only = FALSE")
  }
  
  all_arms <- c(exp_arms, "D")
  mu_vec   <- mu[all_arms]
  n_exp    <- length(exp_arms)
  arm_starts_vec <- c(arm_start[exp_arms], 1L)
  
  # --- counters
  arm_counts   <- integer(length(all_arms))
  local_counts <- integer(length(all_arms))      # per-period counts
  target_exp   <- rep(max_group_size, n_exp)
  
  # --- storage (grown if needed)
  cap          <- max(1024L, as.integer(2.5 * expected_total))
  assign_i     <- integer(cap)
  times_i      <- integer(cap)
  alloc_bias_i <- numeric(cap)
  period_i     <- integer(cap)
  
  # --- trace storage (only if return_detail/log = TRUE) ---
  want_trace <- isTRUE(return_detail) || isTRUE(return_log)
  if (want_trace) {
    expected_code_i      <- integer(cap)
    bias_cat_i           <- character(cap)
    open_set_log         <- vector("list", cap)
    local_counts_log     <- matrix(NA_integer_, nrow = cap, ncol = length(all_arms),
                                   dimnames = list(NULL, all_arms))
    cohort_i             <- integer(cap)           # cohort used for this patient (1=A,2=B, NA otherwise)
    cohort_i[]           <- NA_integer_
    ctrl_by_cohort_log   <- matrix(0L, nrow = cap, ncol = n_exp,
                                   dimnames = list(NULL, exp_arms))
  } else {
    expected_code_i    <- bias_cat_i <- NULL
    open_set_log       <- NULL
    local_counts_log   <- NULL
    cohort_i           <- NULL
    ctrl_by_cohort_log <- NULL
  }
  
  idx <- 0L; t <- 0L
  
  # --- Randomization (1-step)
  current_block   <- integer(0)
  block_pos       <- 0L
  block_signature <- -1L
  .sig_bitmask <- function(open_codes) sum(bitwShiftL(1L, open_codes - 1L))
  build_block   <- function(open_codes) sample(rep(open_codes, each = block_factor))
  
  # Strict big-stick helper: exclude arms that would yield max-min > a
  .bigstick_choice <- function(codes, counts_vec, a) {
    if (!length(codes)) return(NA_integer_)
    counts_open <- counts_vec[codes]
    k <- length(codes)
    allowed <- logical(k)
    for (i in seq_len(k)) {
      tmp <- counts_open
      tmp[i] <- tmp[i] + 1L
      if ((max(tmp) - min(tmp)) <= a) {
        allowed[i] <- TRUE
      }
    }
    cand <- codes[allowed]
    if (!length(cand)) {
      # fallback: choose among arms with minimal current count
      min_c <- min(counts_open)
      cand <- codes[counts_open == min_c]
    }
    k2 <- length(cand)
    u  <- runif(1)
    j  <- 1L + as.integer(floor(u * k2))
    if (j > k2) j <- k2
    cand[j]
  }
  
  next_assignment <- function(open_codes) {
    stopifnot(length(open_codes) > 0L)
    
    # 1) Complete randomization
    if (rand_mode == "complete") {
      choice <- base::sample(open_codes, 1L)
      return(choice)
    }
    
    # 2) Big-stick randomization
    if (rand_mode == "bigstick") {
      counts_open <- local_counts[open_codes]
      if (!length(counts_open)) {
        choice <- base::sample(open_codes, 1L)
      } else {
        choice <- .bigstick_choice(open_codes, local_counts, bigstick_a)
      }
      return(choice)
    }
    
    # 3) Block randomization
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
  
  # Step 2 (within cohort: {A,D} or {B,D}) – separate blocks for A-branch and B-branch
  current_block_within   <- list(
    A = integer(0),
    B = integer(0)
  )
  block_pos_within       <- list(A = 0L, B = 0L)
  block_signature_within <- list(A = -1L, B = -1L)
  
  # Step-2 local pseudo-counts for big-stick within each cohort branch:
  #   exp1, exp2 = the two "copies" of the branch
  #   ctrl       = control
  local_counts_within <- list(
    A = c(exp1 = 0L, exp2 = 0L, ctrl = 0L),
    B = c(exp1 = 0L, exp2 = 0L, ctrl = 0L)
  )
  
  # Cohort-specific control counts (per period)
  # ctrl_by_cohort[1] = controls randomized via cohort A in current period
  # ctrl_by_cohort[2] = controls randomized via cohort B in current period
  ctrl_by_cohort <- integer(n_exp)
  names(ctrl_by_cohort) <- exp_arms
  
  .build_block_simple <- function(codes) sample(rep(codes, each = block_factor))
  
  .next_step1_cohort <- function(open_codes) {
    # Only valid if both A and B are open; otherwise fall back to one-step
    cohort_open <- intersect(open_codes, c(1L, 2L))
    if (length(cohort_open) != 2L) return(NA_integer_)
    
    if (rand_mode == "complete") {
      return(sample(cohort_open, 1L))
    }
    
    if (rand_mode == "bigstick") {
      counts_open <- local_counts[cohort_open]
      if (!length(counts_open)) {
        return(sample(cohort_open, 1L))
      }
      choice <- .bigstick_choice(cohort_open, local_counts, bigstick_a)
      return(choice)
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
  
  # Step-2 allocator (within cohort A or B)
  .next_step2_within <- function(branch, ctrl_code) {
    if (!(branch %in% c(1L, 2L))) {
      return(branch)  # defensive
    }
    key <- if (branch == 1L) "A" else "B"
    
    # COMPLETE RANDOMIZATION (2:1)
    if (rand_mode == "complete") {
      u <- runif(1)
      choice <- if (u < (2/3)) branch else ctrl_code
      return(choice)
    }
    
    # BIG-STICK ON PSEUDO-ARMS
    if (rand_mode == "bigstick") {
      lc <- local_counts_within[[key]]  # named vector: exp1, exp2, ctrl
      codes <- c("exp1", "exp2", "ctrl")
      
      choice_slot <- .bigstick_choice(codes, lc, bigstick_a)
      
      lc[choice_slot] <- lc[choice_slot] + 1L
      local_counts_within[[key]] <<- lc
      
      if (choice_slot %in% c("exp1", "exp2")) {
        return(branch)
      } else {
        return(ctrl_code)
      }
    }
    
    # BLOCK RANDOMIZATION: permuted (branch, branch, ctrl)
    sig      <- .sig_bitmask(c(branch, ctrl_code))
    need_new <- (length(current_block_within[[key]]) == 0L) ||
      (block_pos_within[[key]] >= length(current_block_within[[key]])) ||
      (sig != block_signature_within[[key]])
    
    if (need_new) {
      block <- c(branch, branch, ctrl_code)
      k <- length(block)
      if (k > 1L) {
        for (i in seq_len(k - 1L)) {
          u <- runif(1)
          j <- i + as.integer(floor(u * (k - i + 1L)))
          if (j > k) j <- k
          if (j != i) {
            tmp      <- block[i]
            block[i] <- block[j]
            block[j] <- tmp
          }
        }
      }
      current_block_within[[key]]   <<- block
      block_pos_within[[key]]       <<- 0L
      block_signature_within[[key]] <<- sig
    }
    
    block_pos_within[[key]] <<- block_pos_within[[key]] + 1L
    choice <- current_block_within[[key]][[block_pos_within[[key]]]]
    
    return(choice)
  }
  
  # Two-step wrapper: returns list(code, cohort) or NULL for fallback
  .two_step_next_assignment <- function(open_codes, n_exp, forced_cohort = NA_integer_) {
    exp_open <- intersect(open_codes, seq_len(n_exp))
    # Only implemented when A and B are both open
    if (length(exp_open) != 2L || !all(sort(exp_open) == c(1L, 2L))) {
      return(NULL)  # signal to fall back to one-step
    }
    ctrl_code <- n_exp + 1L
    
    cohort_code <- forced_cohort
    if (is.na(cohort_code)) {
      cohort_code <- .next_step1_cohort(open_codes)
      if (is.na(cohort_code)) return(NULL)
    }
    
    code <- .next_step2_within(cohort_code, ctrl_code)
    list(code = code, cohort = cohort_code)
  }
  
  # --- Arm opening and closing times
  t_open     <- rep(NA_integer_, n_exp)
  close_time <- rep(NA_integer_, n_exp)
  ctrl_n     <- integer(n_exp)
  global_ctrl_total <- 0L
  
  # --- finished flags (speed)
  finished_mask <- rep(FALSE, n_exp)
  
  # ============================
  # period tracking
  # ============================
  period_idx <- 0L
  period_t_start <- 1L
  current_counts <- matrix(0L, nrow = length(all_arms), ncol = 3L,
                           dimnames = list(all_arms, c("pos","neg","neu")))
  period_counts_df <- data.frame(
    period = integer(0), t_start = integer(0), t_end = integer(0),
    arm = factor(character(0), levels = all_arms),
    pos = integer(0), neg = integer(0), neu = integer(0),
    stringsAsFactors = FALSE
  )
  .flush_and_start_new_period <- function(t_now) {
    if (sum(current_counts) > 0L) {
      add <- data.frame(
        period = period_idx,
        t_start = period_t_start,
        t_end = t_now - 1L,
        arm = factor(rownames(current_counts), levels = all_arms),
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
    
    # Reset per-period counts
    local_counts[]            <<- 0L
    local_counts_within$A[]   <<- 0L
    local_counts_within$B[]   <<- 0L
    ctrl_by_cohort[]          <<- 0L   # IMPORTANT: per-period reset
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
    if (any(open_codes %in% seq_len(n_exp)) && !((n_exp+1L) %in% open_codes)) {
      open_codes <- c(open_codes, n_exp+1L)
    }
    
    newly_eligible <- intersect(which(is.na(t_open)), intersect(seq_len(n_exp), open_codes))
    if (length(newly_eligible)) {
      .flush_and_start_new_period(t)
      t_open[newly_eligible] <- t
      # ctrl_n stays 0 for new windows (concurrent controls only)
    }
    
    if (!any(open_codes %in% seq_len(n_exp))) next
    
    if (idx >= cap) {
      grow     <- max(512L, as.integer(0.5 * cap))
      assign_i <- c(assign_i, integer(grow))
      times_i  <- c(times_i,  integer(grow))
      alloc_bias_i <- c(alloc_bias_i, numeric(grow))
      period_i     <- c(period_i,     integer(grow))
      
      if (want_trace) {
        expected_code_i    <- c(expected_code_i, integer(grow))
        bias_cat_i         <- c(bias_cat_i,      character(grow))
        open_set_log       <- c(open_set_log,    vector("list", grow))
        add <- matrix(
          NA_integer_, nrow = grow, ncol = ncol(local_counts_log),
          dimnames = list(NULL, colnames(local_counts_log))
        )
        local_counts_log   <- rbind(local_counts_log, add)
        cohort_i           <- c(cohort_i, rep(NA_integer_, grow))
        add2 <- matrix(
          0L, nrow = grow, ncol = ncol(ctrl_by_cohort_log),
          dimnames = list(NULL, colnames(ctrl_by_cohort_log))
        )
        ctrl_by_cohort_log <- rbind(ctrl_by_cohort_log, add2)
      }
      
      cap <- cap + grow
    }
    # ---------------- Allocation bias policy ----------------
    alloc_bias_next <- 0
    bias_cat        <- "neu"
    expected_code   <- NA_integer_
    current_cohort  <- NA_integer_
    
    if (isTRUE(two_step) && bias_policy %in% c("favor_B_2step","favor_all_exp_2step")) {
      exp_open <- intersect(open_codes, seq_len(n_exp))
      if (length(exp_open) == 2L && all(c(1L,2L) %in% exp_open)) {
        current_cohort <- .next_step1_cohort(open_codes)
        pol <- .decide_bias_by_cohort(
          cohort_code   = current_cohort,
          open_codes    = open_codes,
          local_counts  = local_counts,
          n_exp         = n_exp,
          alloc_bias    = alloc_bias,
          policy        = bias_policy,
          ctrl_by_cohort= ctrl_by_cohort
        )
      } else {
        # fall back to corresponding 1-step policy if cohort-based policy not applicable
        base_policy <- sub("_2step$", "", bias_policy)
        pol <- .decide_bias(open_codes, local_counts, alloc_bias, base_policy, n_exp)
      }
    } else {
      pol <- .decide_bias(open_codes, local_counts, alloc_bias, bias_policy, n_exp)
    }
    
    alloc_bias_next <- pol$bias_val
    bias_cat        <- pol$bias_cat
    expected_code   <- pol$expected
    
    # ---------------- Draw actual assignment ----------------
    assignment_cohort <- NA_integer_
    arm_code <- NA_integer_
    
    if (isTRUE(two_step)) {
      res2 <- .two_step_next_assignment(open_codes, n_exp, forced_cohort = current_cohort)
      if (!is.null(res2)) {
        arm_code          <- res2$code
        assignment_cohort <- res2$cohort
      } else {
        arm_code <- next_assignment(open_codes)
      }
    } else {
      arm_code <- next_assignment(open_codes)
    }
    
    # ---------------- Update indices ----------------
    idx               <- idx + 1L
    assign_i[idx]     <- arm_code
    times_i[idx]      <- t
    alloc_bias_i[idx] <- alloc_bias_next
    period_i[idx]     <- period_idx
    
    # ---------------- TRACE LOGGING (PRE-ASSIGNMENT STATE) ----------------
    # At this point local_counts and ctrl_by_cohort still reflect *previous* patients.
    if (want_trace) {
      expected_code_i[idx]    <- if (is.na(expected_code)) NA_integer_ else as.integer(expected_code)
      bias_cat_i[idx]         <- bias_cat
      open_set_log[[idx]]     <- open_codes
      local_counts_log[idx, ] <- local_counts              # BEFORE increment
      ctrl_by_cohort_log[idx, ] <- ctrl_by_cohort          # BEFORE increment
      cohort_i[idx]           <- assignment_cohort
    }
    
    # ---------------- Now update counts for this assignment ----------------
    arm_counts[arm_code]   <- arm_counts[arm_code]   + 1L
    local_counts[arm_code] <- local_counts[arm_code] + 1L
    
    # update cohort-specific control counts *after* assignment
    ctrl_code <- n_exp + 1L
    if (!is.na(assignment_cohort) &&
        arm_code == ctrl_code &&
        assignment_cohort %in% seq_len(n_exp)) {
      ctrl_by_cohort[assignment_cohort] <- ctrl_by_cohort[assignment_cohort] + 1L
    }
    
    current_counts[all_arms[arm_code], bias_cat] <-
      current_counts[all_arms[arm_code], bias_cat] + 1L
    
    
    # ---------------- Arm closure logic ----------------
    if (arm_code %in% seq_len(n_exp)) {
      k <- arm_code
      if (!finished_mask[k] && (arm_counts[k] + ctrl_n[k]) >= trial_sample_size) {
        close_time[k]    <- t
        finished_mask[k] <- TRUE
        .flush_and_start_new_period(t)
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
      responder_counts  = rc,
      ctrl_by_cohort    = ctrl_by_cohort,
      trace_df          = NULL
    )
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
  
  # trim vectors / matrices
  assign_i     <- assign_i[seq_len(idx)]
  times_i      <- times_i [seq_len(idx)]
  alloc_bias_i <- alloc_bias_i[seq_len(idx)]
  period_i     <- period_i[seq_len(idx)]
  if (want_trace && idx > 0L) {
    expected_code_i      <- expected_code_i[seq_len(idx)]
    bias_cat_i           <- bias_cat_i[seq_len(idx)]
    open_set_log         <- open_set_log[seq_len(idx)]
    local_counts_log     <- local_counts_log[seq_len(idx), , drop = FALSE]
    cohort_i             <- cohort_i[seq_len(idx)]
    ctrl_by_cohort_log   <- ctrl_by_cohort_log[seq_len(idx), , drop = FALSE]
  }
  alloc_bias_i <- alloc_bias_i[seq_len(idx)]
  period_i     <- period_i[seq_len(idx)]          # <<< ADD
  
  # linear vs stepwise vs inverted u vs seasonal chronological bias
  if(chronobias_type == "linear") {
    s_t     <- pmin(times_i, expected_total) / expected_total
    chr_bias <- beta_time * s_t
  } else if (chronobias_type == "stepwise") {
    chr_bias <- chronobias_incr[period_i]
  } else if (chronobias_type == "inv_u") {
    
    mid_point <- expected_total / 2
    times <- pmin(times_i, expected_total)
    chr_bias <- ifelse(times < mid_point,
                       beta_time * times/expected_total,
                       -beta_time * (times - mid_point)/expected_total + 
                         beta_time * (mid_point/expected_total))

  } else if (chronobias_type == "seasonal") {
    chr_bias <- beta_time * sin(4 * pi * times_i / expected_total)
  }
  
  mu_pat  <- mu_vec[assign_i] + chr_bias + alloc_bias_i
  outcomes <- rnorm(idx, mean = mu_pat, sd = 1)
  #print(mu_pat)
  
  responder_cat <- ifelse(alloc_bias_i > 0, "pos",
                          ifelse(alloc_bias_i < 0, "neg", "neu"))
  .df_all <- data.frame(
    t         = times_i,
    code      = assign_i,
    arm_label = all_arms[assign_i],
    responder = factor(responder_cat, levels = c("pos","neg","neu")),
    stringsAsFactors = FALSE
  )
  
  responder_counts <- matrix(
    0L, nrow = 2L * n_exp, ncol = 3L,
    dimnames = list(
      c(exp_arms, paste0("D_for_", exp_arms)),
      c("pos","neg","neu")
    )
  )
  
  for (k in seq_len(n_exp)) {
    arm_name  <- exp_arms[k]
    dfor_name <- paste0("D_for_", arm_name)
    
    t0 <- t_open[k]
    t1 <- close_time[k]
    
    idx_arm <- which(.df_all$code == k & !is.na(t0) & !is.na(t1) &
                       .df_all$t >= t0 & .df_all$t <= t1)
    if (length(idx_arm)) {
      tab_arm <- table(.df_all$responder[idx_arm])
      responder_counts[arm_name, names(tab_arm)] <- as.integer(tab_arm)
    }
    
    idx_ctrl <- which(.df_all$arm_label == "D" & !is.na(t0) & !is.na(t1) &
                        .df_all$t >= t0 & .df_all$t <= t1)
    if (length(idx_ctrl)) {
      tab_ctrl <- table(.df_all$responder[idx_ctrl])
      responder_counts[dfor_name, names(tab_ctrl)] <- as.integer(tab_ctrl)
    }
  }
  
  # per-arm bias/performance metrics (unchanged)
  if (idx > 0L) {
    arm_names <- names(mu_vec)
    n_arms    <- length(arm_names)
    alloc_vec <- alloc_bias_i
    
    mse_col <- sprintf("mse_outcome_vs_%.3f", alpha)
    
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
    
    add_rows <- paste0("D_for_", exp_arms)
    add_mat  <- matrix(NA_real_, nrow = length(add_rows), ncol = 3,
                       dimnames = list(add_rows, colnames(arm_metrics)))
    
    for (k in seq_len(n_exp)) {
      arm_open <- t_open[k]
      t_end    <- close_time[k]
      if (is.na(arm_open) || is.na(t_end) || t_end < arm_open) next
      
      in_win <- times_i >= arm_open & times_i <= t_end
      y_idx <- which(assign_i == (n_exp + 1L) & in_win)
      
      arm_metrics <- matrix(NA_real_, nrow = n_arms, ncol = 3,
                            dimnames = list(arm_names,
                                            c(mse_col,
                                              "mean_allocation_bias",
                                              "mean_chronological_bias")))
      for (code in seq_len(n_arms)) {
        idxs <- which(assign_i == code)
        if (length(idxs) == 0L) next
        
        arm_metrics[code, mse_col]                  <- mean((outcomes[idxs] - alpha)^2)
        arm_metrics[code, "mean_allocation_bias"]   <- mean(alloc_vec[idxs])
        
        if(chronobias_type == "linear") {
          arm_metrics[code, "mean_chronological_bias"] <- mean(beta_time * s_t[idxs])
        } else if(chronobias_type == "stepwise") {
          mean_bias <- mean(chronobias_incr[period_i[idxs]])
          arm_metrics[code, "mean_chronological_bias"] <- mean_bias
        } else if (chronobias_type == "inv_u") {
          mean_bias <- mean(ifelse(times_i[idxs] < mid_point,
                                   beta_time * times_i[idxs]/expected_total,
                                   -beta_time * (times_i[idxs] - mid_point)/expected_total + 
                                     beta_time * (mid_point/expected_total)))
          arm_metrics[code, "mean_chronological_bias"] <- mean_bias
        } else if (chronobias_type == "seasonal") {
          mean_bias <- mean(beta_time * sin(4 * pi * times_i[idxs] / expected_total))
          arm_metrics[code, "mean_chronological_bias"] <- mean_bias
        }
        
      }
    }
    
    arm_metrics <- rbind(arm_metrics, add_mat)
    
  } else {
    mse_col <- sprintf("mse_outcome_vs_%.3f", alpha)
    arm_metrics <- matrix(numeric(0), nrow = 0, ncol = 3,
                          dimnames = list(NULL,
                                          c(mse_col,
                                            "mean_allocation_bias",
                                            "mean_chronological_bias")))
  }
  
  # testing (unchanged)
  reject   <- logical(n_exp)
  realized <- matrix(0L, 2, n_exp, dimnames = list(c("arm","ctrl"), paste0(exp_arms, "_vs_D")))
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
    } else {
      y   <- c(outcomes[x_idx], outcomes[y_idx])
      trt <- c(rep("arm", length(x_idx)), rep("ctrl", length(y_idx)))
      per <- c(period_i[x_idx],            period_i[y_idx])
      
      ans <- .fit_lm_time(
        y, trt, per, alpha = alpha,
        test_side = test_side, alternative = alternative
      )
      reject[k] <- ans$reject
    }
  }
  names(reject) <- paste0(exp_arms, "_vs_D")
  
  .flush_and_start_new_period(t + 1L)
  
  trace_df <- NULL
  if (want_trace && idx > 0L) {
    code_to_lab <- function(code) ifelse(is.na(code), "none", all_arms[code])
    open_labs <- vapply(open_set_log, function(v) paste(all_arms[v], collapse=","), character(1))
    lc_df <- as.data.frame(local_counts_log, stringsAsFactors = FALSE)
    names(lc_df) <- paste0(names(lc_df), "_local")
    
    cohort_lab <- ifelse(is.na(cohort_i), "none", exp_arms[cohort_i])
    
    cbc_df <- as.data.frame(ctrl_by_cohort_log, stringsAsFactors = FALSE)
    names(cbc_df) <- paste0("ctrl_from_", names(cbc_df))
    
    trace_df <- data.frame(
      i        = seq_len(idx),
      t        = times_i,
      period   = period_i,
      open     = open_labs,
      expected = code_to_lab(expected_code_i),
      responder= bias_cat_i,
      bias_val = alloc_bias_i,
      assigned = all_arms[assign_i],
      cohort   = cohort_lab,
      lc_df,
      cbc_df,
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
    responder_counts  = responder_counts,
    ctrl_by_cohort    = ctrl_by_cohort,
    trace_df          = trace_df
  )
  
  res$rand_mode       <- rand_mode
  res$block_factor    <- block_factor
  res$max_group_size  <- max_group_size
  res$concurrent_only <- concurrent_only
  res$alloc_bias      <- alloc_bias
  res$beta_time       <- beta_time
  res$chronobias_incr <- chronobias_incr
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
                           chronobias_incr = 0,
                           chronobias_type = "linear",
                           rand_mode = "block",
                           block_factor = 1,
                           alloc_bias = 0,
                           two_step = FALSE,
                           bigstick_a = 2L) {
  res <- replicate(n_sim,
                   platform_trials_simulation(max_group_size=max_group_size,
                                              mu=mu,
                                              alpha=alpha,
                                              arm_start=arm_start,
                                              concurrent_only=concurrent_only,
                                              expected_total=expected_total,
                                              beta_time=beta_time,
                                              chronobias_incr=chronobias_incr,
                                              chronobias_type = chronobias_type,
                                              rand_mode=rand_mode,
                                              block_factor=block_factor,
                                              alloc_bias=alloc_bias,
                                              two_step=two_step,
                                              bigstick_a=bigstick_a),
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
                               two_step=FALSE,
                               bigstick_a = 2L) {
  if (is.null(n_cores)) {
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
    if (is.na(n_cores) || n_cores < 1) n_cores <- max(1L, parallel::detectCores(logical = FALSE))
  }
  n_cores <- max(1L, n_cores)
  
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
    two_step =two_step,
    bigstick_a = bigstick_a
  )
  
  cl <- parallel::makeCluster(length(chunk_sizes), type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
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

# ============================================================
#  Unified summary for both serial and parallel usage
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
    rand_mode       = c("block","complete","bigstick"),
    block_factor    = 1,
    alloc_bias      = 0,
    exp_arms        = NULL,
    seed            = NULL,
    verbose_every   = 0L,
    test_side       = c("two.sided","one.sided"),
    alternative     = c("greater","less"),
    bias_policy     = c("favor_B","favor_all_exp","average",
                        "favor_B_2step","favor_all_exp_2step"),
    analysis_model  = c("ttest","lm_time"),
    two_step = FALSE,
    bigstick_a = 2L
) {
  rand_mode   <- match.arg(rand_mode)
  test_side   <- match.arg(test_side)
  alternative <- match.arg(alternative)
  bias_policy <- match.arg(bias_policy)
  analysis_model <- match.arg(analysis_model)
  if (!is.null(seed)) set.seed(seed)
  
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
    two_step        = two_step,
    bigstick_a      = bigstick_a
  )
  
  rej_mat   <- matrix(FALSE, nrow = n_sim, ncol = n_exp)
  colnames(rej_mat) <- comp_names
  arm_n_mat  <- matrix(NA_integer_, nrow = n_sim, ncol = n_exp)
  ctrl_n_mat <- matrix(NA_integer_, nrow = n_sim, ncol = n_exp)
  colnames(arm_n_mat)  <- comp_names
  colnames(ctrl_n_mat) <- comp_names
  
  run_serial <- function() {
    for (i in seq_len(n_sim)) {
      if (verbose_every > 0L && (i %% verbose_every == 0L)) {
        message(sprintf("[calc_rejection_summary] %d/%d ...", i, n_sim))
      }
      ans <- do_one(base_args)
      
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
  
  per_cmp <- colMeans(rej_mat)
  fwer <- mean(rowSums(rej_mat) > 0)
  
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
      alloc_bias=alloc_bias, exp_arms=exp_arms, n_cores=n_cores, two_step=two_step,
      bigstick_a=bigstick_a
    )
  )
}


