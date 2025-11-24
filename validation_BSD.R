#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(parallel)
})

# Load your simulation code (with updated big-stick implementation)
source("simulation_functions.R")

# ------------------------------------------------------------
# Helper: recompute "allowed" arms under the strict big-stick rule
# as implemented by .bigstick_choice:
#   - codes: character or numeric vector of arm labels/codes
#   - counts_vec: named vector of counts indexed by labels
#   - a: max allowed imbalance (max - min)
# Returns the subset of `codes` that are allowed.
# ------------------------------------------------------------
.allowed_bigstick_codes <- function(codes, counts_vec, a) {
  if (!length(codes)) return(codes[FALSE])
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
    # Fallback behaviour in .bigstick_choice: restrict to arms with minimal count
    min_c <- min(counts_open)
    cand  <- codes[counts_open == min_c]
  }
  cand
}

# ------------------------------------------------------------
# Validate big-stick rule for a single trial:
#   - For one-step (two_step = FALSE): works on the global open set.
#   - For two-step:
#       * When both A and B are open (two-step path) and final assignment is A or B:
#           validate step 1 in the {A,B} space.
#       * When only one experimental arm or fallback is used:
#           validate as one-step in the open set.
#
# Returns:
#   - n_steps_checked
#   - n_violations
#   - violations_detail: list of violation objects
#   - steps_detail: list of ALL checked steps (for manual inspection)
# ------------------------------------------------------------
validate_bigstick_trial <- function(sim, bigstick_a, two_step = FALSE) {
  tr <- sim$trace_df
  if (is.null(tr) || nrow(tr) == 0L) {
    return(list(
      n_steps_checked   = 0L,
      n_violations      = 0L,
      violations_detail = list(),
      steps_detail      = list()
    ))
  }
  
  # find local-count columns: e.g. A_local, B_local, ...
  arm_cols <- grep("_local$", names(tr), value = TRUE)
  if (!length(arm_cols)) {
    warning("No *_local columns found in trace_df; cannot validate big-stick.")
    return(list(
      n_steps_checked   = 0L,
      n_violations      = 0L,
      violations_detail = list(),
      steps_detail      = list()
    ))
  }
  all_arms <- sub("_local$", "", arm_cols)
  
  # try to infer experimental arms (if available)
  exp_arms <- sim$exp_arms
  if (is.null(exp_arms)) {
    # simple fallback: all arms except "D"
    exp_arms <- setdiff(all_arms, "D")
  }
  
  n_steps_checked <- 0L
  violations      <- list()
  steps_detail    <- list()
  
  for (i in seq_len(nrow(tr))) {
    open_str <- tr$open[i]
    assigned <- as.character(tr$assigned[i])
    
    # parse open arms (labels)
    open_labels <- if (nzchar(open_str)) strsplit(open_str, ",", fixed = TRUE)[[1]] else character(0)
    open_labels <- open_labels[open_labels != ""]
    
    # if nothing open or assigned not in our known arms, skip
    if (!length(open_labels)) next
    if (!(assigned %in% all_arms)) next
    
    # counts BEFORE assignment
    before <- as.integer(tr[i, arm_cols])
    names(before) <- all_arms
    
    # two-step logic: step 1 uses cohorts A/B when both are open & two_step = TRUE
    exp_open <- intersect(open_labels, exp_arms)
    use_two_step_path <- isTRUE(two_step) &&
      length(exp_open) == 2L &&
      all(c("A","B") %in% exp_open)
    
    if (use_two_step_path) {
      # Two-step path:
      #   Step 1 big-stick is applied in {A,B} only, using local_counts on A,B.
      #   We can only check that when final assignment is A or B, since then
      #   we know which cohort was chosen.
      cohort_space <- intersect(c("A","B"), open_labels)
      if (length(cohort_space) == 2L && assigned %in% cohort_space) {
        allowed <- .allowed_bigstick_codes(cohort_space, before, bigstick_a)
        n_steps_checked <- n_steps_checked + 1L
        
        is_violation <- !(assigned %in% allowed)
        
        # store step info (for manual inspection)
        steps_detail[[length(steps_detail) + 1L]] <- list(
          step      = i,
          t         = tr$t[i],
          type      = "step1_two_step",
          open      = open_labels,
          cohort_space = cohort_space,
          assigned  = assigned,
          counts_before = before[cohort_space],
          allowed   = allowed,
          is_violation  = is_violation
        )
        
        if (is_violation) {
          violations[[length(violations) + 1L]] <- steps_detail[[length(steps_detail)]]
        }
      }
      # If final assignment is D (control), we don't know which cohort was chosen,
      # so we cannot validate step 1 for that step. Step 2 uses pseudo-arms,
      # which aren't observable in trace_df, so we don't validate it here.
      next
    }
    
    # One-step path (either two_step = FALSE, or fallback when two-step is not used):
    # Here big-stick is applied directly to open_codes (real arms).
    codes   <- open_labels
    allowed <- .allowed_bigstick_codes(codes, before, bigstick_a)
    
    # Only check steps where the assigned arm is in the big-stick decision space.
    if (!(assigned %in% codes)) next
    
    n_steps_checked <- n_steps_checked + 1L
    is_violation <- !(assigned %in% allowed)
    
    step_obj <- list(
      step      = i,
      t         = tr$t[i],
      type      = "one_step",
      open      = codes,
      cohort_space = NA_character_,
      assigned  = assigned,
      counts_before = before[codes],
      allowed   = allowed,
      is_violation  = is_violation
    )
    steps_detail[[length(steps_detail) + 1L]] <- step_obj
    
    if (is_violation) {
      violations[[length(violations) + 1L]] <- step_obj
    }
  }
  
  list(
    n_steps_checked   = n_steps_checked,
    n_violations      = length(violations),
    violations_detail = violations,
    steps_detail      = steps_detail
  )
}

# ------------------------------------------------------------
# Run many trials and summarize:
#   - For each trial, run validate_bigstick_trial()
#   - Show how many steps were checked and how many violations
#   - Build:
#       * violations_table: only violating steps
#       * steps_table: ALL checked steps (for manual inspection)
# ------------------------------------------------------------
validate_bigstick_many <- function(
    n_sim       = 200,
    bigstick_a  = 2L,
    two_step    = FALSE,
    verbose     = TRUE,
    ...
) {
  set.seed(20251118)
  
  total_steps_checked <- 0L
  total_violations    <- 0L
  per_trial           <- vector("list", n_sim)
  viol_rows           <- list()
  step_rows           <- list()
  
  for (i in seq_len(n_sim)) {
    if (verbose && (i %% 20L == 0L)) {
      cat("[validate_bigstick_many] trial", i, "of", n_sim, "...\n")
    }
    
    sim <- platform_trials_simulation(
      rand_mode     = "bigstick",
      bigstick_a    = bigstick_a,
      return_detail = TRUE,
      two_step      = two_step,
      ...
    )
    
    chk <- validate_bigstick_trial(sim, bigstick_a = bigstick_a, two_step = two_step)
    
    total_steps_checked <- total_steps_checked + chk$n_steps_checked
    total_violations    <- total_violations    + chk$n_violations
    per_trial[[i]]      <- chk
    
    ## ----- flatten ALL checked steps -----
    if (length(chk$steps_detail) > 0L) {
      for (st in chk$steps_detail) {
        counts_before_str <- paste(
          paste0(names(st$counts_before), "=", as.integer(st$counts_before)),
          collapse = ";"
        )
        allowed_str <- paste(st$allowed, collapse = ",")
        open_str    <- paste(st$open, collapse = ",")
        cohort_str  <- if (is.null(st$cohort_space) || all(is.na(st$cohort_space))) {
          NA_character_
        } else {
          paste(st$cohort_space, collapse = ",")
        }
        
        step_rows[[length(step_rows) + 1L]] <- data.frame(
          sim_id            = i,
          step              = st$step,
          t                 = st$t,
          type              = st$type,
          assigned          = st$assigned,
          open              = open_str,
          cohort_space      = cohort_str,
          counts_before_str = counts_before_str,
          allowed_str       = allowed_str,
          is_violation      = st$is_violation,
          stringsAsFactors  = FALSE
        )
      }
    }
  }
  
  steps_table <- if (length(step_rows)) {
    do.call(rbind, step_rows)
  } else {
    data.frame(
      sim_id            = integer(0),
      step              = integer(0),
      t                 = integer(0),
      type              = character(0),
      assigned          = character(0),
      open              = character(0),
      cohort_space      = character(0),
      counts_before_str = character(0),
      allowed_str       = character(0),
      is_violation      = logical(0),
      stringsAsFactors  = FALSE
    )
  }
  
  violations_table <- subset(steps_table, is_violation)
  
  cat("============================================\n")
  cat("Big-stick validation summary\n")
  cat("  two_step          =", two_step, "\n")
  cat("  a                 =", bigstick_a, "\n")
  cat("  n_sim             =", n_sim, "\n")
  cat("  steps checked     =", total_steps_checked, "\n")
  cat("  total violations  =", total_violations, "\n")
  cat("  violation rate    =", if (total_steps_checked > 0)
    sprintf("%.4f", total_violations / total_steps_checked)
    else "NA", "\n")
  cat("============================================\n")
  
  invisible(list(
    total_steps_checked = total_steps_checked,
    total_violations    = total_violations,
    per_trial           = per_trial,
    steps_table         = steps_table,
    violations_table    = violations_table
  ))
}

# ------------------------------------------------------------
# Example calls when run as a script
# ------------------------------------------------------------
if (identical(sys.nframe(), 0L)) {
  ## 1) One-step big-stick with 3 experimental arms always open
  res_one <- validate_bigstick_many(
    n_sim          = 200,
    bigstick_a     = 2L,
    two_step       = FALSE,
    max_group_size = 30,
    mu             = c(A=0, B=0, C=0, D=0),
    arm_start      = c(A=0, B=0, C=0),
    concurrent_only= TRUE,
    expected_total = 200,
    beta_time      = 0,
    alloc_bias     = 0,
    exp_arms       = c("A","B","C")
  )
  
  print(head(res_one$steps_table, 10))
  print(head(res_one$violations_table, 10))
  
  ## 2) Two-step big-stick with A,B as experimental arms
  res_two <- validate_bigstick_many(
    n_sim          = 200,
    bigstick_a     = 2L,
    two_step       = TRUE,
    max_group_size = 30,
    mu             = c(A=0, B=0, C=0, D=0),
    arm_start      = c(A=0, B=0, C=0),
    concurrent_only= TRUE,
    expected_total = 200,
    beta_time      = 0,
    alloc_bias     = 0,
    exp_arms       = c("A","B")
  )
  
  print(head(res_two$steps_table, 10))
  print(head(res_two$violations_table, 10))
}
