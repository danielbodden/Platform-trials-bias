# ----------------------------
# runs per procedure
# ----------------------------
n_runs <- 1000

# ---------- one run: return both Δ and raw A / D_for_A responder counts ----------
one_run_extract <- function(rand_mode, block_factor, rep_seed,
                            B_start_fixed, max_group_size_exp,
                            alloc_bias_val, chrono_beta_val) {
  set.seed(rep_seed)
  res <- platform_trials_simulation(
    exp_arms        = c("A","B"),
    arm_start       = c(A=1L, B=B_start_fixed),
    max_group_size  = max_group_size_exp,
    expected_total  = 96,
    alpha           = 0.025,
    test_side       = "one.sided",
    alternative     = "greater",
    rand_mode       = rand_mode,
    block_factor    = block_factor,
    alloc_bias      = alloc_bias_val,
    beta_time       = chrono_beta_val,
    concurrent_only = TRUE,
    analysis_model  = "ttest",
    bias_policy     = "favor_all_exp"
  )
  bm <- res$bias_metrics
  if (is.null(bm) || !length(bm)) {
    return(list(delta = data.frame(), groupsA = data.frame()))
  }
  
  # --- Δs for plotting (arm − concurrent control) ---
  if (is.data.frame(bm)) bm <- as.matrix(bm)
  need_rows <- c("A","B","D_for_A","D_for_B")
  if (!all(need_rows %in% rownames(bm))) {
    return(list(delta = data.frame(), groupsA = data.frame()))
  }
  keep_cols <- intersect(colnames(bm), c("mean_allocation_bias","mean_chronological_bias"))
  if (!length(keep_cols)) return(list(delta = data.frame(), groupsA = data.frame()))
  
  A  <- bm["A", keep_cols, drop=TRUE]
  B  <- bm["B", keep_cols, drop=TRUE]
  DA <- bm["D_for_A", keep_cols, drop=TRUE]
  DB <- bm["D_for_B", keep_cols, drop=TRUE]
  
  delta_df <- data.frame(
    arm    = rep(c("Arm A − concurrent control","Arm B − concurrent control"), each = length(keep_cols)),
    metric = rep(keep_cols, times = 2L),
    value  = as.numeric(c(A - DA, B - DB)),
    run_id = rep_seed,
    stringsAsFactors = FALSE
  )
  
  # --- RAW responder counts for A and its concurrent control (D_for_A) ---
  rc <- res$responder_counts
  if (is.null(rc) || !all(c("A","D_for_A") %in% rownames(rc))) {
    groupsA_df <- data.frame()
  } else {
    groupsA_df <- data.frame(
      group  = c("Arm A", "Concurrent control of Arm A"),
      pos    = c(rc["A","pos"],        rc["D_for_A","pos"]),
      neg    = c(rc["A","neg"],        rc["D_for_A","neg"]),
      neu    = c(rc["A","neu"],        rc["D_for_A","neu"]),
      run_id = rep_seed,
      stringsAsFactors = FALSE
    )
  }
  
  list(delta = delta_df, groupsA = groupsA_df)
}

# ---------- collect: build both big tables (Δ for plot, raw A/D_for_A for counts) ----------
collect_all <- function() {
  cores_env <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  cores <- if (is.na(cores_env)) 1L else max(1L, cores_env)
  workers <- if (.Platform$OS.type == "windows") {
    max(1L, min(cores, as.integer(Sys.getenv("MAX_WORKERS_WIN", "4"))))
  } else cores
  batch_target <- as.integer(Sys.getenv("BATCH_TARGET", "100"))
  
  tasks <- list()
  for (p_i in seq_along(procedures)) {
    pr <- procedures[[p_i]]
    nbatches <- max(1L, ceiling(n_runs / batch_target))
    for (b in seq_len(nbatches)) {
      start <- (b - 1L) * batch_target + 1L
      stop  <- min(b * batch_target, n_runs)
      if (start > stop) next
      seed_vec <- seed_base + seq.int(start, stop) + 1000L * p_i
      tasks[[length(tasks)+1L]] <- list(
        rand_mode          = pr$rand_mode,
        block_factor       = pr$block_factor,
        seeds              = seed_vec,
        procedure          = pr$key,
        B_start_fixed      = B_start_fixed,
        max_group_size_exp = max_group_size_exp,
        alloc_bias_val     = alloc_bias_val,
        chrono_beta_val    = chrono_beta_val
      )
    }
  }
  
  run_task <- function(task) {
    deltas <- vector("list", length(task$seeds))
    groups <- vector("list", length(task$seeds))
    for (i in seq_along(task$seeds)) {
      ans <- one_run_extract(
        rand_mode          = task$rand_mode,
        block_factor       = task$block_factor,
        rep_seed           = task$seeds[[i]],
        B_start_fixed      = task$B_start_fixed,
        max_group_size_exp = task$max_group_size_exp,
        alloc_bias_val     = task$alloc_bias_val,
        chrono_beta_val    = task$chrono_beta_val
      )
      if (nrow(ans$delta))  { ans$delta$procedure  <- task$procedure;  deltas[[i]] <- ans$delta }
      if (nrow(ans$groupsA)){ ans$groupsA$procedure<- task$procedure;  groups[[i]] <- ans$groupsA }
    }
    list(
      delta  = if (length(Filter(NROW, deltas))) do.call(rbind, deltas) else data.frame(),
      groups = if (length(Filter(NROW, groups))) do.call(rbind, groups) else data.frame()
    )
  }
  
  parts <- if (.Platform$OS.type == "windows" && workers > 1L) {
    cl <- parallel::makeCluster(workers, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, {
      Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1")
      source("simulation_functions.R"); NULL
    })
    parallel::clusterExport(cl,
                            varlist=c("one_run_extract","procedures","n_runs","seed_base",
                                      "B_start_fixed","max_group_size_exp","alloc_bias_val","chrono_beta_val"),
                            envir=environment())
    parallel::parLapply(cl, tasks, run_task)
  } else if (.Platform$OS.type != "windows" && workers > 1L) {
    parallel::mclapply(tasks, run_task, mc.cores=workers)
  } else {
    lapply(tasks, run_task)
  }
  
  # bind both outputs across tasks
  delta_all  <- do.call(rbind, lapply(parts, `[[`, "delta"))
  groups_all <- do.call(rbind, lapply(parts, `[[`, "groups"))
  
  if (!nrow(delta_all))  stop("No Δ rows produced — check bias_metrics for A,B, D_for_A, D_for_B.")
  if (!nrow(groups_all)) stop("No A/D_for_A responder rows produced — check responder_counts.")
  
  # prettify labels in Δ table (for plotting)
  keep_map <- c("mean_allocation_bias"="Mean allocation bias",
                "mean_chronological_bias"="Mean chronological bias")
  delta_all <- subset(delta_all, metric %in% names(keep_map))
  delta_all$metric <- keep_map[delta_all$metric]
  delta_all$arm <- factor(delta_all$arm,
                          levels=c("Arm A − concurrent control","Arm B − concurrent control"))
  delta_all$procedure <- factor(delta_all$procedure,
                                levels=c("Complete randomization",
                                         "Block randomization (1*arms)",
                                         "Block randomization (2*arms)",
                                         "Block randomization (8*arms)"))
  
  # groups table already only concerns A / D_for_A responder counts
  groups_all$group <- factor(groups_all$group,
                             levels=c("Arm A","Concurrent control of Arm A"))
  groups_all$procedure <- factor(groups_all$procedure,
                                 levels=c("Complete randomization",
                                          "Block randomization (1*arms)",
                                          "Block randomization (2*arms)",
                                          "Block randomization (8*arms)"))
  
  list(delta = delta_all, groupsA = groups_all)
}

# ---------- plot & extra responder counts for A and D_for_A ----------
plot_and_save <- function(data_list, file_stub="bias_by_proc_SINGLE") {
  bias_rows  <- data_list$delta     # Δ data for plotting
  groupsA_df <- data_list$groupsA   # raw A / D_for_A responder counts (pos/neg/neu + run_id)
  
  # 1) write raw Δ data
  csv_path <- file.path(out_dir, paste0(file_stub, ".csv"))
  write.csv(bias_rows, csv_path, row.names=FALSE)
  
  # 2) plot Δ as before
  proc_cols <- c(
    "Complete randomization"           = "#1b9e77",
    "Block randomization (1*arms)"     = "#d95f02",
    "Block randomization (2*arms)"     = "#7570b3",
    "Block randomization (8*arms)"     = "#e7298a"
  )
  p <- ggplot2::ggplot(bias_rows, ggplot2::aes(x=arm, y=value, fill=procedure)) +
    ggplot2::geom_violin(trim=FALSE, alpha=0.45, color=NA) +
    ggplot2::geom_boxplot(width=0.20, outlier.shape=16, outlier.size=0.8, alpha=0.9) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", linewidth=0.3, color="black") +
    ggplot2::facet_grid(rows=ggplot2::vars(procedure), cols=ggplot2::vars(metric),
                        scales="fixed", switch="y") +
    ggplot2::scale_fill_manual(values=proc_cols, name="Randomization") +
    ggplot2::labs(
      title="Δ Bias vs Randomization Procedure (arm − concurrent control)",
      x="Contrast", y="Bias value (arm − concurrent control)"
    ) +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(
      strip.text=ggplot2::element_text(face="bold"),
      strip.placement="outside",
      axis.text.x=ggplot2::element_text(angle=15, hjust=1),
      plot.title=ggplot2::element_text(face="bold"),
      legend.position="bottom",
      panel.spacing.y=grid::unit(16,"pt")
    )
  print(p)
  png_path <- file.path(out_dir, paste0(file_stub, ".png"))
  pdf_path <- file.path(out_dir, paste0(file_stub, ".pdf"))
  grDevices::png(filename=png_path, type="cairo", width=14, height=10, units="in", res=150)
  print(p); grDevices::dev.off()
  grDevices::cairo_pdf(file=pdf_path, width=14, height=10)
  print(p); grDevices::dev.off()
  
  # 3) mean responder counts for Arm A and its concurrent control, by procedure
  if (!all(c("group","pos","neg","neu","procedure") %in% names(groupsA_df))) {
    warning("groupsA_df does not contain responder counts; skipping count summary.")
  } else {
    means <- aggregate(cbind(pos,neg,neu) ~ procedure + group,
                       data = groupsA_df, FUN = mean)
    means$total    <- with(means, pos + neg + neu)
    means$pos_prop <- with(means, ifelse(total > 0, pos / total, NA_real_))
    means$neg_prop <- with(means, ifelse(total > 0, neg / total, NA_real_))
    means$neu_prop <- with(means, ifelse(total > 0, neu / total, NA_real_))
    
    cat("\n=== Mean responder counts per trial (Arm A vs concurrent control), by procedure ===\n")
    print(means[order(means$procedure, means$group), ], row.names = FALSE)
    
    means_path <- file.path(out_dir, paste0(file_stub, "_A_vs_DforA_mean_responder_counts.csv"))
    write.csv(means, means_path, row.names = FALSE)
    message("Saved mean responder counts: ", means_path)
  }
  
  # 4) Weighted average of allocation-bias difference (Arm A − concurrent control)
  #    Weight per run = concurrent-window total = N_A + N_DforA
  alloc_bias_df <- subset(
    data_list$delta,
    arm == "Arm A − concurrent control" &
      grepl("Mean allocation bias", metric, fixed = TRUE)
  )
  
  if (nrow(alloc_bias_df) && all(c("procedure","run_id","value") %in% names(alloc_bias_df))) {
    # Build weights from groupsA_df (two rows per run: Arm A and Concurrent control of Arm A)
    weights_raw <- subset(
      groupsA_df,
      group %in% c("Arm A","Concurrent control of Arm A") &
        !is.na(run_id)
    )
    if (nrow(weights_raw)) {
      # concurrent-window N per (procedure, run_id): sum over A and D_for_A rows
      weights_raw$total <- with(weights_raw, pos + neg + neu)
      weights_by_run <- aggregate(total ~ procedure + run_id, data = weights_raw, FUN = sum)
      names(weights_by_run)[names(weights_by_run) == "total"] <- "weight"
      
      # join weights to deltas
      merged <- merge(alloc_bias_df, weights_by_run, by = c("procedure","run_id"), all.x = TRUE)
      
      # drop runs with missing/zero weight
      merged$weight[is.na(merged$weight)] <- 0
      merged <- subset(merged, weight > 0)
      
      if (nrow(merged)) {
        # weighted mean per procedure
        w_sums <- aggregate(cbind(wsum = value * weight, weight = weight) ~ procedure,
                            data = merged, FUN = sum, na.rm = TRUE)
        w_sums$weighted_mean_diff <- with(w_sums, ifelse(weight > 0, wsum / weight, NA_real_))
        
        out_w <- w_sums[, c("procedure","weighted_mean_diff","weight")]
        names(out_w)[names(out_w) == "weight"] <- "total_weight"
        
        cat("\n=== Weighted allocation-bias difference (Arm A − concurrent control), by procedure ===\n")
        print(out_w[order(out_w$procedure), ], row.names = FALSE)
        
        alloc_w_path <- file.path(out_dir, paste0(file_stub, "_A_vs_DforA_alloc_bias_weighted_summary.csv"))
        write.csv(out_w, alloc_w_path, row.names = FALSE)
        message("Saved weighted allocation-bias summary: ", alloc_w_path)
      } else {
        warning("No rows with positive weights after merge; cannot compute weighted mean difference.")
      }
    } else {
      warning("groupsA_df has no rows for Arm A / concurrent control with run_id; cannot weight.")
    }
  } else {
    warning("No allocation-bias Δ rows for Arm A with run_id found.")
  }
  
  message("Wrote: ", png_path, " and ", pdf_path,
          "\nΔ data: ", csv_path)
}

# --- run ---
set.seed(seed_base)
data_list <- collect_all()
plot_and_save(data_list)
