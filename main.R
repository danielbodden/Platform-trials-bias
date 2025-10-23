# ============================================================
# For runs on RWTH High Computing Cluster
# ============================================================
#source("PT_bias/simulation_functions.R")

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[1]) else 120000
set.seed(42)

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
                           alloc_bias = 0) {
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
                                           alloc_bias=alloc_bias),
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
                               show_progress = TRUE) {
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
    alloc_bias=alloc_bias
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


# ============================================================
#  PLOTS: Rejection probability vs timing of start of arm B for
#  (i)  allocation bias only
#  (ii) chronological bias only
# ============================================================

# Sweep helper over B start times for a given setting
sweep_B <- function(b_starts, n_sim=500, concurrent_only=TRUE,
                    beta_time=0, alloc_bias=0,
                    max_group_size=50, expected_total=200,
                    rand_mode="complete", block_factor=1,
                    seed_base=1001) {
  rows <- lapply(seq_along(b_starts), function(i){
    bs <- b_starts[i]
    set.seed(seed_base + i)
    out <- summarize_runs_par(n_sim=n_sim, n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
                              max_group_size=max_group_size,
                          mu=c(A=0,B=0,C=0,D=0),
                          alpha=0.05,
                          arm_start=c(A=0,B=bs,C=0),
                          concurrent_only=concurrent_only,
                          expected_total=expected_total,
                          beta_time=beta_time,
                          rand_mode=rand_mode,
                          block_factor=block_factor,
                          alloc_bias=alloc_bias)
    pr <- out$per_comparison_rejection_rate
    data.frame(B_start=bs,
               arm=c("A","B","C"),
               rate=as.numeric(pr[c("A_vs_D","B_vs_D","C_vs_D")]),
               stringsAsFactors=FALSE)
  })
  do.call(rbind, rows)
}

# Build combined dataset for complete & block
make_plot_data <- function(b_starts=seq(0, 100, by=5), n_sim=500) {
  combos <- expand.grid(
    scenario      = c("Allocation bias only","Chronological bias only"),
    concurrency   = c("Concurrent","Nonconcurrent"),
    rand_mode     = c("complete","block"),
    stringsAsFactors = FALSE
  )
  out_list <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    sc  <- combos$scenario[i]
    cc  <- combos$concurrency[i]
    rm  <- combos$rand_mode[i]
    bt  <- if (sc == "Chronological bias only") 1 else 0.0
    ea  <- if (sc == "Allocation bias only") 0.3 else 0.0
    conc_flag <- (cc == "Concurrent")
    df <- sweep_B(b_starts, n_sim=n_sim, concurrent_only=conc_flag,
                  beta_time=bt, alloc_bias=ea, rand_mode=rm)
    df$scenario    <- sc
    df$concurrency <- cc
    df$rand_mode   <- ifelse(rm=="complete","complete","Block")
    out_list[[i]] <- df
  }
  do.call(rbind, out_list)
}

# Plot and save to file;
plot_and_save <- function(dat, filename, title_sub) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    pdf(paste0(filename, ".pdf"), width=7, height=6)
    on.exit(dev.off(), add=TRUE)
    par(mfrow=c(1,1), mar=c(4,4,3,1))
    plot(NULL,
         xlim = range(dat$B_start),
         ylim = c(0, 1),   # <-- starts at 0
         xlab = "B start (patients)",
         ylab = "Rejection probability",
         main = title_sub)
    cols <- c(A="black", B="blue", C="darkred")
    for (a in c("A","B","C")) {
      df <- dat[dat$arm == a, ]
      lines(df$B_start, df$rate, type="b", pch=19, col=cols[a])
    }
    legend("topleft", legend=c("A","B","C"), col=cols, lty=1, pch=19, bty="n")
  } else {
    library(ggplot2)
    gg <- ggplot(dat, aes(x=B_start, y=rate, color=arm)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(limits = c(0, NA)) +   # <-- starts at 0, auto top
      labs(
        title = "Rejection probability vs timing of addition of arm B",
        subtitle = title_sub,
        x = "B start (patients)",
        y = "Rejection probability",
        color = "Arm"
      ) +
      theme_minimal(base_size = 12)
    ggsave(filename = paste0(filename, ".png"), plot = gg, width=7, height=6, dpi=150)
    ggsave(filename = paste0(filename, ".pdf"), plot = gg, width=7, height=6)
  }
}


# === Build data and write plots ===
set.seed(2025)
b_grid <- seq(0, 100, by=25)     # arm B opens between 0 and 100
plot_dat_all <- make_plot_data(b_starts=b_grid, n_sim=100000)
# === Save combined plot data for reuse ===
write.csv(plot_dat_all, file = "plot_data_all.csv", row.names = FALSE)

# Split and save one file per combination
unique_sc <- unique(plot_dat_all$scenario)
unique_cc <- unique(plot_dat_all$concurrency)
unique_rm <- unique(plot_dat_all$rand_mode)

for (sc in unique_sc) {
  for (cc in unique_cc) {
    for (rm in unique_rm) {
      sub <- subset(plot_dat_all, scenario==sc & concurrency==cc & rand_mode==rm)
      fname <- paste0("rej_vs_B_",
                      gsub(" ","_", tolower(sc)), "_",
                      tolower(cc), "_",
                      tolower(rm))
      subtitle <- paste(sc, "|", cc, "|", rm, "randomization")
      plot_and_save(sub, fname, subtitle)
    }
  }
}
