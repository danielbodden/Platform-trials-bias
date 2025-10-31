#load ggplot2
library(ggplot2)
source("simulation_functions.R")


######################
# --- SERIAL: run many trials (no parallel), collect stacked per-run rows
bias_runs_collect <- function(
    n_sim   = 200,
    sim_fun = platform_trials_simulation,   # or platform_trials_simulation_logged
    sim_args = list(exp_arms = c("A","B"), rand_mode = "block", alloc_bias = 0.3, beta_time = 0.4),
    seed = NULL,
    verbose_every = 0L,   # e.g., 50 for a simple progress print; 0 to disable
    continue_on_error = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  out <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    if (verbose_every > 0L && (i %% verbose_every == 0L)) {
      message(sprintf("[bias_runs_collect] %d/%d ...", i, n_sim))
    }
    
    one <- try(
      {
        res <- do.call(sim_fun, sim_args)
        bias_metrics_to_df(res, sim_id = i)
      },
      silent = TRUE
    )
    
    if (inherits(one, "try-error")) {
      if (!continue_on_error) stop(one)
      # record an empty rowset but keep sim_id for traceability
      out[[i]] <- data.frame(sim_id = i, arm = character(0), metric = character(0), value = numeric(0))
    } else {
      out[[i]] <- one
    }
  }
  
  do.call(rbind, out)
}

# --- turn one trial's bias_metrics into tidy rows
bias_metrics_to_df <- function(res, sim_id = NA_integer_) {
  bm <- res$bias_metrics
  if (is.null(bm) || length(bm) == 0L) {
    return(data.frame(sim_id=sim_id, arm=character(0), metric=character(0), value=numeric(0)))
  }
  data.frame(
    sim_id = sim_id,
    arm    = rep(rownames(bm), times = ncol(bm)),
    metric = rep(colnames(bm), each  = nrow(bm)),
    value  = as.vector(bm),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# --- run many trials (parallel) and collect stacked per-run rows
bias_runs_collect_par <- function(
    n_sim = 500,
    sim_fun = platform_trials_simulation,     # or platform_trials_simulation_logged
    sim_args = list(exp_arms = c("A","B"), rand_mode = "block", alloc_bias = 0.3, beta_time = 0.4),
    seed = NULL,
    cores = max(1L, parallel::detectCores() - 1L)
) {
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_sim)
  
  do_one <- function(i) {
    set.seed(seeds[i])
    res <- do.call(sim_fun, sim_args)
    bias_metrics_to_df(res, sim_id = i)
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("sim_fun","sim_args","bias_metrics_to_df","seeds"), envir = environment())
    parts <- parallel::parLapply(cl, seq_len(n_sim), do_one)
  } else {
    parts <- parallel::mclapply(seq_len(n_sim), do_one, mc.cores = cores)
  }
  do.call(rbind, parts)
}

# --- robust aggregation to mean/SD/SE/CI per arm & metric
aggregate_bias_df <- function(df, conf_level = 0.95) {
  if (is.null(df) || !nrow(df)) return(df)
  ag <- aggregate(value ~ arm + metric, df, function(x) {
    x <- x[is.finite(x)]
    c(n=length(x),
      mean=if (length(x)) mean(x) else NA_real_,
      sd=if (length(x) > 1) sd(x) else NA_real_)
  })
  if (is.matrix(ag$value)) {
    tmp <- as.data.frame(ag$value)
    if (is.null(colnames(tmp)) || any(colnames(tmp) == "")) {
      colnames(tmp) <- c("n","mean","sd")
    }
    ag <- cbind(ag[c("arm","metric")], tmp)
  } else {
    pull <- function(v, nm, idx) if (!is.null(names(v)) && nm %in% names(v)) unname(v[[nm]]) else unname(v[[idx]])
    ag$n    <- vapply(ag$value, pull, numeric(1), nm="n",    idx=1)
    ag$mean <- vapply(ag$value, pull, numeric(1), nm="mean", idx=2)
    ag$sd   <- vapply(ag$value, pull, numeric(1), nm="sd",   idx=3)
    ag$value <- NULL
  }
  ag$se <- ifelse(is.finite(ag$sd) & ag$n > 0, ag$sd / sqrt(ag$n), NA_real_)
  alpha <- 1 - conf_level
  ag$ci_mult <- qt(1 - alpha/2, df = pmax(1, ag$n - 1))
  ag$ci_lo <- ag$mean - ag$ci_mult * ag$se
  ag$ci_hi <- ag$mean + ag$ci_mult * ag$se
  ag[order(ag$metric, ag$arm), c("arm","metric","n","mean","sd","se","ci_lo","ci_hi")]
}

# --- PLOTTING: single panel, grouped by metric on x-axis (shared y-axis)
plot_bias_box_singlepanel <- function(df_runs, title = "Per-run bias (one panel, grouped by metric)") {
  if (is.null(df_runs) || !nrow(df_runs)) stop("Empty df_runs. Collect runs first.")
  df <- df_runs
  df$metric <- factor(
    df$metric,
    levels = c("mse_outcome_vs_0.05","mean_allocation_bias","mean_chronological_bias"),
    labels = c("MSE (y−0.05)","Mean alloc. bias","Mean chrono. bias")
  )
  ggplot(df, aes(x = metric, y = value, fill = arm)) +
    geom_boxplot(position = position_dodge2(width = 0.75, preserve = "single"), outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = arm),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                size = 0.7, alpha = 0.25, show.legend = FALSE) +
    labs(title = title, x = NULL, y = "Value", fill = "Arm") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# --- PLOTTING: faceted by metric (shared y-axis)
plot_bias_box_faceted <- function(df_runs, title = "Per-run bias (faceted, shared y-axis)") {
  if (is.null(df_runs) || !nrow(df_runs)) stop("Empty df_runs. Collect runs first.")
  df <- df_runs
  df$metric <- factor(
    df$metric,
    levels = c("mse_outcome_vs_0.05","mean_allocation_bias","mean_chronological_bias"),
    labels = c("MSE (y−0.05)","Mean alloc. bias","Mean chrono. bias")
  )
  ggplot(df, aes(x = arm, y = value, fill = arm)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.7) +
    facet_wrap(~ metric, nrow = 1, scales = "fixed") +
    labs(title = title, x = "Arm", y = "Value") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}


# 1) Collect per-run bias rows (serial)
set.seed(123)
df_runs <- bias_runs_collect(
  n_sim = 200,
  sim_fun = platform_trials_simulation,   # or platform_trials_simulation_logged
  sim_args = list(
    exp_arms   = c("A","B"),
    rand_mode  = "block",
    alloc_bias = 0.3,
    beta_time  = 0.4,
    expected_total = 200
  ),
  verbose_every = 50
)

# 2) Summarize to mean/SD/SE/CI per arm & metric
#summ <- aggregate_bias_df(df_runs, conf_level = 0.95)
#print(summ)

# 3a) Single-panel plot (x = metric, grouped by arm)
#p1 <- plot_bias_box_singlepanel(df_runs, "Bias per run (one panel)")
#print(p1)

# 3b) Faceted plot (one facet per metric, shared y-axis)
#p2 <- plot_bias_box_faceted(df_runs, "Bias per run (faceted)")
#print(p2)






