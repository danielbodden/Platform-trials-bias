# ============================================================
# For runs on RWTH High Computing Cluster
# ============================================================
#source("PT_bias/simulation_functions.R")
source("simulation_functions.R")

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[1]) else 120000
set.seed(42)


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
    out <- calc_rejection_summary(
      n_sim = n_sim,
      n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      max_group_size = max_group_size,
      mu = c(A=0,B=0,D=0),
      alpha = 0.05,
      arm_start = c(A=0,B=bs,C=0),
      concurrent_only = concurrent_only,
      expected_total = expected_total,
      beta_time = beta_time,
      rand_mode = rand_mode,
      block_factor = block_factor,
      alloc_bias = alloc_bias,
      exp_arms = c("A","B")            # <<< new
      
    )
    pr <- out$per_comparison_rejection_rate
    fwer_val <- out$fwer   # <-- add this line
    
    want <- intersect(c("A_vs_D","B_vs_D"), names(pr))
    df <- data.frame(
      B_start = bs,
      arm     = sub("_vs_.*$", "", want),
      rate    = as.numeric(pr[want]),
      stringsAsFactors = FALSE
    )
    
    # enforce order and convert
    lvl <- intersect(c("A","B"), unique(df$arm))
    df$arm <- factor(df$arm, levels = lvl)
    df <- df[order(df$arm), ]
    df$arm <- as.character(df$arm)
    
    # add FWER row
    df_fwer <- data.frame(
      B_start = bs,
      arm     = "FWER",
      rate    = fwer_val,
      stringsAsFactors = FALSE
    )
    
    # combine
    rbind(df, df_fwer)
    

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
  library(ggplot2)
  
  gg <- ggplot(dat, aes(x = B_start, y = rate, color = arm)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    # red dashed horizontal line at alpha = 0.05
    geom_hline(yintercept = 0.05, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(
      values = c("A" = "#1b9e77", "B" = "#d95f02", "FWER" = "#7570b3"),
      name = "Arm"
    ) +
    labs(
      title = "Rejection probability and FWER vs timing of arm B",
      subtitle = title_sub,
      x = "B start (patients)",
      y = "Rejection probability"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  # save both PNG and PDF
  ggsave(filename = paste0(filename, ".png"), plot = gg, width = 7, height = 6, dpi = 150)
  ggsave(filename = paste0(filename, ".pdf"), plot = gg, width = 7, height = 6)
}



# === Build data and write plots ===
set.seed(2025)
b_grid <- seq(0, 100, by=25)     # arm B opens between 0 and 100
plot_dat_all <- make_plot_data(b_starts=b_grid, n_sim=50)
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




################# Power plots:

# --- Sweep helper for power under H1 (Δ = 0.56 for A and B)
sweep_B_power <- function(b_starts,
                          n_sim = 200,
                          n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1")),
                          max_group_size = 50,
                          expected_total = 200,
                          concurrent_only = TRUE,
                          rand_mode = "block",
                          block_factor = 1,
                          alloc_bias = 0,     # η
                          beta_time = 0,      # β_time
                          delta_A = 0.56,
                          delta_B = 0.56) {
  
  rows <- lapply(seq_along(b_starts), function(i) {
    bs <- b_starts[i]
    out <- calc_rejection_summary(
      n_sim = n_sim,
      n_cores = n_cores,
      max_group_size = max_group_size,
      mu = c(A = delta_A, B = delta_B, D = 0),   # H1: both arms effective by Δ
      alpha = 0.05,
      arm_start = c(A = 0, B = bs),
      concurrent_only = concurrent_only,
      expected_total = expected_total,
      beta_time = beta_time,        # scenario toggle
      rand_mode = rand_mode,
      block_factor = block_factor,
      alloc_bias = alloc_bias,      # scenario toggle
      exp_arms = c("A","B")
    )
    
    pr <- out$per_comparison_rejection_rate  # per-arm power
    any_power <- out$fwer                    # familywise power (ANY)
    
    want <- intersect(c("A_vs_D","B_vs_D"), names(pr))
    df <- data.frame(
      B_start = bs,
      arm     = sub("_vs_.*$", "", want),
      rate    = as.numeric(pr[want]),
      metric  = "per-arm",
      stringsAsFactors = FALSE
    )
    lvl <- intersect(c("A","B"), unique(df$arm))
    df$arm <- factor(df$arm, levels = lvl)
    df <- df[order(df$arm), ]
    df$arm <- as.character(df$arm)
    
    df_any <- data.frame(
      B_start = bs,
      arm     = "ANY",
      rate    = any_power,
      metric  = "any",
      stringsAsFactors = FALSE
    )
    
    rbind(df, df_any)
  })
  
  do.call(rbind, rows)
}


plot_power_all <- function(dat, filename, title_sub, ref_y = 0.80) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 needed")
  library(ggplot2)
  dat <- subset(dat, arm %in% c("A","B","ANY"))
  
  gg <- ggplot(dat, aes(x = B_start, y = rate, color = arm)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    geom_hline(yintercept = ref_y, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = c("A" = "#1b9e77", "B" = "#d95f02", "ANY" = "#7570b3"),
                       name = "Arm") +
    labs(
      title = "Per-arm power and familywise power vs timing of arm B",
      subtitle = title_sub,
      x = "B start (patients)",
      y = "Probability"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(paste0(filename, ".png"), plot = gg, width = 7, height = 6, dpi = 150)
  ggsave(paste0(filename, ".pdf"), plot = gg, width = 7, height = 6)
}
set.seed(2025)
b_grid <- seq(0, 100, by = 25)
n_sim  <- 100000  # or smaller for quick runs

scenarios <- list(
  list(name = "Chronological bias only", beta_time = 1,   alloc_bias = 0),
  list(name = "Allocation bias only",    beta_time = 0,   alloc_bias = 0.3)
)

modes <- c("complete", "block")

for (sc in scenarios) {
  for (rm in modes) {
    dat <- sweep_B_power(
      b_starts = b_grid,
      n_sim = 2000,
      n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1")),
      rand_mode = rm,
      concurrent_only = TRUE,     # change if you also want nonconcurrent runs
      beta_time = sc$beta_time,
      alloc_bias = sc$alloc_bias,
      delta_A = 0.56, delta_B = 0.56
    )
    
    fname <- paste0("power_AB_ANY_",
                    ifelse(rm=="complete","complete","block"), "_",
                    gsub(" ", "_", tolower(sc$name)))
    subtitle <- paste0(sc$name, " | ", ifelse(rm=="complete","Complete","Block"), " randomization")
    
    # Optional: save the data used for the plots
    write.csv(dat, paste0(fname, ".csv"), row.names = FALSE)
    
    plot_power_all(dat, filename = fname, title_sub = subtitle, ref_y = 0.80)
  }
}



############################
# What happens with time trend under ANOVA?

out <- calc_rejection_summary(
  n_sim = 2000,
  exp_arms = c("A","B"),
  rand_mode = "block",
  test_side = "one.sided",
  alternative = "greater",
  bias_policy = "favor_B",
  analysis_model = "anova_period",
  alloc_bias =0,
  beta_time=3,
  alpha=0.025,
  concurrent_only=FALSE
)

print(out)

out <- calc_rejection_summary(
  n_sim = 2000,
  exp_arms = c("A","B"),
  rand_mode = "block",
  test_side = "one.sided",
  alternative = "greater",
  bias_policy = "favor_B",
  analysis_model = "ttest",
  alloc_bias =0,
  beta_time=3,
  alpha=0.025,
  concurrent_only=FALSE
)

print(out)

