#!/usr/bin/env Rscript
# ============================================================
#  Plotting: chronological bias (linear), concurrent controls
#  Three B-start scenarios (1, 25, 49)
#
#  Layout per scenario:
#    - Rows: Type I error (top), Power (bottom)
#    - Cols: Arm A, Arm B, P(Reject ≥ 1 H0)
#    - Curves: Randomization procedures (color + shape)
#
#  Input:
#    PT_bias/results_chronobias_threeScenarios/chronobias_concurrent_threeScenarios_t1e_power.csv
#
#  Output:
#    plots_chronobias_threeBstarts/chronobias_<Bstart>_2x3.pdf
#    plots_chronobias_threeBstarts/chronobias_threeBstarts_2x3_stacked.pdf  (if patchwork available)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) suppressPackageStartupMessages(library(patchwork))

# ---- STYLE SHEET ----
source("PT_bias/plot_style_bimj.R")

# ---------- I/O ----------
infile  <- "PT_bias/results_chronobias_threeScenarios/chronobias_concurrent_threeScenarios_t1e_power.csv"
out_dir <- "plots_chronobias_threeBstarts"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat_raw <- read.csv(infile, check.names = FALSE)

# ---------- factor helpers ----------
measure_factor <- function(x) factor(x, levels = c("Type I error", "Power"))

# plotmath strings so H[0] is guaranteed to render as a subscript
series_factor_parsed <- function(x) {
  recode <- c(
    A    = "Arm~A",
    B    = "Arm~B",
    FWER = "paste('P(Reject ', '\u2265', ' 1 ', H[0], ')')"
  )
  factor(recode[x], levels = recode[c("A", "B", "FWER")])
}

# stable scenario order by numeric b_start
scenario_levels <- sort(unique(dat_raw$scenario_b_start))

dat <- dat_raw %>%
  mutate(
    chrono_strength = as.numeric(chrono_strength),
    measure_fct     = measure_factor(measure),
    series_fct      = series_factor_parsed(series),
    rand_key        = bimj_rand_key(panel),
    scenario_fct    = factor(scenario_b_start, levels = scenario_levels)
  )

# ---------- helper: build 2x3 plot for a single scenario ----------
make_plot_for_scenario <- function(df_scen) {
  b_start <- unique(df_scen$scenario_b_start)
  if (length(b_start) != 1L) b_start <- b_start[1L]
  
  # enforce y-limits via dummy rows:
  dummy_t1e <- expand.grid(
    chrono_strength = NA_real_,
    rate            = c(0, 0.20),
    measure_fct     = factor("Type I error", levels = levels(df_scen$measure_fct)),
    series_fct      = levels(df_scen$series_fct),
    rand_key        = levels(df_scen$rand_key)
  )
  dummy_pow <- expand.grid(
    chrono_strength = NA_real_,
    rate            = c(0, 1.00),
    measure_fct     = factor("Power", levels = levels(df_scen$measure_fct)),
    series_fct      = levels(df_scen$series_fct),
    rand_key        = levels(df_scen$rand_key)
  )
  
  df_plot <- bind_rows(df_scen, dummy_t1e, dummy_pow)
  
  alpha_df <- data.frame(
    measure_fct = factor("Type I error", levels = levels(df_plot$measure_fct)),
    yint        = 0.025
  )
  
  ggplot(
    df_plot,
    aes(
      x     = chrono_strength,
      y     = rate,
      color = rand_key,
      shape = rand_key,
      group = rand_key
    )
  ) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    geom_hline(
      data        = alpha_df,
      aes(yintercept = yint),
      linetype    = "dashed",
      colour      = "grey40",
      inherit.aes = FALSE
    ) +
    facet_grid(
      measure_fct ~ series_fct,
      scales = "free_y",
      labeller = labeller(series_fct = label_parsed)
    ) +
    scale_x_continuous(
      breaks = sort(unique(df_scen$chrono_strength)),
      limits = range(df_scen$chrono_strength, na.rm = TRUE)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    bimj_scale_color_rand(name = "Randomization procedure", include_ref = FALSE) +
    bimj_scale_shape_rand(name = "Randomization procedure", include_ref = FALSE) +
    guides(
      # puts PBR alone on 2nd row (CR + BSD on first row)
      color = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    labs(
      x = expression(paste("Chronological bias strength ", nu)),
      y = NULL,
      title = paste0("Arm B start at ", b_start)
    ) +
    bimj_theme(base_size = 13) +
    theme(
      # override style sheet: legend not stacked vertically
      legend.box       = "horizontal",
      legend.direction = "horizontal"
    )
}

# ---------- build and save per-scenario plots ----------
scenarios <- levels(dat$scenario_fct)
plots <- vector("list", length(scenarios))

for (i in seq_along(scenarios)) {
  scen <- scenarios[i]
  df_scen <- dat %>% filter(scenario_fct == scen)
  
  p_scen <- make_plot_for_scenario(df_scen)
  plots[[i]] <- p_scen
  
  b_start_val <- unique(df_scen$scenario_b_start)[1]
  outfile <- file.path(out_dir, paste0("chronobias_B", b_start_val, "_2x3.pdf"))
  bimj_save_pdf(p_scen, outfile = outfile, width = 11, height = 8)
}

# ---------- combined stacked plot ----------
if (has_patchwork && length(plots) >= 2) {
  p_all <- plots[[1]]
  for (i in 2:length(plots)) p_all <- p_all / plots[[i]]
  
  outfile_all <- file.path(out_dir, "chronobias_threeBstarts_2x3_stacked.pdf")
  bimj_save_pdf(p_all, outfile = outfile_all, width = 11, height = 22)
}

message("Done. Chronobias plots written to: ", normalizePath(out_dir, mustWork = FALSE))
