#!/usr/bin/env Rscript
# ============================================================
#  Plotting: allocation bias (concurrent controls) – POWER
#  Three scenarios × three step schemes
#
#  Input (CSV):
#    PT_bias/results_allocbias_concurrent_threeScenarios/
#      power_vs_alloc_concurrent_threeScenarios_steps.csv
#
#  Output (PDF only):
#    allocbiasPOWER_concurrent_threeScenarios_ALL3_3x3_stacked.pdf
#    allocbiasPOWER_concurrent_threeScenarios_S1_S2_3x3_stacked.pdf
#    allocbiasPOWER_concurrent_threeScenarios_S3_3x3.pdf
#
#  Notes:
#    - y-axis fixed to 0..1
#    - y-label "Power"
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)   # unit()
})

# --- IMPORTANT: source the style sheet (same as your T1E plots) ---
source("PT_bias/plot_style_bimj.R")

# Display titles for the three scenarios (in scenario order)
DISPLAY_STARTS <- c(1, 25, 49)

in_dir  <- "PT_bias/results_allocbias_concurrent_threeScenarios"
out_dir <- "plots_allocbiasPOWER_concurrent_threeScenarios"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

infile <- file.path(in_dir, "power_vs_alloc_concurrent_threeScenarios_steps.csv")
dat_raw <- read.csv(infile, check.names = FALSE)

y_lim <- c(0, 1)

# -----------------------------
# Prepare data (use style helpers)
# -----------------------------
scheme_levels <- c(
  "1-step randomization",
  "2-step randomization\n(total counts only)",
  "2-step randomization\n(total counts by cohort)"
)
scheme_map <- setNames(
  scheme_levels,
  c(
    "1-step randomization",
    "2-step randomization (total counts only)",
    "2-step randomization (total counts by cohort)"
  )
)

scenario_levels <- sort(unique(dat_raw$scenario_b_start))

dat <- dat_raw %>%
  mutate(
    alloc_bias    = as.numeric(alloc_bias),
    
    # factors using style sheet helpers
    series_fct    = bimj_series_factor(series),   # includes "A", "B", "FWER" (disjunctive power)
    rand_key      = bimj_rand_key(rand_label),
    bias_base_fct = bimj_bias_policy_factor(bias_policy_base),
    
    scenario_fct  = factor(scenario_b_start, levels = scenario_levels),
    scheme_fct    = factor(scheme_map[scheme_label], levels = scheme_levels)
  )

y_label_power <- if (exists("BIMJ_Y_POWER")) BIMJ_Y_POWER else "Power"
x_label_ab    <- if (exists("BIMJ_X_ALLOC_LAMBDA")) BIMJ_X_ALLOC_LAMBDA else "Allocation bias"

# -----------------------------
# Plot helper per scenario
# -----------------------------
make_plot_for_scenario <- function(df_scen) {
  scen_idx <- as.integer(df_scen$scenario_fct[1])
  start_disp <- if (!is.na(scen_idx) && scen_idx >= 1 && scen_idx <= length(DISPLAY_STARTS)) {
    DISPLAY_STARTS[scen_idx]
  } else {
    unique(df_scen$scenario_b_start)[1]
  }
  
  ggplot(
    df_scen,
    aes(
      x        = alloc_bias,
      y        = rate,
      color    = rand_key,
      linetype = bias_base_fct,
      group    = interaction(rand_key, bias_base_fct)
    )
  ) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    facet_grid(
      series_fct ~ scheme_fct,
      labeller = labeller(series_fct = label_parsed, scheme_fct = label_value)
    ) +    
    scale_x_continuous(
      breaks = sort(unique(df_scen$alloc_bias)),
      limits = range(df_scen$alloc_bias, na.rm = TRUE)
    ) +
    scale_y_continuous(
      limits = y_lim,
      expand = expansion(mult = c(0, 0))
    ) +
    bimj_scale_color_rand(name = "Randomization procedure", include_ref = FALSE) +
    bimj_scale_linetype_bias_policy(name = "Biasing policy") +
    labs(
      x = x_label_ab,
      y = y_label_power,
      title = paste0("Arm A start at ", start_disp)
    ) +
    bimj_theme(base_size = 13)
}

# -----------------------------
# Stack builder
# -----------------------------
scenarios <- levels(dat$scenario_fct)

build_stack <- function(scenarios_keep) {
  n <- length(scenarios_keep)
  
  plots <- lapply(seq_along(scenarios_keep), function(i) {
    scen <- scenarios_keep[i]
    df_scen <- dat %>% filter(scenario_fct == scen)
    
    make_plot_for_scenario(df_scen) +
      bimj_stack_margin(
        is_top    = (i == 1),
        is_bottom = (i == n)
      )
  })
  
  p <- plots[[1]]
  if (n > 1) for (i in 2:n) p <- p / plots[[i]]
  
  p +
    plot_layout(ncol = 1, guides = "collect") &
    bimj_guides_bias_then_rand(rand_rows = 3) &
    theme(legend.position = "bottom")
}

p_all3 <- build_stack(scenarios)
p_s12  <- build_stack(scenarios[1:2])
p_s3   <- build_stack(scenarios[3])

# -----------------------------
# Save PDFs (same dimensions as your T1E version)
# -----------------------------
bimj_save_pdf(
  p_all3,
  outfile = file.path(out_dir, "allocbiasPOWER_concurrent_threeScenarios_ALL3_3x3_stacked.pdf"),
  width = 12.5, height = 18
)

bimj_save_pdf(
  p_s12,
  outfile = file.path(out_dir, "allocbiasPOWER_concurrent_threeScenarios_S1_S2_3x3_stacked.pdf"),
  width = 12.5, height = 12
)

bimj_save_pdf(
  p_s3,
  outfile = file.path(out_dir, "allocbiasPOWER_concurrent_threeScenarios_S3_3x3.pdf"),
  width = 12.5, height = 6
)

message("Done.")
