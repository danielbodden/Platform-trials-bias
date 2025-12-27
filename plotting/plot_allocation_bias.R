#!/usr/bin/env Rscript
# ============================================================
#  Plotting: allocation bias (concurrent controls)
#  Three scenarios × three step schemes
#
#  Styling/layout aligned with: plot_chrono_bias.R
#   - facet labels parsed (no quotes / no "geq")
#   - robust y-range via coord_cartesian
#   - optional patchwork stacking
#   - color+shape for randomization (combined legend)
#
#  Output (PDF only):
#    allocbias_concurrent_threeScenarios_ALL3_3x3_stacked.pdf
#    allocbias_concurrent_threeScenarios_S1_S2_3x3_stacked.pdf
#    allocbias_concurrent_threeScenarios_S3_3x3.pdf
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) suppressPackageStartupMessages(library(patchwork))

# ---- STYLE SHEET ----
source("PT_bias/plot_style_bimj.R")

# ---- Local readability tweaks (allocation-bias plots only) ----
# Put Permuted block randomization on its own legend row
# by re-ordering the legend entries: CR, BSD, PBR (with 2-row legend).
if (exists("BIMJ_RAND_LEVELS")) {
  BIMJ_RAND_LEVELS <- c("CR", "PBR", "BSD", "PBR0")
}

# Split the FWER facet strip over two lines (parsed plotmath):
#   line 1: P(Reject
#   line 2: ≥ 1 H0)
if (exists("BIMJ_P_REJECT_TEXT")) {
  BIMJ_P_REJECT_TEXT <- "atop('P(Reject', '≥ 1 ' * H[0] * ')')"
}


# Titles you wanted for the three scenarios (in scenario order)
DISPLAY_STARTS <- c(1, 25, 49)

# ---- IO (override via env vars if desired) ----
in_dir  <- Sys.getenv("IN_DIR",  unset = "PT_bias/results_allocbias_concurrent_threeScenarios")
out_dir <- Sys.getenv("OUT_DIR", unset = "plots_allocbias_concurrent_threeScenarios")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

infile <- file.path(in_dir, "t1e_vs_alloc_concurrent_threeScenarios_steps.csv")
dat_raw <- read.csv(infile, check.names = FALSE)

alpha_ref <- 0.025
# Keep the same visible range as the chronobias plots for T1E
ylim_t1e <- c(0, 0.20)

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

# stable scenario order by numeric b_start
scenario_levels <- sort(unique(dat_raw$scenario_b_start))

# build factors
Dat <- dat_raw %>%
  mutate(
    alloc_bias    = as.numeric(alloc_bias),
    rate          = as.numeric(rate),
    
    # facet labels: these are *plotmath strings* in the style sheet
    series_fct    = bimj_series_factor(series),
    
    # randomization mapping (keys used by style scales)
    rand_key      = bimj_rand_key(rand_label),
    
    # bias policy legend labels + linetypes
    bias_base_fct = bimj_bias_policy_factor(bias_policy_base),
    
    scenario_fct  = factor(scenario_b_start, levels = scenario_levels),
    scheme_fct    = factor(scheme_map[scheme_label], levels = scheme_levels)
  )

# ---- sanity checks: fail fast with helpful info ----
if (any(is.na(Dat$rand_key))) {
  bad <- sort(unique(Dat$rand_label[is.na(Dat$rand_key)]))
  stop(
    "Unrecognized rand_label values. Add them to bimj_rand_key() in plot_style_bimj.R.\n",
    "Values: ", paste(bad, collapse = ", ")
  )
}
if (any(is.na(Dat$scheme_fct))) {
  bad <- sort(unique(Dat$scheme_label[is.na(Dat$scheme_fct)]))
  stop(
    "Unrecognized scheme_label values. Update scheme_map in plot_allocation_bias.R.\n",
    "Values: ", paste(bad, collapse = ", ")
  )
}

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
      shape    = rand_key,
      linetype = bias_base_fct,
      group    = interaction(rand_key, bias_base_fct)
    )
  ) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    geom_hline(yintercept = alpha_ref, linetype = "dashed", colour = "grey40") +
    facet_grid(
      series_fct ~ scheme_fct,
      labeller = labeller(series_fct = label_parsed)
    ) +
    scale_x_continuous(
      breaks = sort(unique(df_scen$alloc_bias)),
      limits = range(df_scen$alloc_bias, na.rm = TRUE)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = ylim_t1e) +
    bimj_scale_color_rand(name = "Randomization procedure", include_ref = FALSE) +
    bimj_scale_shape_rand(name = "Randomization procedure", include_ref = FALSE) +
    bimj_scale_linetype_bias_policy(name = "Biasing policy") +
    labs(
      x = BIMJ_X_ALLOC_LAMBDA,
      y = NULL,
      title = paste0("Arm B start at ", start_disp)
    ) +
    bimj_theme(base_size = 13) +
    # make the combined randomization legend compact (CR+BSD on row 1, PBR on row 2)
    guides(
      color = guide_legend(nrow = 3, byrow = TRUE),
      shape = guide_legend(nrow = 3, byrow = TRUE)
    )
}

# -----------------------------
# Stack builder (with spacing between stacked panels)
# -----------------------------
scenarios <- levels(Dat$scenario_fct)

build_stack <- function(scenarios_keep) {
  plots <- lapply(seq_along(scenarios_keep), function(i) {
    scen <- scenarios_keep[i]
    df_scen <- Dat %>% filter(scenario_fct == scen)
    
    p <- make_plot_for_scenario(df_scen)
    
    # spacing between stacked plots (only matters for the patchwork stacks)
    if (has_patchwork) {
      p <- p + bimj_stack_margin(
        is_top    = (i == 1),
        is_bottom = (i == length(scenarios_keep))
      )
    }
    
    p
  })
  
  if (!has_patchwork) {
    stop("patchwork is required for the stacked outputs. Install patchwork or run per-scenario plots.")
  }
  
  p <- plots[[1]]
  if (length(plots) > 1) for (i in 2:length(plots)) p <- p / plots[[i]]
  
  p +
    plot_layout(ncol = 1, guides = "collect") &
    bimj_guides_bias_then_rand(rand_rows = 3) &
    theme(legend.position = "bottom")
}

# -----------------------------
# Build figures
# -----------------------------
p_all3 <- build_stack(scenarios)

if (length(scenarios) >= 2) {
  p_s12 <- build_stack(scenarios[1:2])
} else {
  p_s12 <- NULL
}

if (length(scenarios) >= 3) {
  p_s3 <- build_stack(scenarios[3])
} else {
  p_s3 <- NULL
}

# -----------------------------
# Save PDFs (same dimensions as before)
# -----------------------------
bimj_save_pdf(
  p_all3,
  outfile = file.path(out_dir, "allocbias_concurrent_threeScenarios_ALL3_3x3_stacked.pdf"),
  width = 12.5, height = 18
)

if (!is.null(p_s12)) {
  bimj_save_pdf(
    p_s12,
    outfile = file.path(out_dir, "allocbias_concurrent_threeScenarios_S1_S2_3x3_stacked.pdf"),
    width = 12.5, height = 12
  )
}

if (!is.null(p_s3)) {
  bimj_save_pdf(
    p_s3,
    outfile = file.path(out_dir, "allocbias_concurrent_threeScenarios_S3_3x3.pdf"),
    width = 12.5, height = 6
  )
}

message("Done.")
