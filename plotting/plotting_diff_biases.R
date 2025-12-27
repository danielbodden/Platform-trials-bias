#!/usr/bin/env Rscript
# ============================================================
#  Plotting: nonconcurrent chronobias SHAPES (B-only)
#  Layout:
#    Rows    = Measure (Type I error (arm B), Power (arm B))
#    Columns = Chronological bias shape (linear, stepwise, inverted u, seasonal)
#    Curves  = Randomization procedure (color+shape)
#    Linetype= Analysis model (t-test vs ANCOVA)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("PT_bias/plot_style_bimj.R")

parse_int <- function(x) suppressWarnings(as.integer(as.character(x)))

safe_save_pdf <- function(plot, filename, width, height) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  force(plot)
  
  dev_open <- FALSE
  tryCatch({
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(filename, width = width, height = height, onefile = TRUE)
    } else {
      grDevices::pdf(filename, width = width, height = height, onefile = TRUE)
    }
    dev_open <- TRUE
    print(plot)
  }, finally = {
    if (dev_open && grDevices::dev.cur() != 1) grDevices::dev.off()
  })
  
  message("Written: ", normalizePath(filename, mustWork = FALSE))
  invisible(filename)
}

legend_layout_theme <- function() {
  theme(
    legend.box       = "vertical",
    legend.direction = "horizontal",
    legend.box.just  = "left",
    legend.position  = "bottom"
  )
}

# --- row labels you asked for ---
measure_factor <- function(x) {
  x <- as.character(x)
  x2 <- ifelse(x == "Type I error", "Type I error (arm B)",
               ifelse(x == "Power", "Power (arm B)", x))
  factor(x2, levels = c("Type I error (arm B)", "Power (arm B)"))
}

# ------------------------------------------------------------
# Main plotting function
# ------------------------------------------------------------
plot_chronoshapes_nonconc_Bonly <- function(
    t1e_file,
    power_file,
    out_pdf,
    width = 12,
    height = 8
) {
  if (!file.exists(t1e_file))   stop("Missing: ", t1e_file)
  if (!file.exists(power_file)) stop("Missing: ", power_file)
  
  t1e <- read.csv(t1e_file, check.names = FALSE)
  pow <- read.csv(power_file, check.names = FALSE)
  
  dat <- bind_rows(t1e, pow) %>%
    mutate(
      b_start = parse_int(b_start),
      series  = as.character(series),
      measure = as.character(measure)
    ) %>%
    filter(series == "B") %>%  # safety
    mutate(
      measure_fct = measure_factor(measure),
      rand_key = bimj_rand_key(panel),
      analysis_model_fct = bimj_model_factor(analysis_model)
    ) %>%
    filter(!is.na(rand_key), !is.na(analysis_model_fct), !is.na(rate))
  
  # ---- columns order + labels (NO beta in label) ----
  chrono_levels <- c("linear", "stepwise", "inv_u", "seasonal")
  chrono_label <- function(key) {
    switch(as.character(key),
           linear   = "Linear trend",
           stepwise = "Stepwise trend",
           inv_u    = "Inverted-U trend",
           seasonal = "Seasonal trend",
           as.character(key)
    )
  }
  
  dat$chrono_key <- factor(as.character(dat$chrono_key), levels = chrono_levels)
  dat$chronobias <- factor(
    vapply(dat$chrono_key, chrono_label, character(1)),
    levels = vapply(chrono_levels, chrono_label, character(1))
  )
  
  # ---- x scale (1..49) ----
  BGRID <- as.integer(seq(1, 49, by = 12))
  
  # ---- enforce y-limits per row ----
  blank_t1e <- tidyr::expand_grid(
    b_start     = 1L,
    chronobias  = levels(dat$chronobias),
    measure_fct = factor("Type I error (arm B)", levels = levels(dat$measure_fct)),
    rate        = c(0, 0.2)
  )
  blank_pow <- tidyr::expand_grid(
    b_start     = 1L,
    chronobias  = levels(dat$chronobias),
    measure_fct = factor("Power (arm B)", levels = levels(dat$measure_fct)),
    rate        = c(0, 1.0)
  )
  
  alpha_df <- data.frame(
    measure_fct = factor("Type I error (arm B)", levels = levels(dat$measure_fct)),
    yint = 0.025
  )
  
  X_LAB_TB <- expression(paste("Arm B start ", t[B]))
  
  p <- ggplot(
    dat,
    aes(
      x = b_start, y = rate,
      color = rand_key, shape = rand_key,
      linetype = analysis_model_fct,
      group = interaction(rand_key, analysis_model_fct)
    )
  ) +
    geom_blank(data = blank_t1e, aes(x = b_start, y = rate), inherit.aes = FALSE) +
    geom_blank(data = blank_pow, aes(x = b_start, y = rate), inherit.aes = FALSE) +
    geom_hline(
      data = alpha_df,
      aes(yintercept = yint),
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_point(size = 1.9, na.rm = TRUE) +
    facet_grid(measure_fct ~ chronobias, scales = "free_y") +
    scale_x_continuous(
      breaks = BGRID,
      labels = as.character(BGRID),
      limits = c(1, 49),
      expand = expansion(add = c(2, 2))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = NULL,   # remove title
      x = X_LAB_TB,
      y = NULL
    ) +
    bimj_scale_color_rand(include_ref = FALSE) +
    bimj_scale_shape_rand(include_ref = FALSE) +
    bimj_scale_linetype_model() +
    guides(
      color    = guide_legend(order = 1, nrow = 2, byrow = TRUE),
      shape    = guide_legend(order = 1, nrow = 2, byrow = TRUE),
      linetype = guide_legend(order = 2, nrow = 1, byrow = TRUE)
    ) +
    bimj_theme(base_size = 13) +
    legend_layout_theme()
  
  safe_save_pdf(p, out_pdf, width = width, height = height)
  invisible(p)
}

# ------------------------------------------------------------
# Script entry point
# ------------------------------------------------------------
in_dir  <- Sys.getenv("IN_DIR",  unset = "PT_bias/results_nonconcurrent")
out_dir <- Sys.getenv("OUT_DIR", unset = "plots_nonconc_chronoshapes_Bonly")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

t1e_file   <- file.path(in_dir,  "chronoshapes_nonconc_t1e_Bonly_bothModels.csv")
power_file <- file.path(in_dir,  "chronoshapes_nonconc_power_Bonly_bothModels.csv")
out_pdf    <- file.path(out_dir, "chronoshapes_nonconc_Bonly_2row.pdf")

plot_chronoshapes_nonconc_Bonly(
  t1e_file   = t1e_file,
  power_file = power_file,
  out_pdf    = out_pdf,
  width      = 12,
  height     = 8
)

message("Done. Plot written to: ", normalizePath(out_dir, mustWork = FALSE))
