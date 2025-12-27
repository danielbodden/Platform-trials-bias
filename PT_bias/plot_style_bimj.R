# ============================================================
#  Biometrical Journal plotting style (ggplot2)
#  Save as: PT_bias/plot_style_bimj.R
#
#  Goals:
#   - Consistent colors + labels for randomization procedures:
#       Complete randomization
#       Big stick design(3)
#       Permuted block randomization (a_j)
#   - Bigger legend fonts
#   - P(Reject >= 1 H0) formatted nicely (no "geq" text)
#   - Allocation x-axis: Allocation bias strength λ
#   - More vertical spacing between facet rows
#   - Legend stacking: Biasing policy above Randomization procedure
#   - PDF only (no TIFF)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)   # unit()
})

# ---------- Global strings / expressions ----------

# Plain-text fallback that never prints "geq"
BIMJ_P_REJECT_TEXT <- "'P(Reject ≥ 1 ' * H[0] * ')'"
# Plotmath expression version (for legends/axis titles if you want)
BIMJ_P_REJECT_EXPR <- expression(paste("P(Reject ", "\u2265", " 1 ", H[0], ")"))

BIMJ_X_ALLOC_LAMBDA <- expression(paste("Allocation bias strength ", lambda))
BIMJ_Y_T1E          <- "Type I error"
BIMJ_Y_RATE         <- "Rate"

# ---------- Randomization procedure: keys / colors / labels ----------

# Internal keys (use these in data and aesthetics)
BIMJ_RAND_LEVELS <- c("CR", "BSD", "PBR", "PBR0")  # PBR0 optional reference

# Keep these consistent across *all* plots
BIMJ_RAND_COLORS <- c(
  CR   = "#009E73",  # Complete randomization
  BSD  = "#D55E00",  # Big stick design
  PBR  = "#0072B2",  # Permuted block randomization
  PBR0 = "black"     # optional reference curve
)

# Optional, but nice when you map shape as well
BIMJ_RAND_SHAPES <- c(
  CR   = 16,
  BSD  = 17,
  PBR  = 15,
  PBR0 = 4
)

# Legend labels (expressions allowed)
BIMJ_RAND_LABELS <- list(
  CR  = expression("Complete randomization"),
  BSD = expression("Big stick design(3)"),
  PBR = expression("Permuted block randomization (" * a[j] * ")"),
  PBR0 = expression("PBR (" * eta * " = 0) \u2014 reference")
)

# Map your various raw strings to the keys above
bimj_rand_key <- function(x) {
  out <- rep(NA_character_, length(x))
  
  out[x %in% c("Complete randomization", "CR")] <- "CR"
  
  out[x %in% c(
    "Big stick design (a = 3)",
    "Big stick design(a = 3)",
    "Big stick (a = 3)",
    "BSD"
  )] <- "BSD"
  
  out[x %in% c(
    "Permuted block randomization (factor = 1)",
    "PBR (factor = 1)",
    "Permuted block randomization",
    "PBR"
  )] <- "PBR"
  
  out[x %in% c(
    "PBR (eta=0) — reference",
    "PBR (η = 0) — reference",
    "PBR (eta=0)",
    "PBR0"
  )] <- "PBR0"
  
  factor(out, levels = BIMJ_RAND_LEVELS)
}

# ---------- Common factor helpers ----------

bimj_bias_policy_factor <- function(x) {
  factor(
    x,
    levels = c("favor_all_exp", "favor_B"),
    labels = c("Prefer all experimental arms", "Prefer arm B")
  )
}

# Facet labels for A/B/FWER that do NOT risk "geq" showing up.
# (Use plain text; it’s robust across devices/fonts.)
bimj_series_factor <- function(x) {
  factor(
    x,
    levels = c("A", "B", "FWER"),
    labels = c("'Arm A'", "'Arm B'", BIMJ_P_REJECT_TEXT)
  )
}

# ---------- Scales ----------

bimj_scale_color_rand <- function(name = "Randomization procedure", include_ref = FALSE) {
  lvls <- if (include_ref) BIMJ_RAND_LEVELS else BIMJ_RAND_LEVELS[BIMJ_RAND_LEVELS != "PBR0"]
  
  scale_color_manual(
    name   = name,
    values = BIMJ_RAND_COLORS,
    breaks = lvls,
    labels = BIMJ_RAND_LABELS
  )
}

bimj_scale_shape_rand <- function(name = "Randomization procedure", include_ref = FALSE) {
  lvls <- if (include_ref) BIMJ_RAND_LEVELS else BIMJ_RAND_LEVELS[BIMJ_RAND_LEVELS != "PBR0"]
  
  scale_shape_manual(
    name   = name,
    values = BIMJ_RAND_SHAPES,
    breaks = lvls,
    labels = BIMJ_RAND_LABELS
  )
}

bimj_scale_linetype_bias_policy <- function(name = "Biasing policy") {
  scale_linetype_manual(
    name   = name,
    values = c(
      "Prefer all experimental arms" = "solid",
      "Prefer arm B"                 = "22"
    )
  )
}

# ---------- Analysis model (linetype) ----------

bimj_model_factor <- function(x) {
  factor(
    as.character(x),
    levels = c("ttest", "lm_time"),
    labels = c("t-test", "Linear model")
  )
}

bimj_scale_linetype_model <- function(name = "Analysis model") {
  scale_linetype_manual(
    name   = name,
    values = c("t-test" = "solid", "Linear model" = "22")
  )
}



# ---------- Theme ----------

bimj_theme <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      # Facet spacing (fixes row labels / axis tick overlap in tight grids)
      panel.spacing.y = unit(1.05, "lines"),
      panel.spacing.x = unit(0.90, "lines"),
      
      strip.text = element_text(face = "bold", size = base_size),
      axis.text  = element_text(size = base_size),
      axis.title = element_text(size = base_size + 1),
      
      # Legend sizing (bigger like you requested)
      legend.title      = element_text(size = base_size + 1),
      legend.text       = element_text(size = base_size),
      legend.key.width  = unit(2.8, "cm"),
      legend.key.height = unit(0.9, "cm"),
      
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.direction  = "vertical",
      legend.box.just   = "left",
      
      plot.title        = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor  = element_blank()
    )
}

# Ensures legend order: bias policy first, randomization second (stacked)
bimj_guides_bias_then_rand <- function(rand_rows = 3) {
  guides(
    linetype = guide_legend(order = 1, nrow = 1),
    color    = guide_legend(order = 2, nrow = rand_rows, byrow = TRUE),
    shape    = guide_legend(order = 2, nrow = rand_rows, byrow = TRUE)
  )
}

# ---------- Stacking helper (patchwork-safe) ----------
# Add a bit of extra margin between stacked plots so y-axis ticks don’t collide.
bimj_stack_margin <- function(is_top = FALSE, is_bottom = FALSE,
                              top_pt = 10, bottom_pt = 22, left_pt = 8, right_pt = 8) {
  t <- if (is_top) top_pt else 5
  b <- if (is_bottom) 5 else bottom_pt
  theme(plot.margin = margin(t = t, r = right_pt, b = b, l = left_pt, unit = "pt"))
}

# ---------- PDF writer (no TIFF) ----------
bimj_save_pdf <- function(p, outfile, width, height) {
  grDevices::cairo_pdf(outfile, width = width, height = height)
  print(p)
  grDevices::dev.off()
  message("Written: ", normalizePath(outfile, mustWork = FALSE))
}
