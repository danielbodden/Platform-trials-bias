#!/usr/bin/env Rscript
# Master script to recreate all figures from the pre-computed CSV results.
# Run from the repository root:
#   Rscript make_figures.R

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(getwd())
  normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = FALSE)
}

root <- get_script_dir()
setwd(root)

message("Working directory: ", normalizePath(getwd(), winslash = "/", mustWork = FALSE))

run_plot <- function(script, inputs = character()) {
  missing <- inputs[!file.exists(inputs)]
  if (length(missing) > 0L) {
    message("SKIP  ", script, "  (missing input file(s): ",
            paste(missing, collapse = ", "), ")")
    return(invisible(FALSE))
  }

  message("RUN   ", script)
  # Source in a dedicated environment to avoid cross-script object clashes.
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)
  invisible(TRUE)
}

jobs <- list(
  list(
    script = "plotting/plot_allocation_bias.R",
    inputs = c("PT_bias/results_allocbias_concurrent_threeScenarios/t1e_vs_alloc_concurrent_threeScenarios_steps.csv")
  ),
  list(
    script = "plotting/plot_allocation_bias_power.R",
    inputs = c("PT_bias/results_allocbias_concurrent_threeScenarios/power_vs_alloc_concurrent_threeScenarios_steps.csv")
  ),
  list(
    script = "plotting/plot_chrono_bias.R",
    inputs = c("PT_bias/results_chronobias_threeScenarios/chronobias_concurrent_threeScenarios_t1e_power.csv")
  ),
  list(
    script = "plotting/plot_nonconcurrent.R",
    inputs = c(
      "PT_bias/results_nonconcurrent/chronobias_nonconc_t1e_bothModels.csv",
      "PT_bias/results_nonconcurrent/chronobias_nonconc_power_bothModels.csv",
      "PT_bias/results_nonconcurrent/allocbias_nonconc_t1e_preferB_bothModels.csv",
      "PT_bias/results_nonconcurrent/allocbias_nonconc_power_preferB_bothModels.csv"
    )
  ),
  list(
    script = "plotting/plotting_diff_biases.R",
    inputs = c(
      "PT_bias/results_nonconcurrent/chronoshapes_nonconc_t1e_Bonly_bothModels.csv",
      "PT_bias/results_nonconcurrent/chronoshapes_nonconc_power_Bonly_bothModels.csv"
    )
  ),
  list(
    script = "plotting/plot_chronological_bias_power.R",
    inputs = c("PT_bias/results_chronobiasPOWER_nonconc_threeSchemes/power_vs_chron_nonconc_steps.csv")
  )
)

ok <- vapply(jobs, function(j) run_plot(j$script, j$inputs), logical(1))
message("")
message("Done. Plot scripts executed: ", sum(ok), " / ", length(ok))
message("Output folders typically start with: plots_")
