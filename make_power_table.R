#!/usr/bin/env Rscript

# ============================================================
#  Make Power Table (A_vs_D and B_vs_D separately)
#  - Parallel across scenarios (serial inside workers)
#  - Pretty outputs: CSV + ASCII + Markdown + LaTeX + XLSX (styled)
#  - Column groups: Concurrent, Non-concurrent, Non-concurrent ANOVA
#  - Bias order: Allocation bias -> Chron. Bias -> Both
# ============================================================

suppressWarnings(suppressPackageStartupMessages({
  library(parallel)
  library(data.table)
}))

# ----------------------------
#  User settings (edit here)
# ----------------------------
params <- list(
  sims            = 40000,                                 # simulations per scenario
  outer_cores     = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")),  # workers for scenarios
  simfun          = "simulation_functions.R",             # path to your simulation functions
  out_dir         = "results",                            # output directory
  seed            = 20251031,                             # RNG seed
  max_group       = 24,                                   # per-arm max group size
  expected_total  = 96,                                   # time scale for trend
  arm_start_A     = 1,
  arm_start_B     = 16
)

# ----------------------------
#  Prepare environment
# ----------------------------
set.seed(params$seed)
dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1")

# ----------------------------
#  Fixed experimental setup
# ----------------------------
alpha_one_sided <- 0.025
alternative     <- "greater"
test_side       <- "one.sided"
exp_arms        <- c("A","B")
mu              <- c(A = 0.826, B = 0.826, D = 0)  # delta = 0.56 vs control
arm_start       <- c(A = params$arm_start_A, B = params$arm_start_B)

# Bias components (labels must match exactly)
BIAS_LEVELS <- list(
  `Chron. Bias`     = list(beta_time = 0.16, alloc_bias = 0.00),
  `Allocation bias` = list(beta_time = 0.00, alloc_bias = 0.08),
  `Both`            = list(beta_time = 0.16, alloc_bias = 0.08)
)

# Randomization rows in the output
RAND_ROWS <- list(
  `Block randomization (1*arms)` = list(rand_mode="block",    block_factor=1),
  `Block randomization (2*arms)` = list(rand_mode="block",    block_factor=2),
  `Block randomization (8*arms)` = list(rand_mode="block",    block_factor=8),
  `Complete randomization`       = list(rand_mode="complete", block_factor=1)
)

# Three top-level column groups
GROUPS <- list(
  `Concurrent`             = list(concurrent_only = TRUE,  analysis_model = "ttest"),
  `Non-concurrent`         = list(concurrent_only = FALSE, analysis_model = "ttest"),
  `Non-concurrent ANOVA`   = list(concurrent_only = FALSE, analysis_model = "anova_period")
)

# Subcolumns per bias
comps  <- c("A_vs_D","B_vs_D")
# *** Enforced bias order everywhere ***
biases <- c("Allocation bias","Chron. Bias","Both")

# ----------------------------
#  Build scenario grid
# ----------------------------
scenarios <- list()
for (row_name in names(RAND_ROWS)) {
  rmix <- RAND_ROWS[[row_name]]
  for (gname in names(GROUPS)) {
    gdef <- GROUPS[[gname]]
    for (bc in biases) {                      # enforce order here
      b <- BIAS_LEVELS[[bc]]
      scenarios[[length(scenarios)+1L]] <- list(
        row_name       = row_name,
        rand_mode      = rmix$rand_mode,
        block_factor   = rmix$block_factor,
        group_name     = gname,
        concurrent     = gdef$concurrent_only,
        analysis_model = gdef$analysis_model,
        bias_component = bc,                  # exact label
        beta_time      = b$beta_time,
        alloc_bias     = b$alloc_bias
      )
    }
  }
}

# ----------------------------
#  Worker function
# ----------------------------
run_one <- function(scn, params, mu, arm_start,
                    alpha_one_sided, test_side, alternative, exp_arms) {
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
             OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1")
  if (!file.exists(params$simfun)) {
    stop(sprintf("simulation_functions.R not found at: %s", params$simfun))
  }
  source(params$simfun)
  res <- calc_rejection_summary(
    n_sim           = params$sims,
    n_cores         = 1L,  # serial inside worker to avoid exporting helpers
    max_group_size  = params$max_group,
    mu              = mu,
    alpha           = alpha_one_sided,
    arm_start       = arm_start,
    concurrent_only = scn$concurrent,
    expected_total  = params$expected_total,
    beta_time       = scn$beta_time,
    rand_mode       = scn$rand_mode,
    block_factor    = scn$block_factor,
    alloc_bias      = scn$alloc_bias,
    exp_arms        = exp_arms,
    seed            = params$seed,
    verbose_every   = 0L,
    test_side       = test_side,
    alternative     = alternative,
    bias_policy     = "favor_all_exp",
    analysis_model  = scn$analysis_model
  )
  per <- res$per_comparison_rejection_rate
  data.frame(
    randomization  = scn$row_name,
    group_name     = scn$group_name,
    bias_component = scn$bias_component,
    A_vs_D         = unname(per["A_vs_D"]),
    B_vs_D         = unname(per["B_vs_D"]),
    stringsAsFactors = FALSE
  )
}

# ----------------------------
#  Run scenarios in parallel
# ----------------------------
nworkers <- max(1L, params$outer_cores)
use_parallel <- nworkers > 1L && length(scenarios) > 1L

if (use_parallel) {
  cl <- parallel::makeCluster(nworkers, type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterExport(
    cl,
    varlist = c("params","mu","arm_start","alpha_one_sided",
                "test_side","alternative","exp_arms","run_one"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", BLIS_NUM_THREADS="1"); NULL
  })
  pieces <- parallel::parLapply(cl, scenarios, run_one,
                                params=params, mu=mu, arm_start=arm_start,
                                alpha_one_sided=alpha_one_sided,
                                test_side=test_side, alternative=alternative,
                                exp_arms=exp_arms)
} else {
  pieces <- lapply(scenarios, run_one,
                   params=params, mu=mu, arm_start=arm_start,
                   alpha_one_sided=alpha_one_sided,
                   test_side=test_side, alternative=alternative,
                   exp_arms=exp_arms)
}
tidy_df <- rbindlist(pieces, use.names = TRUE, fill = TRUE)

# ----------------------------
#  Build wide table (data)
# ----------------------------
rows   <- names(RAND_ROWS)
groups <- names(GROUPS)

# Column order: Row + for each group -> for each bias (Allocation, Chron., Both) -> A_vs_D, B_vs_D
col_order <- c(
  "Row",
  as.vector(sapply(groups, function(gg)
    as.vector(sapply(biases, function(bc)
      paste(gg, "|", bc, "|", c("A_vs_D","B_vs_D"))
    )))
  )
)

wide_tbl <- data.table(Row = rows)
for (cn in setdiff(col_order, "Row")) wide_tbl[, (cn) := NA_real_]

for (i in seq_len(nrow(tidy_df))) {
  r   <- tidy_df$randomization[i]
  gg  <- tidy_df$group_name[i]
  bc  <- tidy_df$bias_component[i]
  wide_tbl[Row == r, paste(gg, "|", bc, "|", "A_vs_D") := tidy_df$A_vs_D[i]]
  wide_tbl[Row == r, paste(gg, "|", bc, "|", "B_vs_D") := tidy_df$B_vs_D[i]]
}
wide_tbl <- wide_tbl[, ..col_order]

# ----------------------------
#  Save CSV (numeric) + rounded copy
# ----------------------------
num_tbl <- copy(wide_tbl)
for (j in setdiff(names(num_tbl), "Row")) {
  set(num_tbl, j = j, value = as.numeric(num_tbl[[j]]))
}
csv_path <- file.path(params$out_dir, "power_table_wide_separate_AB.csv")
fwrite(num_tbl, csv_path)

fmt3 <- function(x) ifelse(is.na(x), NA, sprintf("%.3f", x))
disp_tbl <- copy(wide_tbl)
for (j in setdiff(names(disp_tbl), "Row")) {
  set(disp_tbl, j = j, value = fmt3(as.numeric(disp_tbl[[j]])))
}

# ============================================================
#  PRETTY OUTPUTS (ASCII / Markdown / LaTeX)
# ============================================================
parse_key <- function(k) {
  parts <- strsplit(k, "\\|")[[1]]
  parts <- trimws(parts)
  list(group = parts[1], bias = parts[2], comp = parts[3])
}

# ---- ASCII pretty table (generalized to 3 groups) ----
ascii_path <- file.path(params$out_dir, "power_table_wide_separate_AB.txt")
col_width <- function(vec, title) {
  max(nchar(title), max(nchar(ifelse(is.na(vec), "", vec)), na.rm = TRUE))
}
w_Row <- max(28L, nchar("Design"))
w_col <- list()
for (g in groups) for (b in biases) for (c in comps) {
  colname <- paste(g, "|", b, "|", c)
  w_col[[colname]] <- max(7L, nchar(c), max(nchar(disp_tbl[[colname]]), na.rm = TRUE))
}
block_width <- function(g, b) w_col[[paste(g,"|",b,"|","A_vs_D")]] + 3 +
  w_col[[paste(g,"|",b,"|","B_vs_D")]]

group_width <- function(g) {
  # Three bias blocks separated by " | "
  sum(sapply(biases, function(b) block_width(g,b))) + 3*(length(biases)-1)
}

# dynamic top border
hline <- function(char="-") {
  segs <- c(paste0(rep(char, w_Row+2), collapse=""))
  for (g in groups) segs <- c(segs, paste0(rep(char, group_width(g)+2), collapse=""))
  paste0("+", paste(segs, collapse = "+"), "+")
}
center <- function(text, width) {
  pad <- width - nchar(text); left <- floor(pad/2); right <- pad-left
  paste0(strrep(" ", left), text, strrep(" ", right))
}
# Header row 1 (group names)
hdr1 <- paste0(
  "| ", center("Design", w_Row), " | ",
  paste(sapply(groups, function(g) center(g, group_width(g))), collapse = " | "),
  " |"
)
# Header row 2 (bias names per group)
make_bias_line <- function(g) {
  paste(
    sapply(biases, function(b) center(b, block_width(g, b))),
    collapse = " | "
  )
}
hdr2 <- paste0("| ", strrep(" ", w_Row), " | ",
               paste(sapply(groups, make_bias_line), collapse = " | "), " |")

# Header row 3 (A_vs_D | B_vs_D for each bias in each group)
make_comp_line <- function(g) {
  paste(
    sapply(biases, function(b)
      paste0(center("A_vs_D", w_col[[paste(g,"|",b,"|","A_vs_D")]]), " | ",
             center("B_vs_D", w_col[[paste(g,"|",b,"|","B_vs_D")]]))),
    collapse = " | "
  )
}
hdr3 <- paste0("| ", strrep(" ", w_Row), " | ",
               paste(sapply(groups, make_comp_line), collapse = " | "), " |")

make_row <- function(rname) {
  cells <- character(0)
  for (g in groups) for (b in biases) for (c in comps) {
    colname <- paste(g,"|",b,"|",c)
    cells <- c(cells, center(ifelse(is.na(disp_tbl[Row == rname][[colname]]), "",
                                    disp_tbl[Row == rname][[colname]]),
                             w_col[[colname]]))
  }
  body_left  <- center(rname, w_Row)
  body_right <- paste(cells, collapse = " | ")
  paste0("| ", body_left, " | ", body_right, " |")
}
ascii_lines <- c(
  hline("="), hdr1, hline("="), hdr2, hdr3, hline("-"),
  unlist(lapply(rows, make_row)),
  hline("=")
)
writeLines(ascii_lines, ascii_path)

# ---- Markdown (flat header, readable) ----
md_path <- file.path(params$out_dir, "power_table_wide_separate_AB.md")
md_tbl <- copy(disp_tbl); setnames(md_tbl, "Row", "Design")
nice_name <- function(col) {
  if (col == "Design") return(col)
  p <- parse_key(col)
  paste0(p$group, " — ", p$bias, " — ", ifelse(p$comp=="A_vs_D","A vs D","B vs D"))
}
setnames(md_tbl, names(md_tbl), vapply(names(md_tbl), nice_name, ""))
cn <- names(md_tbl)
sep <- paste0("|", paste(rep("---", length(cn)), collapse="|"), "|")
cat(paste0("|", paste(cn, collapse="|"), "|\n", sep, "\n"), file = md_path)
apply(md_tbl, 1, function(r)
  cat(paste0("|", paste(r, collapse="|"), "|\n"), file = md_path, append = TRUE))

# ---- LaTeX (booktabs + dynamic groups) ----
tex_path <- file.path(params$out_dir, "power_table_wide_separate_AB.tex")
to_tex_num <- function(x) ifelse(is.na(x), "", gsub("^NA$", "", sprintf("%.3f", as.numeric(x))))
sink(tex_path)
cat("% Requires: \\usepackage{booktabs} \\usepackage{multirow}\n")
# l + 2 cols per bias * 3 biases per group * (#groups)
ncols <- 1 + length(groups) * length(biases) * length(comps)
cat("\\begin{tabular}{l", paste(rep("c", ncols-1), collapse=""), "}\n", sep="")
cat("\\toprule\n")
# First header line: group multicolumns (each = 6 columns)
cat("\\multirow{3}{*}{Design}")
start_col <- 2
for (g in groups) {
  end_col <- start_col + length(biases)*length(comps) - 1
  cat(" & \\multicolumn{", length(biases)*length(comps), "}{c}{", g, "}", sep="")
  start_col <- end_col + 1
}
cat(" \\\\\n")
# cmidrules per group
start_col <- 2
for (g in groups) {
  end_col <- start_col + length(biases)*length(comps) - 1
  cat("\\cmidrule(lr){", start_col, "-", end_col, "}", sep = "")
  start_col <- end_col + 1
}
cat("\n")
# Second header line: biases (Allocation, Chron., Both)
cat(" ")
for (g in groups) {
  for (b in biases) cat(" & \\multicolumn{2}{c}{", b, "}", sep="")
}
cat(" \\\\\n")
# Third header line: comps
cat(" ")
for (g in groups) for (b in biases) cat(" & A vs D & B vs D")
cat(" \\\\\n\\midrule\n")
# Body
for (r in rows) {
  cat(r)
  for (g in groups) for (b in biases) for (c in comps) {
    colname <- paste(g,"|",b,"|",c)
    cat(" & ", to_tex_num(num_tbl[Row==r][[colname]]), sep="")
  }
  cat(" \\\\\n")
}
cat("\\bottomrule\n\\end{tabular}\n")
sink()

# ============================================================
#  Excel (XLSX) with styling & merged headers (generalized to 3 groups)
# ============================================================
xlsx_path <- file.path(params$out_dir, "power_table_wide_separate_AB.xlsx")

have_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
have_writexl  <- requireNamespace("writexl",  quietly = TRUE)

if (have_openxlsx) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Power (wide)")
  
  # Layout parameters
  n_groups <- length(groups)                   # 3
  cols_per_group <- length(biases)*length(comps)  # 6
  total_data_cols <- n_groups * cols_per_group    # 18
  total_cols <- 1 + total_data_cols               # +1 for Design
  
  # Row 1: group labels merged across their 6 columns
  openxlsx::writeData(wb, "Power (wide)", x = "Design", startRow = 1, startCol = 1)
  for (gi in seq_along(groups)) {
    g <- groups[gi]
    start_col <- 2 + (gi-1)*cols_per_group
    end_col   <- start_col + cols_per_group - 1
    openxlsx::writeData(wb, "Power (wide)", x = g, startRow = 1, startCol = start_col)
    openxlsx::mergeCells(wb, "Power (wide)", rows = 1, cols = start_col:end_col)
  }
  
  # Row 2: bias headers for each group (Allocation bias, Chron. Bias, Both) — write horizontally
  for (gi in seq_along(groups)) {
    start_col <- 2 + (gi-1)*cols_per_group
    openxlsx::writeData(
      wb, "Power (wide)",
      x = matrix(biases, nrow = 1),            # this order controls the header
      startRow = 2, startCol = start_col, colNames = FALSE
    )
    # Merge the two subcolumns per bias, in the same order
    for (bj in 0:(length(biases)-1)) {
      s <- start_col + 2*bj
      openxlsx::mergeCells(wb, "Power (wide)", rows = 2, cols = s:(s+1))
    }
  }
  
  # Row 3: subheaders "A vs D", "B vs D" repeated per bias
  subhdr_row <- matrix(rep(c("A vs D","B vs D"), length(biases)), nrow = 1)
  for (gi in seq_along(groups)) {
    start_col <- 2 + (gi-1)*cols_per_group
    openxlsx::writeData(wb, "Power (wide)", x = subhdr_row,
                        startRow = 3, startCol = start_col, colNames = FALSE)
  }
  
  # Row 4+: data in the same order as headers
  build_block_cols <- function(g) {
    unlist(lapply(biases, function(b)
      paste(g, "|", b, "|", c("A_vs_D","B_vs_D"))), use.names = FALSE)
  }
  ordered_cols <- unlist(lapply(groups, build_block_cols), use.names = FALSE)
  data_block <- as.data.frame(lapply(ordered_cols, function(col) as.numeric(num_tbl[[col]])))
  names(data_block) <- ordered_cols
  wide_xlsx <- cbind(Design = num_tbl$Row, data_block)
  openxlsx::writeData(wb, "Power (wide)", x = wide_xlsx, startRow = 4, startCol = 1, colNames = FALSE)
  
  # Styling
  bold_center <- openxlsx::createStyle(textDecoration = "bold", halign = "center", valign = "center", wrapText = FALSE)
  num_style   <- openxlsx::createStyle(numFmt = "0.000", halign = "center")
  left_style  <- openxlsx::createStyle(halign = "left")
  
  openxlsx::addStyle(wb, "Power (wide)", style = bold_center, rows = 1:3, cols = 1:total_cols, gridExpand = TRUE)
  nrows <- nrow(wide_xlsx)
  openxlsx::addStyle(wb, "Power (wide)", style = left_style, rows = 4:(3+nrows), cols = 1, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Power (wide)", style = num_style,  rows = 4:(3+nrows), cols = 2:total_cols, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Power (wide)", style = openxlsx::createStyle(border="TopBottomLeftRight"),
                     rows = 1:(3+nrows), cols = 1:total_cols, gridExpand = TRUE, stack = TRUE)
  
  openxlsx::freezePane(wb, "Power (wide)", firstActiveRow = 4, firstActiveCol = 2)
  openxlsx::setColWidths(wb, "Power (wide)", cols = 1, widths = 35)
  openxlsx::setColWidths(wb, "Power (wide)", cols = 2:total_cols, widths = 14)
  
  # Tidy sheet
  openxlsx::addWorksheet(wb, "Power (tidy)")
  openxlsx::writeData(wb, "Power (tidy)", tidy_df)
  openxlsx::freezePane(wb, "Power (tidy)", firstActiveRow = 2, firstActiveCol = 2)
  openxlsx::setColWidths(wb, "Power (tidy)", cols = 1:ncol(tidy_df), widths = "auto")
  
  # Notes
  openxlsx::addWorksheet(wb, "Notes")
  notes <- data.frame(
    key   = c("alpha (one-sided)", "alternative", "test_side",
              "delta A", "delta B", "control D",
              "alloc_bias (Allocation bias col)", "beta_time (Chron. Bias col)",
              "sims per scenario", "outer_cores", "max_group", "expected_total",
              "arm_start_A", "arm_start_B", "seed"),
    value = c(alpha_one_sided, alternative, test_side,
              mu[["A"]], mu[["B"]], mu[["D"]],
              "0.08", "0.16",
              params$sims, params$outer_cores, params$max_group, params$expected_total,
              params$arm_start_A, params$arm_start_B, params$seed),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, "Notes", notes)
  openxlsx::setColWidths(wb, "Notes", cols = 1:2, widths = c(28, 22))
  
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  
} else if (have_writexl) {
  # Fallback: flat headers (no merges or styles)
  flat_tbl <- copy(num_tbl); setnames(flat_tbl, "Row", "Design")
  nice_name <- function(col) {
    if (col == "Design") return(col)
    p <- parse_key(col)
    paste0(p$group, " — ", p$bias, " — ", ifelse(p$comp=="A_vs_D","A vs D","B vs D"))
  }
  setnames(flat_tbl, names(flat_tbl), vapply(names(flat_tbl), nice_name, ""))
  writexl::write_xlsx(
    list(
      "Power (wide)" = flat_tbl,
      "Power (tidy)" = tidy_df
    ),
    path = xlsx_path
  )
} else {
  warning("Neither 'openxlsx' nor 'writexl' is installed; skipping XLSX. Install one of them to get an Excel file.")
}

# ----------------------------
#  Console print (ASCII) & path summary
# ----------------------------
cat("\nPower table (A_vs_D and B_vs_D separate) — pretty ASCII:\n")
cat(paste0(ascii_lines, collapse="\n"), "\n")

cat("\nSaved tables:\n")
cat(" CSV      : ", csv_path, "\n", sep="")
cat(" ASCII    : ", ascii_path, "\n", sep="")
cat(" Markdown : ", md_path, "\n", sep="")
cat(" LaTeX    : ", tex_path, "\n", sep="")
cat(" Excel    : ", xlsx_path, "\n\n", sep="")
