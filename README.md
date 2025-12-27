This repository contains the simulation code and results used in the manuscript:

> **“Multiple Treatment Arms, Multiple Biases? Allocation and chronological biases in rare disease platform trials”**

With this supplement you can:
1) **Re-create all manuscript/supplement figures** from the pre-computed CSV results (**fast**), and  
2) **Re-run the full Monte Carlo simulations** to recrate the CSV results (**slow; intended for HPC**).

The simulations were run in **R 4.2.2** on the **NHR4CES** high-performance computing cluster.

**No real patient data are used.** All results are produced via Monte Carlo simulation.

---

## Quick start (recreate figures from the provided CSV results)

From the repository root, run:

```bash
Rscript make_figures.R
```
This will call all plotting scripts and write PDFs into the `plots_*` folders (see “Figure map”).

If you prefer to run a single plot script, you can also do:

```bash
Rscript plotting/plot_allocation_bias.R
```

---

## Software requirements

- R **>= 4.2.2** (tested with 4.2.2)
- R packages used across the scripts:
  - `dplyr`, `ggplot2`, `tidyr`, `patchwork`

---

## Repository structure

### Core simulation library
- `PT_bias/simulation_functions.R`  
  Main library of functions used by all simulation runner scripts (trial simulation, bias policies, analysis models, summaries).

### Simulation runner scripts (generate CSV results)
These scripts are computationally intensive and were run on HPC. They write results (CSV) that are later used for plotting.

- `1_allocation_bias.R`  
  Allocation bias with **concurrent controls** (type I error).  
  **Output folder (default):** `PT_bias/results_allocbias_concurrent_threeScenarios/`

- `2_chronological_bias.R`  
  Chronological bias (linear time trend) with **concurrent controls** (type I error + power).  
  **Output folder (default):** `PT_bias/results_chronobias_threeScenarios/`

- `3_nonconc_generate.R`  
  Analyses **enriched by non-concurrent controls** (allocation bias + chronological bias).  
  **Output folder (default):** `PT_bias/results_nonconcurrent/`

- `4_chronoBstart_generate_Bonly.R`  
  Additional analyses for alternative time-trend “shapes” (supplementary).  
  **Output folder (default):** `PT_bias/results_nonconcurrent/` (or `OUT_DIR`)

### Plotting scripts (turn CSV results into publication-ready PDFs)
All plotting scripts rely on the shared style helper:

- `PT_bias/plot_style_bimj.R`  **(required)**

Plot scripts:

- `plotting/plot_allocation_bias.R`  
- `plotting/plot_allocation_bias_power.R`  
- `plotting/plot_chrono_bias.R`  
- `plotting/plot_nonconcurrent.R`  
- `plotting/plotting_diff_biases.R`  
---


## Figure map (which script creates which figure)

| Paper figure | Plot script | Required CSV inputs | Output folder (PDFs) |
| Main Fig. 3 (allocation bias, concurrent controls) | `plot_allocation_bias.R` | `PT_bias/results_allocbias_concurrent_threeScenarios/t1e_vs_alloc_concurrent_threeScenarios_steps.csv` | `plots_allocbias_concurrent_threeScenarios/` |
| Main Fig. 4 (chronological bias, concurrent controls) | `plot_chrono_bias.R` | `PT_bias/results_chronobias_threeScenarios/chronobias_concurrent_threeScenarios_t1e_power.csv` | `plots_chronobias_threeBstarts/` |
| Main Fig. 5 (chronological bias, enriched with non-concurrent controls) | `plot_nonconcurrent.R` | `PT_bias/results_nonconcurrent/chronobias_nonconc_t1e_bothModels.csv` and `PT_bias/results_nonconcurrent/chronobias_nonconc_power_bothModels.csv` | `plots_nonconc_both/` |
| Main Fig. 6 (allocation bias, enriched with non-concurrent controls) | `plot_nonconcurrent.R` | `PT_bias/results_nonconcurrent/allocbias_nonconc_t1e_preferB_bothModels.csv` and `PT_bias/results_nonconcurrent/allocbias_nonconc_power_preferB_bothModels.csv` | `plots_nonconc_both/` |
| Supplement: alternative time-trend shapes (B-only) | `plotting_diff_biases.R` | `PT_bias/results_nonconcurrent/chronoshapes_nonconc_t1e_Bonly_bothModels.csv` and `PT_bias/results_nonconcurrent/chronoshapes_nonconc_power_Bonly_bothModels.csv` | `plots_nonconc_chronoshapes_Bonly/` |
| Supplement (allocation bias power, concurrent controls) | `plot_allocation_bias_power.R` | `PT_bias/results_allocbias_concurrent_threeScenarios/power_vs_alloc_concurrent_threeScenarios_steps.csv` | `plots_allocbiasPOWER_concurrent_threeScenarios/` |

---

## Re-running the simulations (computationally intensive)

All simulations are started via:

```bash
Rscript <script>.R
```

### Environment variables (simulation runner scripts)
- `OUT_DIR`  
  Override the default output directory for CSV results.
- `N_SIM_PER_POINT`  
  Monte Carlo replicates per grid point.
- `SLURM_CPUS_PER_TASK`  
  Used as the default number of cores for parallel execution on HPC (falls back to 1).
- `SEED_BASE`  
  Base seed.

Example:

```bash
N_SIM_PER_POINT=50 SEED_BASE=1 OUT_DIR=PT_bias/results_test Rscript 1_allocation_bias.R
```

---
## Contact
Daniel Bodden  
Email: daniel.bodden@rwth-aachen.de

