# Simulation Framework

This folder contains a flat simulation layout for the BEAT comparisons.

## Structure

- `helpers.R`: shared setup, common DGP helpers, and plotting utilities
- `metrics.R`: policy evaluation metrics
- `method_cf_fd.R`
- `method_cf_np.R`
- `method_debiased.R`
- `method_beat.R`
- `scenario_1_high_corr_tau_z.R`
- `scenario_2_low_corr_tau_z.R`
- `scenario_3_low_corr_tau_z_x2_same_sign.R`
- `scenario_4_low_corr_tau_z_x2_opposite_sign.R`

Each scenario file is organized in the same order:

1. Define the data generating process.
2. Train and run each method.
3. Calculate metrics.
4. Save a scenario plot.

Each scenario is also runnable on its own, for example:

```powershell
Rscript simulations/scenario_1_high_corr_tau_z.R
```

When called from `scripts/run_simulations.R`, each scenario returns its method-metrics table and the master script aggregates those tables across scenarios.

## Run

From the repository root:

```powershell
Rscript scripts/run_simulations.R
```

Outputs are written to `outputs/simulations/`.
