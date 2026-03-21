# UniMLSeminar

This repository contains a university research project built around the BEAT model from:

Eva Ascarza and Ayelet Israeli, "Eliminating unintended bias in personalized policies using bias-eliminating adapted trees (BEAT)."

The paper PDF and appendix are stored in the repository root. The research proposal and notes are in `main.tex`.

## What is in this repo

- `vendor/beat/`: a vendored copy of the original BEAT R package source imported from `ayeletis/beat`
- `scripts/setup_beat.R`: installs BEAT from the vendored source into a local project library
- `scripts/reinstall_beat.R`: reinstalls only the vendored BEAT package against the existing local project library
- `scripts/run_beat_smoke_test.R`: runs a small balanced causal forest example to verify the package loads and trains
- `scripts/run_simulations.R`: runs the reusable simulation framework for comparing paper-style scenarios and methods
- `simulations/`: flat simulation files with one file per method and one file per scenario
- `scripts/run_btgq_smoke_test.R`: runs a small BTGQ causal forest example and tunes `lambda` to a quota target
- `scripts/run_btgq_max_smoke_test.R`: runs the BTGQ-Max variant with a dynamic budget ceiling
- `main.tex`: project proposal and extension ideas

Vendoring the package instead of linking to it as a submodule means future changes to BEAT are tracked directly in this repository. That is the safer setup for a research codebase where you expect to modify the upstream implementation.

## Local setup

The upstream BEAT README states that the package only works with `R 4.2.3` or earlier. On Windows you also need the matching R build tools.

1. Install `R 4.2.3` for Windows.
2. Install `Rtools42`.
3. From the repository root, run:

```powershell
Rscript scripts/setup_beat.R
```

This installs BEAT into a project-local library at `.r-library/`. The setup script now reuses already-installed local packages when the required pinned versions are present and only reaches CRAN when something is missing or the pinned version is not available locally.

### Setup troubleshooting

If `scripts/setup_beat.R` fails even though the project worked before, first check whether `.r-library/` already contains the required packages. The setup script prefers those local installs, but it still needs CRAN access when a dependency is missing or when the installed version does not match the pinned version required by BEAT, especially for `RcppEigen` and `RcppArmadillo`.

If you only changed BEAT code and the local library is otherwise intact, reinstall the vendored BEAT package directly against the existing local library instead of rebundling dependencies:

```powershell
Rscript scripts/reinstall_beat.R
```

On Windows, restricted shells can still block the `R CMD build` step used by `install_local`. If that happens, rerun the same command from a normal local PowerShell session with sufficient permissions.

## Run a smoke test

After setup succeeds, run:

```powershell
Rscript scripts/run_beat_smoke_test.R
```

If the toolchain is installed correctly, the script trains a small `balanced_causal_forest` model and prints a prediction summary.

## Run the simulation framework

After setup succeeds, you can run the scenario benchmark scaffold with:

```powershell
Rscript scripts/run_simulations.R
```

The default run executes the four BEAT-paper-style scenarios and compares:

- `CF-FD`
- `CF-NP`
- `Debiased`
- `BEAT`

Outputs are written to `outputs/simulations/`. Add or modify scenarios directly in the `simulations/scenario_*.R` files.

You can also run a single scenario directly, for example:

```powershell
Rscript simulations/scenario_1_high_corr_tau_z.R
```

## Run the BTGQ smoke test

To exercise the BTGQ extension, run:

```powershell
Rscript scripts/run_btgq_smoke_test.R
```

This trains a baseline `btgq_causal_forest`, tunes `lambda` with `tune_btgq_causal_forest`, and prints the achieved quota trace.

## Run the BTGQ-Max smoke test

To exercise the BTGQ-Max extension, run:

```powershell
Rscript scripts/run_btgq_max_smoke_test.R
```

This uses the BTGQ local split rule and tunes `lambda` with `tune_btgq_max_causal_forest`, reporting the dynamic viable-pool budget and the achieved quota.

## Where to modify BEAT

For future research extensions, modify the vendored source directly:

- R wrappers: `vendor/beat/R/`
- C++ implementation: `vendor/beat/src/`

After changing BEAT internals, rerun:

```powershell
Rscript scripts/reinstall_beat.R
```

Use `scripts/setup_beat.R` only when you need to bootstrap or repair the local library itself. For normal BEAT development, `scripts/reinstall_beat.R` reinstalls the local package from your modified source tree without touching unrelated dependencies.

## Research direction

The current project plan is to explore BEAT extensions such as `C-BEAT` and `BTGQ`, compare them against original BEAT and other baselines, and evaluate metrics such as efficiency, imbalance, and policy delta under increasingly realistic simulated scenarios.
