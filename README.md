# UniMLSeminar

This repository contains a university research project built around the BEAT model from:

Eva Ascarza and Ayelet Israeli, "Eliminating unintended bias in personalized policies using bias-eliminating adapted trees (BEAT)."

The paper PDF and appendix are stored in the repository root. The research proposal and notes are in `main.tex`.

## What is in this repo

- `vendor/beat/`: a vendored copy of the original BEAT R package source imported from `ayeletis/beat`
- `scripts/setup_beat.R`: installs BEAT from the vendored source into a local project library
- `scripts/run_beat_smoke_test.R`: runs a small balanced causal forest example to verify the package loads and trains
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

This installs BEAT and its pinned dependencies into a project-local library at `.r-library/`.

## Run a smoke test

After setup succeeds, run:

```powershell
Rscript scripts/run_beat_smoke_test.R
```

If the toolchain is installed correctly, the script trains a small `balanced_causal_forest` model and prints a prediction summary.

## Where to modify BEAT

For future research extensions, modify the vendored source directly:

- R wrappers: `vendor/beat/R/`
- C++ implementation: `vendor/beat/src/`

After changing BEAT internals, rerun:

```powershell
Rscript scripts/setup_beat.R
```

That reinstalls the local package from your modified source tree.

## Research direction

The current project plan is to explore BEAT extensions such as `C-BEAT` and `BTGQ`, compare them against original BEAT and other baselines, and evaluate metrics such as efficiency, imbalance, and policy delta under increasingly realistic simulated scenarios.
