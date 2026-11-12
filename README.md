# Porous-Medium Tumor Growth UQ with Multi-Fidelity and RFPS Methods

This repository provides the nutrient-related uncertainty experiment from:

**Huimin Yu, Liu Liu, Yu Feng, and Te Qi.**  
*On multi-fidelity methods for a tumor growth model with uncertainties.*

The paper develops a hierarchical multi-fidelity framework for uncertainty quantification in porous-medium-equation tumor growth models with moving free boundaries. It combines coarse-grid porous-medium simulations, level-set approximations of the Hele–Shaw limit, and fine-grid porous-medium simulations.

Residual-based farthest-point sampling (RFPS) combines projection-residual information with distance-based sample separation to enrich the approximation space. Bi-fidelity (BF) and tri-fidelity (TRF) reconstructions use this hierarchy to approximate high-fidelity solutions with a limited high-fidelity sampling budget.

## Experiment: Test 2(a) and Test 2(b)

The repository contains the nutrient-related uncertainty experiment in Section 5.2 of the paper, corresponding to the original code entry named `Test1a`.

Five independent random variables, each uniformly distributed on $[-1,1]$, perturb nutrient consumption, boundary nutrient concentration, the growth coefficient, and the initial tumor radius. The paper defines

$$
\lambda=50\left(1+\sum_{k=1}^{5}\frac{z_k}{2k}\right),
\qquad
c_B=20\left(1+\sum_{k=1}^{5}\frac{z_k}{2k}\right),
$$

$$
G_0=0.5(1+z_1),
\qquad
R_0=0.45(1+z_1).
$$

The initial density is 0.95 inside the disk of radius $R_0$ and zero outside. The spatial domain is $[-2.5,2.5]^2$, with final time $T=1$ and pressure-law exponents $m=8$ and $m=64$.

- **Test 2(a):** Compare pivoted-Cholesky BF reconstruction with RFPS-enriched BF reconstruction, using the Hele–Shaw level-set model and the fine-grid porous-medium model.
- **Test 2(b):** Use the coarse-PME–level-set–fine-PME hierarchy for TRF reconstruction.

The quantities of interest include density means, standard deviations, spatial slices, and reconstruction errors at matched high-fidelity budgets.

## Running the Code

In MATLAB, set the current folder to the repository root. The level-set solver requires the Image Processing Toolbox.

```matlab
run_experiment                     % Default experiment: fresh computation for m = 8
run_experiment('fresh',64)          % Run the same experiment for m = 64
run_experiment('fresh_reference',8) % Also recompute all 100 HF reference samples
run_experiment('cached',8)         % Recompute selection and reconstruction from saved snapshots
run_experiment('pde_check',8)       % Run one full-time PDE sample
plot_saved                         % Plot saved BF errors and TRF statistics
```

The `fresh` mode regenerates the low-fidelity candidates, coarse-grid snapshots, and selected high-fidelity basis. It reuses the supplied high-fidelity reference samples, including values used for enrichment. The `fresh_reference` mode also recomputes those reference samples.

Computational inputs are stored in `data/bf_compute_m08.mat` and `data/bf_compute_m64.mat`. Generated outputs are written to `results/` without overwriting the supplied datasets.

## Files and Dependencies

| Path | Purpose |
|---|---|
| `run_experiment.m` | Main BF, RFPS, and TRF computation |
| `plot_saved.m` | Entry point for saved-data plotting |
| `plot_error.m` | Reconstruction-error plots |
| `plot_fields.m` | Statistical fields and spatial slices |
| `solver/solve_pme.m` | Fine-grid porous-medium solver |
| `solver/solve_heleshaw.m` | Hele–Shaw level-set solver |
| `solver/solve_cmpme.m` | Coarse-grid porous-medium solver |
| `solver/select_basis.m` | Initial basis selection |
| `solver/enrich_basis.m` | RFPS basis enrichment |
| `data/bf_compute_m08.mat` | Computational inputs and snapshots for $m=8$ |
| `data/bf_compute_m64.mat` | Computational inputs and snapshots for $m=64$ |
| `data/bf_enrichment_paper.mat` | Saved BF enrichment-error data |
| `data/trifidelity_paper_fields.mat` | Saved TRF plotting data |

Nutrient, pressure, flux, and distance-function calculations are retained within the corresponding solvers.

Keep `solver/` and `data/` in their existing locations and retain the dataset filenames. Entry points resolve paths relative to the repository. The generated `results/` directory is excluded by `.gitignore`.
