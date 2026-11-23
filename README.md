# Porous-Medium Tumor Growth UQ with Multi-Fidelity and RFPS Methods

This repository accompanies Huimin Yu, Liu Liu, Yu Feng, and Te Qi, *On multi-fidelity methods for a tumor growth model with uncertainties*. It contains one complete nutrient-related circular experiment with **$m=8$ and $m=64$**, corresponding to manuscript **Test 2(a)/(b)** (original code label `Test1a`) and defense slides **14–15**.

Fine-grid porous-medium equations (PME), a level-set Hele–Shaw approximation (LS), and coarse-grid PME support bi-fidelity (BF), residual-based farthest-point sampling (RFPS), and tri-fidelity (TRF) reconstruction. All parameters and initial fields are generated in code. No input datasets are required.

## Run

Set the MATLAB current folder to this repository root:

```matlab
run_experiment('fresh',8)   % Build bases and evaluate on 500 independent parameters
run_reference(8)            % Mean and standard deviation: 1200 independent parameters

run_experiment('fresh',64)
run_reference(64)
```

`run_reference` reuses a basis with matching settings, or constructs it if needed. It computes both HF reference statistics and BF/RFPS/TRF statistics on the same 1200 parameters. Checkpoints are written every 10 samples; rerunning the same command resumes a matching checkpoint. Statistics are accumulated online, without retaining 1200 full HF fields.

For a short calculation through selection, enrichment, evaluation, statistics, and plotting:

```matlab
opts = struct('T',0.031,'n_train',24,'n_pool',24,'n_eval',8,'n_stats',12);
run_experiment('fresh',8,opts)
run_reference(8,opts)
% Replace 8 with 64 to check the other pressure exponent.
```

These explicit reductions do not change the defaults. Full statistical runs require many PDE solves. MATLAB or GNU Octave can be used; no additional MATLAB toolbox is needed. Octave needs an available graphics toolkit. All figures are displayed and saved by local plotting functions in `run_experiment.m`.

## Model and settings

For five independent inputs $z_k\sim\mathcal U[-1,1]$, define

$$
s=\max\left(1+\sum_{k=1}^{5}\frac{z_k}{2k},10^{-6}\right),\qquad
a_1=\max(1+z_1,10^{-6}),
$$

$$
\lambda=50s,\quad c_B=20s,\quad G_0=0.5a_1,\quad R_0=0.45a_1,\qquad
\rho_0(x,y)=0.95\,\mathbf 1_{\{\sqrt{x^2+y^2}<R_0\}}.
$$

Every solver receives `growth_multiplier = a1`. PME and LS each convert it to the physical coefficient with `G0 = 0.5 * growth_multiplier`; coarse PME calls the same PME solver. The multiplier must not be halved before this call.

| Quantity | Setting |
|---|---|
| Domain / final time | $[-2.5,2.5]^2$ / $T=1$ |
| Fine-grid spacing | $\Delta x=\Delta y=0.05$ |
| HF PME time step | 0.0015 for $m=8$; 0.0005 for $m=64$ |
| LS time step | 0.03 |
| Coarse PME spacing | 0.1 |
| Coarse PME time step | 0.006 for $m=8$; 0.002 for $m=64$ |
| Coarse PME exponent | The same $m$ as the HF model |
| LF training candidates | 500 |
| Supplementary candidate pool | 100 |
| Independent error samples | **500** |
| Independent mean/std samples | **1200**, in `run_reference.m` |
| Initial basis / maximum RFPS additions | 15 / 3 |
| RFPS residual weight / regularization | 0.4 / $10^{-8}$ |
| Spatial RFPS weights | $\exp[-(x^2+y^2)/(2\cdot0.15^2)]$ |
| Random seeds | 20261123 + 0, 1, 2, 3 for training, pool, evaluation, statistics |

The four parameter sets are separate. All methods share the same evaluation and statistical parameters. The final time step is shortened to reach $T=1$ exactly. Both PME levels use the same predictor–corrector discretization with limited conservative fluxes and a semi-implicit growth update. LS uses elliptic nutrient/pressure solves, fast-marching normal-velocity extension, forward Euler with Godunov gradients, and distance reinitialization every 20 steps, following the presentation's Appendix B workflow.

## Reconstruction and error definition

BF selects up to 18 LF solutions by pivoted Cholesky. RFPS shares the first 15 BF parameters and accepts up to three supplementary parameters. TRF selects its initial parameters using coarse PME, forms the paired LS/HF basis, and applies the same enrichment procedure. Each reconstruction transfers LS projection coefficients to the paired HF basis. The supplemental RFPS/TRF comparisons use the same circular test as the BF presentation example.

The primary error follows slide 15, with its evaluation count set to **500**:

$$
E=\frac{1}{N_{\mathrm{eval}}}\sum_{i=1}^{N_{\mathrm{eval}}}
\frac{\|U^H(T,z_i)-\widehat U^H(T,z_i)\|_2}{N_XN_Y}.
$$

Error curves show all available basis sizes. `common_budget` gives the largest size available to all methods for equal-basis-count comparison. Rank limits or rejected enrichment candidates can yield fewer than 18 basis functions; actual counts and selected parameters are saved. Standard deviations use population normalization by $N_s=1200$. Aggregate RMSE and relative matrix errors are also saved as auxiliary metrics.

## Files and outputs

| Path | Purpose |
|---|---|
| `experiment_config.m` | All experiment settings |
| `run_experiment.m` | Generate samples, solve models, build bases, evaluate errors, and plot |
| `run_reference.m` | Independent 1200-sample statistics entry point |
| `solver/model_parameters.m`, `sample_parameters.m`, `solve_sample.m` | Shared input and solver dispatch |
| `solver/solve_pme.m`, `solve_cmpme.m`, `solve_heleshaw.m`, `fast_marching.m` | PDE solvers and interface propagation |
| `solver/select_basis.m`, `enrich_basis.m`, `projection_matrix.m`, `evaluate_basis.m` | Selection, reconstruction, and errors |

The generated `results/` directory contains basis files, error curves, reference statistics, checkpoints, and PNG figures, labeled by `m08` or `m64`. Files record the settings and selected parameters. A new calculation with different settings replaces the same-exponent output files. `results/` is excluded from Git. Keep the repository's nested `solver/` directory when uploading.
