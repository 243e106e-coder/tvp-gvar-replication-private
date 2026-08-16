# TVP-GVAR stability update

Replace only these two repository files:

1. `8.12/gpr_structural_irf.R`
2. `8.12/run_gpr_tvp_gvar_structural.R`

Do NOT replace:
- `8.12/patch_BVAR_gpr_lag.R`
- `.github/workflows/run-tvp-gvar_STRUCTURAL_repo_ready.yml`
- `R/BVAR_ttvp.r`
- `8.12/prepare_data_structural.R`

## New diagnostics

The run now calculates, for every historical date and retained posterior draw, the
spectral radius of the global VAR transition/companion matrix.

Hard stability rule:
- stable: rho < 1
- unstable: rho >= 1

Warning only:
- near-unit: 0.98 <= rho < 1

Near-unit draws are NOT deleted. Only rho >= 1 draws are excluded from the
separate stable-only IRF summaries.

New output files include:

- `results/stability_diagnostic.csv`
- `results/stability_summary_by_date.csv`
- `results/stability_summary_selected_dates.csv`
- `results/irf_gpr_structural_stable_only.rda`
- `results/gpr_structural_irf_summary_stable_only.csv`
- `results/gpr_structural_cumulative_gdp_stable_only.csv`
- `results/gdp_sign_diagnostic_stable_only.csv`
- stable-only PNG figures prefixed with `GPR_STRUCT_STABLE_`

The original all-draw outputs are still produced.

## Historical comparison dates

The figure/summary dates are expanded to:

- 2003Q1
- 2008Q3
- 2014Q3
- 2020Q1
- 2022Q1
- 2023Q4

This lets the paper distinguish high-GPR episodes from financial-crisis and
pandemic states without changing the estimation sample.
