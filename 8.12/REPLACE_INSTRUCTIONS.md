# GPR current-only update: replacement map

## Replace these existing files

1. `8.12/gpr_structural_irf.R`
   - Replace with the new `8.12/gpr_structural_irf.R`.
   - It supports both 5-column non-US lag blocks (`current_only`) and the old 6-column block (`current_and_lag`).

2. `8.12/run_gpr_tvp_gvar_structural.R`
   - Replace with the new `8.12/run_gpr_tvp_gvar_structural.R`.
   - It sources `R/BVAR_ttvp_gprlag.r`, reads `TVPGVAR_GPR_LAG_MODE`, and writes a regressor-structure diagnostic.

3. `.github/workflows/run-tvp-gvar_STRUCTURAL_repo_ready.yml`
   - Replace with the new workflow file.
   - GitHub Actions gains a `gpr_lag_mode` selector.

## Add this new file

4. `8.12/patch_BVAR_gpr_lag.R`
   - NEW file. Do not replace `R/BVAR_ttvp.r`.
   - During each GitHub Actions run it reads the original `R/BVAR_ttvp.r` and creates `R/BVAR_ttvp_gprlag.r`.
   - Thus the original replication estimator is preserved unchanged.

## Do NOT replace

- `R/BVAR_ttvp.r` (leave original untouched)
- `8.12/prepare_data_structural.R`
- `8.12/02_build_14e_input_structural.R`
- `8.12/gpr_quarterly_processed.csv`
- macro workbook
- trade-weight workbook

## First test run

Use:
- saves = 100
- burns = 100
- thin = 0.5
- horizon = 12
- shock_pct = 10
- gpr_column = LN_GPR_QMEAN
- gpr_lag_mode = current_only

Expected `results/gpr_lag_spec_diagnostic.csv` under `current_only`:
- every country: `Wexlag1 = 5`
- non-US: `Wex = 6`, `Ylag1 = 5`
- US: `Wex = 5`, `Ylag1 = 6`

For robustness, rerun exactly the same workflow with:
- gpr_lag_mode = current_and_lag

Expected under `current_and_lag`:
- non-US: `Wexlag1 = 6`
- US: `Wexlag1 = 5`
