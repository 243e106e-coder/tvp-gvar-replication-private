# 14-Economy TVP-GVAR: Structural Global GPR Version

## What changed

This package fixes the identification and data-handling issues in the previous endogenous-GPR/GIRF smoke test.

1. The source variable is **Global GPR** (`GPR` in `data_gpr_export.xls`), not country-specific `GPRC_USA`.
2. GPR processing is deterministic: **monthly GPR -> quarterly arithmetic mean -> natural log**.
3. No automatic log guessing, forward-fill, interpolation, winsorization, or sign reversal is allowed.
4. Macro ordering is standardized to `y, dp, r, de, deq`.
5. The US endogenous block is ordered **GPR -> y -> dp -> r -> de -> deq**.
6. The baseline IRF is a **recursive structural GPR shock**, not a generalized innovation.
7. The baseline shock size is fixed at **+10% GPR on impact** at every historical date, which makes time-varying transmission comparable across dates.
8. GDP output includes both the ordinary growth-rate IRF and the **cumulative GDP IRF**.

## Existing repository files still required

Keep the original model engine files in `R/`:

- `BVAR_ttvp.r`
- `Datahandling.r`
- `auxilliary_functions_tvp.r`

Also keep the `threshtvp_source/` directory used by the GitHub Actions patch/install step.

## New files to copy into the repository

- `R/01_process_gpr.R`
- `R/02_build_14e_input_structural.R`
- `R/prepare_data_structural.R`
- `R/gpr_structural_irf.R`
- `run_gpr_tvp_gvar_structural.R`
- `.github/workflows/run-tvp-gvar-structural.yml`

## Input files

The workflow searches automatically for:

- `data_gpr_export*.xls`
- `TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成*.xlsx`
- `Trade_Weights_14_Economies_2000_2012*.xlsx`

The macro workbook must contain `MODEL_COMPLETE_14C`; the trade workbook must contain `Trade_Weights`.

## First run: structural smoke test

Use the workflow defaults:

- saves = 100
- burns = 100
- thin = 0.5
- horizon = 12
- shock_pct = 10
- gpr_column = `LN_GPR_QMEAN`

This run is for checking identification and signs, not final inference.

Check these files first:

1. `results/gpr_shock_validation.csv` — `US_gpr` impact must equal `log(1.10)`.
2. `results/gdp_sign_diagnostic.csv` — inspect GDP signs at horizons 0, 1, 4, 8, 12.
3. `results/GPR_STRUCT_2003Q1_y.png`, `GPR_STRUCT_2008Q3_y.png`, `GPR_STRUCT_2020Q1_y.png`, etc. — GDP-growth IRFs by historical regime.
4. `results/GPR_STRUCT_2003Q1_CUM_GDP.png`, etc. — cumulative GDP-level implications.

Do **not** force negative GDP responses. If a properly identified structural shock produces positive responses for some countries or dates, treat that as a result. If all 14 countries remain positive at nearly every date, continue diagnosing the model before interpretation.

## Paper run

Only after the smoke test is structurally sensible, increase the MCMC substantially (for example saves 2000+, burns 3000+; assess posterior stability/convergence rather than relying on a magic number).

## Recommended robustness exercises

- Existing GIRF specification as robustness, not baseline.
- Alternative financial ordering: `GPR -> y -> dp -> r -> deq -> de`.
- Set workflow `gpr_column=LN_GPR_QMAX` for the quarterly-maximum robustness.
- Set `gpr_column=LN_GPRT_QMEAN` or `LN_GPRA_QMEAN` to separate threats and acts.
- 13-country sample excluding Türkiye for the longer balanced macro sample.
