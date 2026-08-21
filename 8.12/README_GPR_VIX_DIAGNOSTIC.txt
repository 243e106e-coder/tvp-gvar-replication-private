GPR + VIX / P2-STRONG DIAGNOSTIC
================================

Purpose
-------
This is the NEXT diagnostic after the existing P2-current-only runs.

It does NOT try to force GDP negative.

It asks:
  After controlling for global financial risk (VIX), does the response
  of GDP to a structural GPR shock change materially?

Baseline held fixed
-------------------
Domestic lag:
    p = 2

Foreign lag:
    q = 1

Shrinkage:
    B_1 = 8
    B_2 = 0.01
    kappa0 = -0.005

GPR:
    LN_GPR_QMEAN
    +10% normalized GPR impact shock
    non-US direct lagged GPR excluded (current_only)

Trade weights:
    current 2000-2014 weights

Seed:
    20260816

VIX specification
-----------------
Official series:
    FRED VIXCLS

Original source:
    Chicago Board Options Exchange (CBOE)

Raw frequency:
    Daily close

Quarterly processing:
    Arithmetic mean of available daily closes

Model transform:
    log(quarterly mean)

Internal name:
    US_vix

Recursive identification
------------------------
US dominant block:

    GPR -> VIX -> y -> dp -> r -> de -> deq

Interpretation:
    The GPR shock is ordered before VIX, so contemporaneous financial-market
    response through VIX is allowed to be part of the total GPR transmission.

Non-US block
------------
Current foreign/global regressors:
    foreign y
    foreign dp
    foreign r
    foreign de
    foreign deq
    current global GPR
    current global VIX

Lagged foreign/global regressors under current_only:
    lagged foreign y
    lagged foreign dp
    lagged foreign r
    lagged foreign de
    lagged foreign deq
    lagged global VIX

Direct non-US GPR(t-1) is omitted.

Files to upload
---------------
Copy these files into the repository preserving their folders:

    8.12/03_add_vix_to_structural_input.R
    8.12/patch_run_gpr_vix.R
    R/prepare_data_structural_gpr_vix.R
    R/gpr_vix_overrides.R
    .github/workflows/run-p2-strong-gpr-vix-diagnostic.yml

Do not overwrite the current P2 baseline files.

How to run
----------
GitHub -> Actions -> "Run P2 Strong GPR + VIX Diagnostic" -> Run workflow

For the first test use the defaults:
    saves = 500
    burns = 500
    thin = 0.5
    horizon = 12
    shock_pct = 10

What to send back
-----------------
Send the new Actions run URL.

The most important outputs will be:

    results/p2_design_preflight.csv
    results/stability_summary_selected_dates.csv
    results/stability_source_summary.csv
    results/local_stability_summary_by_country.csv
    results/residual_serial_correlation_diagnostic.csv
    results/us_local_coefficient_summary.csv
    results/gdp_sign_diagnostic.csv
    results/gdp_sign_diagnostic_stable_only.csv
    results/gpr_structural_irf_summary.csv
    results/gpr_structural_irf_summary_stable_only.csv
    results/gpr_structural_cumulative_gdp.csv
    results/gpr_structural_cumulative_gdp_stable_only.csv
    results/vix_input_summary.txt
    results/gpr_vix_run_parameters.csv

Decision rule
-------------
Compare against the existing P2-strong current-GPR baseline.

A. Stability similar/better + GDP response becomes smaller/negative:
   VIX/global financial-risk channel is economically important.

B. Stability similar/better + GDP remains positive:
   Do not force it negative. Move to country/crisis heterogeneity and mechanisms.

C. Stability collapses:
   Adding VIX increases dimensional/stability pressure. Do not use this as the
   main specification until shrinkage/stability is re-tuned in a separate grid.

D. GDP sign changes only in unstable/all-draw summaries:
   Ignore that apparent sign change. Stable posterior inference remains primary.

Important
---------
Do NOT:
  - change p, weights, GPR definition, VIX ordering, and priors simultaneously;
  - clip unstable coefficients;
  - discard most posterior draws just to obtain a desired GDP sign;
  - interpret a positive GDP response as automatically wrong.
