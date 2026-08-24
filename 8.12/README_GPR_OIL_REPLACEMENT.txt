GPR + BRENT OIL REPLACEMENT
============================

Purpose
-------
Replace VIX with IMF Brent Oil in the CURRENT GDP log-level dominant-unit TVP-GVAR,
while keeping the rest of the specification unchanged.

Files to add
------------
1) 8.12/13_add_oil_to_structural_input.R
2) .github/workflows/run-gdp-loglevel-dominant-gpr-oil.yml

Oil data already in repository
------------------------------
8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv

Oil transformation
------------------
GL_oil = ln(Brent USD/barrel quarterly price level)

No:
- differencing
- growth rate
- standardization
- deflation
- new aggregation

Dominant block
--------------
[GL_gpr, GL_oil]

Baseline recursive identification
---------------------------------
GPR -> Oil

Clean comparison settings
-------------------------
Use p=2 for this first run.
Keep:
- same q logic
- same trade weights
- same country variables
- same GPR definition
- same GPR +10% normalization
- same six crisis dates
- same MCMC inputs
- same stable-only rule

Why generated templates are used
--------------------------------
The workflow mechanically converts the current validated GPR+VIX mapping/IRF/run
scripts into GPR+Oil versions. This minimizes accidental changes unrelated to the
single intended replacement VIX -> Oil.

Expected posterior
------------------
results/predDens_dominant_gpr_oil.rda

Expected IRF draw file
----------------------
results/irf_dominant_gpr_oil.rda
