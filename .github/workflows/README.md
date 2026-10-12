# TVP-GVAR research repository — historical specifications

This repository contains **multiple experimental, non-equivalent specifications**, not one validated baseline. Before interpreting any IRF, identify the GitHub Actions run ID, commit SHA, workflow and actual input transforms. See [historical audit](../../8.12/AUDIT_HISTORICAL_MODELS_20261012.md).

## Important: independent four-shock experiment

The successful run [37872702188](https://github.com/243e106e-coder/tvp-gvar-replication-private/actions/runs/37872702188), historical commit `1ff13aa298148620728452f5f97d33f5b9443a0b`, used `.github/workflows/compare-independent-shocks.yml`.

- **Four separately estimated systems**, each containing exactly one global dominant-unit variable `GL_shock`: log quarterly GPR, log quarterly VIX, log GEPU or log quarterly Brent. This is **not** a simultaneous four-risk conditional model.
- GDP `y` was **reconstructed as an indexed log level from GDP_DLOG**, with first in-sample index 100 and scale divisor 100; **not** the original GDP growth rate. CPI `dp`, REER `de`, equities `deq` retained their growth/inflation/return forms; `r` remained a rate level.
- The dominant-unit shock is injected through a scalar variance, propagated via `solve(G,u)`, and **normalized to a +10% change in the corresponding logged global index** (exact log impulse `log(1.10)`). It is **not** the one-SD generalized IRF claimed in some older documentation.
- The experiment used nominal settings `p=1`, `saves=500`, `burns=500`, `thin=0.5` unless manually overridden. Check the actual run inputs and posterior diagnostics before publication.
- IRF summaries were computed **only from posterior draws with global spectral radius below 1**. Report stable fractions and assess selection sensitivity before substantive claims.

## Trade weights — unresolved provenance warning

The four-shock workflow explicitly read `8.12/Trade_Weights_14_Economies_2000_2014.csv`. The associated builder contains a comment and method description saying some **2013–2014 missing bilateral observations had previously been mirror-completed**. The final 14x14 CSV is a normalized matrix: its numbers alone do **not** identify which underlying bilateral observations were missing or mirror-filled. This is an **unresolved audit finding**, not verified pair-level lineage and not evidence of deliberate fabrication. Do not present this matrix as strictly observed/unimputed until original source tables and construction logs are reconciled. No source data have been altered by this documentation update.

## Other variants

- `8.12/README_RUN.md` describes an **earlier structural GPR growth-rate workflow**; its cumulative GDP-growth IRFs do **not** describe run 37872702188, where `y` is already a log level and **must not be cumulated**.
- `8.12/README_GPR_VIX_DIAGNOSTIC.txt` describes a **two-variable GPR+VIX dominant unit**, distinct from the single-`GL_shock` four-run experiment.
- `8.12/README_GPR_OIL_REPLACEMENT.txt` describes the GPR+Oil two-factor replacement experiment, also distinct.
- The new `8.12/TVPGVAR_14C_5VAR_LOG_LEVELS.xlsx` candidate has all five macro series in levels (except the rate, already a level); it is **not yet a validated cointegrated TVP-GVAR**.

See the individual workflow YAML, `provenance.txt` artifacts, variable mappings, and dated audit before referring to a model as “the baseline.” The presence of a successful GitHub Actions check does not establish model validity, stationarity, shock exogeneity, or Bayesian convergence.
