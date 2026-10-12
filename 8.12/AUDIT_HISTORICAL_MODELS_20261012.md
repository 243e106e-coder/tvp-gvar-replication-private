# Historical model and provenance audit — 2026-10-12

**Scope:** read-only inspection of committed files and the run metadata for GitHub Actions [37872702188](https://github.com/243e106e-coder/tvp-gvar-replication-private/actions/runs/37872702188). Reference commit: `1ff13aa298148620728452f5f97d33f5b9443a0b`. **Not** an independent reconstruction of all raw data or re-estimation of posterior draws. Claims below are tagged by evidence strength.

## Verified directly from historical workflow and code

| Area | Observed implementation | Consequence |
|---|---|---|
| Workflow | `.github/workflows/compare-independent-shocks.yml`; matrix GPR/VIX/EPU/OIL | Four independent estimations, not a joint four-shock regression |
| Macro source | `8.12/TVP_GVAR_14经济体_5变量_2000Q1_2026Q2_处理完成.xlsx`, `MODEL_COMPLETE_14C` read by `02_build_14e_input_structural.R` | Original macro source workbook, not the later all-levels workbook |
| GDP transformation | `02b_build_14e_input_gdp_loglevel.R`; GDP_DLOG divided by `TVPGVAR_GDP_DLOG_SCALE=100`, cumulated from the second sample period, anchor index=100 | `y` is an *indexed log GDP level*, **not** delta-log GDP |
| Other variables | `*_CPI_DLOG -> dp`, `*_RATE_LEVEL -> r`, `*_REER_DLOG -> de`, `*_EQ_RETURN -> deq` | Mixed-transform system: inflation, rate level, REER changes, equity returns |
| GPR | `8.12/gpr_quarterly_processed.csv`, `LN_GPR_QMEAN` | Global quarterly log-mean GPR, subject to upstream source verification |
| VIX | `8.12/VIXCLS_daily_fred.csv` processed through `03_add_vix_to_structural_input.R` | Daily-origin series; inspect processor for handling of missing days |
| EPU | `8.12/Global_EPU_Quarterly_1997Q1_2026Q3.csv`, `GEPU_current` logged in `compare_four_shocks.R` | Quarterly data and month-count check (if available) |
| Oil | `8.12/IMF_Brent_quarterly_log_2000Q1_2026Q2.csv`, `brent_usd_per_barrel` logged in `compare_four_shocks.R` | Brent price series; verify upstream observation provenance |
| Trade W | `8.12/Trade_Weights_14_Economies_2000_2014.csv`, 14 reporter-by-partner rows and columns | Final matrix numeric/nonnegative/zero-diagonal; does not prove raw bilateral observations unfilled |
| Shock | Scalar dominant `GL_shock`; `solve(G,u)`, normalized to `log1p(shock_pct/100)` | Structural transmission under model assumptions, **not** generalized 1-SD GIRF |
| Posterior | Default `saves=500`, `burns=500`, `thin=0.5`; IRF kept only for draws with spectral radius `<1` | Check retained draws, stable shares, and convergence; workflow success is not sufficient |

## Known documented concern: trade weights

`8.12/02_build_14e_input_structural.R` (historical lines around 177–181, 233–237) says **2013–2014 missing bilateral observations were mirror-completed prior to assembling the row-normalized matrix**. Both the provenance comment and the method string repeat this claim. The final CSV has no observation-by-observation missingness/source flags. It is not possible from the normalized final W alone to establish which country pairs and source-year observations were affected. The builder's check `any(!is.finite(W_raw))` applies to the **already aggregated matrix** and does not establish that the upstream source contained no gaps.

**Status: unresolved lineage; do not falsely claim zero imputation.** Further action: locate original official bilateral trade exports, exact construction script, annual raw matrices, and missing/mirror-completion log; reconcile to the 14x14 W. If no authentic provenance exists, build a replacement W using only documented observed bilateral cells (and explicit coverage limitations); never replace missing with 0, mirror values, or interpolate.

## Documentation errors and experiment distinctions

1. Root `README.md` and `.github/workflows/README.md` formerly characterized the whole repository as `US_gpr` and **generalized 1-SD IRFs**. This applied, at most, to an older workflow and is false as a universal description of the independent four-shock experiment.
2. `8.12/README_RUN.md` documents a **GDP-growth structural GPR** workflow including cumulative GDP IRFs; these instructions must not be applied to 37872702188, whose GDP is already a log level.
3. `8.12/README_GPR_VIX_DIAGNOSTIC.txt` describes a **GPR+VIX joint dominant unit** with ordering `GPR -> VIX`, whereas run 37872702188 has **one independent scalar global risk per run**.
4. The `no missing-data imputation` statement in four-shock `provenance.txt` is too broad when weights have a documented upstream mirror-completion concern. It can be used only for the **final estimation input alignment**, not upstream historical weight lineage.

## What remains unverified

- Original macro values, transformations, quarterly aggregation and revisions against named OECD/BIS/FRED/IMF national data records.
- Raw 2000–2014 partner-year bilateral trade cells, exact mirrored cells, and weight denominator/coverage.
- Numerical consistency of all posterior IRFs and credible intervals with stored draws, convergence/ESS and stable draw fractions.
- Sensitivity of identification to shock ordering and confounding; country-specific IRF sign explanations.
- Formal validity of an all-log-level cointegrated GVAR (not established in this audit).

## Safe next decision

**First** reconstruct the lineage of historical trade W and check all four historical `stability.csv` artifacts before any interpretation or rerun. **Then** verify global-shock scaling, GDP response units, stability filtering and posterior adequacy. Only once these checks pass should we rerun a controlled original-versus-alternative specification. Changes to README are documentation corrections **only**, no empirical result is overwritten.
