Dominant-unit p=2/q=1 IRF recovery fix

Root cause
----------
The country model is GVAR(p=2,q=1). At lag 2, Theta_2 exists but Wexlag_2 does
not. During global reconstruction the missing foreign/global lag-2 coefficient
block must be represented as a structural zero matrix.

For a 5-variable country:
  W_i rows = 5 own + 7 foreign/global = 12
  lag 1 B  = [Theta_1  Lambda_1] = 5 x 12
  lag 2 B  = [Theta_2  0_(5x7)] = 5 x 12

The failed code built lag 2 as [Theta_2] only, producing 5 columns and the
observed "B/W mismatch for AU lag 2: 5 vs 12".

Upload
------
8.12/dominant_gpr_vix_irf.R
8.12/06_recover_dominant_irf_from_saved_posterior.R
.github/workflows/recover-dominant-irf.yml

Then run:
Actions -> Recover Dominant-Unit IRF Without MCMC

Default source run:
32576318744

This downloads the already-completed posterior from the failed run and resumes
IRF/stability calculations without re-estimating the 15 TVP-BVAR units.
