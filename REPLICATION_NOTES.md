TVP-GVAR replication

Successful date:
2026-06-17

R version:
4.3.3

Working settings:
CPU=1
saves=50
burns=50
thin=0.1

Generated files:
ttvp_gvar_ssr_sp_level.rda
irf_ttvp_gvar_ssr.rda

Important patches:
- svsample2 -> svsample
- add get_sv_para()
- add get_sv_latent()
- replace sfLapply with lapply
- KF dimension mismatch fixed

Runtime:
about 32 minutes
