# TVP-GVAR with endogenous global GPR

The global GPR series is included as `US_gpr` in the dominant (US) block. It is
therefore endogenous in the stacked global system and weakly exogenous only in
the individual non-US country models. The series is the natural logarithm of
the quarterly GPR index. Responses are generalized impulse responses to a
one-standard-deviation positive GPR innovation.

## Run

```r
install.packages("threshtvp_0.2.tar.gz", repos = NULL, type = "source")
Sys.setenv(TVPGVAR_SAVES=50, TVPGVAR_BURNS=50, TVPGVAR_THIN=0.1)
source("run_gpr_tvp_gvar.R")
```

After the test succeeds, use at least `SAVES=10000`, `BURNS=30000`. The PNG
files and posterior objects are written to `results/`.
