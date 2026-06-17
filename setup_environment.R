install.packages(c(
  "snowfall", "zoo", "mvtnorm", "MCMCpack",
  "bayesm", "psych", "coda", "GIGrvg",
  "dlm", "rARPACK", "stochvol", "remotes"
))

# local patched threshtvp
install.packages("Rcpp")
install.packages("RcppArmadillo")
system("R CMD INSTALL --preclean threshtvp_source")
