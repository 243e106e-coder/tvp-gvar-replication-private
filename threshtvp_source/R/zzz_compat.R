get_sv_para <- function(z) {
  if (inherits(z, "svdraws")) {
    out <- as.numeric(tail(as.matrix(z$para[[1]]), 1)[1, 1:3])
    names(out) <- c("mu", "phi", "sigma")
    return(out)
  } else {
    out <- as.numeric(z$para[1:3])
    names(out) <- c("mu", "phi", "sigma")
    return(out)
  }
}

get_sv_latent <- function(z) {
  if (inherits(z, "svdraws")) {
    return(as.numeric(tail(as.matrix(z$latent[[1]]), 1)))
  } else {
    return(as.numeric(z$latent))
  }
}
