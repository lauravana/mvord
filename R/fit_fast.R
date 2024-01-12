###### fast ######
build_error_struct_fast <- function(eobj, tpar)
{
  ## takes the transformed parameters and builds initializes some attributes of cor_general eobjs
  ndim <- attr(eobj, "ndim")
  covar <- attr(eobj, "covariate")
  nlev <- NCOL(covar)
  npar1 <- attr(eobj, "npar")/nlev
  corr_pars <- sapply(seq_len(nlev), function(l) {
    nu <- tpar[(l - 1) * npar1 + seq_len(npar1)]
    angles <- pi * exp(nu)/(1 + exp(nu))
    cosmat <- diag(ndim)
    cosmat[lower.tri(cosmat)] <- cos(angles)
    S1 <- matrix(0, nrow = ndim, ncol = ndim)
    S1[, 1L] <- 1
    S1[lower.tri(S1, diag = T)][-(1:ndim)] <- sin(angles)
    #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
    tLmat <- sapply(seq_len(ndim),
                    function(j) cosmat[j, ] * cumprod(S1[j, ]))
    sigma <- crossprod(tLmat)
    sigma[lower.tri(sigma)]
  })
  if (npar1 == 1) dim(corr_pars) <- c(1, nlev)
  rVec <- corr_pars
  sd <- rep(1, ndim)
  return(list(rVec = rVec, sdVec = sd))
}



transf_par_fast <- function(par, rho) {
  par_sigma <- par[rho[["npar.thetas"]] + rho[["npar.betas"]] +
                     seq_len(attr(rho[["error.structure"]], "npar"))]
  sigmas <- build_error_struct_fast(rho[["error.structure"]], par_sigma)
  par_beta <- par[rho[["npar.thetas"]] + seq_len(rho[["npar.betas"]])]
  betatilde <- rho[["constraints_mat"]] %*% par_beta
  betatilde_fast <- par_beta
  par_theta <- rho[["transf_thresholds"]](par[seq_len(rho[["npar.thetas"]])], rho,
                                          betatilde)

  thetatilde <- lapply(seq_len(rho[["ndim"]]), function(j)
    par_theta[[j]] + rho[["thold_correction"]][[j]](betatilde, k = j, rho = rho))

  #   <- sapply(1:rho[["mult.obs, function(j) rho[["x[[j]] %*% beta[[j]])
  pred.fixed <- lapply(seq_len(rho[["ndim"]]), function(j) as.double(rho[["xfast"]] %*% betatilde_fast[rho[["indjbeta_mat"]][j,]]))
  #pred.fixed <- lapply(rho[["indjbeta_fast, function(j) as.double(crossprod(rho[["xfast, betatilde_fast[j])))

  # pred.upper  <- lapply(seq_len(rho[["ndim), function(j) {
  #   th_u <- c(thetatilde[[j]], rho[["inf.value)[rho[["y[, j]]
  #       (th_u - pred.fixed[[j]] - rho[["offset[[j]])/sigmas[["sdVec
  # })#/sigmas[["sdVec
  pred.upper  <- vapply(seq_len(rho[["ndim"]]), function(j) {
    th_u <- c(thetatilde[[j]], rho[["inf.value"]])[rho[["y"]][, j]]
    th_u - pred.fixed[[j]] - rho[["offset"]][[j]]
  }, FUN.VALUE =  double(rho$n))/sigmas[["sdVec"]]
  pred.lower  <- vapply(seq_len(rho[["ndim"]]), function(j) {
    th_l <- c(-rho[["inf.value"]], thetatilde[[j]])[rho[["y"]][, j]]
    th_l - pred.fixed[[j]] - rho[["offset"]][[j]]
  }, FUN.VALUE =  double(rho$n))/sigmas[["sdVec"]]

  predu <- do.call("rbind",lapply(rho[["combis_fast"]], function(h){
    pred.upper[h[["ind_i"]], h[["combis"]], drop = F]
  }))

  predl <- do.call("rbind",lapply(rho[["combis_fast"]], function(h){
    pred.lower[h[["ind_i"]], h[["combis"]], drop = F]
  }))
  predr <- unlist(lapply(rho[["combis_fast"]], function(h){
    sigmas$rVec[h[["r"]]]
  }))

  predu_univ <- pred.upper[rho[["ind_univ"]]]
  predl_univ <- pred.lower[rho[["ind_univ"]]]

  list(U = predu, L = predl, U_univ = predu_univ, L_univ = predl_univ,
       corr_par = predr)
}



PLfun_fast <- function(par, rho){
  tmp <- transf_par_fast(par, rho)
  #r_mat <- tmp[["corr_par"]]#[, rho$dummy_pl_lag == 1, drop = F]
  logp <- double(rho[["n"]])
  ## check for q_i = 1
  pr <- rho[["link"]][["F_uni"]](tmp[["U_univ"]]) -
    rho[["link"]][["F_uni"]](tmp[["L_univ"]])
  pr[pr < .Machine$double.eps] <- .Machine$double.eps
  logp[rho[["ind_univ"]][,1]] <- log(pr)
  ## iterate over bivariate pairs
  prh <- rho[["link"]][["F_biv_rect"]](
    U = tmp[["U"]],
    L = tmp[["L"]],
    r = tmp[["corr_par"]])
  prh[prh < .Machine$double.eps] <- .Machine$double.eps
  -sum(rho[["weights_fast"]] * log(c(prh, pr)))
}
