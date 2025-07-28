##########################################################
###### PL Standard Errors - no parallel    ###############
##########################################################
PL_se <- function(rho){
  # if (is.null(rho$link$deriv.fun)) {
  # par <- rho$optpar
  # ## construct Jacobian
  # J <- make_jac_matrix(rho)
  # J.inv <- solve(J)
  # if(rho$control$trace != 0) cat("Computing variability and hessian matrix numerically ... \n")
  # Vi_num <- matrix(0, ncol = length(par), nrow = rho$n)
  # for (i in seq_len(rho$n)) {
  #   if (i %% 100 == 0)  cat('Computed gradient for', i, 'out of', rho$n,'subjects\n')
  #   Vi_num[i, ] <- grad(function(par) neg_logPL_comp_i(par, rho, i), par, method = "Richardson")
  # }
  # cat("\n")
  # rho$Vi_num <-  Vi_num %*% J.inv
  # rho$V <- rho$n/(rho$n - length(par)) * crossprod(rho$Vi_num) # original variability matrix
  # if(rho$control$trace != 0) cat("\nComputing Hessian numerically ... \n")
  # Ht <-  hessian(function(par) PLfun(par, rho), par,
  #                method = "Richardson",
  #                method.args=list(eps=1e-6)) # Fisher matrix H(Gamma transf)
  #   rho$H.inv <-   J %*%  solve(Ht) %*% t(J)
  # } else {
  par <- rho$optpar
  if (rho$evalOnly)   {
    rho$H.inv <- matrix(nrow = 0L, ncol = 0L)
    rho$V <- matrix(nrow = 0L, ncol = 0L)
    rho$varGamma <- matrix(nrow = 0L, ncol = 0L)
    rho$df <- 0L
  } else {
    if ((rho$ndim) == 2 && rho$n <= 1000) {
      if(rho$control$trace != 0) cat("Computing hessian matrix numerically ... \n")
      J <- make_jac_matrix(rho)
      if (length(par) > 0) {
        Ht <-  hessian(function(par) PLfun(par, rho), par,
                       method = "Richardson") # Fisher matrix H(Gamma transf)
        Ht.inv <- tryCatch(chol2inv(chol(Ht)))
        rho$H.inv <- J %*% tcrossprod(Ht.inv, J)
      } else {
        rho$H.inv <- matrix(nrow = 0L, ncol = 0L)
      }
      rho$varGamma <- rho$H.inv ## inverse observed fisher information matrix
      rho$claic <- 2 * rho$objective + 2 * ncol(rho$H.inv)
      rho$clbic <- 2 * rho$objective + log(rho$n) * ncol(rho$H.inv)
      rho$df <- ncol(rho$H.inv)
    } else {
      if(rho$control$trace != 0) cat("Computing variability and hessian matrix analytically ... \n")
      derivs_for_se <- derivs_ana(rho)
      rho$g_list <- derivs_for_se$g_list
      rho$V <- rho$n/(rho$n - NCOL(derivs_for_se$V)) * derivs_for_se$V  ## correct for degrees of freedom
      rho$H.inv <- tryCatch(chol2inv(chol(derivs_for_se$H)))
      rho$varGamma <- rho$H.inv %*% rho$V %*% rho$H.inv ## inverse godambe
      rho$claic <- 2 * rho$objective + 2 * sum(diag(rho$V %*% rho$H.inv))
      rho$clbic <- 2 * rho$objective + log(rho$n) * sum(diag(rho$V %*% rho$H.inv))
      rho$df <-  sum(diag(rho$V %*% rho$H.inv))
    }
  }
  rho$seGamma <- sqrt(diag(rho$varGamma))
  if(rho$control$trace != 0) cat("Done computing the standard errors!\n")
  rho
}
#############################################################################
# #' @title Derivative of rectangle probability wrt correlation parameter
# #' @description This function computes the derivative of the rectangle probability wrt correlation parameter
# #' @param Uk upper predictor for the k-th response
# #' @param Ul upper predictor for the l-th response
# #' @param Lk lower predictor for the k-th response
# #' @param Ll lower predictor for the l-th response
# #' @param r value or vector of the correlation parameters
# #' @param fun function computing the derivative dF(x,y,r)/dr
d_corr_rect <- function(Uk, Ul, Lk, Ll, r, fun) {
  - fun(Uk, Ul, r) + fun(Uk, Ll, r) + fun(Lk, Ul, r) - fun(Lk, Ll, r)
}
# #' @title Derivative of rectangle probability wrt threshold or regression parameter
# #' @description This function computes the derivative of the rectangle probability wrt threshold or regression parameter
# #' @param Uk upper predictor for the k-th response
# #' @param Ul upper predictor for the l-th response
# #' @param Lk lower predictor for the k-th response
# #' @param Ll lower predictor for the l-th response
# #' @param r value or vector of the correlation parameters
# #' @param dUkmat d Uk/d theta or d Uk/d beta
# #' @param dLkmat d Lk/d theta or d Lk/d beta
# #' @param d_biv_fun function computing the derivative dF(x,y,r)/dx of the bivariate pdf (depends on link)
d_rect <- function(Uk, Ul, Lk, Ll, r,
                   dUkmat, dLkmat, d_biv_fun) {
  UU <- d_biv_fun(Uk, Ul, r)
  UL <- d_biv_fun(Uk, Ll, r)
  LU <- d_biv_fun(Lk, Ul, r)
  LL <- d_biv_fun(Lk, Ll, r)
  ((UU - UL) * dUkmat  - (LU - LL) * dLkmat)
}
# #' @title Derivative of rectangle probability wrt standard deviation parameters
# #' @description This function computes the derivative of the rectangle probability wrt standard deviation parameters
# #' @param Uk upper predictor for the k-th response
# #' @param Ul upper predictor for the l-th response
# #' @param Lk lower predictor for the k-th response
# #' @param Ll lower predictor for the l-th response
# #' @param r value or vector of the correlation parameters
# #' @param d_biv_fun function computing the derivative dF(x,y,r)/dx of the bivariate pdf (depends on link)
d_sd_rect <- function(Uk, Ul, Lk, Ll, sd_k, r, d_biv_fun){
  UU <- d_biv_fun(Uk, Ul, r)
  UL <- d_biv_fun(Uk, Ll, r)
  LU <- d_biv_fun(Lk, Ul, r)
  LL <- d_biv_fun(Lk, Ll, r)
  1/sd_k * (UU * Uk - UL * Uk - LU * Lk + LL * Lk)
}

d_theta_rect_kl <- function(k, l, indkl, rho) {
  FOO <- d_rect(Uk = rho$U[indkl, k],
                Ul = rho$U[indkl, l],
                Lk = rho$L[indkl, k], Ll = rho$L[indkl, l],
                r = rho$r,
                dUkmat = rho$B1[[k]][indkl, ]/ rho$std.dev.mat[indkl, k],
                dLkmat = rho$B2[[k]][indkl,] /rho$std.dev.mat[indkl, k] ,
                d_biv_fun = rho$link$deriv.fun$dF2dx)
  return(FOO)
}

d_psi_rect_kl <- function(k, l, indkl, rho, r, sd_k, sd_l) {
  Ck <- rho$C[rho$indjbeta2[[k]], , drop = F]
  Cl <- rho$C[rho$indjbeta2[[l]], , drop = F]
  XmatUk <- rho$XU[[k]] %*% Ck
  XmatUl <- rho$XU[[l]] %*% Cl
  XmatLk <- rho$XL[[k]] %*% Ck
  XmatLl <- rho$XL[[l]] %*% Cl
  - as.matrix(d_rect(Uk = rho$U[indkl, k], Ul = rho$U[indkl, l],
                     Lk = rho$L[indkl, k], Ll = rho$L[indkl, l],
                     r = r, dUkmat =  (XmatUk)[indkl, ]/sd_k,
                     dLkmat = (XmatLk)[indkl, ]/sd_k,
                     d_biv_fun = rho$link$deriv.fun$dF2dx)
              +
                d_rect(Uk = rho$U[indkl, l], Ul = rho$U[indkl, k],
                       Lk = rho$L[indkl, l], Ll = rho$L[indkl, k], r = r,
                       dUkmat =  (XmatUl)[indkl, ]/sd_l,
                       dLkmat = (XmatLl)[indkl, ]/sd_l,
                       d_biv_fun = rho$link$deriv.fun$dF2dx))
}


nll_grad_univ <- function(U, L, XU, XL, sd_j, dhdsigma, wts, dfun, pfun) {
  ## Univariate Probability
  pr <- pfun(U) - pfun(L)
  pr[pr < .Machine$double.eps] <- .Machine$double.eps
  #############################################
  ## Gradient components in univariate case: ##
  #############################################
  dpsi <- - 1/sd_j * (dfun(U) * XU -  dfun(L) * XL)
  dsd <-  1/sd_j * (dfun(U) * U - dfun(L) * L) * dhdsigma
  l <- list(dpsi = wts * 1/pr * dpsi,
            dsd  = wts * 1/pr * dsd)
  return(l)
}

## derivative functions dr/dalpha for the correlation part of error structures
dr_dalpha_general <- function(S, lag = NULL){
  S
}

dr_dalpha_equi <- function(alpha, S, lag = NULL){
  salpha <- c(S %*% alpha)
  exp(2 * salpha)/(exp(2 * salpha) + 1)^2  * 4 * S
}

dr_dalpha_ar1 <- function(alpha, S, lag){
  salpha <- c(S %*% alpha)
  r <- (exp(2 * salpha) - 1)/(exp(2 * salpha) + 1)
  const <-  exp(2 * salpha)/(exp(2 * salpha) + 1)^2  * 4

  lag * r^(lag - 1) * const * S
}
## derivative functions dh/dgamma for the variance part of error structures
dh_dgamma_general <- function(S){
  S
}

constraints_theta <- function(rho) {
  pick.col.theta <- lapply(seq_len(rho$ndim), function(j)
    switch(rho$threshold,
           flexible      = seq_len(rho$ntheta[j]),
           fix1first     = seq_len(rho$ntheta[j])[-1],
           fix2first     = seq_len(rho$ntheta[j])[-c(1,2)],
           fix2firstlast = seq_len(rho$ntheta[j] - 1)[-1]))
  C_theta <- bdiag(lapply(unique(rho$threshold.constraints), function(j) {
    idq <- which(rho$threshold.constraints == j)
    do.call("rbind", lapply(idq, function(k){
      Ik <- diag(rho$ntheta[k])
      Ik[, pick.col.theta[[k]], drop = F]
    }))}))

  id_rows <- NULL
  indtheta <- cumsum(c(0, rho$ntheta[seq_len(rho$ndim)])) + 1

  for (j in unique(rho$threshold.constraints)) {
    idq <- which(rho$threshold.constraints == j)
    id_rows <- c(id_rows,
                 c(sapply(idq, function(i) indtheta[i]: (indtheta[i + 1] - 1))))
  }
  C_theta[id_rows, ] <- C_theta
  # C_theta <- bdiag(lapply(unique(rho$threshold.constraints), function(j) {
  #   idq <- which(rho$threshold.constraints == j)
  # do.call("rbind", lapply(idq, function(k){
  #   Ik <- diag(rho$ntheta[k])
  #   Ik[, pick.col.theta[[k]], drop = F]
  # }))}))
  C_theta
}

set_offset_threshold_u <- function(rho) {
  lapply(seq_len(rho$ndim), function(j) {
    x <- rho$threshold.values[[j]]
    x[is.na(x)] <- 0
    ou <- c(double(rho$ntheta[j]), rho$inf.value)[rho$y[, j]]
    rho$XU[[j]][, seq_along(x), drop = F] %*% x + ou - rho$offset[[j]]
  })
}

set_offset_threshold_l <- function(rho) {
  lapply(seq_len(rho$ndim), function(j) {
    x <- rho$threshold.values[[j]]
    x[is.na(x)] <- 0
    ol <- c(-rho$inf.value, double(rho$ntheta[j]))[rho$y[, j]]
    rho$XL[[j]][, seq_along(x), drop = F] %*% x  + ol - rho$offset[[j]]
  })
}
derivs_ana <- function(rho){
  ############################################
  ## function for analytic gradient and hessian
  #############################################
  par <- rho$optpar
  #######################################
  ## Upper and lower matrices
  #######################################
  ## making theta a vector
  ## ideally theta would be the raw parameters in par
  theta <- lapply(seq_len(rho$ndim), function(j)
    switch(rho$threshold,
           flexible      = rho$theta[[j]],
           fix1first     = rho$theta[[j]][-1],
           fix2first     = rho$theta[[j]][-c(1,2)],
           fix2firstlast = rho$theta[[j]][seq_len(rho$ntheta[j] - 1)[-1]]))
  #theta <- unlist(theta[unique(rho$threshold.constraints)])
  theta <- unlist(theta[!duplicated(rho$threshold.constraints)])

  C_theta <- constraints_theta(rho)
  if (length(rho$constraints) > 0) {
    rho$C <- as.matrix(bdiag(list(C_theta, bdiag(rho$constraints))))
  } else {
    rho$C <- as.matrix(C_theta)
  }
  rho$indjbeta2 <- lapply(seq_len(rho$ndim), function(j) {
    c(outer(rho$inds.cat[[j]], (seq_len(rho$p + 1) - 1) * rho$nthetas, "+"))
  })
  rho$XU <- rho$XL <- vector("list", rho$ndim)
  for (j in seq_len(rho$ndim)) {
    XcatU <- XcatL <- integer(0)
    ncat <- rho$ntheta[j] + 1
    if (rho$p > 0) {
      mm <- model.matrix(~ -1 + rho$y[,j] : rho$x[[j]],
                         model.frame(~ -1 + rho$y[,j] : rho$x[[j]],
                                     na.action = function(x) x))
      XcatL <- mm[,-(ncat * (seq_len(rho$p) - 1) + 1), drop = FALSE]
      XcatU <- mm[,-(ncat * seq_len(rho$p)), drop = FALSE]
    }
    B2 <- (col(matrix(0, rho$n, ncat)) == c(unclass(rho$y[, j])))
    BU <- B2[,-ncat, drop = FALSE]
    BL <- B2[,-1, drop = FALSE]
    rho$XU[[j]] <- cbind(BU, -XcatU)
    rho$XL[[j]] <- cbind(BL, -XcatL)
  }

  psi <- rho$C %*% c(theta, rho$beta)
  ##
  npar.err <-  attr(rho$error.structure, "npar")
  #npar.err.start <-  attr(rho$error.structure, "npar_start")
  npar.cor <-  attr(rho$error.structure, "npar.cor")
  npar.sd  <-  attr(rho$error.structure, "npar.sd")


  if (npar.err > 0) {
    par_sigma <-  par[rho$npar.thetas + rho$npar.betas + seq_len(npar.err)]
    alpha <- par_sigma[seq_len(npar.cor)]
    gamma <- par_sigma[-seq_len(npar.cor)]
  }

  S <- attr(rho$error.structure, "covariate")

  sigmas <- rho$build_error_struct(rho$error.structure,
                                   par_sigma,
                                   rveclen = length(rho$combis)) ## TODO: not working nicely with specified combis
  ###############################################
  rho$offsetu <- set_offset_threshold_u(rho)
  rho$offsetl <- set_offset_threshold_l(rho)

  rho$U  <- sapply(seq_len(rho$ndim), function(j) {
    rho$XU[[j]] %*% psi[rho$indjbeta2[[j]]] + rho$offsetu[[j]]
  })/sigmas$sdVec

  rho$L  <- sapply(seq_len(rho$ndim), function(j) {
    rho$XL[[j]] %*% psi[rho$indjbeta2[[j]]] + rho$offsetl[[j]]
  })/sigmas$sdVec
  ########################################################
  rho$std.dev.mat <- sigmas$sdVec ## matrix of standard deviations
  if (is.null(dim(rho$std.dev.mat))){
    rho$std.dev.mat <- matrix(1, ncol = rho$ndim, nrow = rho$n)
  }
  ######################################################
  ## First the univariate case (q_i = 1)
  ######################################################
  g_list <- NULL
  if (NROW(rho$ind_univ) > 0){
    ## univariate functions
    dfun <- rho$link$deriv.fun$dF1dx
    pfun <- rho$link$F_uni
    ## univariate observations:
    g_list <- lapply(unique(rho$ind_univ[, 2]), function(j) {
      idj <- rho$ind_univ[rho$ind_univ[, 2] == j, , drop = F]
      subj <- idj[, 1] # subject ids for j-th response
      ## for each j
      U <- rho$U[idj]
      L <- rho$L[idj]
      sd_j <- rho$std.dev.mat[idj]
      wts  <- rho$weights[subj]
      XU <- rho$XU[[j]][subj, ] %*% rho$C[rho$indjbeta2[[j]], ]
      XL <- rho$XL[[j]][subj, ] %*% rho$C[rho$indjbeta2[[j]], ]
      Sj <- S[subj, , drop = F]
      gr <- nll_grad_univ(U, L, XU, XL, sd_j, Sj, wts, dfun, pfun)
      ## fill into an (n x length(par)) gradient matrix
      if (npar.sd > 0) {
        colpos_sd <- rho$npar.thetas + rho$npar.betas + npar.cor +
          (seq_len(NCOL(Sj)) - 1) * rho$ndim + j
      } else {
        colpos_sd <-  0
      }
      gr_mat <- matrix(0, nrow = rho$n, length(par))
      gr_mat[subj, seq_len(rho$npar.thetas + rho$npar.betas)] <- gr$dpsi
      gr_mat[subj, colpos_sd]  <- gr$dsd
      gr_mat
    })
  }
  #####################################
  ## take each possible pair (k, l)
  ######################################
  r_mat <- sigmas$rVec[, rho$dummy_pl_lag == 1, drop = F]
  it0 <- length(g_list)
  #it <- 1
  for (it in (it0 + seq_along(rho$combis))) {
    comb <- rho$combis[[it - it0]]
    ## for each pair make an indicator for each subject where the pair applies
    indkl <- rho$ind_kl[[it - it0]]
    k <- comb[1]
    l <- comb[2]
    Uk <- rho$U[indkl, k] ## upper predictor k-th response
    Ul <- rho$U[indkl, l] ## upper predictor l-th response
    Lk <- rho$L[indkl, k] ## lower predictor k-th response
    Ll <- rho$L[indkl, l] ## lower predictor l-th response
    sd_k <- rho$std.dev.mat[indkl, k] ## standard deviation of k-th response
    sd_l <- rho$std.dev.mat[indkl, l] ## standard deviation of l-th response
    Skl <- S[indkl, , drop = F]
    ## correlation
    rkl <- rho$r <- r_mat[indkl, (it - it0)]
    ## pr_{kl}
    pr <- rep(1, rho$n)
    pr[indkl] <- rho$link$F_biv_rect(
      U = cbind(Uk, Ul), L = cbind(Lk, Ll), r = rkl)
    pr[pr < .Machine$double.eps] <- .Machine$double.eps
    ## vector h_kl will contain the gradient for all d -log p_{kl}/d pars
    dpsi <- dcorr <- dstddev <- NULL

    ##################
    ## dpsi
    ##################
    if (rho$npar.betas + rho$npar.thetas > 0) {
      dpsi <- matrix(0,
                     ncol = rho$npar.betas + rho$npar.thetas,
                     nrow = rho$n)

      dpsi[indkl, ] <-  d_psi_rect_kl(k, l, indkl, rho, rkl, sd_k, sd_l)
    }
    ##################
    ## dcorr
    ##################
    if (npar.cor > 0){
      dcorr <- matrix(0, ncol = npar.cor, nrow = rho$n)
      dLdr <- d_corr_rect(Uk, Ul, Lk, Ll, r = rkl,
                          fun = rho$link$deriv.fun$dF2dr)
      drdalpha <- switch(rho$error.structure$name,
                         cor_general = dr_dalpha_general(Skl),
                         cov_general = dr_dalpha_general(Skl),
                         cor_equi    = dr_dalpha_equi(alpha, Skl),
                         cor_ar1     = dr_dalpha_ar1(alpha, Skl, lag = l - k))
      colpos <- switch(rho$error.structure$name,
                       cor_general = (seq_len(NCOL(Skl)) - 1) * rho$ndim * (rho$ndim - 1)/2 + (it - it0),
                       cov_general = (seq_len(NCOL(Skl)) - 1) * rho$ndim * (rho$ndim - 1)/2 + (it - it0),
                       cor_equi    = seq_len(npar.cor),
                       cor_ar1     = seq_len(npar.cor))

      dcorr[indkl, colpos] <- dLdr * drdalpha
    }
    ##################
    ## dstddev
    ##################
    if (npar.sd > 0){
      dstddev <- matrix(0, nrow = rho$n, ncol = npar.sd)
      poslevk <- (seq_len(NCOL(Skl)) - 1) * rho$ndim + k
      poslevl <- (seq_len(NCOL(Skl)) - 1) * rho$ndim + l
      dhdsigma <- dh_dgamma_general(Skl)
      dstddev[indkl, poslevk] <-
        d_sd_rect(Uk, Ul, Lk, Ll, sd_k, rkl,
                  d_biv_fun = rho$link$deriv.fun$dF2dx) * dhdsigma
      dstddev[indkl, poslevl] <-
        d_sd_rect(Ul, Uk, Ll, Lk, sd_l, rkl,
                  d_biv_fun = rho$link$deriv.fun$dF2dx) * dhdsigma
    }
    g_list[[it]] <- rho$weights * 1/pr * cbind(dpsi, dcorr, dstddev)
  }
  ## matrix containing the gradients for each subject
  Vi <- Reduce("+", g_list)

  ## Variability matrix
  V <- crossprod(Vi)
  ## Hessian matrix - expected Fisher information matrix.
  H <- Reduce("+", lapply(g_list, crossprod))
  list(g_list = g_list, V = V, H = H)
}#
## This will help general purpose optimizers
#### make Jacobian
# make_jacobian <- function(rho) {
#  diag(length(rho$optpar))
#}
#############################################################
###### neg loglikelihood component for each subject i #######
######                 for numeric gradient           #######
# #############################################################
# transf_par_i <- function(par, rho, i) {
#   par_sigma <- par[rho$npar.thetas + rho$npar.betas + seq_len(attr(rho$error.structure, "npar_start"))]
#   sigmas <- build_error_struct(rho$error.structure, par_sigma)
#   if (is.null(dim(sigmas$sdVec))){
#     sdi <- sigmas$sdVec
#   } else {
#     sdi <- sigmas$sdVec[i, ]
#   }
#   par_beta <- par[rho$npar.thetas + seq_len(rho$npar.betas)]
#   betatilde <- rho$constraints_mat %*% par_beta
#   par_theta <- rho$transf_thresholds(par[seq_len(rho$npar.thetas)], rho, betatilde)
#   thetatilde <- lapply(seq_len(rho$ndim), function(j)
#     par_theta[[j]] + rho$thold_correction[[j]](betatilde, k = j, rho = rho))
#
#   pred.upper  <- sapply(seq_len(rho$ndim), function(j) {
#    th_u <- c(thetatilde[[j]], rho$inf.value)[rho$y[i, j]]
#    xbeta_u <- sum(rho$XcatU[[j]][i, ] * betatilde[rho$indjbeta[[j]]])
#    th_u - xbeta_u - rho$offset[[j]]
#   })/sdi
#   pred.lower  <- sapply(seq_len(rho$ndim), function(j) {
#     th_l <- c(-rho$inf.value, thetatilde[[j]])[rho$y[i, j]]
#     xbeta_l <- sum(rho$XcatL[[j]][i, ] * betatilde[rho$indjbeta[[j]]])
#     th_l - xbeta_l - rho$offset[[j]]
#   })/sdi
#   list(U = pred.upper, L = pred.lower,
#        corr_par = sigmas$rVec[i, , drop=F])
# }
# ##############
# neg_logPL_comp_i <- function(par, rho, i) {
#   # transform parameters and get upper and lower bounds
#   tmp <- transf_par_i(par, rho, i)
#   U <- tmp$U
#   L <- tmp$L
#   r_mat <- tmp$corr_par
#   q <- which(!is.na(rho$y[i, ]))
#   if (length(q) == 1){
#     pr <- rho$link$F_uni(U[q]) - rho$link$F_uni(L[q])
#     logPLi <- rho$weights[i] * log(max(pr, .Machine$double.eps))
#   } else {
#     combis <- combn(q, 2)
#     combis <- combis[,which((combis[2,] - combis[1,])  <= rho$PL.lag), drop = FALSE]
#     logPLi <- 0
#     dim(U) <- dim(L) <- c(1, length(U))
#     for (h in seq_len(ncol(combis))) {
#       r <- r_mat[, h]
#       pr <- rho$link$F_biv_rect(U = U[,combis[, h], drop = F],
#                         L = L[,combis[, h], drop = F],
#                         r)
#       logPLi <- logPLi + rho$weights[i] * log(max(pr, .Machine$double.eps))
#   }
# }
#   -logPLi
# }
#
# transf_thresholds_fix2_firstlast_jac <- function(rho, j, gamma_j, i){
#   recursive.theta <- function(i) {
#     if (i == 0) 0
#     else return ((exp(gamma_j[i]) + recursive.theta(i - 1))/(1 + exp(gamma_j[i])))
#   }
#     theta <- sapply(seq_along(gamma_j), function(i)
#       recursive.theta(i))
#     theta[i]
# }
#
#
# transf_thresholds_fix2_first_jac <- function(rho,j,gamma_j,i){
#       c(0, cumsum(c(1 ,exp(gamma_j))))[i+2]
# }
#
make_jac_matrix <- function(rho) {
  par <- rho$optpar
  parthetabeta <- unlist(par[seq_len(rho$npar.thetas + rho$npar.betas)])
  jac <- list()
  if (length(parthetabeta)>0) {
    jac <- list(numDeriv::jacobian(function(x) transf_psi(x, rho), x = parthetabeta))
  }
  # first.ind.theta <- sapply(rho$ind.thresholds, "[", 1)
  # transf_thresholds_jac <- switch(rho$threshold,
  #                                 fix2firstlast = transf_thresholds_fix2_firstlast_jac,
  #                                 fix2first = transf_thresholds_fix2_first_jac)
  #
  # gamma <- par[seq_len(rho$npar.thetas)]
  # if (rho$threshold == "flexible") {
  #   jac <- lapply((seq_len(rho$ndim))[which(rho$npar.theta.opt > 0)], function(j){ #rho$npar.theta
  #     emat <- diag(rho$ntheta[j])
  #     if (ncol(emat) >= 2) {
  #       emat[,1] <- 1
  #       for (k in 2:ncol(emat))
  #         emat[(k:nrow(emat)), k] <-
  #           exp(gamma[(first.ind.theta[j]) + seq_len(rho$ntheta[j]-1)])[k - 1]
  #     }
  #     emat
  #   })
  # } else {
  #   if (rho$threshold == "fix1first") {
  #     jac <- lapply((seq_len(rho$ndim))[which(rho$npar.theta.opt > 0)], function(j){ #rho$npar.theta
  #       emat <- diag(rho$ntheta[j])
  #       if (ncol(emat) >= 2) {
  #         emat[,1] <- 1
  #         for (k in 2:ncol(emat))
  #           emat[(k:nrow(emat)), k] <-
  #             exp(gamma[(first.ind.theta[j]) + seq_len(rho$npar.theta.opt[j])-1])[k - 1] #rho$npar.theta
  #       }
  #       emat[-1,-1]
  #     })
  #   } else {
  #     jac <- lapply((seq_len(rho$ndim))[which(rho$npar.theta.opt > 0)], function(j){ #rho$npar.theta
  #       gamma_j <- gamma[first.ind.theta[j] + seq_len(rho$npar.theta.opt[j]) - 1] #rho$npar.theta
  #       t(sapply(1:length(gamma_j),
  #                function(i) grad(function(x) transf_thresholds_jac(rho, j, x, i), x=gamma_j)))
  #     })
  #   }
  # }
  ## Jacobian for BETAS: no transform
  # jac[sum(rho$npar.theta.opt > 0) + seq_len(rho$npar.betas)] <- 1/rho$fi_scale
  ## Jacobian for ERROR STRUCTURE
  parsigma <- par[rho$npar.thetas + rho$npar.betas + seq_len(attr(rho$error.structure, "npar"))]
  if (length(parsigma) > 0) {
    corr.jac <- corr_jac(rho$error.structure, parsigma)
    # jac[sum(rho$npar.theta.opt > 0) + rho$npar.betas + seq_along(corr.jac)] <- corr.jac #rho$npar.theta
    jac <-  c(jac, corr.jac)
  }
  ## make Jacobian matrix
  J <- as.matrix(bdiag(jac))
  return(J)
}
