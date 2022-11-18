mvord.fit <- function(rho){
  ## number of thresholds per outcome
  rho$ntheta <- sapply(seq_len(rho$ndim), function(j) nlevels(rho$y[, j]) - 1) # no of categories - 1

  if (is.null(rho$threshold.values)) {
    if ((rho$error.structure$type == "correlation") && (rho$intercept.type == "fixed")) {
      rho$threshold.values <- lapply(seq_len(rho$ndim), function(j) rep.int(NA, rho$ntheta[j]))
    } else if((rho$error.structure$type %in% c("correlation")) && (rho$intercept.type == "flexible")) {
      rho$threshold.values <- lapply(seq_len(rho$ndim), function(j) if(rho$ntheta[j] == 1) 0 else c(0,rep.int(NA,max(rho$ntheta[j]-1,0))))
      cat("Note: First threshold for each response is fixed to 0 in order to ensure identifiability!\n")
    } else if ((rho$error.structure$type %in% c("covariance")) && (rho$intercept.type == "flexible")) {
      rho$threshold.values <- lapply(seq_len(rho$ndim), function(j) if(rho$ntheta[j] == 1) 0 else c(0,1,rep.int(NA,max(rho$ntheta[j]-2,0))))
      cat("Note: First two thresholds for each response are fixed to 0 and 1 in order to ensure identifiability!\n")
    } else if ((rho$error.structure$type %in% c("covariance")) && (rho$intercept.type == "fixed")) {
      rho$threshold.values <- lapply(seq_len(rho$ndim), function(j) if(rho$ntheta[j] == 1) 0 else c(0,rep.int(NA,max(rho$ntheta[j]-1,0))))
      cat("Note: First threshold for each response is fixed to 0 in order to ensure identifiability!\n")
    }
  } else {
    if (length(rho$threshold.values) != rho$ndim) stop("Length of threshold values does not match number of outcomes")
  }
  #values of fixed (non-NA) threshold parameters
  rho$threshold.values.fixed <- lapply(rho$threshold.values, function(x) x[!is.na(x)])
  ## number of non-fixed thresholds per rater
  rho$npar.theta <- sapply(seq_len(rho$ndim), function(j) sum(is.na(rho$threshold.values[[j]])))

  #checks if binary outcome is present
  rho$binary <- any(sapply(seq_len(rho$ndim), function(j) length(rho$threshold.values[[j]]) == 1))

  #rho$threshold.type
  rho$threshold <- set_threshold_type(rho)

  rho$ncat <- rho$ntheta + 1

  rho$nthetas <- sum(rho$ntheta)
  rho$ncat.first.ind <- cumsum(c(1,rho$ntheta))[- (rho$ndim + 1)]
  #number of flexible threshold parameters (in optimizer)
  rho$npar.theta.opt <- rho$npar.theta
  rho$npar.theta.opt[duplicated(rho$threshold.constraints)] <- 0
  rho$npar.thetas <- sum(rho$npar.theta.opt)

  #set offsets from coef.values and updates
  rho <- set_offset_up(rho)
  #for(j in 1L:rho$ndim){
  #  check_rank(j, rho)
  #}

  rho$coef.names <- attributes(rho$x[[1]])$dimnames[[2]]

  rho$constraints <- get_constraints(rho)
  # vector of number of parameters for each coefficient
  rho$npar.beta <- 0
  if (length(rho$constraints) > 0) rho$npar.beta <- sapply(rho$constraints, NCOL)
  #number of total coefficients
  rho$npar.betas <- sum(rho$npar.beta)

  check_args_constraints(rho)

  ##INCLUDE CHECKS here
  check_args_thresholds(rho)

  ##############################################################################################
  rho$ind.thresholds <- get_ind_thresholds(rho$threshold.constraints, rho)

  rho$inf.value <- 10000

  rho$evalOnly <- ifelse((rho$npar.thetas + rho$npar.betas + attr(rho$error.structure, "npar") == 0), TRUE, FALSE)

  ###############################
  # starting values
  ###############################
  if (is.null(rho$start.values)) {
    rho$start <- c(get_start_values(rho), start_values(rho$error.structure))
  } else {
    if(length(unlist(rho$start.values)) !=  (rho$npar.thetas + rho$npar.betas)){
      cat(paste0("length should be ", rho$npar.thetas + rho$npar.betas))
      stop("start.values (theta + beta) has false length", call. = FALSE)
    }
    #transform theta starting values
    theta.start <- rho$threshold.values
    rho$start.values$theta <- lapply(seq_len(rho$ndim), function(j) {
      rho$start.values$theta[[rho$threshold.constraints[j]]]
      #rho$start.values$theta[[1]]
    })

    gamma.start <- lapply(seq_len(rho$ndim), function(j) {
      theta.start[[j]][is.na(theta.start[[j]])] <-  rho$start.values$theta[[j]]
      theta2gamma(theta.start[[j]])
    })
    start.val.gamma <- lapply(seq_len(rho$ndim),
                              function(j) gamma.start[[j]][is.na(rho$threshold.values[[j]])]) #gammas
    start.val.gamma <- unlist(start.val.gamma[!duplicated(rho$threshold.constraints)])
    #transform error.struct parameters
    rho$start <- c(start.val.gamma,
                   (rho$start.values$beta),#betas
                   start_values(rho$error.structure))
  }
  #source("mvord/R/utilities.R")
  rho$transf_thresholds <- switch(rho$threshold,
                                  flexible      = transf_thresholds_flexible,
                                  fix1first     = transf_thresholds_fix1_first,
                                  fix2first     = transf_thresholds_fix2_first,
                                  fix2firstlast = transf_thresholds_fix2_firstlast,
                                  fixall        = transf_thresholds_fixall)
  rho$build_error_struct <- ifelse(attr(rho$error.structure, "npar") == 0, build_error_struct_fixed,  build_error_struct)

  #################################
  ## make covariate matrices
  #################################
  rho$p <- NCOL(rho$x[[1]])
  rho$inds.cat <- lapply(seq_len(rho$ndim), function(j)
    seq_len(rho$ntheta[j]) +  rho$ncat.first.ind[j] - 1)
  rho$indjbeta <- lapply(seq_len(rho$ndim), function(j) {
    c(outer(rho$inds.cat[[j]], (seq_len(rho$p) - 1) * rho$nthetas, "+"))
  })

  dim_ind <- rep.int(seq_len(rho$ndim), rho$ntheta)
  x_new <- rho$x
  # x_new[[4]][416,]
  # rho$XcatU[[4]][416,]
  # x_new[[4]][2,]
  # rho$XcatU[[4]][1,]
  # rho$XcatU[[2]][5,]
  # rho$y[416,]
  # rho$y[1:10,]
  if(rho$p > 0) {
    x_center <- vector("list", rho$ndim)
    x_scale <- vector("list", rho$ndim)
    for (j in seq_len(rho$ndim)) {
      ## standardize numeric variables
      mu <- double(rho$p)
      sc <- rep.int(1, rho$p)
      ind_int <- attr(x_new[[j]], "assign") != 0
      tmp <- scale_mvord(x_new[[j]][, ind_int, drop = FALSE])
      x_new[[j]][, ind_int] <- tmp$x
      mu[ind_int] <- tmp$mu
      sc[ind_int] <- tmp$sc
      x_center[[j]] <- mu
      x_scale[[j]] <- sc
    }
    constraints_scaled <- lapply(seq_along(rho$constraints), function(p) {
      sxp <- sapply(x_scale, "[", p)
      x <- rho$constraints[[p]]
      #first index of 1
      fi <- apply(x, 2, function(y) which(y == 1)[1])
      fi[is.na(fi)] <- 1
      sweep(x * sxp[dim_ind], 2, sxp[dim_ind[fi]], "/")
    })
    if (length(constraints_scaled) > 0) {
      rho$constraints_mat <- bdiag(constraints_scaled)
    } else {
      rho$constraints_mat <- numeric(0)
    }
    ## scale and center for par_beta
    tmp <- (lapply(seq_along(rho$constraints), function(p) {
      sxp <- sapply(x_scale, "[", p)
      mxp <- sapply(x_center, "[", p)
      fi <- apply(rho$constraints[[p]], 2, function(y) which(y == 1)[1])
      fi[is.na(fi)] <- 1
      list(scale = sxp[dim_ind[fi]], center = mxp[dim_ind[fi]])
    }))
    rho$fi_scale  <- unlist(sapply(tmp, "[", "scale"))
    rho$fi_center <- unlist(sapply(tmp, "[", "center"))
  } else rho$constraints_mat <- integer()
  ## Upper and lower matrices
  rho$XcatU <- vector("list", rho$ndim)
  rho$XcatL <- vector("list", rho$ndim)
  if (rho$p > 0) {
    for (j in seq_len(rho$ndim)) {
      ncat <- rho$ntheta[j] + 1
      mm <- model.matrix(~ - 1 + rho$y[,j] : x_new[[j]],
                         model.frame(~ - 1 + rho$y[,j] : x_new[[j]],
                                     na.action = function(x) x))
      rho$XcatL[[j]] <- mm[,-(ncat * (seq_len(rho$p) - 1) + 1), drop = F]
      rho$XcatU[[j]] <- mm[,-(ncat * seq_len(rho$p)), drop = F]
    }
  } else {
    rho$XcatU <- lapply(seq_len(rho$ndim), function(x) integer()) #creates warning, but OK
    rho$XcatL <- lapply(seq_len(rho$ndim), function(x) integer()) #In th_u - xbeta_u : Recycling array of length 1 in vector-array arithmetic is deprecated
  }


  ####################################################
  ## build contrasts for the thetas
  rho$contr_theta <- do.call("rbind", rep.int(list(diag(rho$nthetas)), rho$p))
  ## make function which compute correction factor for thresholds:
  if(rho$p > 0){
    rho$mat_center_scale <- bdiag(rho$constraints) %*% c(rho$fi_center/rho$fi_scale)
  } else rho$mat_center_scale <- integer()
  rho$thold_correction <- vector("list", rho$ndim)
  is.dup <- duplicated(rho$threshold.constraints)
  if(rho$p > 0){
    rho$thold_correction <-lapply(seq_len(rho$ndim), build_correction_thold_fun0, rho = rho)
    if (any(is.dup)) {
      tmp <- lapply(which(is.dup), build_correction_thold_fun, rho = rho)
      rho$thold_correction[which(is.dup)] <- tmp
    }
  } else {
    rho$thold_correction <-lapply(seq_len(rho$ndim), build_correction_thold_fun0, rho = rho)
  }
  #############################################################################
  ## help variables to save computation in the likelihood function
  if (is.null(rho$combis)) {
    rho$combis <- combn(rho$ndim, 2, simplify = FALSE)
  } else {
    check_combis(rho)
  }

  rho$dummy_pl_lag <- sapply(rho$combis, function(x)
    diff(x) <= rho$PL.lag)
  rho$combis <- rho$combis[rho$dummy_pl_lag]

  ## for which subjects is q_i = 1
  ind_i <- rowSums(!is.na(rho$y)) == 1
  rho$ind_univ <- which(!is.na(rho$y) & ind_i, arr.ind=T)
  rho$n_univ <- NROW(rho$ind_univ)
  ## index for subjects containing pair c(k,l)
  rho$ind_kl <- lapply(rho$combis, function(kl)
    rowSums(!is.na(rho$y[, kl])) == 2)

  rho$ind_i <- lapply(seq_along(rho$combis), function(h)  which(rho$ind_kl[[h]]))
  rho$combis_fast <- lapply(seq_along(rho$combis),function(h){
    list("combis" = rho$combis[[h]],
         "ind_i" = rho$ind_i[[h]],
         "r" = rep(h,length(rho$ind_i[[h]])))
  })

  rho$n_biv <- sum(sapply(rho$ind_i, length))

  rho$weights_fast <- rho$weights[c(unlist(rho$ind_i), rho$ind_univ[,1])]

  ### fast ###
  if(("cor_general" %in% class(rho$error.structure)) & is.null(rho$coef.values_input) & (rho$coef.constraints_VGAM == FALSE) &
     attr(rho$error.structure, "formula")[[2]] == 1 & !isTRUE(rho$error.structure$fixed) &
     !anyNA(rho$coef.constraints_input) &
     all(sapply(seq_len(length(rho$x)-1), function(j){
       if((NCOL(rho$x[[j]]) == 0) & (NCOL(rho$x[[j + 1]]) == 0)) {TRUE
       }else if(any(sapply(rho$x, anyNA))) FALSE else{
         rho$x[[j]] == rho$x[[j + 1]]
       }
     }))){
    rho$fast_fit <- TRUE
  } else rho$fast_fit <- FALSE
  #rho$fast_fit <- FALSE #TODO at the moment always use PL_fun

  if(rho$fast_fit){
    #TODO for constraints
    rho$indjbeta_mat <-  do.call("cbind",lapply(rho$constraints, function(x){x[rho$ncat.first.ind,] == 1
    }))
    rho$indjbeta_fast <- lapply(seq_len(rho$ndim), function(j) {
      j + (seq_len(rho$p)-1) * rho$npar.beta
    })
    rho$fi_scale_fast <- rho$fi_scale[1 + (seq_len(rho$p)-1) * rho$ndim]
    rho$fi_center_fast <- rho$fi_center[1 + (seq_len(rho$p)-1) * rho$ndim]
    rho$xfast <- x_new[[1]]
  }
  ##############################################
  ## OPTIMIZE PAIRWISE LIKELIHOOD ##############
  ##############################################
  if (rho$evalOnly) {
    rho$optpar <- numeric(0)
    rho$objective <- c(value = PLfun_old(rho$start , rho))
  } else {
    if (is.character(rho$solver)) {
      if(rho$fast_fit){
        rho$optRes <- suppressWarnings(optimx(rho$start, function(x) PLfun_fast(x, rho),
                                              method = rho$solver,
                                              hessian = FALSE,
                                              control = rho$control))
      } else{
        rho$optRes <- suppressWarnings(optimx(rho$start, function(x) PLfun(x, rho),
                                              method = rho$solver,
                                              hessian = FALSE,
                                              control = rho$control))
      }
      rho$optpar <- unlist(rho$optRes[1:length(rho$start)])
      rho$objective <- unlist(rho$optRes["value"])
    }
    if (is.function(rho$solver)){
      rho$optRes <- rho$solver(rho$start, function(x) PLfun(x, rho))
      if (is.null(rho$optRes$optpar)|is.null(rho$optRes$objvalue)|is.null(rho$optRes$convcode)) stop("Solver function does not return the required objects.")
      rho$optpar <- rho$optRes$optpar
      rho$objective <- rho$optRes$objvalue
      rho$message <- rho$optRes$message

    }
    if (rho$optRes$convcode != 0){
      stop("NO/FALSE CONVERGENCE - choose a different optimizer or different starting values.")
    }
  }
  ## Construct original parameters
  par_beta <- rho$optpar[rho$npar.thetas + seq_len(rho$npar.betas)]
  if(rho$p > 0){
    betatilde <- rho$constraints_mat %*% par_beta
    par_theta <- rho$transf_thresholds(rho$optpar[seq_len(rho$npar.thetas)], rho, betatilde)
    thetatilde <- lapply(seq_len(rho$ndim), function(j)
      par_theta[[j]] + rho$thold_correction[[j]](betatilde, k = j, rho = rho))
    betatildemu <-  betatilde * rho$mat_center_scale
    br <- drop(crossprod(rho$contr_theta, betatildemu[,1]))
    correction <- lapply(rho$inds.cat, function(k) br[k])
  } else {
    correction <- rep.int(list(0), rho$ndim)
    par_theta <- rho$transf_thresholds(rho$optpar[seq_len(rho$npar.thetas)], rho, 0)
    thetatilde <- lapply(seq_len(rho$ndim), function(j) par_theta[[j]])
  }
  ## Thresholds
  rho$theta <- lapply(seq_len(rho$ndim), function(j) thetatilde[[j]] + correction[[j]])
  ## Regression coefficients
  rho$beta <- par_beta/rho$fi_scale
  ############################################
  ##############################################
  ## Compute Standard Errors
  #############################################
  if (rho$se) {
    rho <- PL_se(rho)
  }
  ##############################################
  res <- mvord_finalize(rho)
  if (rho$se) {
    rownames(rho$varGamma) <- colnames(rho$varGamma) <- c(names(unlist(res$theta))[is.na(unlist(rho$threshold.values))][!duplicated(unlist(rho$ind.thresholds))],                                                        names(res$beta),                                                     attr(res$error.struct, "parnames"))
  }
  ## clean up
  rho$XcatU <- NULL
  rho$XcatL <- NULL
  rho$transf_thresholds <- NULL
  rho$get_ind_thresholds <- NULL
  #rho$ind.thresholds <- NULL
  rho$y.NA.ind <- NULL
  rho$binary <- NULL
  rho$intercept.type <- NULL
  rho$intercept <- NULL
  rho$threshold.values.fixed <- NULL
  rho$ncat.first.ind <- NULL
  rho$constraints_mat <- NULL
  rho$ind_univ <- NULL
  rho$ind_kl <- NULL
  # rho$combis <- NULL
  rho$fi_scale <- NULL
  rho$fi_center <- NULL
  rho$mat_center_scale <- NULL
  rho$dummy_pl_lag <- NULL
  rho$inds.cat <- NULL
  rho$thold_correction <- NULL
  rho$contr_theta <- NULL
  rho$error.structure <- NULL
  rho$beta <- NULL
  #rho$theta <- NULL
  #rho$nthetas <- NULL
  #rho$npar.theta.opt <- NULL
  #rho$npar.theta <- NULL
  #rho$npar.thetas <- NULL
  #rho$npar.betas <- NULL
  rho$coef.names <- NULL
  ##
  rho$link$deriv.fun <- NULL
  rho$link$F_biv_rect <- NULL
  rho$link$F_biv <- NULL
  ##
  res$rho <- rho
  res$rho$formula <- rho$formula.input
  res$rho$formula.input <- NULL
  attr(res, "contrasts") <- rho$contrasts
  res$rho$contrasts <- NULL
  rho$timestamp2 <- proc.time()
  res$rho$runtime <- rho$timestamp2 - rho$timestamp1
  res$rho$timestamp1 <- NULL

  class(res) <- "mvord"

  return(res)
}
##########################################################
######                      PL             ###############
##########################################################
# k<-1
# h<-1
# par <- rho$start
# PLfun(par, rho)
PLfun <- function(par, rho){
  tmp <- transf_par(par, rho)
  ## check for q_i = 1
  pr <- double(rho[["n_univ"]])
  pr <- rho[["link"]][["F_uni"]](tmp[["U_univ"]]) -
    rho[["link"]][["F_uni"]](tmp[["L_univ"]])
  pr[pr < .Machine$double.eps] <- .Machine$double.eps
  ## iterate over bivariate pairs
  prh <- double(rho[["n_biv"]])
  prh <- rho[["link"]][["F_biv_rect"]](
    U = tmp[["U"]],
    L = tmp[["L"]],
    r = tmp[["corr_par"]])
  prh[prh < .Machine$double.eps] <- .Machine$double.eps
  # return(-(sum(log(prh)) + sum(log(pr))))
  # -sum(log(c(prh, pr)))
  -sum(rho[["weights_fast"]] * log(c(prh, pr)))
  # - sum(rho$weights * logp)
}
# PLfun <- function(par, rho){
#   tmp <- transf_par(par, rho)
#   #r_mat <- tmp[["corr_par"]]#[, rho$dummy_pl_lag == 1, drop = F]
#   logp <- double(rho[["n"]])
#   ## check for q_i = 1
#   pr <- rho[["link"]][["F_uni"]](tmp[["U_univ"]]) -
#     rho[["link"]][["F_uni"]](tmp[["L_univ"]])
#   pr[pr < .Machine$double.eps] <- .Machine$double.eps
#   logp[rho[["ind_univ"]][,1]] <- log(pr)
#   ## iterate over bivariate pairs
#   prh <- rho[["link"]][["F_biv_rect"]](
#     U = tmp[["U"]],
#     L = tmp[["L"]],
#     r = tmp[["corr_par"]])
#   prh[prh < .Machine$double.eps] <- .Machine$double.eps
#   #logp[ind_i] <- logp[ind_i] +  log(prh)
#   -sum(log(c(prh, pr)))
#   #  - sum(rho$weights * logp)
# }
PLfun_old <- function(par, rho){
  tmp <- transf_par_old(par, rho)
  pred.upper <- tmp$U
  pred.lower <- tmp$L
  r_mat <- tmp$corr_par[, rho$dummy_pl_lag == 1, drop = F]
  logp <- double(rho$n)
  ## check for q_i = 1
  pr <- rho$link$F_uni(pred.upper[rho$ind_univ]) -
    rho$link$F_uni(pred.lower[rho$ind_univ])
  pr[pr < .Machine$double.eps] <- .Machine$double.eps
  logp[rho$ind_univ[,1]] <- log(pr)
  ## iterate over bivariate pairs
  for (h in seq_along(rho$combis)){
    ind_i <- which(rho$ind_kl[[h]])
    r <- r_mat[ind_i, h]
    prh <- rho$link$F_biv_rect(
      U = pred.upper[ind_i, rho$combis[[h]], drop = F],
      L = pred.lower[ind_i, rho$combis[[h]], drop = F],
      r = r)
    prh[prh < .Machine$double.eps] <- .Machine$double.eps
    logp[ind_i] <- logp[ind_i] +  log(prh)
  }
  - sum(rho$weights * logp)
}

.onLoad <- function(library, pkg)
{
  library.dynam("mvord", pkg, library)
  invisible()
}
