#' @title Error Structures in mvord
#' @description Different \code{error.structures} are available in \pkg{mvord}:
#' \itemize{
#' \item general correlation structure (default) \code{cor_general(~ 1)},
#' \item general covariance structure \code{cov_general(~ 1)},
#' \item factor dependent correlation structure \code{cor_general(~ f)},
#' \item factor dependent covariance structure \code{cov_general(~ f)},
#' \item equicorrelation structure \code{cor_equi(~ 1)},
#' \item covariate dependent equicorrelation structure \code{cor_equi(~ S)},
#' \item AR(1) correlation structure \code{cor_ar1(~ 1)}, or
#' \item covariate dependent AR(1) correlation structure \code{cor_ar1(~ S)}.
#' }
#' For more details see vignette.
#' @param formula \code{\link{formula}} object
#' @param value specifies values of the correlation (and variance) parameters.
#' For \code{cor_equi()} and \code{cor_ar1()} it can be either a vector of correlations (in (-1,1)) of length one
#' (same correlation for all subjects) or a vector of length equal to the number of subjects.
#' For \code{cor_general()}, it can be either a vector of the lower triangular correlation matrix elements
#' (same structure for all subjects) or a matrix with number of rows equal to the number of subjects.
#' For \code{cov_general()}, it can be either a vector of the lower triangular covariance matrix elements (including the diagonal)
#' (same structure for all subjects) or a matrix with number of rows equal to the number of subjects.
#' Default is \code{value = numeric(0)} object.
#' In this case the correlation parameters are initialized with zero (and variance parameters with 1 for \code{cov_general})
#' @param fixed logical specifying whether the parameters of the error structure should not be optimized in the procedure, but will be \cr
#' fixed to the values specified in the argument \code{value}. Defaults to \code{fixed = FALSE}.
#' @export
#' @name error_struct
cov_general <-
  ## Constructor for the cov_general class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    eobj <- list(name = "cov_general",
                  formula = formula,
                  type = "covariance", value = value, fixed = fixed)
    attr(eobj, "formula") <- formula
    class(eobj) <- c("cov_general", "error_struct")
    eobj
  }
#' @rdname error_struct
#' @export
cor_general <-
  ## Constructor for the cor_general class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    obj <- list(name = "cor_general", formula = formula,
                  type = "correlation", value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_general", "error_struct")
    obj
  }
#' @rdname error_struct
#' @export
cor_ar1 <-
  ## Constructor for the cor_ar1 class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    obj <- list(name = "cor_ar1",
                  formula = formula,
                  type = "correlation", value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_ar1", "error_struct")
    obj
  }
#' @rdname error_struct
#' @export
cor_equi <-
  ## Constructor for the cor_equi class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    obj <- list(name = "cor_equi",
                  formula = formula,
                  type = "correlation", value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_equi", "error_struct")
    obj
  }
########################################

build_error_struct <-
  ## extractor for correlation matrix
  function(eobj, ...) UseMethod("build_error_struct")

build_error_struct_fixed <-
  ## extractor for correlation matrix
  function(eobj, ...) UseMethod("build_error_struct_fixed")

start_values <-
  ## extractor for correlation matrix
  function(eobj, ...) UseMethod("start_values")

initialize <-
  ## initializes the structures
  function(eobj, ...) UseMethod("initialize")

finalize_fun <-
  ## finalizes the structures
  function(eobj, ...) UseMethod("finalize_fun")

finalize <-
  ## initializes the structures
  function(eobj, ...) UseMethod("finalize")

get_covariate <-
  ## initializes the structures
  function(eobj, ...) UseMethod("get_covariate")

init_fun <-
  ## initializes the structures
  function(eobj, ...) UseMethod("init_fun")

#################
##   Methods for error_struct
#################
formula.error_struct <-
  ## Accessor for the covariate formula
  function(x, ...) eval(attr(x, "formula"))

get_covariate.error_struct <- function(eobj, data.x, contrasts) {
  covar_mat <- lapply(data.x, function(x)
    suppressWarnings(model.matrix(formula(eobj),
                                  model.frame(formula(eobj), x, na.action = function(x) x),
                                  contrasts.arg = contrasts)))
  ## check if covariate matrices are equal
  if (!all(sapply(1:(length(covar_mat) - 1), function(i)
    all(covar_mat[[i]] == covar_mat[[i+1]], na.rm = T)))) {
    stop("Covariates in error structure must be
         constant across outcomes!")
  }
  # make one matrix
  covar_mat1 <- sapply(1:ncol(covar_mat[[1]]), function(k){
      xtcol <- do.call(cbind,lapply(covar_mat, `[`,, k))
      xtcol_final <- apply(xtcol,1,function(i) unique(i[!is.na(i)]))
      xtcol_final
  })
  attributes(covar_mat1) <- attributes(covar_mat[[1]])[1:2]
  covar_mat1
}

initialize.error_struct <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, data, contrasts)
  {
    attr(eobj, "ynames") <- colnames(data$y)
    #attr(eobj, "subjnames") <- rownames(data$y)
    attr(eobj, "ndim") <- length(data$x)
    attr(eobj, "nobs") <- nrow(data$y)
    attr(eobj, "covariate") <-
      get_covariate(eobj, data.x = data$x, contrasts = contrasts)
    eobj
  }

update.error_struct <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, data, contrasts)
  {
    attr(eobj, "ynames") <- colnames(data$y)
    #attr(eobj, "subjnames") <- rownames(data$y)
    attr(eobj, "ndim") <- length(data$x)
    attr(eobj, "nobs") <- nrow(data$y)
    attr(eobj, "covariate") <-
      get_covariate(eobj, data.x = data$x, contrasts = contrasts)
    eobj
  }

finalize.error_struct <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, tpar)
  {
    eobj <- finalize_fun(eobj, tpar)

    #eobj$value_tmp <- NULL #comment 26.02.2021

    #  attr(eobj, "subjnames") <- NULL
    #attr(eobj, "ynames") <- NULL
    #  attr(eobj, "ndim") <- NULL
   # attr(eobj, "nobs") <- NULL
    # attr(eobj, "covariate") <- NULL
    attr(eobj, "npar.cor") <- NULL
    attr(eobj, "npar.sd") <- NULL
    eobj
  }


###############################
### Methods for cov_general ###
###############################
start_values.cov_general <- function(eobj) {
  ## builds starting values for the correlation structure
  tmp <- double(attr(eobj, "npar"))
  ## TODO for the given values
  tmp
}

init_fun.cov_general <-
  function(eobj,  data, contrasts)
  {
    form <- formula(eobj)
    if (length(all.vars(form)) > 1)
      stop("Only one factor is supported in cov_general.")
    ## if intercept included rewrite formula without
    if (length(all.vars(form)) == 1 & attr(terms(form), "intercept") == 1)
      attr(eobj, "formula") <- as.formula(sprintf("~ 0 + %s", all.vars(form)))
    eobj <- initialize(eobj, data, contrasts)
    n    <- attr(eobj, "nobs")
    ndim <- attr(eobj, "ndim")
    r <- eobj$value
    ## check value
    if (length(r) == 0) {
      ind_sd <- cumsum(c(1, ndim - seq_len(ndim - 1) + 1))
      r <- numeric(ndim *  (ndim - 1)/2 + ndim)
      r[ind_sd] <- 1
    } # if default set to zero
    if (length(r) ==  (ndim *  (ndim - 1)/2  + ndim))  r <- matrix(rep.int(r, n), ncol = ndim * (ndim - 1)/2 + ndim, byrow= T) # if only one vector of length ndim *  (ndim - 1)/2 default set to zero
    if (nrow(r) != n) stop("Number of rows of argument value in cor_general() is not equal to number of subjects.")
    ## TODO - check positive semi-definiteness??
    ## end check
    eobj$value_tmp <- r
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, ndim * (ndim - 1)/2 * NCOL(attr(eobj, "covariate")))
    attr(eobj, "npar.sd") <-  ifelse(eobj$fixed, 0, ndim *  NCOL(attr(eobj, "covariate")))
    attr(eobj, "npar") <-   attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    if(length(all.vars(form)) == 1 && !is.factor(data$x[[1]][, all.vars(form)]))
      stop("For cov_general covariate must be factor!")

    eobj
  }

## ind_sd <- cumsum(c(1, ndim - seq_len(ndim - 1) + 1))
build_error_struct_fixed.cov_general <-
  ## builds the correlation matrix when fixed = T of cor_general objects
  function(eobj, tpar = NULL, rveclen = NULL)
  {
    ## takes the transformed parameters and builds initializes some attributes of cor_general objects
    ndim <- attr(eobj, "ndim")
    ind_sd <- cumsum(c(1, ndim - seq_len(ndim - 1) + 1))
    smat <- sqrt(eobj$value_tmp[, ind_sd])
    tmp <- sapply(combn(ndim, 2, simplify = F), function(x) smat[,x[1]] * smat[,x[2]])
    rmat <- eobj$value_tmp[, - ind_sd] / tmp
    return(list(rVec = rmat, sdVec = smat))
  }


build_error_struct.cov_general <-
  function(eobj, tpar, rveclen = NULL)
  {
    ## takes the transformed parameters and builds/initializes some attributes of
    ## cor_general objects
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    nlev <- NCOL(covar)
    npar.cor <- attr(eobj, "npar.cor")/nlev
    corr_pars <- sapply(1:nlev, function(l) {
      nu <- tpar[(l - 1) * npar.cor + seq_len(npar.cor)]
      angles <- pi * exp(nu)/(1 + exp(nu))
      cosmat <- diag(ndim)
      cosmat[lower.tri(cosmat)] <- cos(angles)
      S1 <- matrix(0, nrow = ndim, ncol = ndim)
      S1[, 1L] <- 1
      S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
      tLmat <- sapply(1:ndim,
                      function(j) cosmat[j, ] * cumprod(S1[j, ]))
      sigma <- crossprod(tLmat)
      sigma[lower.tri(sigma)]
    })
    if (is.null(ncol(corr_pars))) dim(corr_pars) <- c(1, nlev)
    rVec <- covar %*% t(corr_pars)
    pos <- npar.cor * nlev
    sd_lev <- sapply(1:nlev,
                     function(l) exp(tpar[pos + (l - 1) * ndim + seq_len(ndim)]))
    sdVec <- tcrossprod(covar, sd_lev)
    return(list(rVec = rVec, sdVec = sdVec))
  }

finalize_fun.cov_general <-
  function(eobj, tpar)
  {
    ## takes the transformed parameters and finalizez cor_general eobjs
    if (eobj$fixed){
      attr(eobj, "par") <- numeric(0)
    } else {
      ndim <- attr(eobj, "ndim")
      covar <- attr(eobj, "covariate")
      nlev <- NCOL(covar)
      npar.cor <- attr(eobj, "npar")/nlev - ndim
      corr_mat <- lapply(1:nlev, function(l) {
      nu <- tpar[(l - 1) * npar.cor + seq_len(npar.cor)]
      angles <- pi * exp(nu)/(1 + exp(nu))
      cosmat <- diag(ndim)
      cosmat[lower.tri(cosmat)] <- cos(angles)
      S1 <- matrix(0, nrow = ndim, ncol = ndim)
      S1[, 1L] <- 1
      S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
      tLmat <- sapply(1:ndim,
                      function(j) cosmat[j, ] * cumprod(S1[j, ]))
      sigma <- crossprod(tLmat)
      sigma[lower.tri(sigma)]
    })
    pos <- npar.cor * nlev
    cov_vec <- c(unlist(corr_mat),
                 exp(tpar[(pos + 1) : attr(eobj, "npar")]))
    ## names
    ynames <- attr(eobj, "ynames")
    ind <- combn(ndim,2)
    fnames <- colnames(covar)
    if (eobj$formula == ~1) {
      ## correlation names
      names.corr <-
        sapply(seq_len(NCOL(ind)), function(j)
          sprintf("corr %s %s", attr(eobj, "ynames")[ind[1,j]],
                  attr(eobj, "ynames")[ind[2,j]]))

      ## std deviation names
      names.sigma <- paste("sigma", ynames)
    } else { ## if factor dependent
      ## correlation names
      names.corr.pair <- apply(ind, 2, function(i)
        paste(ynames[i], collapse = " "))
      names.corr <- paste("corr", rep(fnames, each = NCOL(ind)),
                          rep(names.corr.pair, nlev))
      ## std deviation names
      names.sigma <- paste("sigma",
                           rep(fnames, each = ndim),
                           rep(ynames, nlev))
    }
    attr(eobj, "par") <- cov_vec
    attr(eobj, "parnames") <- c(names.corr, names.sigma)
    }
    eobj
  }

#################################
#### Methods for cor_general ####
#################################
start_values.cor_general <- function(eobj) {
  ## builds starting values for the correlation structure
  tmp <- double(attr(eobj, "npar"))
  ## TODO for the given values
  tmp
}

init_fun.cor_general <-
  function(eobj,  data, contrasts)
  {
    ## initializes some attributes of cor_general eobjs
    form <- formula(eobj)
    if (length(all.vars(form)) > 1)
      stop("Only one factor is supported in cor_general.")
    ## if intercept included rewrite formula without
    if (length(all.vars(form)) == 1 & attr(terms(form), "intercept") == 1)
      attr(eobj, "formula") <- as.formula(sprintf("~ 0 + %s", all.vars(form)))
    eobj <- initialize(eobj, data, contrasts)

    n    <- attr(eobj, "nobs")
    ndim <- attr(eobj, "ndim")

    npar <-  ndim *  (ndim - 1)/2 * NCOL(attr(eobj, "covariate"))

    r <- eobj$value
    ## check value
    if (length(r) == 0) r <- numeric(ndim *  (ndim - 1)/2) # if default set to zero
    if (length(r) == n) r <- cbind(r)
    if (length(r) == ndim *  (ndim - 1)/2)  r <- matrix(rep.int(r, n), ncol = ndim * (ndim - 1)/2, byrow= T) # if only one vector of length ndim *  (ndim - 1)/2 default set to zero

    if (nrow(r) != n) stop("Number of rows of argument value in cor_general() is not equal to number of subjects.")

    ## TODO - check positive semi-definiteness??
    ## end check
    eobj$value_tmp <- r
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, npar)
    attr(eobj, "npar.sd") <- 0
    attr(eobj, "npar") <-   attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    if(length(all.vars(form)) == 1 && !is.factor(data$x[[1]][, all.vars(form)]))
      stop("For cor_general covariate must be factor!")
    eobj
  }


build_error_struct_fixed.cor_general <-
  ## builds the correlation matrix when fixed = T of cor_general objects
  function(eobj, tpar = NULL, rveclen = NULL)
  {
    ## takes the transformed parameters and builds initializes some attributes of cor_general objects
    sd <- rep.int(1, attr(eobj, "ndim"))
    return(list(rVec = eobj$value_tmp, sdVec = sd))
  }


build_error_struct.cor_general <-
  function(eobj, tpar, rveclen = NULL)
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
    rVec <- tcrossprod(covar, corr_pars)
    sd <- rep(1, ndim)
    return(list(rVec = rVec, sdVec = sd))
  }

finalize_fun.cor_general <-
  function(eobj, tpar)
  {
  if (eobj$fixed){
    attr(eobj, "par") <- numeric(0)
  } else {
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    nlev <- NCOL(covar)
    npar1 <- attr(eobj, "npar")/nlev
    corr_mat <- lapply(seq_len(nlev), function(l) {
      nu <- tpar[(l - 1) * npar1 + seq_len(npar1)]
      angles <- pi * exp(nu)/(1 + exp(nu))
      cosmat <- diag(ndim)
      cosmat[lower.tri(cosmat)] <- cos(angles)
      S1 <- matrix(0, nrow = ndim, ncol = ndim)
      S1[, 1L] <- 1
      S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
      tLmat <- sapply(seq_len(ndim),
                      function(j) cosmat[j, ] * cumprod(S1[j, ]))
      sigma <- crossprod(tLmat)
      sigma[lower.tri(sigma)]
    })
    corr_vec <- unlist(corr_mat)
    ## names
    ynames <- attr(eobj, "ynames")
    ind <- combn(ndim,2)
    if (eobj$formula == ~1) {
      ## correlation names
      names.corr <-
        sapply(seq_len(NCOL(ind)), function(j)
          sprintf("corr %s %s", attr(eobj, "ynames")[ind[1,j]],
                  attr(eobj, "ynames")[ind[2,j]]))
    } else { ## if factor dependent
      ## correlation names
      names.corr.pair <- apply(ind, 2, function(i)
        paste(ynames[i], collapse = " "))
      names.corr <- paste("corr", rep(colnames(covar), each = NCOL(ind)),
                          rep(names.corr.pair, ncol(covar)))
    }

    attr(eobj, "par") <- corr_vec
    attr(eobj, "parnames") <- names.corr[seq_along(tpar)]
  }
    eobj
  }

#############################
#### Methods for cor_ar1 ####
#############################
start_values.cor_ar1 <- function(eobj) {
  ## builds starting values for the correlation structure
  if (eobj$fixed) {
    tmp <- double(0)
  } else {
    n <- attr(eobj, "nobs")
    S <- attr(eobj, "covariate")
    r <- eobj$value_tmp
    tr <- 1/2 * (log(1 + r) - log(1 - r)) # transformed r
    StS <- crossprod(S)
    Chol <- chol(StS)
    Cholinv <- backsolve(Chol, diag(ncol(Chol)))
    StSinv <- chol2inv(Chol, size = ncol(StS))
    Sty <- crossprod(S, tr)
    alpha <- drop(crossprod(StSinv, Sty))
    tmp <- alpha
  }
  tmp
}

init_fun.cor_ar1 <-
  function(eobj,  data, contrasts)
  {
    ## initializes some attributes of cor_equi eobjs
    #form <- formula(eobj)
    eobj <- initialize(eobj, data, contrasts)
    n <- attr(eobj, "nobs")
    r <- eobj$value
    ## check value
    if (any(abs(r) > 1)) stop("The argument value in cor_ar1() is not a valid correlation parameter.")
    if (length(r) == 0) r <- 0
    if (length(r) == 1) {
      r <- rep.int(r, n)
    } else {
      if (length(r) != n) stop("Length of argument value in cor_ar1() is not equal to number of subjects.")
    }
    ## end check
    eobj$value_tmp <- r
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, NCOL(attr(eobj, "covariate")))
    attr(eobj, "npar.sd") <- 0
    attr(eobj, "npar") <- attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    eobj
  }


#eobj <- rho[["error.structure"]]
build_error_struct_fixed.cor_ar1 <-
  ## builds the correlation matrix when fixed = T of cor_ar1 objects
  function(eobj, tpar = NULL, rveclen = NULL){
    r <- eobj$value_tmp
    ndim <- attr(eobj, "ndim")
    corr_pars <- do.call(rbind, lapply(seq_along(r), function(i){
      sigma <- diag(ndim)
      sigma[lower.tri(sigma)]  <- r[i]^sequence((ndim-1):1)
      sigma <- sigma + t(sigma) - diag(ndim)
      sigma[lower.tri(sigma)]
    }))
    if (is.null(ncol(corr_pars))) dim(corr_pars) <- c(length(corr_pars), 1)
    sdVec <- rep(1, ndim)
    list(rVec = corr_pars, sdVec=sdVec)
  }

build_error_struct.cor_ar1 <-
  function(eobj, tpar, rveclen = NULL)
  {
    ## takes the transformed parameters and builds initializes some attributes of cor_general eobjs
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    z <- covar %*% tpar
    r <- z2r(z)
    corr_pars <- sapply(seq_along(r), function(i){
      sigma <- diag(ndim)
      sigma[lower.tri(sigma)]  <- r[i]^sequence((ndim-1):1)
      sigma <- sigma + t(sigma) - diag(ndim)
      sigma[lower.tri(sigma)]
    })
    if (is.null(ncol(corr_pars)))
      dim(corr_pars) <- c(1, length(corr_pars))
    sdVec <- rep(1, ndim)
    list(rVec = t(corr_pars), sdVec=sdVec)
  }

finalize_fun.cor_ar1 <-
  function(eobj, tpar)
  {
    covar <- attr(eobj, "covariate")
    attr(eobj, "par") <- tpar
    attr(eobj, "parnames") <- colnames(covar)[seq_along(tpar)]
    eobj
  }
##############################
#### Methods for cor_equi ####
##############################
start_values.cor_equi <- function(eobj) {
  ## builds starting values for the correlation structure
  if (eobj$fixed) {
    tmp <- double(0)
  } else {
    n <- attr(eobj, "nobs")
    S <- attr(eobj, "covariate")
    r <- eobj$value_tmp
    tr <- 1/2 * (log(1 + r) - log(1 - r)) # transformed r
    StS <- crossprod(S)
    Chol <- chol(StS)
    Cholinv <- backsolve(Chol, diag(ncol(Chol)))
    StSinv <- chol2inv(Chol, size = ncol(StS))
    Sty <- crossprod(S, tr)
    alpha <- drop(crossprod(StSinv, Sty))
    if (any(abs(tr - S %*% alpha) > 1E-10)) cat("Supplied starting values in argument value might be inconsistent with the formula of the error.structure.\n")
    tmp <- alpha
  }
  tmp
}

init_fun.cor_equi <-
  function(eobj,  data, contrasts)
  {
    ## initializes some attributes of cor_equi eobjs
    #form <- formula(eobj)
    eobj <- initialize(eobj, data, contrasts)
    n <- attr(eobj, "nobs")
    r <- eobj$value
    if (length(r) == 0) r <- 0
    ## check
    if (any(abs(r) > 1)) stop("The argument value in cor_equi() is not a valid correlation parameter.")
    if (length(r) == 1) {
      r <- rep.int(r, n)
    } else {
      if (length(r) != n) stop("Length of argument value in cor_equi() is not equal to number of subjects.")
    }
    ## end checks
    eobj$value_tmp <- r
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, NCOL(attr(eobj, "covariate")))
    attr(eobj, "npar.sd") <- 0
    attr(eobj, "npar") <- attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    eobj
  }

build_error_struct_fixed.cor_equi <-
  function(eobj, tpar = NULL, rveclen = NULL)
  {
    ## tpar argument: transformed parameters (from optimizer)
    ## builds the correlation and standard deviation parameters for cor_equi eobjs
    r <- eobj$value_tmp
    ndim <- attr(eobj, "ndim")
    npar1 <- rveclen #ndim * (ndim - 1)/2
    corr_pars <- matrix(rep(r, npar1), ncol = npar1)
    if (is.null(ncol(corr_pars)))
      dim(corr_pars) <- c(length(corr_pars), 1)
    sdVec <- rep(1, ndim)
    list(rVec = corr_pars, sdVec=sdVec)
  }

build_error_struct.cor_equi <-
  function(eobj, tpar, rveclen = NULL)
  {
    ## tpar argument: transformed parameters (from optimizer)
    ## builds the correlation and standard deviation parameters for cor_equi eobjs
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    z <- covar %*% tpar
    r <- z2r(z)
    npar1 <- rveclen#ndim * (ndim - 1)/2
    corr_pars <- matrix(rep(r, npar1), ncol = npar1)
    if (is.null(ncol(corr_pars)))
      dim(corr_pars) <- c(length(corr_pars), 1)
    sdVec <- rep(1, ndim)
    list(rVec = corr_pars, sdVec=sdVec)
  }

finalize_fun.cor_equi <-
  function(eobj, tpar)
  {
    ## finalizes some attributes of cor_equi eobjs
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    attr(eobj, "par") <- tpar
    attr(eobj, "parnames") <- colnames(covar)[seq_along(tpar)]
    eobj
  }

########################################################
########################################################
#' @title Extracts Error Structure of Multivariate Ordinal Regression Models.
#' @description
#' A generic function which extracts for each subject the estimated
#' error structure parameters from objects of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param type choose type \code{"sigmas"}, \code{"alpha"}, \code{"corr"}, or \code{"z"}.
#' @param ... further arguments passed to or from other methods.
#' @details \itemize{
#' \item{\code{sigmas}} {extracts the correlation/covariance matrices corresponding to each subject.
#'             Applicable in line with \code{cor_general, cov_general, cor_equi, cor_ar1}.}
#' \item{\code{alpha}} {extracts the parameters of the covariate dependent error structure.
#' Applicable in line with \code{cor_equi, cor_ar1}.}
#' \item{\code{corr}} {extracts the subject-specific correlation parameters. Applicable in
#' line with \code{cor_equi}, \code{cor_ar1}.}
#' \item{\code{z}} {extracts the subject-specific Fisher-z score. Applicable in line
#' with \code{cor_equi, cor_ar1}.}}
#' @export
error_structure <- function(object, type, ...) UseMethod("error_structure")
#' @rdname error_structure
#' @export
error_structure.mvord <- function(object, type = NULL, ...)  {
  val <- error_structure(object$error.struct, type = type , ...)
  val
}
error_structure.cor_general <- function(eobj, type, ...){
  par <- attr(eobj, "par")
  npar <- attr(eobj, "npar")# npar <- length(par)
  ndim <- attr(eobj, "ndim")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  nlev <- NCOL(covar)
  if(npar == 0){
    par <- eobj$value
    npar <- length(par)
  }
  npar.cor <- npar/nlev
  corr_lev <- lapply(seq_len(nlev), function(l) {
    sigma <- diag(ndim)
    sigma[lower.tri(sigma)] <- par[(l - 1) * npar.cor + seq_len(npar.cor)]
    s <- sigma + t(sigma) - diag(ndim)
    colnames(s) <- rownames(s) <- ynames
    s
  })
  indlev <- apply(covar, 1, function(x) which(x == 1))
  corr_n <- corr_lev[indlev]
  names(corr_n) <- rownames(attr(eobj, "covariate"))
  corr_n
}

error_structure.cov_general <- function(eobj, type, ...){
  npar <- attr(eobj, "npar")
  par <- attr(eobj, "par")
  ndim <- attr(eobj, "ndim")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  nlev <- NCOL(covar)
  #NEW
  if(npar == 0){
    par <- eobj$value
    npar <- length(par)
#TODO level
    #NEW
    ndim <- attr(eobj, "ndim")
    ind_sd <- cumsum(c(1, ndim - seq_len(ndim - 1) + 1))
    smat <- sqrt(eobj$value_tmp[, ind_sd])
    tmp <- sapply(combn(ndim, 2, simplify = F), function(x) smat[,x[1]] * smat[,x[2]])
    rmat <- eobj$value_tmp[, - ind_sd] / tmp

    cov_n <- lapply(seq_len(NROW(smat)), function(i){
      R <- diag(ndim) * smat[i,]
      R[lower.tri(R)] <- rmat[i,]
      R
    })

    ###
  } else{
  npar.cor <- npar/nlev - ndim
  cov_lev <- lapply(seq_len(nlev), function(l) {
    R <- diag(ndim)
    R[lower.tri(R)] <- par[(l - 1) * npar.cor + seq_len(npar.cor)]
    R <- R + t(R) - diag(ndim)
    s <- par[nlev * npar.cor + (l - 1) * ndim + seq_len(ndim)]
    sigma <- t(s * R) * s
    colnames(sigma) <- rownames(sigma) <- ynames
    sigma
  })
  indlev <- apply(covar, 1, function(x) which(x == 1))
  cov_n <- cov_lev[indlev]
  }
  names(cov_n) <- rownames(attr(eobj, "covariate"))
  cov_n
}

error_structure.cor_equi <- function(eobj, type, ...){
  npar <- attr(eobj, "npar")
  par <- attr(eobj, "par")
  ndim <- attr(eobj, "ndim")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  covar <- attr(eobj, "covariate")
  if(length(par) == 0) par <- eobj$value
  z <- covar %*% par
  colnames(z) <- "Fisher-z Score"
  r <- z2r(z)
  colnames(r) <- "Correlation"
  sigmas <- lapply(seq_along(r), function(i) {
    tmp <- matrix(r[i], nrow = ndim, ncol = ndim)
    diag(tmp) <- 1
    rownames(tmp) <- colnames(tmp) <- ynames
    tmp
  })
  names(sigmas)  <-  rownames(attr(eobj, "covariate"))

  if (!is.null(type)){
    par <- switch(type,
                  alpha = par,
                  sigmas =  sigmas,
                  corr = r,
                  z = z)
  }
  return(par)
}

error_structure.cor_ar1 <- function(eobj, type, ...){
  npar <- attr(eobj, "npar")
  par <- attr(eobj, "par")
  ndim <- attr(eobj, "ndim")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  if(length(par) == 0) par <- eobj$value
  z <- covar %*% par
  colnames(z) <- "Fisher-z Score"
  r <- z2r(z)
  colnames(r) <- "Correlation"
  sigmas <-  lapply(seq_along(r), function(i){
    tmp <- diag(ndim)
    tmp[lower.tri(tmp)]  <- r[i]^sequence((ndim-1):1)
    tmp <- tmp + t(tmp) - diag(ndim)
    rownames(tmp) <- colnames(tmp) <- ynames
    tmp
  })
  names(sigmas)  <-  rownames(attr(eobj, "covariate"))

  if (!is.null(type)){
    par <- switch(type,
                  alpha = par,
                  sigmas =  sigmas,
                  corr = r,
                  z = z)
  }
  return(par)
}

# error_structure.cor_rel_var <- function(eobj, type, ...){
#   npar <- attr(eobj, "npar")
#   par <- attr(eobj, "par")
#   ndim <- attr(eobj, "ndim")
#   covar <- attr(eobj, "covariate")
#   ynames <- attr(eobj, "ynames")
#   nobs <- attr(eobj, "nobs")
#   npar.cor <- attr(eobj, "npar.cor")
#   R <- diag(ndim)
#   R[lower.tri(R)] <- par[seq_len(npar.cor)]
#   R <- R + t(R) - diag(ndim)
#   s <- c(1, par[npar.cor + seq_len(attr(eobj, "npar.sd"))])
#   sigma <- t(s * R) * s
#   colnames(sigma) <- rownames(sigma) <- ynames
#   sigmas <- rep(list(sigma), nobs)
#   names(sigmas) <- rownames(attr(eobj, "covariate"))
#   sigmas
# }

##########################################
## IF numeric SE should be computed ... ##
##########################################
# corr_jac <- function(eobj, tpar, ...) UseMethod("corr_jac")
# corr_jac.cor_general <- function(eobj, tpar){
#   nlev <- NCOL(attr(eobj, "covariate"))
#   ndim <- attr(eobj, "ndim")
#   npar.cor <- ndim * (ndim - 1)/2
#   lapply(seq_len(nlev), function(l)
#     t(sapply(seq_len(npar.cor), function(i) grad(function(x) corr_jac_num_fct(ndim, x, i),
#                                           x=tpar[(l - 1) * npar.cor + seq_len(npar.cor)]))))
# }
# corr_jac.cov_general <- function(eobj, tpar){
#   nlev <- NCOL(attr(eobj, "covariate"))
#   ndim <- attr(eobj, "ndim")
#   npar.cor <- ndim * (ndim - 1)/2
#   l <- lapply(1:nlev, function(l)
#     t(sapply(1:npar.cor, function(i) grad(function(x) corr_jac_num_fct(ndim, x, i),
#                                           x=tpar[(l - 1) * npar.cor + seq_len(npar.cor)]))))
#   l[length(l) + seq_len(nlev * ndim)] <-
#     exp(tpar[npar.cor * nlev + seq_len(ndim * nlev)])
#   l
# }
#
# corr_jac.cor_ar1 <- function(eobj, tpar){
#   list(diag(attr(eobj, "npar")))
# }
# corr_jac.cor_equi<- function(eobj, tpar){
#   list(diag(attr(eobj, "npar")))
# }
# corr_jac.cor_rel_var <- function(eobj, tpar){
#   ndim <- attr(eobj, "ndim")
#   npar.cor <- attr(eobj, "npar.cor")
#   npar.sd <- attr(eobj, "npar.sd")
#   l <- list(sapply(1:npar.cor, function(i) grad(function(x) corr_jac_num_fct(ndim, x, i),
#                                                 x=tpar[seq_len(npar.cor)])))
#   l[1 + seq_len(npar.sd)] <- exp(tpar[npar.cor + seq_len(npar.sd)])
#   l
# }
#
# corr_jac_num_fct <- function(ndim, nu, i){
#   # i is the ith correlation parameter
#   L <- diag(ndim)
#   angles <- pi * exp(nu)/(1 + exp(nu))
#   L[lower.tri(L)] <- cos(angles)
#   S <-  matrix(0, nrow = ndim - 1, ncol = ndim - 1)
#   S[lower.tri(S,diag=T)] <- sin(angles)
#   S <- apply(cbind(1, rbind(0, S)), 1, cumprod)
#   L <- L * t(S)
#   sigma <- tcrossprod(L)
#   sigma[lower.tri(sigma)][i]
# }


z2r <- function (z) {
  ifelse(z > 354, 1, (exp(2 * z) - 1)/(1 + exp(2 * z)))
}

#' @title Print Method for class error_struc.
#' @description Prints error structure of class \code{\link{error_struct}}.
#' @param x object of class \code{\link{error_struct}}
#' @param ... further arguments passed to or from other methods.
#' @method print error_struct
#' @export
print.error_struct <- function(x, ...){
  cat("Parameters of the error structure:\n")
  print(attr(x, "par"), ...)
}
