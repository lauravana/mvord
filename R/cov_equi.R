cov_equi <-
  ## Constructor for the cor_equi class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    obj <- list(name = "cov_equi",
                formula = formula,
                type = "covariance", value = value, fixed = fixed)
    attr(obj, "formula") <- formula
    class(obj) <- c("cov_equi", "error_struct")
    obj
  }


###############################
### Methods for cov_equi ###
###############################
# eobj <- cov_equi(~1)
# error.structure <- init_fun(eobj, data.mvord, contrasts)
# error.structure <- init_fun.cov_equi(eobj, data.mvord, contrasts)
# attr(error.structure, "npar")

start_values.cov_equi <- function(eobj) {
  ## builds starting values for the correlation structure
  tmp <- double(attr(eobj, "npar"))
  ## TODO for the given values
  tmp
}


init_fun.cov_equi <- function(eobj,  data, contrasts){
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
    attr(eobj, "npar.cor") <- 1
    attr(eobj, "npar.sd") <-  ifelse(eobj$fixed, 0, ndim *  NCOL(attr(eobj, "covariate")))
    attr(eobj, "npar") <-   attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    if(length(all.vars(form)) == 1 && !is.factor(data$x[[1]][, all.vars(form)]))
      stop("For cov_general covariate must be factor!")

    eobj
  }


build_error_struct.cov_equi <-function(eobj, tpar){
    ## tpar argument: transformed parameters (from optimizer)
    ## builds the correlation and standard deviation parameters for cor_equi eobjs
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    nobs <- length(covar)
    #TODO no levels supported
    z <- covar %*% tpar[1]
    #r <- z2r(z)
    npar1 <- ndim * (ndim - 1)/2



    cov_pars <- matrix(rep(z, npar1), ncol = npar1)
    if (is.null(ncol(cov_pars)))
      dim(cov_pars) <- c(length(cov_pars), 1)
    s <- exp(tpar[1 + seq_len(attr(eobj, "npar.sd"))])
    cov_tmp <- diag(s)
    cov_tmp[lower.tri(cov_tmp)] <- cov_tmp[upper.tri(cov_tmp)] <- z[1]
    cor_tmp <- cov2cor(cov_tmp)
    sdVec <- sqrt(s)
    corr_pars <- matrix(rep(cor_tmp[lower.tri(cor_tmp)], nobs), ncol = npar1, byrow = TRUE)

    list(rVec = corr_pars, sdVec=sdVec)
  }

finalize_fun.cov_equi <- function(eobj, tpar){
    ## finalizes some attributes of cor_equi eobjs
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")
    attr(eobj, "par") <- tpar

    ###
    tau <- covar %*% tpar[1]
    #r <- z2r(z)
    # npar1 <- ndim * (ndim - 1)/2



    # cov_pars <- matrix(rep(z, npar1), ncol = npar1)
    # if (is.null(ncol(cov_pars)))
    #   dim(cov_pars) <- c(length(cov_pars), 1)
    s <- sqrt(exp(tpar[1 + seq_len(attr(eobj, "npar.sd"))]))

    ###

    ynames <- attr(eobj, "ynames")
    ind <- combn(ndim,2)
    fnames <- colnames(covar)
    # if (eobj$formula == ~1) {
      ## correlation names
      names.corr <-
        sapply(seq_len(NCOL(ind)), function(j)
          sprintf("corr %s %s", attr(eobj, "ynames")[ind[1,j]],
                  attr(eobj, "ynames")[ind[2,j]]))

      ## std deviation names
      names.sigma <- paste("sigma", ynames)

    attr(eobj, "parnames") <- c("tau", names.sigma)
    attr(eobj, "par") <- c(tau[1], s)
    eobj
  }
