#' @title Marginal Predictions for Multivariate Ordinal Regression Models.
#'
#' @description
#' Obtains marginal predictions/fitted measures for objects of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param type types \code{"prob"}, \code{"class"}, \code{"linpred"}, \code{"cum.prob"},  \code{"all.prob"} are available.
#' @param newdata (optional) data frame of new covariates and new responses.
#' The names of the variables should correspond to the names of the
#'  variables used to fit the model. By default the data on which the model
#'  was estimated is considered.
#' @param subjectID (optional) vector specifying for which subjectIDs the predictions\cr or fitted values should be computed.
#' @param newoffset (optional) list of length equal to the number of outcomes, each element containing a vector of offsets to be considered.
#' @param ... further arguments passed to or from other methods.
#' @details The following types can be chosen in \code{marginal_predict}:
#' \tabular{ll}{
#'   \code{type} \tab description\cr
#'   \code{"prob"} \tab  predicted marginal probabilities for the observed response categories. Used as default if newdata contains column(s) for response variable(s).\cr
#'   \code{"class"} \tab predicted marginal classes of the observed responses.\cr
#'   \code{"linpred"} \tab predicted linear predictor \cr
#'   \code{"cum.prob"} \tab predicted marginal cumulative probabilities for the observed response categories.\cr
#'   \code{"all.prob"} \tab predicted marginal probabilities for all ordered classes of each response. Used as default if newdata contains no column(s) for response variable(s).
#'   }
#'
##' If provided, the row names of the output correspond to the subjectIDs, otherwise they correspond to the row id of the observations.
#' @seealso \code{\link{predict.mvord}}, \code{\link{joint_probabilities}}
#' @export
marginal_predict <- function(object, newdata = NULL,
                             type = NULL,
                             subjectID = NULL,
                             newoffset = NULL, ...){
  typearg <- type
  if (is.null(typearg)) type <- "prob"

  if(!(type %in% c("prob", "linpred", "class", "cum.prob", "all.prob")))
    stop("Invalid type chosen. Only types 'prob','linpred', 'class', 'cum.prob' and 'all.prob' are available.")

  if(is.null(newdata)){
    x <- object$rho$x
    y <- object$rho$y
    offset <- object$rho$offset
  } else {
    newdata <- as.data.frame(newdata)
    tmp <- prepare_newdata(object, newdata, newoffset)
    x <- tmp$x
    y <- tmp$y
    if (all(is.na(y))) {
      ## If no response is available in newdata:
      if (is.null(typearg)) {
        message("Setting type default to all.prob")
        type <- "all.prob"
      }
      ## Also, prob and linpred cannot be computed, as they need a response.
      if (type %in% c("prob", "linpred")) stop("newdata contains no response columns. Only types 'class', 'cum.prob' and 'all.prob' are supported.")
    }
    object$error.struct <- tmp$error.struct
    offset <- tmp$offset
  }
  if (is.null(subjectID)) ind <- seq_len(NROW(y)) else {
    if(!all(subjectID %in% rownames(y))) stop("Not all subjectIDs in data!")
    ind <- match(subjectID, rownames(y))
  }
  sigma <- error_structure(object, type = "sigmas")
  stddevs <- sqrt(t(sapply(sigma, diag)))

  marg_predictions <- switch(type,
                             prob     = marg_pred_prob_fun(object, y, x, offset, stddevs, ind),
                             linpred  = marg_pred_linpred_fun(object, y, x, offset, stddevs, ind),
                             all.prob = marg_pred_allprob_fun(object, y, x, offset, stddevs, ind),
                             cum.prob = marg_pred_cumprob_fun(object, y, x, offset, stddevs, ind),
                             class    = marg_pred_class_fun(object, y, x, offset, stddevs, ind))
  return(marg_predictions)
}


#' @title Predict method for Multivariate Ordinal Regression Models.
#' @description Obtains predicted or fitted values for objects of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param type types \code{"class"}, \code{"prob"} and \code{"cum.prob"} are available.
#' @param newdata (optional) data frame of new covariates and new responses.
#' @param subjectID (optional) vector specifying for which subjectIDs the predictions\cr or fitted values should be computed.
#' @param newoffset (optional) list of length equal to the number of outcomes, each element containing a vector of offsets to be considered.
#' @param ... further arguments passed to or from other methods.
#' @details
#' \tabular{ll}{
#'   \code{type} \tab description\cr
#'   \code{"class"} \tab combination of response categories with the highest probability. Used as default if newdata contains no column(s) for response variable(s) or all responses are NA.\cr
#'   \code{"prob"} \tab fitted joint probability for the observed response categories \cr
#'   \tab or the categories provided in the response column(s) in \code{newdata}. Used as default if newdata contains column(s) for response variable(s).\cr
#'   \tab If response column(s) in \code{newdata} contain only NAs or is missing, this type is not supported.\cr
#'   \code{"cum.prob"} \tab fitted joint cumulative probability for the observed response\cr
#'   \tab categories or the categories provided in the response column(s) in \code{newdata}.\cr
#'   \tab If response column(s) in \code{newdata} contain only NAs or is missing, this type is not supported.\cr
#'   }
#' If provided, the (row) names of the output correspond to the subjectIDs,
#' otherwise they correspond to the row id of the observations.
#' @seealso \code{\link{marginal_predict}}, \code{\link{joint_probabilities}}
#' @method predict mvord
#' @export
predict.mvord <- function(object, newdata = NULL, type = NULL,
                          subjectID = NULL, newoffset = NULL, ...){
  typearg <- type
  if (is.null(typearg)) type <- "prob"

  # checks
  if (is.null(object$rho$link$F_multi)) stop("Multivariate probabilities cannot be computed! Try marginal_predict()!")
  if(!(type %in% c("prob", "class", "cum.prob")))  stop("Invalid type chosen. Only types 'prob', 'class' and 'cum.prob' are available.")

  if(is.null(newdata)){
    x <- object$rho$x
    y <- object$rho$y
    offset <- object$rho$offset
  } else {
    newdata <- as.data.frame(newdata)
    tmp <- prepare_newdata(object, newdata, newoffset)
    x <- tmp$x
    y <- tmp$y
    if (all(is.na(y))) {
      ## If no response is available in newdata:
      if (is.null(typearg)) {
        message("Setting type default to class")
        type <- "class"
      }
      ## Also, prob and linpred cannot be computed, as they need a response.
      if (type %in% c("prob", "cum.prob")) stop("newdata contains no response columns. Only type 'class' is supported.")
    }
    object$error.struct <- tmp$error.struct
    offset <- tmp$offset
  }

  if(is.null(subjectID)) ind <- seq_len(nrow(y)) else {
    if(!all(subjectID %in% rownames(y))) stop("Not all subjectIDs in data!")
    ind <- match(subjectID, rownames(y))
  }
  ## get correlation/covariance matrices
  sigma <- error_structure(object, type ="sigmas")
  stddevs <- sqrt(t(sapply(sigma, diag)))

  predictions <- switch(type,
                        prob     = pred_prob_fun(object, y, x, offset, stddevs, sigma)[ind],
                        cum.prob = pred_cumprob_fun(object, y, x, offset, stddevs, sigma)[ind],
                        class    = pred_class_fun(object, y, x, offset, stddevs, sigma)[ind, ])
  return(predictions)
}

#' @title Extracts fitted Probabilities for Multivariate Ordinal Regression Models.
#'
#' @description
#' Extracts fitted probabilities for given response categories from a fitted model of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param response.cat vector or matrix with response categories (for each subject one row of length equal to the number of multiple measurements).
#' @param newdata (optional) data frame of new covariates and new responses. The names of the variables should correspond to the names of the
#'  variables used to fit the model. By default the data on which the model was estimated is considered.
#' @param type \code{"prob"} for joint probabilities and \code{"cum.prob"} for joint cumulative probabilities.
#' @param subjectID (optional) vector specifying for which subjectIDs the predictions\cr or fitted values should be computed.
#' @param newoffset (optional) list of length equal to the number of outcomes, each element containing a vector of offsets to be considered.
#' @param ... further arguments passed to or from other methods.
#' @details
# #' \code{newdata} has to be in the same data format as in the fitted object of class \code{'mvord'}.
#'
##' The current implementation supports only in-sample predictions.
#' The row names of the output correspond to the subjectIDs.
#' @seealso \code{\link{predict.mvord}}, \code{\link{marginal_predict}}
#' @export
joint_probabilities <- function(object, response.cat, newdata = NULL,
                                type = "prob", subjectID = NULL, newoffset = NULL,...) {
  #checks
  if (is.null(object$rho$link$F_multi)) stop("Multivariate probabilities cannot be computed! Try marginal_predict()!")
  # args <- list(...)
  # exist <- "newdata" %in% names(args)
  if(!type %in% c("prob", "cum.prob")) stop("Invalid type chosen. Only types prob and cum.prob are available.")
  # if(!exist) newdata <- NULL
  # if (!is.null(newdata)) stop("newdata is not supported at the moment!")
  ndim <- object$rho$ndim
  if (is.null(newdata)){
    nobs <- nobs(object)
    x <- object$rho$x
    y <- object$rho$y
    offset <- object$rho$offset
  } else {
    newdata <- as.data.frame(newdata)
    ## check response.cat
    if (is.vector(response.cat) && length(response.cat) != ndim)
      stop("response.cat must have length equal to the number of outcomes.")
    tmp <- prepare_newdata(object, newdata, newoffset)
    x <- tmp$x
    y <- tmp$y
    if (!is.vector(response.cat) && nrow(response.cat) != nrow(y))
      stop("Number of rows of response.cat does not correspond to number of subjects in newdata.")

    object$error.struct <- tmp$error.struct
    offset <- tmp$offset

  }

  if (is.null(subjectID)) {
    ind <- seq_len(nrow(y))
  } else {
    if(!all(subjectID %in% rownames(y))) stop("Not all subjectIDs in data!")
    ind <- match(subjectID, rownames(y))
  }
  ## get correlation/covariance matrices
  sigma <- error_structure(object, type ="sigmas")
  stddevs <- sqrt(t(sapply(sigma, diag)))

  if(is.vector(response.cat)) {
    response.cat <- matrix(response.cat, ncol = length(response.cat), nrow = nrow(y), byrow = TRUE)
  }

  ytmp <- cbind.data.frame(lapply(seq_len(ndim), function(j){
    if (!all(response.cat[,j] %in% c(NA,levels(y[,j])))) {
      stop("response.cat are different from the categories in the original data set")
    } else {
      ordered(response.cat[,j], levels = levels(y[,j]))
    }
  }))
  Xcat <- make_Xcat(object, ytmp, x)
  betatilde <- bdiag(object$rho$constraints) %*% object$beta

  pred.upper  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_u <- c(object$theta[[j]], object$rho$inf.value)[ytmp[, j]]
    xbeta_u <- as.double(Xcat$U[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_u - xbeta_u - offset[[j]]
  })/stddevs
  pred.lower  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_l <- c(-object$rho$inf.value, object$theta[[j]])[ytmp[, j]]
    xbeta_l <- as.double(Xcat$L[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_l - xbeta_l - offset[[j]]
  })/stddevs
  pred.lower[is.na(pred.lower)] <- -object$rho$inf.value
  pred.upper[is.na(pred.upper)] <- object$rho$inf.value

  if(type == "cum.prob") pred.lower <- matrix(-object$rho$inf.value,
                                              ncol = ndim,
                                              nrow = nrow(pred.upper))
  prob <- object$rho$link$F_multi(U = pred.upper, L = pred.lower,
                                  list_R = lapply(sigma, cov2cor))[ind]
  names(prob) <- rownames(y)[ind]
  return(prob)
}


prepare_newdata <- function(object, newdata, newoffset) {
  if (object$rho$function.name == "mvord") {
    if (!all(object$rho$index %in% colnames(newdata)))
      stop("Subject and outcome index do not appear in column names of newdata.")
    if (!(object$rho$response.name %in% colnames(newdata))) {
      newdata[, object$rho$response.name] <- NA_integer_
    }

    data.mvord <- mvord_data(newdata, object$rho$index,
                             object$rho$response.name,
                             unique(c(object$rho$x.names, object$rho$weights.name)),
                             y.levels = object$rho$levels,
                             response.names = object$rho$response.names)

    y <- data.mvord$y
    x <- lapply(seq_len(object$rho$ndim), function(j) {
      if(all(is.na(data.mvord$x[[j]]))){
        tmp <- data.mvord$x[[j]]
      } else {
        # rhs.form <- object$rho$formula #RH260820
        # rhs.form[[2]] <- NULL #RH260820
        rhs.form <- as.formula(paste0("~", paste(c(0,colnames(data.mvord$x[[j]])), collapse = " + "))) #RH260820 TODO 0 or 1
        rhs.form <- as.formula(paste(as.character(object$rho$formula[-2]), collapse = " "))
        new.rhs.form <- update(rhs.form, ~ . + 1)
        tmp <-  suppressWarnings(model.matrix(new.rhs.form,
                                              model.frame(new.rhs.form,  data.mvord$x[[j]], na.action = function(x) x),
                                              contrasts.arg = attr(object, "contrasts")))
        attribute <- attr(tmp, "assign")
        intercept <- ifelse(attr(terms.formula(object$rho$formula), "intercept") == 1, TRUE, FALSE)
        if (intercept == FALSE){
          attribute <- attribute[-1]
          tmp <- tmp[,-1, drop = FALSE]
        }
        tmp <- tmp[match(rownames(data.mvord$x[[j]]),rownames(tmp)), , drop=FALSE]
        rownames(tmp) <- rownames(data.mvord$x[[j]])
        attr(tmp, "assign") <- attribute
      }
      tmp
    })
    error.struct <- init_fun(object$error.struct, data.mvord, attr(object, "contrasts"))
  } else {
    # if (is.null(object$rho$index)) stop("Estimated model uses MMO2, newdata is long format.")
    if (!all(object$rho$y.names %in% colnames(newdata))) {
      # stop("Response names in newdata do not match with the outcome names in the estimated model.")
      y <- matrix(NA, nrow = nrow(newdata), ncol = object$rho$ndim)
    } else {
      y <- newdata[,object$rho$y.names]
      y <- do.call("cbind.data.frame", lapply(seq_len(ncol(y)), function(j)
        ordered(y[,j], levels = object$rho$levels[[j]])))
    }
    colnames(y) <- object$rho$y.names
    x <- lapply(seq_len(object$rho$ndim), function(j) {
      rhs.form <- as.formula(paste(as.character(object$rho$formula[-2]), collapse = " "))
      new.rhs.form <- update(rhs.form, ~ . + 1)
      tmp <-  suppressWarnings(model.matrix(new.rhs.form,
                                            model.frame(new.rhs.form,  newdata, na.action = function(x) x),
                                            contrasts.arg = attr(object, "contrasts")))
      attribute <- attr(tmp, "assign")
      intercept <- ifelse(attr(terms.formula(object$rho$formula), "intercept") == 1, TRUE, FALSE)
      if (intercept == FALSE){
        attribute <- attribute[-1]
        tmp <- tmp[,-1, drop = FALSE]
      }
      attr(tmp, "assign") <- attribute
      as.data.frame(tmp)
    })
    data.mvord <-  list(y = y, x = x)
    error.struct <- init_fun(object$error.struct, data.mvord,
                             attr(object, "contrasts"))
  }
  if (is.null(newoffset)) {
    newoffset <- lapply(seq_len(object$rho$ndim), function(j) {
      rhs.form <- object$rho$formula
      rhs.form[[2]] <- NULL
      newdata_tmp <- switch(object$rho$function.name,
                            "mvord" = data.mvord$x[[j]],
                            "mvord2" = newdata)
      mf <- model.frame(rhs.form, newdata_tmp,
                        na.action = function(x) x)
      mf[is.na(mf)] <- 0
      if (is.null(model.offset(mf))) {
        ofs <- double(NROW(y))
      } else {
        ofs <- model.offset(mf)
      }
      ofs
    })
  }
  return(list(error.struct = error.struct, y = y, x = x, offset = newoffset))
}

make_Xcat <- function(object, y, x) {
  ndim <- NCOL(y)
  ncat <- NULL
  for (j in seq_len(NCOL(y))) {
    ncat <- c(ncat, nlevels(y[, j]))
  }
  XcatU <- lapply(seq_len(ndim), function(x) integer())
  XcatL <- lapply(seq_len(ndim), function(x) integer())
  if (object$rho$p > 0) {
    for (j in seq_len(ndim)) {
      B2 <- 1 * (col(matrix(0, nrow(y), ncat[j])) ==
                   c(unclass(y[, j])))
      mf <- do.call("cbind", lapply(as.data.frame(x[[j]]),
                                    function(x) B2 * x))
      XcatL[[j]] <- mf[,-(ncat[j] * (seq_len(object$rho$p) - 1) + 1), drop = FALSE]
      XcatU[[j]] <- mf[,-(ncat[j] * seq_len(object$rho$p)), drop = FALSE]
    }
  }
  list(U = XcatU, L = XcatL)
}

## marginal_predict: type == "linpred"
marg_pred_linpred_fun <- function(object, y, x, offset, stddevs, ind) {
  Xcat <- make_Xcat(object, y, x)
  betatilde <- bdiag(object$rho$constraints) %*% object$beta
  pred.upper  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_u <- c(object$theta[[j]], object$rho$inf.value)[y[, j]]
    xbeta_u <- as.double(Xcat$U[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_u - xbeta_u - offset[[j]]
  })/stddevs
  pred.lower  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_l <- c(-object$rho$inf.value, object$theta[[j]])[y[, j]]
    xbeta_l <- as.double(Xcat$L[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_l - xbeta_l - offset[[j]]
  })/stddevs

  colnames(pred.lower) <- colnames(pred.upper) <- object$rho$y.names
  rownames(pred.lower) <- rownames(pred.upper) <- rownames(y)[ind]
  return(list(U = pred.upper[ind, ], L = pred.lower[ind, ]))
}
## marginal_predict: type == "prob"
marg_pred_prob_fun <- function(object, y, x, offset, stddevs, ind) {
  Xcat <- make_Xcat(object, y, x)
  betatilde <- bdiag(object$rho$constraints) %*% object$beta
  pred.upper  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_u <- c(object$theta[[j]], object$rho$inf.value)[y[, j]]
    xbeta_u <- as.double(Xcat$U[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_u - xbeta_u - offset[[j]]
  })/stddevs
  pred.lower  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_l <- c(-object$rho$inf.value, object$theta[[j]])[y[, j]]
    xbeta_l <- as.double(Xcat$L[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_l - xbeta_l - offset[[j]]
  })/stddevs
  prob <- object$rho$link$F_uni(pred.upper) -
    object$rho$link$F_uni(pred.lower)
  prob <- prob[ind, ]
  colnames(prob) <- object$rho$y.names
  rownames(prob) <- rownames(y)[ind]
  return(prob)
}
## marginal_predict: type == "all.prob"
marg_pred_allprob_fun <- function(object, y, x, offset, stddevs, ind) {
  betatilde <- bdiag(object$rho$constraints) %*% object$beta
  probs <- lapply(seq_len(object$rho$ndim), function(j){
    pr <- sapply(seq_len(object$rho$ncat[j]), function(k){
      ytmp <- y[,j]
      ytmp[seq_along(ytmp)] <- object$rho$levels[[j]][k]
      ytmp <- ordered(ytmp, levels = object$rho$levels[[j]])
      dim(ytmp)  <- c(length(ytmp),1)
      Xcat <- make_Xcat(object, ytmp, x[j])
      th_u <- c(object$theta[[j]], object$rho$inf.value)[ytmp[, 1]]
      xbeta_u <- as.double(Xcat$U[[1]] %*% betatilde[object$rho$indjbeta[[j]]])
      pred.upper <- (th_u - xbeta_u - offset[[j]])/stddevs[,j]
      th_l <- c(-object$rho$inf.value, object$theta[[j]])[ytmp[, 1]]
      xbeta_l <- as.double(Xcat$L[[1]] %*% betatilde[object$rho$indjbeta[[j]]])
      pred.lower <- (th_l - xbeta_l - offset[[j]])/stddevs[, j]
      object$rho$link$F_uni(pred.upper) - object$rho$link$F_uni(pred.lower)
    })[ind, ]
    colnames(pr) <- levels(object$rho$y[, j])
    rownames(pr) <- rownames(y)[ind]
    pr
  })
  names(probs) <- object$rho$y.names
  return(probs)
}
marg_pred_cumprob_fun <-  function(object, y, x, offset, stddevs, ind) {
  probs <- marg_pred_allprob_fun(object, y, x, offset, stddevs, ind)
  cum.probs <- lapply(probs, function(x) t(apply(x, 1, cumsum)))
  return(cum.probs)
}

marg_pred_class_fun <-  function(object, y, x, offset, stddevs, ind) {
  probs <- marg_pred_allprob_fun(object, y, x, offset, stddevs, ind)
  y.ord <- as.data.frame(sapply(seq_along(probs), function(j){
    apply(probs[[j]], 1, function(i) {
      class.max <- object$rho$levels[[j]][which.max(i)]
      ifelse(length(class.max)==0, NA, class.max)
    })
  }))
  for (j in seq_along(object$rho$levels))
    y.ord[,j] <- ordered(y.ord[, j], levels = object$rho$levels[[j]])
  colnames(y.ord) <- object$rho$y.names
  return(y.ord[ind, ])
}

pred_prob_fun <- function(object, y, x, offset, stddevs, sigma) {
  Xcat <- make_Xcat(object, y, x)
  betatilde <- bdiag(object$rho$constraints) %*% object$beta
  pred.upper  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_u <- c(object$theta[[j]], object$rho$inf.value)[y[, j]]
    xbeta_u <- as.double(Xcat$U[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_u - xbeta_u - offset[[j]]
  })/stddevs
  pred.lower  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_l <- c(-object$rho$inf.value, object$theta[[j]])[y[, j]]
    xbeta_l <- as.double(Xcat$L[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_l - xbeta_l - offset[[j]]
  })/stddevs
  pred.upper[is.na(pred.upper)] <- object$rho$inf.value
  pred.lower[is.na(pred.lower)] <- -object$rho$inf.value
  prob <- object$rho$link$F_multi(
    U = pred.upper, L = pred.lower,
    list_R = lapply(sigma, cov2cor))
  names(prob) <- rownames(y)
  return(prob)
}

pred_cumprob_fun <-  function(object, y, x, offset, stddevs, sigma) {
  Xcat <- make_Xcat(object, y, x)
  betatilde <- bdiag(object$rho$constraints) %*% object$beta
  pred.upper  <- sapply(seq_len(object$rho$ndim), function(j) {
    th_u <- c(object$theta[[j]], object$rho$inf.value)[y[, j]]
    xbeta_u <- as.double(Xcat$U[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
    th_u - xbeta_u - offset[[j]]
  })/stddevs
  pred.upper[is.na(pred.upper)] <- object$rho$inf.value
  pred.lower <- matrix(-object$rho$inf.value, ncol = object$rho$ndim,
                       nrow = nrow(pred.upper))
  cum.prob <- object$rho$link$F_multi(U = pred.upper, L = pred.lower,
                                      list_R = lapply(sigma, cov2cor))
  names(cum.prob) <- rownames(y)
  return(cum.prob)
}

pred_class_fun <-  function(object, y, x, offset, stddevs, sigma) {
  ndim <- object$rho$ndim
  betatilde <- bdiag(object$rho$constraints) %*% object$beta

  if (prod(object$rho$ncat) > 1e6) {
    stop("Number of class combinations over 1000000. Try joint_probabilities() for desired class combinations.")
  } else {
    cmbn <- expand.grid(lapply(object$rho$ncat, seq_len))
    cmbn.labels <- expand.grid(object$rho$levels)
    probs <- sapply(seq_len(nrow(cmbn)), function(i){
      print(i)
      if (i %% 100 == 0)  cat('Computed probabilities for', i, 'out of',
                              nrow(cmbn),'combinations\n')
      ytmp <- sapply(seq_len(ndim),
                     function(j) object$rho$levels[[j]][cmbn[i,j]])
      ytmp <- matrix(ytmp, ncol = ndim, nrow = nrow(y), byrow = TRUE)
      # ytmp <-cbind.data.frame(lapply(seq_len(object$rho$ndim), function(j){
      #   if (!all(ytmp[,j] %in% c(NA, levels(y[,j]))))  stop("response.cat are different from the categories in the original data set")
      #  else ordered(ytmp[,j], levels = object$rho$levels[[j]])
      # }))
      ytmp <- cbind.data.frame(
        lapply(seq_len(object$rho$ndim), function(j) {
          ordered(ytmp[,j], levels = object$rho$levels[[j]])
      }))
      Xcat <- make_Xcat(object, ytmp, x)
      pred.upper  <- sapply(seq_len(ndim), function(j) {
        th_u <- c(object$theta[[j]], object$rho$inf.value)[ytmp[, j]]
        xbeta_u <- as.double(Xcat$U[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
        th_u - xbeta_u - offset[[j]]
      })/stddevs
      pred.lower  <- sapply(seq_len(object$rho$ndim), function(j) {
        th_l <- c(-object$rho$inf.value, object$theta[[j]])[ytmp[, j]]
        xbeta_l <- as.double(Xcat$L[[j]] %*% betatilde[object$rho$indjbeta[[j]]])
        th_l - xbeta_l - offset[[j]]
      })/stddevs
      pred.lower[is.na(pred.lower)] <- -object$rho$inf.value
      pred.upper[is.na(pred.upper)] <- object$rho$inf.value
      object$rho$link$F_multi(U = pred.upper, L = pred.lower,
                              list_R = lapply(sigma, cov2cor))
    })
    ind.max <- apply(probs,1,which.max)
    class <- cmbn.labels[ind.max,]
    rownames(class) <- rownames(y)
    colnames(class) <- object$rho$y.names
    return(class)
  }
}
