initialize_MMO <- function(rho, formula, data, error.structure, contrasts){
  rho$function.name <- "mvord"
  rho$formula <- as.formula(paste(formula[[2]][[2]], paste(as.character(formula[-2]), collapse = " ")))
  if (length(formula[[2]]) == 2){
    rho$index <- colnames(data)[1:2]
  } else rho$index <- c(as.character(formula[[2]][[3]]), as.character(formula[[2]][[4]]))
  rho$response.name <- all.vars(rho$formula[[2]])
  rho$response.names <- levels(ordered(data[, rho$index[2]]))
  ## checks specific to mvord
  if(length(rho$response.name) > 1) stop("only one response needed", call.=FALSE)
  if (any(is.na(data[, rho$index[1]]))) stop("Missing values are not allowed in the subject index.", call.=FALSE)
  if (any(is.na(data[, rho$index[2]]))) stop("Missing values are not allowed in the measurement index.", call.=FALSE)
  if (!is.null(rho$response.levels) & length(rho$response.levels) != length(rho$response.names))
      stop("Length of response levels must be equal to the number of responses.", call.=FALSE)

  check_args_input1(rho, data) ## used also for mvord2

  if(!is.null(rho$weights.name)){
    if(any(is.na(data[,rho$weights.name]))) {
      data[,rho$weights.name][is.na(data[,rho$weights.name])] <- 0
      warning("Weights with values of NA are set to 0.")
    }
  }
  ## check if more than one response --

  rho$intercept <- ifelse(attr(terms.formula(rho$formula),
                               "intercept") == 1, TRUE, FALSE)
  rho$x.names <- c(all.vars(rho$formula[[3]]),
                   all.vars(formula(error.structure)[[2]]))
  if (any(is.na(data[, rho$x.names]))) stop("Missing values in the covariates are not allowed.")

  names.rest <- colnames(data)
  data.mvord <- mvord_data(data, index = rho$index, y.names = rho$response.name, x.names = names.rest,
                           y.levels = rho$response.levels,
                           response.names = rho$response.names)

  rho$levels <- data.mvord$ylevels

  rho$y <- data.mvord$y
  rho$y.names <- colnames(rho$y)

  rho$ndim <- ncol(rho$y)
  rho$n <- nrow(rho$y)
  rho$x <- lapply(seq_len(rho$ndim), function(j) {
    rhs.form <- rho$formula
    rhs.form[[2]] <- NULL
    new.rhs.form <- update(rhs.form, ~  . + 1)

    tmp <- suppressWarnings(model.matrix(new.rhs.form,
              model.frame(new.rhs.form,  data.mvord$x[[j]], na.action = function(x) x),
              contrasts.arg = contrasts))
    attribute <- attr(tmp, "assign")
    if(rho$intercept == FALSE){
      attribute <- attribute[-1]
      tmp <- tmp[,-1, drop = F]
    }
    attr(tmp, "assign") <- attribute
    tmp
  })

  if (is.null(rho$weights.name)) {
    rho$weights <- rep(1, nrow(rho$y))
  } else {
    tmp <- sapply(data.mvord$x, function(j) as.numeric(j[,rho$weights.name]))
    rho$weights <- apply(tmp,1,function(x) unique(x[!is.na(x)]))
    if(is.list(rho$weights)) stop("Weights need to be constant across multiple measurements", call.=FALSE)
    if(any(rho$weights < 0)) stop("Weights must be non-negative", call.=FALSE)
  }
  ## initialize error structure
  rho$error.structure <- init_fun(error.structure, data.mvord, contrasts)
  ## set offset
  if (is.null(rho$offset)) {
    rho$offset <- lapply(seq_len(rho$ndim), function(j) {
      rhs.form <- rho$formula
      rhs.form[[2]] <- NULL
      mf <- model.frame(rhs.form, data.mvord$x[[j]],
                        na.action = function(x) x)
      #mf[is.na(mf)] <- 0
      if(NCOL(mf) > 0){
        for(k in seq_along(NCOL(mf))){
          if(is.numeric(mf[,k]))
            mf[is.na(mf[,k]),k] <- 0
        }
      }
      model.offset(mf)
    })
  }
  check_args_input2(rho, data)
  rho
}

initialize_MMO2 <- function(rho, formula, data, error.structure, contrasts){
  rho$function.name <- "mvord2"
  rho$formula <- as.formula(paste0("cbind(", paste(formula[[2]][-1], collapse = ", "), ") ", paste(as.character(formula[-2]), collapse = " ")))
  rho$y.names <- all.vars(rho$formula[[2]])
  rho$y <- data[,rho$y.names]
  rho$ndim <- ncol(rho$y)
  rho$n <- nrow(rho$y)

  check_args_input1(rho, data)
  if(!is.null(rho$weights.name)){
    if(any(is.na(data[,rho$weights.name]))) {
      data[,rho$weights.name][is.na(data[,rho$weights.name])] <- 0
      warning("Weights with values of NA are set to 0.")
    }
  }

  rho$x.names <- c(all.vars(rho$formula[[3]]), all.vars(formula(error.structure)[[2]]))
  if (any(is.na(data[, rho$x.names]))) stop("Missing values in the covariates are not allowed.")

  if(!all(sapply(rho$y, is.ordered))){ warning("Responses are unordered. Natural ordering is used.", call.=FALSE)
    rho$y <- cbind.data.frame(lapply(rho$y, ordered))
    colnames(rho$y) <- rho$y.names
  }

  rho$levels <- lapply(rho$y, levels)
  for (j in seq_len(rho$ndim)){
    if (!all(rho$levels[[j]] %in% unique(rho$y[, j])))
       warning(sprintf("For response %i, not all response
          levels are observed. Model might be non-identifiable if
          the thresholds for this response are not restricted.", j),
        call.=FALSE)
  }

  #rho$rownames.constraints <- unlist(rho$levels)
  rho$intercept <- ifelse(attr(terms.formula(rho$formula), "intercept") == 1, TRUE, FALSE)
  rho$x <- lapply(seq_len(rho$ndim), function(j) {
    rhs.form <- as.formula(paste(as.character(formula[-2]), collapse = " "))
    new.rhs.form <- update(rhs.form, ~ . + 1)
    tmp <-  suppressWarnings(model.matrix(new.rhs.form,
                  model.frame(new.rhs.form,  data, na.action = function(x) x),
                  contrasts.arg = contrasts))
    attribute <- attr(tmp, "assign")
    if(rho$intercept == FALSE){
      attribute <- attribute[-1]
      tmp <- tmp[,-1, drop = F]
    }
    attr(tmp, "assign") <- attribute
    tmp
  })

  if (is.null(rho$offset)) {
    rho$offset <- lapply(seq_len(rho$ndim), function(j) {
      mf <- model.frame(rho$formula, data, na.action = function(x) x)
      model.offset(mf)
    })
  }

  if (is.null(rho$weights.name)) {
    rho$weights <- rep(1, nrow(rho$y))
  } else {
    rho$weights <- as.vector(data[, rho$weights.name])
    if(any(rho$weights < 0)) stop("Weights must be non-negative", call.=FALSE)
  }
  data.mvord2 <- list(y = data[, rho$y.names],
                      x = rep(list(data[, rho$x.names, drop = F]), rho$ndim))

  rho$error.structure <- init_fun(error.structure, data.mvord2, contrasts)
  check_args_input2(rho, data)
  rho
}
