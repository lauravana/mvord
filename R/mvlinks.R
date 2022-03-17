#' @title Multivariate link functions in mvord
#' @description Different \code{link} functions are available in \pkg{mvord}:
#' @details We allow for two different link functions, the multivariate probit
#' link and the multivariate logit link.
#' For the multivariate probit
#' link a multivariate normal distribution for the errors is applied. The
#' normal bivariate probabilities which enter the pairwise log-likelihood
#' are computed with the package \pkg{pbivnorm}.
#'
#' For the multivariate logit link a \eqn{t} copula based multivariate
#' distribution with logistic margins is used.
#'   The \code{mvlogit()} function has an optional integer valued argument
#' \code{df} which specifies the degrees of freedom to be used for the
#' \eqn{t} copula.  The default value of the degrees of freedom parameter is
#' 8. We restrict the degrees of freedom to be integer valued because the
#' most efficient routines for computing bivariate \eqn{t} probabilities do
#' not support non-integer degrees of freedom. For further details see vignette.
#' @param df integer specifying the degrees of freedom of the t copula
#' @return The functions \code{mvlogit()} and \code{mvprobit()} returns an object
#' of \code{\link{class}} \code{'mvlink'}.
#' An object of \code{\link{class}} \code{'mvlink'} is a list containing the following components:
#'
#' \itemize{
#'  \item{\code{name}}{
#'
#'  name of the multivariate link function}
#'   \item{\code{df}}{
#'
#'   degrees of freedom of the t copula; returned only for \code{mvlogit()}
#'   }
#'  \item{\code{F_uni}}{
#'
#'   a function corresponding to the univariate margins of the
#'   multivariate distribution \eqn{F}  of the subject errors; the function returns \eqn{Pr(X \leq x) = F_1(x)}
#'   }
#'   \item{\code{F_biv}}{
#'
#'   a function corresponding to the bivariate distribution of the
#'   multivariate distribution \eqn{F} of the subject errors \eqn{Pr(X \leq x, Y\leq y|r) = F_2(x, y, r)};
#'   }
#'   \item{\code{F_biv_rect}}{
#'
#'   the function computes the rectangle probabilities from based on \code{F_biv};
#'   the function has the matrices \code{U} (upper bounds) and \code{L} (lower bounds)
#'   as well as vector \code{r} containing the correlation coefficients
#'   corresponding to the bivariate distribution as arguments; the matrices
#'   \code{U} and \code{L} both have two columns, first corresponding to the bounds of x,
#'   second to the bounds of y; the number of rows corresponds to the number of observations;
#'   the rectangle probabilities are defined as
#'   \eqn{Pr(L[,1]\leq X\leq U[,1], L[,2]\leq Y \leq U[,2]|r) = F_2(U[,1], U[,2],r) - F_2(U[,1], L[,2],r)- F_2(L[,1], U[,2],r) + F_2(L[,1], L[,2],r)}
#'   }
#'   \item{\code{F_multi}}{
#'
#'   the function computes the multivariate probabilities for distribution function \eqn{F};
#'   the function has the matrices \code{U} (upper bounds) and \code{L} (lower bounds)
#'   as well as the list \code{list_R} containing for each observation the correlation matrix;
#'   F is needed for the computation of the fitted/predicted joint probabilities.  If NULL only marginal probabilities can be computed.
#'   }
#'   \item{\code{deriv.fun}}{
#'
#'   (needed for computation of analytic standard errors) a list containing the following gradient functions:
#'         \itemize{
#'            \item{\code{dF1dx}}{ derivative \eqn{dF_1(x)/dx} function,}
#'            \item{\code{dF2dx}}{ derivative \eqn{dF_2(x,y,r)/dx} function,}
#'            \item{\code{dF2dr}}{ derivative \eqn{dF_2(x,y,r)/dr } function.}
#'          }
#'   If \code{deriv.fun = NULL} numeric standard errors will be computed.
#'   }
#' }
#' @name mvlinks
#' @export
mvprobit <- function() {
  mvlink <- list(name = "mvprobit",
                 F_uni = pnorm,  # F_1(x)
                 F_biv = pbivnorm, # F_2(x, y, r)
                 F_biv_rect = rectbiv_norm_prob,  # Pr(L_x <= X <= U_x, L_y <= Y <= U_y|r)
                 F_multi = function(U, L, list_R) {
                   sapply(1:nrow(U), function(i) sadmvn(lower = L[i, ],
                                                        upper = U[i, ],
                                                        mean = rep.int(0, NCOL(U)),
                                                        varcov = list_R[[i]]))
                 },
                 deriv.fun = list(
                   dF1dx = dnorm, # dF_1(x)/dx
                   dF2dx = function(x, y, r) dnorm(x) * pnorm((y - r * x)/sqrt(1 - r^2)),
                   dF2dr = function(x, y, r){
                     1/(2 * pi * sqrt(1 - r^2)) *
                       exp(-(x^2 - 2 * r * x * y + y^2)/(2 * (1 - r^2)))})
  )
  class(mvlink) <- "mvlink"
  return(mvlink)
}
#' @export
#' @rdname mvlinks
mvlogit <- function(df = 8L){
  if(!is.integer(df)) {
    warning("Degrees of freedom in mvlogit() need to be integer.
            Rounding to the closest integer.")
    df <- as.integer(df)
  }
  mvlink <- list(name = "mvlogit",
                 df = df,
                 sqrv = pi/sqrt(3),
                 F_uni = plogis, # F_1(x)
                 F_biv = function(x, y, r) {
                  #  x <- qt(plogis(x), df = df)
                  #  y <- qt(plogis(y), df = df)
                   unlist(lapply(seq_along(x), function(i)
                     biv_nt_prob2(nu = df,
                                  lower = c(-10000, -10000),
                                  upper = c(x[i], y[i]),
                                  r     = r[i])))},
                 F_biv_rect = function(U, L, r) {
                   # U <- qt(plogis(U), df = df)
                   # L <- qt(plogis(L), df = df)
                   #L[is.infinite(L)] <- -10000
                   #U[is.infinite(U)] <- 10000
                   unlist(lapply(seq_len(nrow(U)), function(i)
                     biv_nt_prob2(nu = df,
                                  lower = L[i, ],
                                  upper = U[i, ],
                                  r     = r[i])))},
                 F_multi = function(U, L, list_R) {
                  # U <- qt(plogis(U), df = df)
                  # L <- qt(plogis(L), df = df)
                   unlist(lapply(seq_len(nrow(U)), function(i) sadmvt(df = df,
                                                               lower = L[i, ],
                                                               upper = U[i, ],
                                                               mean = rep.int(0, NCOL(U)),
                                                               S = list_R[[i]])))},

                 deriv.fun =  list(
                   dF1dx = dlogis,
                   dF2dx = function(x, y, r) deriv_biv_t_copula(x, y, r, df = df) ,
                   dF2dr = function(x, y, r) deriv_corr_t_copula(x, y, r, df = df))
  )
  class(mvlink) <- "mvlink"
  return(mvlink)
}
# #' @title Print Method for class \code{"mvlink"}
# #' @description Prints name of the multivariate link
# #' @param object object of class \code{"mvlink"}
# #' @noRd
# # #' @rdname mvord
# #' @method print mvlink
# #' @export
#print.mvlink <- function(object) sprintf("Multivariate link: %s", object$name)

rectbiv_norm_prob <- function(U, L, r) {
  # computes the rectangle probabilities for biv.normal-distribution
  p1 <- pbivnorm(U[, 1], U[, 2], r)
  p2 <- pbivnorm(L[, 1], U[, 2], r)
  p3 <- pbivnorm(U[, 1], L[, 2], r)
  p4 <- pbivnorm(L[, 1], L[, 2], r)
  ## replace NaN
  p1[is.nan(p1)] <- 0
  p2[is.nan(p2)] <- 0
  p3[is.nan(p3)] <- 0
  p4[is.nan(p4)] <- 0
  pr <- p1 - p2 - p3 + p4
  return(pr)
}

biv_nt_prob2 <-function (nu, lower, upper, mean = c(0, 0), r) {
  # computes the rectangle probabilities for biv.t-distribution
  #nu <- df
  rho <- as.double(r)
  infin <- c(2, 2)
  infin <- as.integer(infin)
  prob <- as.double(0)
  a <- .Fortran("smvbvt", prob, nu, lower, upper, infin, rho,
                PACKAGE = "mvord")
  return(a[[1]])
}

deriv_biv_t_copula <- function(x, y, r, df){
  ## add transformation
  newx <- qt(plogis(x), df = df)
  newy <- qt(plogis(y), df = df)
  inf.value <- 10000#sqrt(.Machine$double.xmax)/2
  newx <- replace(newx, newx == Inf, inf.value)
  newx <- replace(newx, newx == - Inf, -inf.value)
  newy <- replace(newy, newy == Inf, inf.value)
  newy <- replace(newy, newy == - Inf, -inf.value)
  ## conditional parameters
  mu_c <- r * newx
  sigma_c <- sqrt((df + newx^2)/(df + 1) * (1 - r^2))
  df_c <- df + 1
  dlogis(x) * pt((newy - mu_c)/sigma_c, df = df_c)
}

deriv_corr_t_copula <- function(x, y, r, df){
  ## add transformation
  newx <- qt(plogis(x), df = df)
  newy <- qt(plogis(y), df = df)
  inf.value <- 10000#sqrt(.Machine$double.xmax)/2
  newx <- replace(newx, newx == Inf, inf.value)
  newx <- replace(newx, newx == - Inf, -inf.value)
  newy <- replace(newy, newy == Inf, inf.value)
  newy <- replace(newy, newy == - Inf, -inf.value)
  1/(2 * pi * sqrt(1 - r^2)) *
    (1 + (newx^2 - 2 * r * newx * newy + newy^2)/(df * (1 - r^2)))^(- df/2)
}
