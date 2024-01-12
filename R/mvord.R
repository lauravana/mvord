#' @rdname mvord-package
#' @useDynLib mvord
#' @title Multivariate Ordinal Regression Models in R.
#' @description  The R package mvord implements composite likelihood
#' estimation in the class of multivariate ordinal regression models with probit and logit link.
#' A flexible modeling framework for multiple ordinal measurements on
#' the same subject is set up, which takes into consideration the
#' dependence among the multiple observations by employing different
#' error structures.
#' Heterogeneity in the error structure across the subjects can be
#' accounted for by the package, which allows for covariate dependent
#' error structures.
#' In addition, regression coefficients and threshold parameters are
#' varying across the multiple response dimensions in the default
#' implementation. However, constraints can be defined by the user if a
#' reduction of the parameter space is desired.
#' @details see \code{\link{mvord}}
#' @name mvord-package
#' @aliases mvord-package
#' @references Hirk R, Hornik K, Vana L (2020). “\pkg{mvord}: An \strong{R} Package for Fitting Multivariate Ordinal Regression Models.” \emph{Journal of Statistical Software}, \strong{93}(4), 1–41, \doi{10.18637/jss.v093.i04}.
NULL

#' @title Simulated panel of credit ratings
#'
#' @description A data set containing simulated credit ratings assigned by one rater and simulated perfomance measures for firms in different years.
#'
#' \itemize{
#'   \item \code{rating} credit ratings
#'   \item \code{firm_id} firm index
#'   \item \code{year} year index
# #'   \item \code{ICR} interest coverage ratio, which measures how well the interest expenses can be
# #' covered from the free operating cash-flow of a company
#'   \item \code{LR} liquidity ratio, relating the cash held by a company to the current liabilities
#'   \item \code{LEV} leverage ratio relating debt to earnings before interest and taxes
# #'   \item \code{LEV2} leverage ratio measuring the percentage of debt in the long-term capital of a firm
#'   \item \code{PR} profitability ratio of retained earnings to assets
#'   \item \code{RSIZE} log of relative size of the company in the market
#'   \item \code{BETA} a measure of systematic risk
#'   \item \code{BSEC} business sector of a firm (factor with 8 levels)
#' }
#'
#' @name data_cr_panel
#' @docType data
# #' @keywords datasets
#' @usage data("data_cr_panel", package = "mvord")
#' @format A data frame with 11320 rows and 9 variables
NULL


#' Simulated credit ratings
#'
#' A data set containing simulated credit ratings and simulated perfomance measures from four raters.
#'
#' \itemize{
#'   \item \code{rater1} credit ratings assigned by rater 1
#'   \item \code{rater2} credit ratings assigned by rater 2
#'   \item \code{rater3} credit ratings assigned by rater 3
#'   \item \code{rater4} credit ratings assigned by rater 4
#'   \item \code{firm_id} firm index
##'   \item \code{rater_id} rater index
# #'   \item \code{ICR} interest coverage ratio, which measures how well the interest expenses can be
#' covered from the free operating cash-flow of a company
#'   \item \code{LR} liquidity ratio, relating the cash held by a company to the current liabilities
#'   \item \code{LEV} leverage ratio relating debt to earnings before interest and taxes
# #'   \item \code{LEV2} leverage ratio measuring the percentage of debt in the long-term capital of a firm
#'   \item \code{PR} profitability ratio of retained earnoings to assets
#'   \item \code{RSIZE} log of relative size of the company in the market
#'   \item \code{BETA} a measure of systematic risk
# #'   \item \code{BSEC} business sector of a firm (factor with 8 levels)
#' }
#'
#' @name data_cr
#' @docType data
# #' @keywords datasets
#' @usage data("data_cr", package = "mvord")
#' @format A data frame with 690 rows and 11 columns
NULL

#' Simulated credit ratings
#'
#' A simulated data set where three different raters (\code{rater1, rater2} and \code{rater3})
#'  assign ordinal ratings on different firms. \code{rater3} uses a different rating scale
#'   compared to \code{rater1} and \code{rater2}. The IDs for each subject \eqn{i} of the \eqn{n = 1000}
#' firms are stored in the column \code{firm_id}.
#'
#' \itemize{
#'   \item \code{firm_id} firm index
#'   \item \code{rater1}  ordinal rating outcome of rater 1
#'   \item \code{rater2}  ordinal rating outcome of rater 2
#'   \item \code{rater3}  ordinal rating outcome of rater 3
#'   \item \code{X1} covariate X1
#'   \item \code{X2} covariate X2
#'   \item \code{X3} covariate X3
#'   \item \code{X4} covariate X4
#'   \item \code{X5} covariate X5
#'   \item \code{X6} covariate X6 (factor)
#' }
#' @name data_mvord2
#' @docType data
# #' @keywords datasets
#' @usage data("data_mvord2", package = "mvord")
#' @format A data frame with 1000 rows and 10 variables
NULL


#' Simulated panel of credit ratings
#'
#' A simulated data set where one rater assigns ratings over the years \eqn{2001} to \eqn{2010} for a set of firms.
#' The IDs for each subject \eqn{i} of the \eqn{n = 1000} firms are stored in the column \code{firm_id}.
#' The year of the rating observation is stored in the column \code{year}.
#' The ordinal ratings are provided in the column \code{rating} and all the covariates in the remaining columns.
#'
#' \itemize{
#'   \item \code{firm_id} firm index
#'   \item \code{year}  year index (2001 - 2010)
#'   \item \code{rating} ordinal credit ratings
#'   \item \code{X1} covariate X1
#'   \item \code{X2} covariate X2
#'   \item \code{X3} covariate X3
#'   \item \code{X4} covariate X4
#'   \item \code{X5} covariate X5
#'   \item \code{X6} covariate X6 (factor)
#' }
#' @name data_mvord_panel
#' @docType data
# #' @keywords datasets
#' @usage data("data_mvord_panel", package = "mvord")
#' @format A data frame with 10000 rows and 9 variables
NULL


#' Simulated credit ratings
#'
#' A simulated data set where three different raters (\code{rater1, rater2} and \code{rater3})
#' assign ordinal ratings on different firms. \code{rater3} uses a different rating scale
#' compared to \code{rater1} and \code{rater2}, i.e., the number of threshold categories is different.
#' For each firm we simulate five different covariates \code{X1, ..., X5} from a standard
#' normal distribution. Additionally, each firm is randomly assigned to a business sector (sector \code{X}, \code{Y} or \code{Z}), captured by the covariate \code{X6}. Furthermore, we simulate
#'  multivariate normally distributed errors. For a given set of parameters we obtain the three rating variables for
#'  each firm by slotting the latent scores according to the corresponding threshold parameters.
#' The IDs for each subject \eqn{i} of the \eqn{n = 1000} firms are stored in the column \code{firm_id}. The IDs of the raters are stored
#' in the column \code{rater_id}. The ordinal ratings are provided in the column \code{rating} and all the covariates in the remaining columns.
#' Overall, the data set has 3000 rows, for each of the \eqn{n = 1000} firms it has three rating observations.
#'
#' \itemize{
#'   \item \code{firm_id} firm index
#'   \item \code{rater_id}  rater index
#'   \item \code{rating} ordinal credit ratings
#'   \item \code{X1} covariate X1
#'   \item \code{X2} covariate X2
#'   \item \code{X3} covariate X3
#'   \item \code{X4} covariate X4
#'   \item \code{X5} covariate X5
#'   \item \code{X6} covariate X6 (factor)
#' }
#'
#' @name data_mvord
#' @docType data
# #' @keywords datasets
#' @usage data("data_mvord", package = "mvord")
#' @format A data frame with 3000 rows and 9 variables
NULL

#' Data set toy example
#'
#' A data set containing two simulated ordinal responses with three categories,
#' two quantitative covariates \code{X1} and \code{X2} and two categorical covariates
#' \code{f1} and \code{f2}.
#'
#' \itemize{
#'   \item \code{Y1} ordinal outcome \code{Y1} (three categories)
#'   \item \code{Y2} ordinal outcome \code{Y2} (three categories)
#'   \item \code{X1} covariate \code{X1}
#'   \item \code{X2} covariate \code{X2}
#'   \item \code{f1} categorical covariate \code{f1}
#'   \item \code{f2} categorical covariate \code{f2}
#' }
#'
#' @name data_mvord_toy
#' @docType data
# #' @keywords datasets
#' @usage data("data_mvord_toy", package = "mvord")
#' @format A data frame with 100 rows and 6 variables
NULL

#' Essay data
#'
#' The multirater agreement data set is taken from Chapter 5 in "Ordinal Data Modeling",
#' from Johnson, Valen E and Albert, J. The data consists of grades assigned to
#' 198 essays by 5 experts, each of whom rated all essays on a 10-point
#' scale. A score of 10 indicates an excellent essay. In addition, the
#' average word length is also available as an essay characteristic.
#'
#' \itemize{
#'   \item \code{Judge1} ordinal outcome:  grades assigned by expert 1
#'   \item \code{Judge2} ordinal outcome:  grades assigned by expert 2
#'   \item \code{Judge3} ordinal outcome:  grades assigned by expert 3
#'   \item \code{Judge4} ordinal outcome:  grades assigned by expert 4
#'   \item \code{Judge5} ordinal outcome:  grades assigned by expert 5
#'   \item \code{wl} covariate: word length
#' }
#'
#' @name essay_data
#' @docType data
# #' @keywords datasets
#' @usage data("essay_data", package = "mvord")
#' @format A data frame with 198 rows and 6 variables
NULL

##IMPORTS
#' @importFrom optimx optimx
#' @importFrom stats formula update model.frame model.matrix coef as.formula cov2cor dnorm dt pnorm pt terms.formula plogis dlogis qt nobs predict terms printCoefmat model.offset na.omit logLik binomial glm.fit sd na.pass
#' @importFrom pbivnorm pbivnorm
#' @importFrom MASS polr
#' @importFrom utils combn data write.table
#' @importFrom mnormt sadmvn
#' @importFrom mnormt sadmvt
#' @importFrom mvtnorm pmvt
#' @importFrom Matrix bdiag
#' @importFrom numDeriv grad hessian
#' @import minqa
#' @import BB
#' @import dfoptim
#' @import ucminf


#############################################################################################
#' @title Multivariate Ordinal Regression Models.
#'
#' @description
#' Multivariate ordinal regression models in the R package  \code{mvord} can be fitted using the function
#' \code{mvord()}. Two different data structures can be passed on to \code{mvord()} through
#' the use of two different multiple measurement objects \code{MMO} and \code{MMO2} in the left-hand side of
#' the model formula. \code{MMO} uses a long data format, which has the advantage that it allows for
#' varying covariates across multiple measurements. This flexibility requires the specification a
#' subject index as well as a multiple measurement index. In contrast to \code{MMO}, the function \code{MMO2}
#' has a simplified data structure, but is only applicable in settings where the covariates do not
#' vary between the multiple measurements. In this case, the multiple ordinal observations as
#' well as the covariates are stored in different columns of a \code{\link{data.frame}}. We refer to this data
#' structure as wide data format.
#' @details
#' \describe{
#' \item{Implementation \code{MMO}:}{
#'   \describe{
#'     \item{\code{data}:}{
#' In \code{MMO} we use a long format for the input of data, where each row contains a subject index
#' (\code{i}), a multiple measurement index (\code{j}), an ordinal
#' observation (Y) and all the covariates (X1 to Xp). This long format data structure is
#' internally transformed to a matrix of responses which contains NA in the case of missing
#' entries and a list of covariate matrices. This is performed by the multiple measurement object
#' \code{MMO(Y, i, j)}
#' specifying the column names of the subject index and the multiple measurement index in data.
#' The column containing the ordinal observations can contain integer or character values or can
#' be of class (ordered) 'factor'. When using the long data structure, this column is basically
#' a concatenated vector of each of the multiple ordinal responses. Internally, this vector is
#' then split according to the measurement index. Then the ordinal variable corresponding to
#' each measurement index is transformed into an ordered factor. For an integer or a character
#' vector the natural ordering is used (ascending, or alphabetical). If for character vectors the
#' alphabetical order does not correspond to the ordering of the categories, the optional argument
#' response.levels allows to specify the levels for each response explicitly. This is performed
#' by a list of length q, where each element contains the names of the levels of the ordered
#' categories in ascending (or if desired descending) order. If all the multiple measurements use
#' the same number of classes and same labelling of the classes, the column Y can be stored as
#' an ordered 'factor' (as it is often the case in longitudinal studies).
#' The order of the multiple measurements is needed when specifying constraints on the thresh-
#' old or regression parameters. This order is based on the type of the
#' multiple measurement index column in data. For 'integer', 'character' or 'factor' the
#' natural ordering is used (ascending, or alphabetical). If a different order of the multiple responses is desired,
#' the multiple measurement index column should be an ordered factor with
#' a corresponding ordering of the levels.
#'
#' If the categories differ across multiple measurements (either the number of categories or the category labels)
#' one needs to specify the \code{response.levels} explicitly. This is performed by a list
#' of length \eqn{J} (number of multiple measurements), where each element contains
#' the names of the levels of the ordered categories in ascending or descending order.}
#' \preformatted{response.levels = list(c("G","F","E", "D", "C", "B", "A"),
#'                        c("G","F","E", "D", "C", "B", "A"),
#'                        c("O","N","M","L", "K", "J", "I", "H"))}
#'
#' \item{\code{formula}}{
#' The ordinal responses (e.g., \code{rating}) are passed by a \code{formula} object.
#' Intercepts can be included or excluded in the model depending on the model paramterization:
#' \describe{
#' \item{Model without intercept:}{ If the intercept should be removed the \code{formula} for a given response (\code{rating})
#' and covariates (\code{X1} to \code{Xp}) has the following form:
#'
#'      \code{formula = MMO(rating, firm_id, rater_id) ~ 0 + X1 + ... + Xp}.
#' }
#' \item{Model with intercept:}{ If one wants to include an intercept in the model, there are two equivalent possibilities
#' to set the model \code{formula}. Either one includes the intercept explicitly by:
#'
#'     \code{formula = MMO(rating, firm_id, rater_id) ~ 1 + X1 + ... + Xp},
#'
#' or by
#'
#'   \code{formula = MMO(rating, firm_id, rater_id) ~ X1 + ... + Xp}.
#' }
#' }
#' }
#' }
#' }
#' \item{Implementation \code{MMO2}:}{
#'   \describe{
#'     \item{\code{data}:}{The data structure applied by \code{MMO2} is slightly simplified, where the multiple ordinal
#' observations as well as the covariates are stored as columns in a \code{\link{data.frame}}. Each subject \eqn{i}
#' corresponds to one row of the data frame, where all outcomes (with missing
#' observations set to NA) and all the covariates are stored in different columns.
#' Ideally each outcome column is of type ordered factor. For column types like 'integer',
#' 'character' or 'factor' a warning is given and the natural ordering is used (ascending, or
#' alphabetical).}
#'
#' \item{\code{formula}}{
#' The ordinal responses (e.g., \code{rating}) are passed by a \code{formula} object.
#' Intercepts can be included or excluded in the model depending on the model parameterization:
#'
#'   \code{formula = MMO2(rater1, rater2, rater3) ~ X1 + ... + Xp}.
#' }
#' }
#' }
#'   \item{\code{error.structure}}{
#'  We allow for different error structures depending on the model parameterization:
#'\itemize{
#'   \item {Correlation:}
#'   \itemize{
#'   \item \code{cor_general}
#' The most common parameterization is the general correlation matrix.
#'
#'  \code{error.structure = cor_general(~ 1)}
#'
#' This parameterization can be extended by allowing a factor dependent
#' correlation structure, where the correlation of each subject \eqn{i} depends
#' on a given subject-specific factor \code{f}. This factor \code{f} is not allowed to vary
#' across multiple measurements \eqn{j} for the same subject \eqn{i} and due to numerical
#' constraints only up to maximum 30 levels are allowed.
#'
#'       \code{error.structure = cor_general(~ f)}
#'
#'   \item \code{cor_equi}
#' A covariate dependent equicorrelation structure, where the correlations
#' are equal across all \eqn{J} dimensions and depend on subject-specific covariates \code{S1, ..., Sm}.
#' It has to be noted that these covariates \code{S1, ..., Sm} are not allowed to vary across
#'  multiple measurements \eqn{j} for the same subject \eqn{i}.
#'
#'          \code{error.structure = cor_equi(~ S1 + ... + Sm)}
#'
#'   \item \code{cor_ar1}
#' In order to account for some heterogeneity the \eqn{AR(1)} error structure
#' is allowed to depend on covariates \code{X1, ..., Xp} that are constant
#' over time for each subject \eqn{i}.
#'
#'       \code{error.structure = cor_ar1(~ S1 + ... + Sm)}
#'}
#'
#'\item {Covariance:}
#'\itemize{
#'\item \code{cov_general}
#'
#' In case of a full variance-covariance parameterization the standard parameterization
#'  with a full variance-covariance is obtained by:
#'
#'  \code{error.structure = cov_general(~ 1)}
#'
#'  This parameterization can be extended to the factor dependent covariance structure,
#'   where the covariance of each subject depends on a given factor \code{f}:
#'
#'  \code{error.structure = cov_general(~ f)}
#'   }
#'   }
#'   }
#'
#'   \item{\code{coef.constraints}}{
#'   The package supports
#'   constraints on the regression coefficients. Firstly, the
#'   user can specify whether the regression coefficients should be equal
#'   across some or all response dimensions. Secondly, the values of some
#'   of the regression coefficients can be fixed.
#'
#'   As there is no unanimous way to specify such constraints, we offer
#'   two options. The first option is similar to the specification of constraints on the thresholds.
#'    The constraints can be specified in this case as a vector or matrix of integers,
#'     where coefficients getting same integer value are set equal.
#'   Values of the regression coefficients can be fixed through a matrix.
#'   Alternatively constraints on the regression coefficients can be specified
#'   by using the design employed by the \pkg{VGAM} package.
#'   The constraints in this setting are set through a named list,
#'   where each element of the list contains a matrix full-column rank.
#'   If the values of some regression coefficients should be fixed, offsets can be used.
#'   This design has the advantage that it supports
#'   constraints on outcome-specific as well as category-specific
#'   regression coefficients. While the first option has the advantage of requiring a more concise input,
#'    it does not support category-specific coefficients.
#'   The second option offers a more flexible design in this respect. For further information
#'   on the second option we refer to the vignette and to the documentation of \code{\link[VGAM]{vglm}}.
#'
#' Using the first option, constraints can be specified by a vector or a matrix \cr
#'    \code{coef.constraints}.
#'     First, a simple and less flexible way by specifying a vector \cr
#'     \code{coef.constraints}
#'      of dimension \eqn{J}.
#'      This vector is allocated in the following way:
#' The first element of the vector \code{coef.constraints} gets a value of 1. If the coefficients
#'  of the multiple measurement \eqn{j = 2} should be equal to the coefficients of the first dimension (\eqn{j=1}) again
#'   a value of 1 is set. If the coefficients should be different to the coefficients of the first dimension
#'   a value of 2 is set. In analogy, if the coefficients of dimensions two and three
#'    should be the same one sets both values to 2 and if they should be different,
#'     a value of 3 is set. Constraints on the regression coefficients of the remaining multiple measurements are set analogously.
#'
#'  \code{coef.constraints <- c(1,1,2,3)}
#'
#'  This vector \code{coef.constraints} sets the coefficients of the first two raters equal
#'  \deqn{\beta_{1\cdot} = \beta_{2\cdot}}
#'  A more flexible way to specify constraints on the regression coefficients is a matrix with \eqn{J} rows and \eqn{p} columns,
#'   where each column specifies constraints on one of the \eqn{p} coefficients in the same way as above.
#'    In addition, a value of \code{NA} excludes a corresponding coefficient (meaning it should be fixed to zero).
#'
#'    \preformatted{coef.constraints <- cbind(c(1,2,3,4), c(1,1,1,2), c(NA,NA,NA,1),
#'                           c(1,1,1,NA), c(1,2,3,4), c(1,2,3,4))}
#'
#'        This matrix \code{coef.constraints} gives the following constraints:
#'\itemize{
#'  \item \eqn{\beta_{12} = \beta_{22} = \beta_{32}}
#'    \item \eqn{\beta_{13} = 0}
#'    \item \eqn{\beta_{23} = 0}
#'    \item \eqn{\beta_{33} = 0}
#'    \item \eqn{\beta_{44} = 0}
#'    \item \eqn{\beta_{14} = \beta_{24} = \beta_{34}}
#'}
#'}
#'
#'
#'   \item{\code{coef.values}}{
#'   In addition, specific values on regression coefficients can be set in the matrix \cr
#'   \code{coef.values}.
#'    Parameters are removed if the value is set to zero (default for \code{NA}'s in \cr
#'    \code{coef.constraints})
#'     or to some fixed value. If constraints on parameters are set, these dimensions need to have
#'      the same value in \code{coef.values}. Again each column corresponds to one regression coefficient.
#'
#'  Together with the \code{coef.constraints} from above we impose:
#'
#'    \preformatted{coef.constraints <- cbind(c(1,2,2), c(1,1,2), c(NA,1,2),
#'                           c(NA,NA,NA), c(1,1,2))}
#'
#'  \preformatted{coef.values <- cbind(c(NA,NA,NA), c(NA,NA,NA), c(0,NA,NA),
#'                      c(1,1,1), c(NA,NA,NA))}
#' Interaction terms: When constraints on the regression coefficient should be specified in models with interaction terms,
#' the \code{coef.constraints} matrix has to be expanded manually. In case of interaction terms
#' (specified either by \code{X1 + X2 + X1:X2} or equivalently by \code{X1*X2}), one additional
#' column at the end of \code{coef.constraints} for the interaction term has to be specified for
#' numerical variables. For interaction terms including factor variables suitably more columns have
#' to be added to the \code{coef.constraints} matrix.
#' }
#'
#'
#'   \item{\code{threshold.constraints}}{
#'   Similarly, constraints on the threshold parameters can be imposed by a vector of positive integers,
#'    where dimensions with equal threshold parameters get the same integer. When restricting the thresholds of two
#'     outcome dimensions to be the same, one has to be careful that the number of categories in
#'      the two outcome dimensions must be the same. In our example with \eqn{J=4} different outcomes we impose:
#'
#'  \code{threshold.constraints <- c(1,1,2)}
#'
#'    gives the following restrictions:
#'  \itemize{
#'  \item \eqn{\bm\theta_{1} = \bm\theta_{2}}
#'  \item \eqn{\bm\theta_{3}} arbitrary.
#' }
#' }
#'

#'   \item{\code{threshold.values}}{
#'   In addition, threshold parameter values can be specified by \code{threshold.values}
#'    in accordance with identifiability constraints. For this purpose we use a \code{list}
#'     with \eqn{J} elements, where each element specifies the constraints of the particular
#'      dimension by a vector of length of the number of threshold parameters (number of categories - 1).
#'      A number specifies a threshold parameter to a specific value and \code{NA} leaves the parameter flexible.
#'       For \code{\link{data_mvord}} we have

#' \preformatted{threshold.constraints <- NULL}
#'
#' \preformatted{threshold.values <- list(c(-4,NA,NA,NA,NA,4.5),
#'                          c(-4,NA,NA,NA,NA,4.5),
#'                          c(-5,NA,NA,NA,NA,NA,4.5))}
#' }
#' }
#'
#' @return The function \code{mvord} returns an object of \code{\link{class}} \code{"mvord"}.
#'
#' The functions \code{summary} and \code{print} are used to display the results.
#' The function \code{coef} extracts the regression coefficients, a function \code{thresholds} the threshold coefficients
#' and the function \cr
#' \code{error_structure} returns the estimated parameters of the corresponding error structure.
#'
#' An object of \code{\link{class}} \code{"mvord"} is a list containing the following components:
#'
#' \describe{
#'  \item{\code{beta}}{
#'
#'  a named \code{\link{matrix}} of regression coefficients}
#'  \item{\code{theta}}{
#'
#'  a named \code{\link{list}} of threshold parameters}
#'   \item{\code{error.struct}}{
#'
#'   an object of class \code{\link{error_struct}} containing the parameters of the error
#'   structure}
#'   \item{\code{sebeta}}{
#'
#'     a named \code{\link{matrix}} of the standard errors of the regression coefficients}
#'   \item{\code{setheta}}{
#'
#'     a named \code{\link{list}} of the standard errors of the threshold parameters}
#'   \item{\code{seerror.struct}}{
#'
#'   a \code{vector} of standard errors for the parameters of the error structure}
#'   \item{\code{rho}}{
#'
#'     a \code{\link{list}} of all objects that are used in \code{mvord()}}
#' }
#'
#' @seealso %\code{\link{predict.mvord}},
#' \code{\link{print.mvord}}, \code{\link{summary.mvord}}, \code{\link{coef.mvord}},
#'  \code{\link{thresholds.mvord}}, \code{\link{error_structure.mvord}}, \cr
#'  \code{\link{mvord.control}}, \code{\link{data_cr_panel}},\code{\link{data_cr}},
#'  \code{\link{data_mvord_panel}},\code{\link{data_mvord}}, \code{\link{data_mvord2}}
#'
#' @references Hirk R, Hornik K, Vana L (2020). “\pkg{mvord}: An \strong{R} Package for Fitting Multivariate Ordinal Regression Models.” \emph{Journal of Statistical Software}, \strong{93}(4), 1–41, \doi{10.18637/jss.v093.i04}.
#'
#' @param formula an object of class \code{\link{formula}} of the form \code{y ~ X1 + ... + Xp}.
#' @param data \code{\link{data.frame}} containing a subject index, an index for the multiple measurements,
#' an ordinal response \code{y} and covariates \code{X1, ..., Xp}.
# #' @param index (optional) argument to specify the column names of the subject index and the multiple measurement index
# #' by a vector \cr
# #' \code{c("subject", "multiple_measurement")} in \code{data}.
# #' The default value of \code{index} is \code{NULL} assuming that the first column of \code{data} contains
# #' the subject index and the second column the multiple measurement index.
#' @param response.levels (optional) \code{\link{list}} of length equal to the number of multiple measurements to specify the category labels
#' in case of varying categories across multiple measurements
# #' @param response.names (optional) \code{\link{vector}} of the labels of the multiple measurement index in order to
# #' specify the ordering of the responses which is essential when setting constraints on the model parameters.
# #'  The default value of \code{response.names} is \code{NULL} giving the natural ordering of the levels of the factor variable
# #'  of multiple measurements.
#' @param link specifies the link function by \code{mvprobit()} (multivariate normally distributed errors - default)
#' or \code{mvlogit(df = 8)} (multivariate logistically distributed errors), where \code{df} specifies the degrees of freedom of the t copula.
#' @param error.structure different \code{error.structures}: general correlation structure (default)\cr
#' \code{cor_general(~1)},
#' general covariance structure \code{cov_general(~ 1)}, factor dependent correlation structure \code{cov_general(~ f)},
#' factor dependent covariance structure \code{cov_general(~ f)}, a constant  \code{cor_equi(~ 1)} or a covariate dependent equicorrelation structure \cr
#' \code{cor_equi(~ S)},
#' AR(1) correlation structure \code{cor_ar1(~ 1)} or a covariate dependent
#' AR(1) correlation structure \code{cor_ar1(~ S)}.
#' See \code{\link{error_struct}} or 'Details'.
#' @param weights.name (optional) character string with the column name of subject-specific weights in \code{data} which need to be
#' constant across multiple measurements. Negative weights are not allowed.
#' @param offset (optional) this can be used to specify an a priori known component to be included in the linear predictor during fitting.
#'  This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included
#'  in the formula instead or as well, and if more than one is specified their sum is used. See model.offset.
# #' @param scale If \code{scale = TRUE}, then for each response the corresponding continuous covariates are standardized before fitting,
# #'  i.e., by substracting the mean and dividing by the standard deviation.
#' @param coef.constraints (optional) \code{\link{vector}} or \code{\link{matrix}} of constraints on the regression coefficients. See 'Details'.
#' @param coef.values (optional) \code{\link{matrix}} setting fixed values on the regression coefficients. See 'Details'.
#' @param threshold.constraints (optional) \code{\link{vector}} of constraints on the threshold parameters. See 'Details'.
#' @param threshold.values (optional) \code{\link{list}} of (optional) fixed values for the threshold parameters. See 'Details'.
# #' @param se logical, if \code{TRUE} standard errors are computed.
# #' @param start.values vector of (optional) starting values.
# #' @param solver character string containing the name of the applicable solver of \code{\link{optimx}} (default is \code{"newuoa"})
# #'  or wrapper function for user defined solver.
#' @param PL.lag (optional) specifies the time lag of the pairs in the pairwise likelihood approach to be optimized (can be used with \code{cor_ar1}).
#' @param contrasts (optional) an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
#' @param control (optional) a list of parameters for controlling the fitting process. See \code{\link{mvord.control}} for details.
# #' @param control a list of control arguments. See \code{\link{optimx}}.
# #' Only applicable with \code{error.structure = cor_ar1}.
#'
#' @examples
#' library(mvord)
#'
#' #toy example
#' data(data_mvord_toy)
#'
#' #wide data format with MMO2
#' res <- mvord(formula = MMO2(Y1, Y2) ~ 0 + X1 + X2,
#'              data = data_mvord_toy)
#' print(res)
#' summary(res)
#' thresholds(res)
#' coefficients(res)
#' head(error_structure(res))
#'
#' # convert data_mvord_toy into long format
#' df <- cbind.data.frame("i" = rep(1:100,2), "j" = rep(1:2,each = 100),
#'                        "Y" = c(data_mvord_toy$Y1,data_mvord_toy$Y2),
#'                        "X1" = rep(data_mvord_toy$X1,2),
#'                        "X2" = rep(data_mvord_toy$X2,2))
#'
#' #for long format data, use MMO instead of MMO2
#' res <- mvord(formula = MMO(Y, i, j) ~ 0 + X1 + X2, #or formula = MMO(Y) ~ 0 + X1 + X2
#'                data = df)
#' print(res)
#' summary(res)
#' thresholds(res)
#' coefficients(res)
#' head(error_structure(res))
#' \donttest{
#' res2 <- mvord(formula = MMO(Y) ~ 0 + X1 + X2,
#'                data = df,
#'                control = mvord.control(solver = "BFGS"),
#'                threshold.constraints = c(1,1),
#'                coef.constraints = c(1,1))
#' print(res2)
#' summary(res2)
#' thresholds(res2)
#' coefficients(res2)
#' head(error_structure(res2))
#'
#' ## examples
#' #load data
#' data(data_mvord)
#' head(data_mvord)
#'
#' #-------------
#' # cor_general
#' #-------------
#' # approx 2 min
#' res_cor <- mvord(formula = MMO(rating) ~ 0 + X1 + X2 + X3 + X4 + X5,
#'                  data = data_mvord,
#'                  coef.constraints = cbind(c(1,2,2),
#'                                           c(1,1,2),
#'                                           c(NA,1,2),
#'                                           c(NA,NA,NA),
#'                                           c(1,1,2)),
#'                  coef.values = cbind(c(NA,NA,NA),
#'                                      c(NA,NA,NA),
#'                                      c(0,NA,NA),
#'                                      c(1,1,1),
#'                                      c(NA,NA,NA)),
#'                  threshold.constraints = c(1,1,2),
#'                  control = mvord.control(solver = "newuoa"))
#' print(res_cor)
#' summary(res_cor)
#' thresholds(res_cor)
#' coefficients(res_cor)
#' head(error_structure(res_cor))
#'
#' #-------------
#' # cov_general
#' #-------------
#' #approx 4 min
#' res_cov <- mvord(formula = MMO(rating) ~ 1 + X1 + X2 + X3 + X4 + X5,
#'                  data = data_mvord,
#'                  error.structure = cov_general(~1),
#'                  threshold.values = list(c(-4,NA,NA,NA,NA,4.5),
#'                                          c(-4,NA,NA,NA,NA,4),
#'                                          c(-5,NA,NA,NA,NA,NA,4.5))
#' ) #does not converge with BFGS
#'
#' print(res_cov)
#' summary(res_cov)
#' thresholds(res_cov)
#' coefficients(res_cov)
#' head(error_structure(res_cov))
#'
#'
#' #-------------
#' # cor_ar1
#' #-------------
#' #approx 4min
#' data(data_mvord_panel)
#' head(data_mvord_panel)
#'
#' #select subset of data
#' subset_dat <- data_mvord_panel$year %in% c("year3", "year4", "year5", "year6", "year7")
#' data_mvord_panel <- data_mvord_panel[subset_dat,]
#'
#' mult.obs <- 5
#' res_AR1 <- mvord(formula = MMO(rating) ~ 0 + X1 + X2 + X3 + X4 + X5,
#'                  data = data_mvord_panel,
#'                  error.structure = cor_ar1(~1),
#'                  threshold.constraints = c(1,1,1,2,2),
#'                  coef.constraints = c(1,1,1,2,2),
#'                  control = mvord.control(solver = "BFGS"))
#' print(res_AR1)
#' summary(res_AR1)
#' thresholds(res_AR1)
#' coefficients(res_AR1)
#' head(error_structure(res_AR1))
#' head(error_structure(res_AR1, type = "corr"))
#'
#' data(data_mvord2)
#' # approx 2 min
#' res_cor <- mvord(formula = MMO2(rater1, rater2, rater3) ~ 0 + X1 + X2 + X3 + X4 + X5,
#'                  data = data_mvord2,
#'                  coef.constraints = cbind(c(1,2,2),
#'                                           c(1,1,2),
#'                                           c(NA,1,2),
#'                                           c(NA,NA,NA),
#'                                           c(1,1,2)),
#'                  coef.values = cbind(c(NA,NA,NA),
#'                                      c(NA,NA,NA),
#'                                      c(0,NA,NA),
#'                                      c(1,1,1),
#'                                      c(NA,NA,NA)),
#'                  threshold.constraints = c(1,1,2),
#'                  control = mvord.control(solver = "newuoa"))
#' print(res_cor)
#' summary(res_cor)
#' thresholds(res_cor)
#' coefficients(res_cor)
#' head(error_structure(res_cor))
#' }
#'
#' @name mvord
#' @export
mvord <- function(formula,
                  data,
                  error.structure = cor_general(~1),
                  link = mvprobit(),
                  response.levels = NULL,
                  coef.constraints = NULL,
                  coef.values = NULL,
                  threshold.constraints = NULL,
                  threshold.values = NULL,
                  weights.name = NULL,
                  offset = NULL,
                  PL.lag = NULL,
                  contrasts = NULL,
                  control = mvord.control()
){
  #check arguments
  rho <- list()
  rho$timestamp1 <- proc.time()
  rho$link <- link
  rho$se <- control$se
  rho$start.values <- control$start.values
  rho$threshold.constraints <- threshold.constraints
  rho$threshold.values <- threshold.values
  rho$coef.constraints_input <- coef.constraints
  rho$coef.values_input <- coef.values
  rho$formula.input <- formula
  rho$PL.lag <- PL.lag
  rho$solver <- control$solver
  rho$combis <- control$combis
  rho$control <- control$solver.optimx.control
  check_args_optimizer(rho)
  check_args_error.structure(error.structure, data)
  rho$mc <- match.call(expand.dots = FALSE)
  rho$weights.name <- weights.name
  rho$response.levels <- response.levels
  rho$offset <- offset
  #rho$fast_fit <- fast
  if (!all(names(contrasts) %in% c(labels(terms(rho$formula.input)), labels(terms(error.structure$formula))))) {
    Terms <- c(labels(terms(rho$formula.input)), labels(terms(error.structure$formula)))
    id_tmp <- which(!(names(contrasts) %in% Terms))
    str_tmp <- paste(names(contrasts)[id_tmp], collapse = " and ")
    warning(ifelse(length(id_tmp) == 1, sprintf("Variable %s is absent, the contrasts will be ignored.", str_tmp),
                 sprintf("Variables %s are absent, the contrasts will be ignored.", str_tmp)))
  }
  rho$contrasts <- contrasts
  nm <- names(as.list(rho$mc))
  if(!"data" %in% nm) stop("Model needs formula and data.", call.=FALSE)
  #if(!"link" %in% nm) stop("Model needs formula, data and link.", call.=FALSE)
  if(!"formula" %in% nm) stop("Model needs formula and data.", call.=FALSE)
  #formula
  if(!inherits(formula, "formula")) stop("formula has to be of class formula.", call. = FALSE)
  ##  check if data is a data.frame
  if(!is.data.frame(data)) {
    warning("data has to be of type data.frame. Automatically applied as.data.frame() to data.")
    data <- as.data.frame(data)
  }
  ### MMO ###
  if(formula[[2]][[1]] == "MMO"){
    rho <- initialize_MMO(rho, formula, data, error.structure, contrasts)
  ### mvord2 ###
  } else if(formula[[2]][[1]] == "MMO2"){
    rho <- initialize_MMO2(rho, formula, data, error.structure, contrasts)
  }
##############
  rho <- set_args_other(rho)
  check_response_missings(rho)
  mvord.fit(rho)
}
#-----------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------
mvord_finalize <- function(rho){
  est <- list()
  ## Thresholds
  est$theta <- rho$theta
  for(j in seq_len(rho$ndim)){
    names(est$theta[[j]]) <- get_labels_theta(rho,j)
  }
  names(est$theta) <- rho$y.names
  tmp <- rho$threshold.values
  ## Regression coefficients
  est$beta <- rho$beta
  names(est$beta) <- unlist(lapply(seq_along(rho$constraints),
  	function(p) colnames(rho$constraints[[p]])))

  ## Error Structure
  sigma_par <- rho$optpar[rho$npar.thetas + rho$npar.betas +  seq_len(attr(rho$error.structure, "npar"))]
  est$error.struct <- finalize(rho$error.structure, sigma_par)
  ## If standard errors should be computed
  if (rho$se) {
	  est$setheta <- lapply(seq_len(rho$ndim), function(j){
	      if(rho$threshold.constraints[j] %in% rho$threshold.constraints[seq_len(j-1)]){#check threshold.constraints
	        k <- match(rho$threshold.constraints[j], rho$threshold.constraints)
	        tmp[[j]][!is.na(tmp[[k]])] <- 0
	        tmp[[j]][is.na(tmp[[k]])] <- rho$seGamma[rho$ind.thresholds[[k]]]
	      } else {
	        tmp[[j]][!is.na(tmp[[j]])] <- 0
	        tmp[[j]][is.na(tmp[[j]])] <- rho$seGamma[rho$ind.thresholds[[j]]]
	      }
	      tmp[[j]]
	    })
	  names(est$setheta) <- rho$y.names
	  for (j in seq_len(rho$ndim)) names(est$setheta[[j]]) <- get_labels_theta(rho,j)
   	  if (rho$npar.betas > 0) {
   	  	est$sebeta <- rho$seGamma[rho$npar.thetas + seq_len(rho$npar.betas)]
   	  	names(est$sebeta) <- unlist(lapply(rho$constraints, function(x) colnames(x)))
   	  } else {
   	  	est$sebeta <- numeric(0)
   	  }
      est$seerror.struct <- rho$seGamma[rho$npar.thetas + rho$npar.betas +  seq_len(attr(rho$error.structure, "npar"))]
  }
  est
}

#' @title Print Method for Multivariate Ordinal Regression Models.
#' @description Prints thresholds, regression coefficients
#'  and parameters of the error structure of class \code{'mvord'}.
#' @param x object of class \code{'mvord'}
#' @param call displays function call if \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @method print mvord
#' @export
print.mvord <- function(x, call = TRUE, ...){
  if(call){
  cat("\nCall:\n",
      paste(deparse(x$rho$mc), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }
  cat("Thresholds:\n")
  print(x$theta, ...)

  cat("Coefficients:\n")
  print(x$beta, ...)
  cat("\n")
  cat("Parameters of the error structure:\n")
  print(attr(x$error.struct, "par"), ...)
  cat("\n")
  invisible(x)
}

#' @title Summary method for Multivariate Ordinal Regression Models.
#' @description Summary of thresholds, regression coefficients
#' and parameters of the error structure of class \code{'mvord'}.
#' @param object object of class \code{'mvord'}
#' @param short if \code{TRUE} short summary, otherwise extended summary
#' @param call displays function call if \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @method summary mvord
#' @export
summary.mvord <- function(object, short = TRUE, call = TRUE, ...){
  ntotal <-  sum(object$rho$ntheta) + length(object$beta) +
    attr(object$error.struct, "npar")

  names.theta <- unlist(lapply(1:object$rho$ndim, function(j) paste(names(object$theta)[j], names(object$theta[[j]]))))
  names.beta <- names(object$beta)
  names.sigma <- attr(object$error.struct, "parnames")

  coef <- as.matrix(c(unlist(object$theta), as.vector(object$beta),
  	attr(object$error.struct, "par")), ncol = 1)
  rownames(coef) <- c(names.theta, names.beta, names.sigma)
  colnames(coef) <- c("Estimate")
  if(object$rho$se){
    tmp <- matrix(0, ntotal, 4,
                  dimnames = list(c(names.theta, names.beta, as.vector(names.sigma)),
                                  c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
    tmp[, 1] <- coef
    tmp[, 2] <- c(unlist(object$setheta), as.vector(object$sebeta), object$seerror.struct)
    tmp[, 3] <- tmp[, 1]/tmp[, 2]
    tmp[, 4] <- 2 * pnorm(abs(tmp[, 3]), lower.tail=FALSE)
    coef <- cbind.data.frame(tmp)
    coef[tmp[, 2] == 0, 3] <- NA
    coef[tmp[, 2] == 0, 4] <- NA
  }
  mat <- cbind.data.frame(c("link",object$rho$link$name), c("threshold",object$rho$threshold),
                          c("nsubjects", object$rho$n), c("ndim", object$rho$ndim),
                          c("logPL", round(-object$rho$objective,2)),c("CLAIC", ifelse(object$rho$se,round(object$rho$claic,2),NA)),
                          c("CLBIC", ifelse(object$rho$se,round(object$rho$clbic,2),NA)),
                          c("fevals", if(is.null(object$rho$optRes$fevals)) NA else object$rho$optRes$fevals))
  summary.output <- list()
  if(call){
    summary.output$call <-  object$rho$mc
    cat("\nCall: ",
        paste(deparse(object$rho$mc, width.cutoff = getOption("deparse.cutoff")), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }
  summary.output$formula <- object$rho$formula
  cat("Formula: ",
      paste(deparse(object$rho$formula, width.cutoff = getOption("deparse.cutoff")), sep = "\n", collapse = "\n"), "\n", sep = "")
  # summary.output$formula <- print(object$rho$formula)

  cat("\n")
  summary.output$info <- format(mat, justify="right")
  write.table(summary.output$info, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("\n")
  if(short){
    tmp.duplicated <- duplicated(object$rho$threshold.constraints)#duplicated(lapply(object$theta, function(x) unname(round()))
    len.thresh <- object$rho$ntheta
    lower.ind <- cumsum(c(1, len.thresh[-length(len.thresh)]))
    upper.ind <- cumsum(len.thresh)
    tmp.ind <- lapply(seq_along(len.thresh), function(j) if(tmp.duplicated[j] == FALSE) seq(lower.ind[j], upper.ind[j]) else c())
    cat("Thresholds:\n")
    if(object$rho$se) {
      summary.output$thresholds <- printCoefmat(coef[unlist(tmp.ind),])
    } else summary.output$thresholds <- print(coef[unlist(tmp.ind),], ...)

    cat("\nCoefficients:\n")
    tmp.ind <- length(names.theta) + seq_along(names.beta)
    if(object$rho$se) {
      summary.output$coefficients <- printCoefmat(coef[unlist(tmp.ind),])
    } else summary.output$coefficients <- print(coef[unlist(tmp.ind),], ...)
  } else{
    cat("Thresholds:\n")
    if(object$rho$se) {
      summary.output$thresholds <- printCoefmat(coef[seq_along(names.theta),])
    } else summary.output$thresholds <- print(coef[seq_along(names.theta),], ...)

    cat("\nCoefficients:\n")
    if(object$rho$se) {
      summary.output$coefficients <- printCoefmat(coef[length(names.theta) + seq_along(names.beta),])
    } else summary.output$coefficients <- print(coef[length(names.theta) + seq_along(names.beta),], ...)
  }
  cat("\nError Structure:\n")
    if(object$rho$se) {
    summary.output$error.structure <- printCoefmat(coef[length(names.theta) + length(names.beta) + seq_along(names.sigma),])
  } else summary.output$error.structure <- print(coef[length(names.theta) + length(names.beta) + seq_along(names.sigma),], ...)
  class(summary.output) <- "summary.mvord"
  invisible(summary.output)
}

print.summary.mvord <- function(x, ...){
  if(!is.null(x)){
    cat("\nCall: ",
        x$call, "\n\n", sep = "")
  }
  cat("Formula: ")
  print(x$formula, ...)
  cat("\n")
  write.table(x$info, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("\n")
    cat("Thresholds:\n")
    if(ncol(x$thresholds) > 1) printCoefmat(x$thresholds) else print(x$thresholds, ...)
    cat("\nCoefficients:\n")
    if(ncol(x$coefficients) > 1) print(x$coefficients) else print(x$coefficients, ...)
  cat("\nError Structure:\n")
  if(ncol(x$error.structure) > 1) printCoefmat(x$error.structure) else print(x$error.structure, ...)

}

#' @title Coefficients of Multivariate Ordinal Regression Models.
#' @description \code{coef} is a generic function which extracts
#'  regression coefficients from objects of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method coef mvord
#' @export
coef.mvord <- function(object, ...) object$beta

#' @title Thresholds of Multivariate Ordinal Regression Models.
#' @description
#' \code{thresholds} is a generic function which extracts threshold coefficients from objects of class \cr
#' \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
# #' @rdname mvord
#' @export
thresholds <- function(object, ...) UseMethod("thresholds")
#' @rdname thresholds
#' @export
thresholds.mvord <- function(object, ...) object$theta

# #' @title CLAIC of Multivariate Ordinal Regression Models.
# #' @description
# #' \code{claic} is a generic function which extracts the composite likelihood AIC from objects of class \cr
# #' \code{"mvord"}.
# #' @param object object of class \code{'mvord'}
# #' @export
claic <- function(object) UseMethod("claic")
# #' @rdname claic
# #' @export
claic.mvord <- function(object) ifelse(object$rho$se, object$rho$claic, NA)


# #' @title CLBIC of Multivariate Ordinal Regression Models.
# #' @description
# #' \code{clbic} is a generic function which extracts the composite likelihood BIC from objects of class \cr
# #' \code{"mvord"}.
# #' @param object object of class \code{'mvord'}
# #' @export
clbic <- function(object) UseMethod("clbic")
# #' @rdname clbic
# #' @export
clbic.mvord <- function(object) ifelse(object$rho$se, object$rho$clbic, NA)


# #' @title logPL of Multivariate Ordinal Regression Models.
# #' @description
# #' \code{logPL} is a generic function which extracts the log pairwise likelihood from objects of class \cr
# #' \code{"mvord"}.
# #' @param object object of class \code{'mvord'}
# #' @export
logPL <- function(object) UseMethod("logPL")
# #' @rdname logPL
# #' @export
logPL.mvord <- function(object) as.numeric(-object$rho$objective)


#' @title nobs of Multivariate Ordinal Regression Models.
#' @description
#' \code{nobs} is a generic function which extracts the number of observations from objects of class \cr
#' \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method nobs mvord
#' @export
nobs.mvord <- function(object, ...) object$rho$n

#' @title vcov of Multivariate Ordinal Regression Models.
#' @description
#' \code{vcov} is a generic function which extracts the Godambe information matrix from objects of class \cr
#' \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method vcov mvord
#' @export
vcov.mvord <- function(object, ...) object$rho$varGamma

#' @title terms of Multivariate Ordinal Regression Models.
#' @description
#' \code{terms} is a generic function which can be used to extract terms from objects of class \cr
#' \code{'mvord'}.
#' @param x an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method terms mvord
#' @export
terms.mvord <- function(x, ...) terms(x$rho$formula)

#' @title model.matrix of Multivariate Ordinal Regression Models.
#' @description
#' \code{model.matrix} is a generic function which extracts the model matrix
#'  from objects of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method model.matrix mvord
#' @export
model.matrix.mvord <- function(object, ...) object$rho$x

#' @title Fitted Probabilities of Multivariate Ordinal Regression Models.
#' @description
#' A generic function which extracts fitted probabilities for the observed categories from objects of class
#' \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method fitted mvord
#' @export
fitted.mvord <- function(object, ...) predict(object, ...)

#' @title Pairwise Log-Likelihood of Multivariate Ordinal Regression Models.
#' @description
#' \code{logLik} is a generic function which extracts the pairwise log-likelihood from objects of class \cr
#' \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @param ... further arguments passed to or from other methods.
#' @method logLik mvord
#' @export
logLik.mvord <- function(object, ...) structure(-object$rho$objective,
                                                df = sum(diag(object$rho$V %*% object$rho$H.inv)), class = "logLik")

#' @title Constraints on the Regression Coefficients of Multivariate Ordinal Regression Models.
#' @description
#' An extractor function for the constraint matrices corresponding to the regression
#' coefficients from objects of class \code{'mvord'}.
#' @param object an object of class \code{'mvord'}.
#' @export
constraints <- function(object) UseMethod("constraints")
#' @rdname constraints
#' @export
constraints.mvord <- function(object) object$rho$constraints


#' @title Names of regression coefficient constraints in mvord
#' @description
#' An extractor function for the names of the regression coefficient constraints based on the model \code{formula} and \code{data}.
#' @param formula model formula
#' @param data a given data set.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
#' @export
names_constraints <- function(formula, data, contrasts = NULL){
  intercept <- ifelse(attr(terms.formula(formula), "intercept") == 1, TRUE, FALSE)
  nam <-  colnames(model.matrix(update(as.formula(paste(as.character(formula[-2]), collapse = " ")), ~ . + 1), data, contrasts.arg = contrasts))
  if(intercept) nam else nam[-1]
}
