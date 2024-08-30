library(mvord)
library(MASS)


#data(data_toy_example)
tolerance <- 1e-6

mult.obs <- 2


sigma <- matrix(c(1,0.8,0.8,1), ncol = 2)

thresholds <- list(c(-1,1),c(-1,1))
nobs <- 1000
suppressWarnings(RNGversion("3.5.0"))
set.seed(2017)
errors <-  mvrnorm(n = nobs, mu = rep(0, mult.obs), Sigma = sigma)

p <- 10
pa <- 5 # active betas
betas <- rep(list(c(sample(-10:10, pa)/10, rep(0, p-pa))),
                  mult.obs)


pred <- matrix(rnorm(nobs * p), nrow = nobs)

y <- sapply(1:mult.obs, function(j) pred %*% betas[[j]] + errors[, j], simplify = "array")
y.ord <- sapply(1:mult.obs, function(j) cut(y[, , j], c(min(y[, , j]) - 1,
                                                        c(thresholds[[j]]), max(y[, , j]) + 1),
                                            labels = FALSE), simplify = "array")

predictors.fixed <- lapply(1:mult.obs, function(j) pred)
y <- as.data.frame(y.ord)

for(i in 1:mult.obs){
  y[, i] <- factor(y[, i], levels = sort(unique(y[, i])),
                   ordered = TRUE)
}




data_toy_example <- cbind.data.frame(y, predictors.fixed[[1]])
colnames(data_toy_example) <- c("Y1","Y2", paste0("X", 1:p))
w <- c(rep(1/20, 20), rep(1/30, 30), rep(1/20, 20), rep(1/30, 30))

formula <- paste0("MMO2(Y1, Y2) ~ 0 + ", paste0("X", 1:p, collapse = " + "))
res_mvord <- mvord:::mvord(formula = as.formula(formula),
                              data = data_toy_example,
                              link = mvprobit(),
                              error.structure = cor_general(~1),
                              threshold.constraints = c(1,1),
                              coef.constraints = c(1,1),
                              control= mvord.control(solver="BFGS",se=TRUE))

res_mvord$beta

compute_lambda_max <- function(y, X, alpha) {
  if (NCOL(y) > 1) {
    yvec <- c(sapply(y, function(yy) as.numeric(yy)))
    X <- do.call("rbind", rep(list(X), mult.obs))
  } else {
    yvec <- as.numeric(yvec)
  }

  ybar <- mean(yvec, na.rm = TRUE)
  out <- apply(X, 2, function(xj){
    1/ifelse(alpha == 0, 0.01, alpha) *
      mean((yvec - ybar * (1 - ybar)) * rep(xj, ncol(y)))
  })
  max(abs(out))
}

## TODO determine via coordinate descent
lambda_max <- 1000 # compute_lambda_max(y, scale(pred), alpha = 1)
lambda_min_ratio <- 0.0001
lambda_min <- lambda_max * lambda_min_ratio
lambda_grid <- exp(seq(log(lambda_max),
                       log(lambda_min), length.out = nlambda))
claic_vec_ridge <- NULL
beta_mat_ridge <- matrix(nrow=p,ncol=0)
for (i in seq_along(lambda_grid)) {
  ## warm starts don't work well as they cause false convergence
  lambda <- lambda_grid[i]
  # if (i == 1) {
  start_vals <- NULL
  # } else {
  #   start_vals <- list(theta=res_ridge$theta[[1]],
  #                      beta=res_ridge$beta)
  # }
  res_ridge <- mvord:::mvordnet(formula = as.formula(formula),
                                data = data_toy_example,
                                link = mvprobit(),
                                alpha = 0,
                                lambda = lambda,
                                error.structure = cor_general(~1),
                                threshold.constraints = c(1,1),
                                coef.constraints = c(1,1),
                                control= mvord.control(
                                  start.values = start_vals,
                                  solver="Nelder-Mead",se=TRUE))
  claic_vec_ridge <- c(claic_vec_ridge, AIC(res_ridge))
  beta_mat_ridge <- cbind(beta_mat_ridge, coef(res_ridge))
}

claic_vec_lasso <- NULL
beta_mat_lasso <- matrix(nrow=p,ncol=0)
for (i in seq_along(lambda_grid)) {
  ## warm starts don't work well as they cause false convergence
  lambda <- lambda_grid[i]
  # if (i == 1) {
  start_vals <- NULL
  # } else {
  #   start_vals <- list(theta=res_ridge$theta[[1]],
  #                      beta=res_ridge$beta)
  # }
  res_lasso <- mvord:::mvordnet(
    formula = as.formula(formula),
    data = data_toy_example,
    link = mvprobit(),
    alpha = 1L,
    lambda = lambda,
    error.structure = cor_general(~1),
    threshold.constraints = c(1,1),
    coef.constraints = c(1,1),
    control= mvord.control(solver="Nelder-Mead",
                           solver.optimx.control = list(trace=1),
                           se=TRUE))
  claic_vec_lasso <- c(claic_vec_lasso, AIC(res_lasso))
  beta_mat_lasso <- cbind(beta_mat_lasso, coef(res_lasso))
}
plot(claic_vec_lasso)
beta_mat_lasso
