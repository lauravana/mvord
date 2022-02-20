library(MASS)
source("work/simulate.R")

## simulate
mult.obs <- 2
sigma <- matrix(c(1,0.8,0.8,1), ncol = 2)
betas <- list(c(1.5,-3),
              c(1.5,-3))
thresholds <- list(c(-2,2),c(-2,2))
nobs <- 100
seed <- 3
6
seed <- 5
set.seed(seed)
errors <-  mvrnorm(n = nobs, mu = rep(0, mult.obs), Sigma = sigma)

X1 <- rbeta(nobs, shape1 = 3, shape2 = 5) * 6
X2 <- rbeta(nobs, shape1 = 2, shape2 = 5) * 4
X1 <- rnorm(nobs, 1, 1)
X2 <- rnorm(nobs, 0.5, 1)



pred <- cbind(X1, X2)

y <- sapply(1:mult.obs, function(j) pred %*% betas[[j]] + errors[, j], simplify = "array")
plot(density(y))
summary(y)
y.ord <- sapply(1:mult.obs, function(j) cut(y[, , j], c(min(y[, , j]) - 1, c(thresholds[[j]]), max(y[, , j]) + 1),
                                            labels = FALSE), simplify = "array")

table(y.ord[,1])
table(y.ord[,2])
predictors.fixed <- lapply(1:mult.obs, function(j) pred)


y <- as.data.frame(y.ord)
table(y)

for(i in 1:mult.obs){
  y[, i] <- factor(y[, i], levels = sort(unique(y[, i])),
                   ordered = TRUE)
}

data_mvord_toy <- cbind.data.frame(y, predictors.fixed[[1]])
colnames(data_mvord_toy) <- c("Y1","Y2", "X1", "X2")
save(data_mvord_toy, file = "data/data_mvord_toy.rda")


df <- cbind.data.frame("i" = rep(1:100,2), "j" = rep(1:2,each = 100),
                       "Y" = c(data_toy_example_new$Y1,data_toy_example_new$Y2),
                       "X1" = rep(data_toy_example_new$X1,2),
                       "X2" = rep(data_toy_example_new$X2,2))

res <- mvord(formula = Y ~ 0 + X1 + X2,
             data = df)
res


mvord(formula = Y ~ 0 + X1 + X2,
      data = df,
      threshold.constraints = c(1,1),
      coef.constraints = c(1,1),
      control = mvord.control(solver = "BFGS"))

table(marginal.predict(res, type = "class")[,1], res$rho$y[,1])
table(predict(res, type = "class.max")[,1], res$rho$y[,1])

table(marginal.predict(res, type = "class")[,2], res$rho$y[,2])
table(predict(res, type = "class.max")[,2], res$rho$y[,2])
