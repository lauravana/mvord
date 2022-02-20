set.seed(20170818)

###########################
#####  Load real data #####
###########################
#load("work/data_ordinal_failInd.rda")
setwd("~/svn/baRpkgs/trunk/mvord/sim_data/")
#load("work/data_ordinal_win_ratios_impute_nfu_sectors_20002013_7rc.rda")
#load("work/data_ordinal_failInd.rda")
load("data_ordinal_noNA_win_ratios_sectors_without40_55_20002013_7rc_weights.rda")
#data_ordinal <- data_ordinal_fail

str(data_ordinal)

## pick only 10 years
dat <- data_ordinal[data_ordinal$fyear %in% 2005:2012, ]
#dat <- data_ordinal#[data_ordinal$fyear %in% 2001:2010, ]
with(dat, table(SPR))

# ## make ratings nice
levels(dat$SPR)[1:2] <- "Atmp"



with(dat, table(fyear))
with(dat, table(SPR))


dat$R8[dat$R8 == 0] <- .Machine$double.eps
dat$R24[dat$R24 == 0] <- .Machine$double.eps
dat$R31[dat$R31 == 0] <- .Machine$double.eps
#dat$RSIZE <- exp(dat$RSIZE)

## "anonimize" the gvkey
dat$gvkey <- factor(dat$gvkey, labels = 1:length(unique(dat$gvkey)))

## "anonimize" the gsector
dat$gsector <- factor(dat$gsector, labels = paste0("BSEC", 1:length(unique(dat$gsector))))

gsec_per_firm <- sapply(unique(dat$gvkey), function(gv) unique(dat[dat$gvkey==gv, "gsector"]))
#############################
## simulate the covariates ##
#############################
## fit a distribution to the original covariates

dat <- dat[!is.na(dat$SPR),]





ind <- match(unique(dat$gvkey), dat$gvkey)
weights <- 1/table(dat$SPR[ind]) * 100
weights_ind <- weights[match(dat$SPR[ind], names(weights))]
ind_rev <- match(dat$gvkey, unique(dat$gvkey))

dat$weights <- weights_ind[ind_rev]


BSEC <- dat$gsector[ind]

n <- length(ind)


library(mvord)

#2005-2008
res_probit_ar1__real_dat_05_08 <- mvord(
  formula = MMO(SPR, gvkey, fyear) ~ 0 + R5 + R13 + R20 + RSIZE + BETA,
  data = dat[dat$fyear %in% 2005:2008, ],
  weights = "weights",
  # response.levels = rep(list(rev(paste0("RC", 1:6, "_R1"))), 10),
  # response.names = c(2002:2011),
  coef.constraints = c(1,1,1,1),
  error.structure = cor_ar1(~ gsector),
  threshold.constraints = rep(1, 4),
  #coef.constraints = c(rep(1, 5), rep(2, 5)),
  link = mvprobit(),
  control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))
summary(res_probit_ar1__real_dat_05_08)
res_probit_ar1__real_dat_05_08

table(res_probit_ar1__real_dat_05_08$rho$y[,1], marginal_predict(res_probit_ar1__real_dat_05_08, type = "class")[,1])
table(res_probit_ar1__real_dat_05_08$rho$y[,2], marginal_predict(res_probit_ar1__real_dat_05_08, type = "class")[,2])
table(res_probit_ar1__real_dat_05_08$rho$y[,3], marginal_predict(res_probit_ar1__real_dat_05_08, type = "class")[,3])
table(res_probit_ar1__real_dat_05_08$rho$y[,4], marginal_predict(res_probit_ar1__real_dat_05_08, type = "class")[,4])
res_probit_ar1__real_dat_05_08$rho$runtime

res_probit_ar1__real_dat_09_12 <- mvord(
  formula = MMO(SPR, gvkey, fyear) ~ 0 + R5 + R13 + R20 + RSIZE + BETA,
  data = dat[dat$fyear %in% 2009:2012, ],
  weights = "weights",
  # response.levels = rep(list(rev(paste0("RC", 1:6, "_R1"))), 10),
  # response.names = c(2002:2011),
  coef.constraints = c(1,1,1,1),
  error.structure = cor_ar1(~ gsector),
  threshold.constraints = rep(1, 4),
  threshold.values = thresholds(res_probit_ar1__real_dat_05_08),
  #coef.constraints = c(rep(1, 5), rep(2, 5)),
  link = mvprobit(),
  control = mvord.control(solver = "BFGS", solver.optimx.control = list(trace = TRUE)))
summary(res_probit_ar1__real_dat_09_12)

table(res_probit_ar1__real_dat_09_12$rho$y[,1], marginal_predict(res_probit_ar1__real_dat_09_12, type = "class")[,1])
table(res_probit_ar1__real_dat_09_12$rho$y[,2], marginal_predict(res_probit_ar1__real_dat_09_12, type = "class")[,2])
table(res_probit_ar1__real_dat_09_12$rho$y[,3], marginal_predict(res_probit_ar1__real_dat_09_12, type = "class")[,3])
table(res_probit_ar1__real_dat_09_12$rho$y[,4], marginal_predict(res_probit_ar1__real_dat_09_12, type = "class")[,4])
res_probit_ar1__real_dat_09_12$rho$runtime


source("../R/utilities.R")
xs <- mvord_data(dat, index = c("gvkey", "fyear"), y.names = "SPR", x.names = c("R5", "R13", "R18", "R20", "RSIZE", "BETA"),
                 y.levels = NULL, response.names = 2005:2012)[[3]]


winsorize <- function(x, win_lev){
  x[x>quantile(x,1-win_lev/2)] <- quantile(x,1-win_lev/2)
  x[x<quantile(x,win_lev/2)] <- quantile(x,win_lev/2)
  x
}
win_lev <- 0.01

xdiff <- lapply(1:7, function(i){
  sapply(1:6, function(p){
x <- xs[[i+1]][,p] - xs[[i]][,p]
x <- x[!is.na(x)]
#if(p == 6) rep(0,n) else{
est <- MASS::fitdistr(winsorize(x, win_lev), densfun = "t")$estimate
rt(n, df = est[3]) * est[2] + est[1]
#}
  })
})


## unbounded ratios
library(MASS)
dat$R5[dat$R5 == 0] <- .Machine$double.eps
dat$R13[dat$R13 == 0] <- .Machine$double.eps
dat$R18[dat$R18 == 0] <- .Machine$double.eps

distr <- c("cauchy",   "logistic", "normal",  "t")
df_R3    <- distr[which.min(sapply(distr, function(f) AIC(MASS::fitdistr(dat$R3, densfun = f))))]
df_R20 <- distr[which.min(sapply(distr, function(f) AIC(MASS::fitdistr(dat$R20, densfun = f))))]
df_RSIZE <- distr[which.min(sapply(distr, function(f) AIC(MASS::fitdistr(dat$RSIZE, densfun = f))))]

## positive ratios
distr_pos <- c("cauchy",  "exponential", "geometric",
               "log-normal", "lognormal", "logistic", "normal",
               "Poisson", "t", "weibull")

df_R5 <- distr_pos[which.min(sapply(distr_pos, function(f) AIC(fitdistr(dat$R5, densfun = f))))]
df_R13 <- distr_pos[which.min(sapply(distr_pos, function(f) AIC(fitdistr(dat$R13, densfun = f))))]
df_R18 <- distr_pos[which.min(sapply(distr_pos, function(f) AIC(fitdistr(dat$R18, densfun = f))))]
df_SIGMA <- distr_pos[which.min(sapply(distr_pos, function(f) AIC(fitdistr(dat$SIGMA, densfun = f))))]
df_BETA <- distr_pos[which.min(sapply(distr_pos, function(f) AIC(fitdistr(dat$BETA, densfun = f))))]


est_tab <- list(
  fitdistr(dat$R3,  densfun =  df_R3)$estimate,
  fitdistr(dat$R5,  densfun =  df_R5)$estimate,
  fitdistr(dat$R13, densfun = df_R13)$estimate,
  fitdistr(dat$R18, densfun = df_R18)$estimate,
  fitdistr(dat$R20, densfun = df_R20)$estimate,
  fitdistr(dat$RSIZE, densfun = df_RSIZE)$estimate,
  fitdistr(dat$SIGMA, densfun = df_SIGMA)$estimate,
  fitdistr(dat$BETA, densfun = df_BETA)$estimate
)


#dat <- dat[dat$fyear %in% 2010:2012, ]
n <- nrow(dat[dat$fyear %in% 2010:2013, ])
n <- 1415




sim_dat <- data.frame(firm_id = 1:n)
#sim_dat$R3    <- rt(n, df = est_tab[[1]][3]) * est_tab[[1]][2] + est_tab[[1]][1]
sim_dat$R5    <- rlnorm(n, meanlog = est_tab[[2]][1], sdlog = est_tab[[2]][2])
sim_dat$R13   <- rweibull(n, shape = est_tab[[3]][1], scale = est_tab[[3]][2])
#sim_dat$R13   <- rt(n, df = est_tab[[3]][3]) * est_tab[[3]][2] + est_tab[[3]][1]
#sim_dat$R18   <- rweibull(n, shape = est_tab[[4]][1], scale = est_tab[[4]][2])
#sim_dat$R18   <- rt(n, df = est_tab[[4]][3]) * est_tab[[4]][2] + est_tab[[4]][1]
sim_dat$R20   <- rt(n, df = est_tab[[5]][3]) * est_tab[[5]][2] + est_tab[[5]][1]
sim_dat$RSIZE   <- rnorm(n, mean = est_tab[[3]][1], sd = est_tab[[3]][2])
#sim_dat$SIGMA <- rlnorm(n, meanlog = est_tab[[7]][1], sdlog = est_tab[[7]][2])
sim_dat$BETA   <- rweibull(n, shape = est_tab[[8]][1], scale = est_tab[[8]][2])
#sim_dat$BETA  <- rt(n, df = est_tab[[8]][3]) * est_tab[[8]][2] + est_tab[[8]][1]
sim_dat$R13[sim_dat$R13 < 0] <- 0
sim_dat$R18[sim_dat$R18 < 0] <- 0
#sim_dat[,-1] <- apply(sim_dat[,-1], 2, winsorize, win_lev=0.00)
summary(sim_dat)



xall <- c(list(sim_dat[,2:6]), xdiff)

str(sim_dat[,2:6])
str(xall)
str(xdiff)

sapply(xall, nrow)
j=8


xsim <- lapply(1:8, function(j){
  Reduce('+', xall[1:j])
})


summary(res_probit_ar1__real_dat_05_08)
coef(res_probit_ar1__real_dat_05_08)

summary(res_probit_ar1__real_dat_09_12)
coef(res_probit_ar1__real_dat_09_12)




beta_tab <- cbind(#R3 = c(0,0,0,0),
  R5 = c(rep(0.02,4),rep(0.033,4)),
  R13 = c(rep(0.045,4),rep(0.02,4)),
  #R18 = c(rep(0.4,4),rep(0.48,4)),
  R20 = c(rep(-1.15,4),rep(-0.8,4)),
  RSIZE = c(rep(-0.38,4),rep(-0.38,4)),
  BETA = c(rep(0.18,4),rep(0.16,4))
)

#beta_tab <- beta_tab *2

thresholds(res_probit_ar1__real_dat_05_08)
thresholds(res_probit_ar1__real_dat_09_12)

theta<- c(-1.5,-0.5,1,2.5)

theta_labels <- LETTERS[1:5]

r <- c(0.9, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.6)

r.corr <- r[as.numeric(BSEC)]
length(r.corr)

library(MASS)
errors <- do.call(rbind, lapply(1:n, function(f){
  #print(f)
mult.obs <- 8
#num.covariates <- 7
sigma <- diag(mult.obs)
for (i in 1:mult.obs){
  for (j in 1:mult.obs){
    sigma[i,j] <- r.corr[f]^abs(i-j)
  }
}
  mvrnorm(n = 1, mu = rep(0, mult.obs), Sigma = sigma)
}))

summary(errors)





ysim <- sapply(1:8, function(j) {
  pred <- as.matrix(xsim[[j]]) %*% beta_tab[j,] + errors[,j]
  cut(pred, breaks = c(-Inf,theta, Inf), labels = theta_labels)
})

pred <- sapply(1:8, function(j) {
  pred <- as.matrix(xsim[[j]]) %*% beta_tab[j,] #+ errors[,j]
  #cut(pred, breaks = c(-Inf,theta, Inf), labels = theta_labels)
})

summary(pred)
plot(density(pred))

table(ysim[,1])
table(ysim[,2])
table(ysim[,3])
table(ysim[,4])
table(ysim[,5])


yt <- sapply(1:4, function(j) {
  pred <- as.matrix(xsim[[j]]) %*% beta_tab[j,] #+ errors[,j]
  cut(pred, breaks = c(-Inf,theta, Inf), labels = theta_labels)
})

table(yt[,1])
table(yt[,2])
table(yt[,3])
table(yt[,4])

summary(pred)

ysim[1:2,]

as.vector(ysim)

str(xsim)


data_sim_long <- cbind.data.frame(rating = ordered(as.vector(ysim)),
                             firm = rep(1:n, 8),
                  year = paste0("year", rep(1:8, each = n)),
                  do.call(rbind,xsim), BSEC = rep(BSEC,8))

colnames(data_sim_long) <- c("rating", "firm_id", "year", "LR", "LEV", "PR",
                        "RSIZE", "BETA", "BSEC")

str(data_sim_long)
data_cr_panel <- data_sim_long

save(data_cr_panel,
     file = "../data/data_cr_panel.rda", compress='xz')



load("../data/data_cr_panel.rda")
library(mvord)
res_probit_data_sim_long <- mvord(
  formula = MMO(rating, firm_id, year) ~ 0 +  LR + LEV + PR + RSIZE +  BETA,
  data = data_cr_panel,
  error.structure = cor_ar1(~ BSEC),
  coef.constraints = c(rep(1,4), rep(2,4)),
  threshold.constraints = c(rep(1,8)),
  link = mvprobit(),
  control = mvord.control(solver = "BFGS",  solver.optimx.control = list(trace = TRUE)))
summary(res_probit_data_sim_long)
res_probit_data_sim_long$rho$runtime

get_error_struct(res_probit_data_sim_long, type = "corr")

save(res_probit_data_sim_long, file = "res_probit_data_sim_long.rda")
table(res_probit_data_sim_long$rho$y[,1], marginal_predict(res_probit_data_sim_long, type = "class")[,1])
table(res_probit_data_sim_long$rho$y[,2], marginal_predict(res_probit_data_sim_long, type = "class")[,2])
table(res_probit_data_sim_long$rho$y[,3], marginal_predict(res_probit_data_sim_long, type = "class")[,3])
table(res_probit_data_sim$rho$y[,4], marginal_predict(res_probit_data_sim, type = "class")[,4])
res_probit_data_sim_long$rho$runtime
