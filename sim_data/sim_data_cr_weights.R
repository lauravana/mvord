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
dat <- data_ordinal#[data_ordinal$fyear %in% 2000:2013, ]
table(dat$SPR)
aggregate(SPR ~ fyear, data = dat, FUN = table)
aggregate(Moodys ~ fyear, data = dat, FUN = table)
aggregate(Fitch ~ fyear, data = dat, FUN = table)
#dat <- data_ordinal#[data_ordinal$fyear %in% 2001:2010, ]

# ## make ratings nice
levels(dat$SPR)
levels(dat$SPR)[1:2] <- "AAA/AA"
levels(dat$SPR)[5:6] <- "C"

table(dat$SPR)

# levels(dat$SPR)[2:4] <- "A"
# levels(dat$SPR)[3:5] <- "BBB"
# levels(dat$SPR)[4:6] <- "BB"
# levels(dat$SPR)[5:12] <- "B"
#levels(dat$SPR)[6:10] <- "CCC/C"

levels(dat$Fitch)
levels(dat$Fitch)[1:2] <- "AAA/AA"
levels(dat$Fitch)[5:6] <- "C"
table(dat$Fitch)
# levels(dat$Fitch)[2:4] <- "A"
# levels(dat$Fitch)[3:5] <- "BBB"
# levels(dat$Fitch)[4:6] <- "BB"
# levels(dat$Fitch)[5:12] <- "B"
#levels(dat$Fitch)[6:10] <- "CCC/C"

levels(dat$Moodys)
levels(dat$Moodys)[1:2] <- "Aaa/Aa"
levels(dat$Moodys)[6:7] <- "C"

table(dat$Moodys)

# levels(dat$Moodys)[2:4] <- "A"
# levels(dat$Moodys)[3:5] <- "Baa"
# levels(dat$Moodys)[4:6] <- "Ba"
# levels(dat$Moodys)[5:7] <- "B"
#levels(dat$Moodys)[6:7] <- "C"
#levels(dat$Moodys)[6:8] <- "Caa"
#at$Moodys[dat$Moodys == "C"] <- NA

# levels(dat$SPR)[1:2] <- "AAA/AA"
# levels(dat$Fitch)[1:2] <- "AAA/AA"
# levels(dat$Moodys)[1:2] <- "Aaa/Aa"
# str(dat)
#
# levels(dat$SPR)[5:6] <- "B/C"
# levels(dat$Fitch)[5:6] <- "B/C"
# levels(dat$Moodys)[6:7] <- "C"
str(dat)

dat$SPR <- factor(dat$SPR, labels = c("A","B","C","D","E"))#paste0("RC", 1:6, "_R1"))
dat$Fitch <- factor(dat$Fitch, labels = c("A","B","C","D","E"))#paste0("RC", 1:6, "_R2"))
dat$Moodys <- factor(dat$Moodys, labels = c("F", "G","H","I","J","K"))#paste0("RC", 1:7, "_R3"))
yyy <- cbind(dat$SPR, dat$Fitch, dat$Moodys)

yyy <- cbind(dat$SPR, dat$Fitch, dat$Moodys)
head(cbind.data.frame(dat$SPR, dat$Fitch, dat$Moodys))
rater4 <- rowSums(yyy>= 3, na.rm = T) >= 2 | (rowSums(yyy>= 3, na.rm =
                                                        T) == 1 &  rowSums(!is.na(yyy)) == 1)
head(rater4)
#FALSE investment
#TRUE speculative
dat$rater4 <- factor(rater4, labels = c("L", "M"), ordered = T)
head(dat$rater4)

dat <- dat[!is.na(dat$SPR) & !is.na(dat$Fitch) & !is.na(dat$Moodys),]
nrow(dat)

table(dat$SPR)
table(dat$Moodys)
table(dat$Fitch)
table(dat$rater4)


#dat$R8[dat$R8 == 0] <- .Machine$double.eps
#dat$R24[dat$R24 == 0] <- .Machine$double.eps
#dat$R31[dat$R31 == 0] <- .Machine$double.eps
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

weightsSPR <- 1/table(dat$SPR) * 100
weightsMoodys <- 1/table(dat$Moodys) * 100
weightsFitch <- 1/table(dat$Fitch) * 100
weightsrater4 <- 1/table(dat$rater4) * 100

dat$weights <- (weightsSPR[match(dat$SPR, names(weightsSPR))] +
                  weightsMoodys[match(dat$Moodys, names(weightsMoodys))] +
                  weightsFitch[match(dat$Fitch, names(weightsFitch))]+
                  weightsrater4[match(dat$rater4, names(weightsrater4))])/4

dat$weights

summary(dat)

fit <- FALSE
if(fit == TRUE) {
  library(mvord)
  # res_probit_real_dat_weights <- mvord(
  #   formula = MMO2(SPR, Fitch, Moodys, rater4) ~ 0 + R3 + R5 + R13 + R18 + R20 + RSIZE + SIGMA + BETA,
  #   #0 + R4 + R6 + R13 + R18 + R23 + RSIZE + SIGMA + BETA,
  #   data = dat,
  #   weights = "weights",
  #   # response.levels = rep(list(rev(paste0("RC", 1:6, "_R1"))), 10),
  #   # response.names = c(2002:2011),
  #   error.structure = cor_general(~ 1),
  #   #threshold.constraints = rep(1, 10),
  #   #coef.constraints = c(rep(1, 5), rep(2, 5)),
  #   link = mvprobit(),
  #   control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))
  #
  # summary(res_probit_real_dat_weights)
  # res_probit_real_dat_weights$rho$runtime
  # table(res_probit_real_dat_weights$rho$y[,1], marginal.predict(res_probit_real_dat_weights, type = "class")[,1])
  # table(res_probit_real_dat_weights$rho$y[,2], marginal.predict(res_probit_real_dat_weights, type = "class")[,2])
  # table(res_probit_real_dat_weights$rho$y[,3], marginal.predict(res_probit_real_dat_weights, type = "class")[,3])
  # table(res_probit_real_dat_weights$rho$y[,4], marginal.predict(res_probit_real_dat_weights, type = "class")[,4])


  ####
  res_probit_real_dat_weights_constr <- mvord(
    formula = MMO2(SPR, Fitch, Moodys, rater4) ~ 0 + R5 + R13 + R20 + RSIZE + BETA,
    #0 + R4 + R6 + R13 + R18 + R23 + RSIZE + SIGMA + BETA,
    data = dat,
    weights = "weights",
    # response.levels = rep(list(rev(paste0("RC", 1:6, "_R1"))), 10),
    # response.names = c(2002:2011),
    error.structure = cor_general(~ 1),
    #threshold.constraints = rep(1, 10),
    #coef.constraints = c(rep(1, 5), rep(2, 5)),
    link = mvprobit(),
    coef.constraints = cbind(#c(1,1,1,1), #R3
      c(1,1,1,1), #R5
      c(1,2,3,4), #R13
      #c(1,1,2,1), #R18
      c(1,1,1,1), #R20
      c(1,1,1,2), #RSIZE
      #c(1,2,1,1), #SIGMA
      c(1,1,2,3) #BETA
    ),
    threshold.constraints = c(1,1,2,3),
    control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))

  summary(res_probit_real_dat_weights_constr)
  res_probit_real_dat_weights_constr$rho$runtime



  # res_logit_real_dat_weights <- mvord(
  #   formula = MMO2(SPR, Fitch, Moodys, rater4) ~ 0 + R5 + R13 + R20 + RSIZE + BETA,
  #   #0 + R4 + R6 + R13 + R18 + R23 + RSIZE + SIGMA + BETA,
  #   data = dat,
  #   weights = "weights",
  #   # response.levels = rep(list(rev(paste0("RC", 1:6, "_R1"))), 10),
  #   # response.names = c(2002:2011),
  #   error.structure = cor_general(~ 1),
  #   #threshold.constraints = rep(1, 10),
  #   #coef.constraints = c(rep(1, 5), rep(2, 5)),
  #   link = mvlogit(),
  #   coef.constraints = cbind(#c(1,1,1,1), #R3
  #                            c(1,1,1,1), #R5
  #                            c(1,2,3,4), #R13
  #                            #c(1,1,2,1), #R18
  #                            c(1,1,2,3), #R20
  #                            c(1,1,1,2), #RSIZE
  #                            #c(1,2,1,1), #SIGMA
  #                            c(1,1,2,1) #BETA
  #   ),
  #   threshold.constraints = c(1,1,2,3),
  #   control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))
  # summary(res_logit_real_dat_weights)
  # save(res_logit_real_dat_weights, file = "res_logit_real_dat_weights_v4.rda")
  # load("res_logit_real_dat_weights_v4.rda")
  #
  # table(res_logit_real_dat_weights$rho$y[,1], marginal.predict(res_logit_real_dat_weights, type = "class")[,1])
  # table(res_logit_real_dat_weights$rho$y[,2], marginal.predict(res_logit_real_dat_weights, type = "class")[,2])
  # table(res_logit_real_dat_weights$rho$y[,3], marginal.predict(res_logit_real_dat_weights, type = "class")[,3])
  # table(res_logit_real_dat_weights$rho$y[,4], marginal.predict(res_logit_real_dat_weights, type = "class")[,4])
  # res_logit_real_dat_weights$rho$runtime

}

#fit with constraints

#####################################################
win_lev <- 0.01
win_lev <- 0.0

winsorize <- function(x, win_lev){
  x[x>quantile(x,1-win_lev/2)] <- quantile(x,1-win_lev/2)
  x[x<quantile(x,win_lev/2)] <- quantile(x,win_lev/2)
  x
}



# #R3
# x <- dat$R3 #t
# x <- dat$R5 #lognormal
# x <- dat$R13 #lognormal
# x <- dat$R18 #t
# x <- dat$R20 #t
# x <- dat$RSIZE #
# x <- dat$BETA #t
# x <- dat$SIGMA #lognormal
#
#
# winsorize <- function(x, win_lev){
#   x[x>quantile(x,1-win_lev/2)] <- quantile(x,1-win_lev/2)
#   x[x<quantile(x,win_lev/2)] <- quantile(x,win_lev/2)
#   x
# }

# plot(density(x))
# any(x<0)
# win_lev <- 0.01
# x[x>quantile(x,1-win_lev/2)] <- quantile(x,1-win_lev/2)
# x[x<quantile(x,win_lev/2)] <- quantile(x,win_lev/2)
# plot(density(x))
#
#
# ## positive  R5, R13, R18, SIGMA,
# ## unbounded R3, R20, RSIZE, BETA
#
# #x <- dat$R3
# rm(x1)
# x1 <- MASS::fitdistr(x, densfun = "beta")
# x1$loglik
# rm(x2)
# x2 <- MASS::fitdistr(x, densfun = "cauchy")
# x2$loglik
# rm(x3)
# x3 <- MASS::fitdistr(x, densfun = "chi-squared")
# x3$loglik
# rm(x4)
# x4 <- MASS::fitdistr(x, densfun = "exponential")
# x4$loglik
# rm(x5)
# x5 <- MASS::fitdistr(x, densfun = "f")
# x5$loglik
# rm(x6)
# x6 <- MASS::fitdistr(x, densfun = "gamma")
# x6$loglik
# #rm(x7)
# # x7 <- MASS::fitdistr(x, densfun = "geometric")
# # x7$loglik
# rm(x8)
# x8 <- MASS::fitdistr(x, densfun = "log-normal")
# x8$loglik
# rm(x9)
# x9 <- MASS::fitdistr(x, densfun = "lognormal")
# x9$loglik
# rm(x10)
# x10 <- MASS::fitdistr(x, densfun = "logistic")
# x10$loglik
# #rm(x11)
# # x11 <- MASS::fitdistr(x, densfun = "negative binomial")
# # x11$loglik
# rm(x12)
# x12 <- MASS::fitdistr(x, densfun = "normal")
# x12$loglik
# #rm(x13)
# # x13 <- MASS::fitdistr(x, densfun = "Poisson")
# # x13$loglik
# rm(x14)
# x14 <- MASS::fitdistr(x, densfun = "t")
# x14$loglik
# rm(x15)
# x15 <- MASS::fitdistr(x, densfun = "weibull")
# x15$loglik
#
# which.max(c(x1$loglik, x2$loglik, x3$loglik, x4$loglik, x5$loglik, x6$loglik, x8$loglik, x9$loglik,
#             x10$loglik, x12$loglik, x14$loglik, x15$loglik))
#
#
# est_tab <- rbind(MASS::fitdistr(winsorize(dat$R5, win_lev), densfun = "t")$estimate,
#                  c(MASS::fitdistr(winsorize(dat$R8, win_lev), densfun = "weibull")$estimate, NA),
#                  MASS::fitdistr(winsorize(dat$R24, win_lev), densfun = "t")$estimate,
#                  MASS::fitdistr(winsorize(dat$R31, win_lev), densfun = "t")$estimate,
#                  c(MASS::fitdistr(winsorize(dat$R39, win_lev), densfun = "logistic")$estimate, NA),
#                  c(MASS::fitdistr(winsorize(dat$RSIZE, win_lev), densfun = "normal")$estimate, NA),
#                  c(MASS::fitdistr(winsorize(dat$SIGMA, win_lev), densfun = "lognormal")$estimate, NA)
# )
#
#
# n <- nrow(dat)
#
# sim_dat <- data.frame(firm_id = 1:n)
# sim_dat$R5 <- rt(n, df = est_tab[1,3])* est_tab[1,2] + est_tab[1,1]
# sim_dat$R8 <- rweibull(n, shape = est_tab[2,1], scale = est_tab[2,2])
# sim_dat$R24 <- rt(n, df = est_tab[3,3])* est_tab[3,2] + est_tab[3,1]
# sim_dat$R31 <- rt(n, df = est_tab[4,3])* est_tab[4,2] + est_tab[4,1]
# sim_dat$R39 <- rlogis(n, location = est_tab[5,1], scale = est_tab[5,2])
# sim_dat$RSIZE <- rnorm(n, mean = est_tab[6,1], sd = est_tab[6,2])
# sim_dat$SIGMA <- rlnorm(n, meanlog = est_tab[7,1], sdlog = est_tab[7,2])
# #sim_dat$SIGMA <- rlnorm(n, meanlog = est_tab[10,1], sdlog = est_tab[10,2])
# #sim_dat$MB <- rt(n, df = est_tab[11,3])* est_tab[11,2] + est_tab[11,1]
#

##################################################
################### NEW CODE #####################
##################################################
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
n




sim_dat <- data.frame(firm_id = 1:n)
#sim_dat$R3    <- rt(n, df = est_tab[[1]][3]) * est_tab[[1]][2] + est_tab[[1]][1]
sim_dat$R5    <- rlnorm(n, meanlog = est_tab[[2]][1], sdlog = est_tab[[2]][2])
sim_dat$R13   <- rt(n, df = est_tab[[3]][3]) * est_tab[[3]][2] + est_tab[[3]][1]
#sim_dat$R18   <- rt(n, df = est_tab[[4]][3]) * est_tab[[4]][2] + est_tab[[4]][1]
sim_dat$R20   <- rt(n, df = est_tab[[5]][3]) * est_tab[[5]][2] + est_tab[[5]][1]
sim_dat$RSIZE <- rt(n, df = est_tab[[6]][3]) *est_tab[[6]][2] + est_tab[[6]][1]
#sim_dat$SIGMA <- rlnorm(n, meanlog = est_tab[[7]][1], sdlog = est_tab[[7]][2])
sim_dat$BETA  <- rt(n, df = est_tab[[8]][3]) * est_tab[[8]][2] + est_tab[[8]][1]
sim_dat$R13[sim_dat$R13 < 0] <- 0
sim_dat$R13[sim_dat$R13 == 0] <- runif(length(sim_dat$R13[sim_dat$R13 == 0]), 0,4)
#sim_dat$R18[sim_dat$R18 < 0] <- 0
#sim_dat[,-1] <- apply(sim_dat[,-1], 2, winsorize, win_lev=0.00)
summary(sim_dat)

coef(res_logit_real_dat_weights)




# beta_tab <- cbind(#c(-0.025, -0.025, -0.025, -0.025),#R3
#                   c(0.15, 0.15, 0.15, 0.15),#R5
#                   c(0.28, 0.3, 0.24, 0.32),#R13
#                   c(2, 2, 1.8, 2),#R18
#                   c(-1.2, -1.2, -1.3, -1.7),#R20
#                   c(-0.75, -0.75, -0.75, -0.85),#RSIZE
#                   #c(45, 30, 45, 45),#SIGMA
#                   c(1.3, 1.3, 1, 1.3)#BETA
# )

beta_tab <- cbind(#c(-0.025, -0.025, -0.025, -0.025),#R3
  c(0.05, 0.05, 0.05, 0.05),#R5
  c(0.19, 0.19, 0.16, 0.26),#R13
  #c(2, 2, 1.8, 2),#R18
  c(-1, -1, -1, -1),#R20
  c(-0.45, -0.45, -0.45, -0.52),#RSIZE
  #c(45, 30, 45, 45),#SIGMA
  c(0.72, 0.72, 0.66, 0.75)#BETA
)

beta_tab <- beta_tab*4

# beta_tab <- cbind(#R3 = c(0,0,0,0),
#   R5 = c(-0.05,-0.03,-0.04,-0.045),
#   R8 = c(0.5, 0.5, 0.5, 0.5),
#   R24 = c(0.2, 0.16, 0.08, 0.23),
#   R31 = c(4, 4, 4, 4),
#   R39 = c(-4, -5.4, -5.2, -4.2),
#   RSIZE = c(-0.66, -0.6, -0.7, -0.75),
#   SIGMA = c(30, 30, 27, 30)
# )


# beta_tab <- cbind(#R3 = c(0,0,0,0),
#   R9 = c(0.2,0,0,0),
#   R12 = c(0.8, 0.9, 0.9, 1.1),
#   R18 = c(1, 1.1, 0.9, 0.7),
#   R20 = c(-0.9, -0.7, -0.8, -1),
#   R23 = c(-0.3, - 0.8, -0.6, -0.5),
#   R24 = c(-1, -1, -1, -1),
#   BETA = c(0.7, 0.6, 0.6, 0.6),
#   MB = c(-0.4, -0.4, -0.4, -0.4)
# )

#beta_tab <- beta_tab * 2#pi/sqrt(3)

# theta_tab <- list(R1 = c(5,7.5,10,12),
#                   R1 = c(5,7.5,10,12),
#                   R3 = c(5,7,9,11,13),
#                   R4 = c(8.5))
thresholds(res_logit_real_dat_weights)
thresholds(res_probit_real_dat_weights_constr)

# theta_tab <- list(R1 = c(2.8,3.8,5,6.2),
#                   R2 = c(2.8,3.8,5,6.2),
#                   R3 = c(2.6,3.6,4.6,5.5,6.5),
#                   R4 = c(4.5))

theta_tab <- list(R1 = c(13,15.5,18.5,22),
                  R2 = c(13,15.5,18.5,22),
                  R3 = c(12.5,15,18,20,22),
                  R4 = c(18.5))


theta_labels <- list(R1 = LETTERS[1:5],
                     R2 = LETTERS[1:5],
                     R3 = LETTERS[6:11],
                     R4 = LETTERS[12:13])

library(mvord)
error_structure(res_probit_real_dat_weights_constr)[[1]]
sigma <- matrix(c(1,   0.89, 0.9, 0.93,
                  0.89,   1, 0.86, 0.92,
                  0.9, 0.86,   1, 0.88,
                  0.93, 0.92, 0.88,   1), nrow = 4)

# #is.positive.semi.definite(sigma)
# sigma <- matrix(c(1,   0.9, 0.9, 0.9,
#                   0.9,   1, 0.85, 0.9,
#                   0.9, 0.85,   1, 0.9,
#                   0.9, 0.9, 0.9,   1), nrow = 4)

library(copula)
sigma_cor <- cov2cor(sigma)
sd_sigma <-  sqrt(diag(sigma))
myCop <- tCopula(sigma_cor[lower.tri(sigma_cor)],
                 dim = 4, df = 8,
                 df.fixed = TRUE,
                 dispstr = "un")
getSigma(myCop)
ui <- rCopula(n, myCop)
errors <-  t(t(qlogis(ui)) * sd_sigma)#/(pi/sqrt(3))
mean(errors)
sd(errors)
summary(errors)
plot(density(errors[,1]))
#sd(errors[,1]/(pi/sqrt(3)))

ysim <- sapply(1:4, function(j) {
  pred <- as.matrix(sim_dat[, 1+seq_len(5)]) %*% beta_tab[j,] + errors[,j]
  cut(pred, breaks = c(-Inf,theta_tab[[j]], Inf), labels = theta_labels[[j]])
})

summary(ysim)

pred <- sapply(1:4, function(j) {
  as.matrix(sim_dat[, 1+seq_len(5)]) %*% beta_tab[j,]
})
summary(pred)



yt <- sapply(1:4, function(j) {
  pred <- as.matrix(sim_dat[, 1+seq_len(5)]) %*% beta_tab[j,] #+ errors[,j]
  cut(pred, breaks = c(-Inf,theta_tab[[j]], Inf), labels = theta_labels[[j]])
})

summary(yt)

###missings (makes code faster)
sum(!is.na(data_ordinal$SPR))/nrow(data_ordinal)
sum(!is.na(data_ordinal$Fitch))/nrow(data_ordinal)
sum(!is.na(data_ordinal$Moodys))/nrow(data_ordinal)

ysim[sample(1:n, round(n * (1-0.95))),1] <- NA
ysim[sample(1:n, round(n * (1-0.7))),2] <- NA
ysim[sample(1:n, round(n * (1-0.5))),3] <- NA

summary(ysim)


data_cr <- cbind(rater1 = ordered(ysim[,1]),
                 rater2 = ordered(ysim[,2]),
                 rater3 = ordered(ysim[,3]),
                 rater4 = ordered(ysim[,4]),
                 sim_dat)

colnames(data_cr) <- c("rater1", "rater2", "rater3", "rater4", "firm_id",  "LR", "LEV","PR",
                       "RSIZE", "BETA")

str(data_cr)
save(data_cr,
    file = "../data/data_cr.rda", compress='xz')




library(mvord)
res_probit_data_sim <- mvord(
  formula = MMO2(rater1, rater2, rater3, rater4) ~ 0 + LR + LEV +
    PR  + RSIZE + BETA,
  data = data_cr,
  error.structure = cor_general(~ 1),
  coef.constraints = cbind(#c(1,1,1,1), #R3
    c(1,1,1,1), #R5
    c(1,2,3,4), #R13
    #c(1,1,2,1), #R18
    c(1,1,1,1), #R20
    c(1,1,1,2), #RSIZE
    #c(1,2,1,1), #SIGMA
    c(1,1,2,3) #BETA
  ),
  threshold.constraints = c(1,1,2,3),
  link = mvprobit(),
  control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))

summary(res_probit_data_sim)
res_probit_data_sim$rho$runtime

res_probit_data_sim_all <- mvord(
  formula = MMO2(rater1, rater2, rater3, rater4) ~ 0 + LR + LEV +
    PR  + RSIZE + BETA,
  data = data_cr,
  link = mvprobit(),
  control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))
res_probit_data_sim_all$rho$runtime
table(res_probit_data_sim_all$rho$y[,1], marginal_predict(res_probit_data_sim_all, type = "class")[,1])
table(res_probit_data_sim_all$rho$y[,2], marginal_predict(res_probit_data_sim_all, type = "class")[,2])
table(res_probit_data_sim_all$rho$y[,3], marginal_predict(res_probit_data_sim_all, type = "class")[,3])
table(res_probit_data_sim_all$rho$y[,4], marginal_predict(res_probit_data_sim_all, type = "class")[,4])




res_logit_data_sim <- mvord(
  formula = MMO2(rater1, rater2, rater3, rater4) ~0 + LR + LEV +
    PR  + RSIZE + BETA,
  data = data_cr,
  error.structure = cor_general(~ 1),
  coef.constraints = cbind(#c(1,1,1,1), #R3
    c(1,1,1,1), #R5
    c(1,2,3,4), #R13
    #c(1,1,2,1), #R18
    c(1,1,1,1), #R20
    c(1,1,1,2), #RSIZE
    #c(1,2,1,1), #SIGMA
    c(1,1,2,3) #BETA
  ),
  threshold.constraints = c(1,1,2,3),
  link = mvlogit(),
  control = mvord.control(solver = "newuoa", solver.optimx.control = list(trace = TRUE)))
summary(res_logit_data_sim)
save(res_logit_data_sim, file = "res_logit_data_sim.rda")
table(res_logit_data_sim$rho$y[,1], marginal_predict(res_logit_data_sim, type = "class")[,1])
table(res_logit_data_sim$rho$y[,2], marginal_predict(res_logit_data_sim, type = "class")[,2])
table(res_logit_data_sim$rho$y[,3], marginal_predict(res_logit_data_sim, type = "class")[,3])
table(res_logit_data_sim$rho$y[,4], marginal_predict(res_logit_data_sim, type = "class")[,4])
res_logit_data_sim$rho$runtime
