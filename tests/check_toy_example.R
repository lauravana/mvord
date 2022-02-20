library(mvord)
library(MASS)


#data(data_toy_example)
tolerance <- 1e-6

mult.obs <- 2
sigma <- matrix(c(1,0.8,0.8,1), ncol = 2)
betas <- list(c(0.8,-0.5),
              c(0.8,-0.5))
thresholds <- list(c(-1,1),c(-1,1))
nobs <- 100
suppressWarnings(RNGversion("3.5.0"))
set.seed(2017)
errors <-  mvrnorm(n = nobs, mu = rep(0, mult.obs), Sigma = sigma)

X1 <- rnorm(nobs, 0, 1)
X2 <- rnorm(nobs, 0, 1)

pred <- cbind(X1, X2)

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
colnames(data_toy_example) <- c("Y1","Y2", "X1", "X2")
w <- c(rep(1/20, 20), rep(1/30, 30), rep(1/20, 20), rep(1/30, 30))



# convert data_toy_example into long format
df <- cbind.data.frame("i" = rep(1:100,2),
                       "j" = rep(1:2,each = 100),
                       "Y" = c(data_toy_example$Y1, data_toy_example$Y2),
                       "X1" = rep(data_toy_example$X1,2),
                       "X2" = rep(data_toy_example$X2,2),
                       "f1" = factor(sample(rep(data_toy_example$Y2,2)),
                                     ordered = FALSE),
                       "f2" = factor(rep(data_toy_example$Y1,2), ordered = FALSE),
                       w    = rep(w,2))
df$X3 <- cut(df$X2, c(-Inf, -0.2, 0.2, Inf))



# library(ROI)
# ROI_solver <- function(starting.values, objFun, control){
#     n <- length(starting.values)
#     op <- OP(objective = F_objective(objFun, n = n),
#              bounds = V_bound(li = seq_len(n), lb = rep.int(-Inf, n)))
#     optRes <- ROI_solve(op, solver = "nlminb",
#         control = c(list(start = starting.values),
#                     control))
#     list(optpar    = optRes$solution,
#          objective = optRes$objval) # a vector of length equal to number of parameters to optimize
# }
#
#
#
res <- mvord::mvord(formula = MMO(Y) ~  0 + X1 + X2,
                     data = df,
                     link = mvprobit(),
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,1),
                     coef.constraints = c(1,1),
                     contrasts = list(f1 = function(x)
                       contr.treatment(nlevels(df$f1), base = 1),
                       f2 = "contr.sum"),
                     control= mvord.control(solver="BFGS",se=TRUE))
options(digits = 22)

res.summary <- summary(res, short = FALSE)
# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.96257386663519672876,  1.03347036873223707687, -0.96257386663519672876,  1.03347036873223707687), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.63801035062404309883, -0.42672524816263474046), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.85426672822122684536), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.16613952738605441972, 0.15004482537617935822, 0.16613952738605441972, 0.15004482537617935822), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.13576747301357894315, 0.13643665802446261481), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.062440889990360744222), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -134.90867383086322207, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 280.3436634512001433, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 294.05508548271637892, tolerance = tolerance))


res2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,1),
                     coef.constraints = list(matrix(rep(1,4), ncol = 1), matrix(rep(1,4), ncol = 1)))

mvord:::check(all.equal(res$beta, res2$beta, tolerance = tolerance))
mvord:::check(all.equal(res$sebeta, res2$sebeta, tolerance = tolerance))
########################################################################
## No coefficients
res <- mvord::mvord(formula = MMO(Y) ~ -1,
                      data = df,
                      link = mvprobit(),
                      error.structure = cor_general(~1),
                      threshold.constraints = c(1,1),
                      control= mvord.control(solver="BFGS",se=TRUE))
res.summary <- summary(res, short = FALSE)
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.75479706538110091785,  0.86086304364935783973, -0.75479706538110091785,  0.86086304364935783973), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(numeric(0)), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.90579517144642240911), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.13105603927553965971, 0.13360642334434202905, 0.13105603927553965971, 0.13360642334434202905), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(numeric(0)), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.038855435949034601573), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -153.66397119528727444, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 313.51350940088383368, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 321.57073678022845797, tolerance = tolerance))

#polychor
res <- mvord::mvord(formula = MMO(Y) ~ 1,
                     data = df,
                     link = mvprobit(),
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,1),
                     control= mvord.control(solver="BFGS",se=TRUE))

res.summary <- summary(res, short = FALSE)
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(0.0000000000000000000, 1.6157607546978449697, 0.0000000000000000000, 1.6157607546978449697), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.73786225358709867095, 0.77255634837212872057), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.90637393266264265623), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.00000000000000000000, 0.15607113684416776267, 0.00000000000000000000, 0.15607113684416776267), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.13500906208906049422, 0.14063945843173727979), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.038916052675185497378), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -153.5640565887688922, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 315.46144651087109878, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 326.31632228582151356, tolerance = tolerance))
#######################################################################
## cor_general(~factor)
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_general(~f2),
                     threshold.constraints = c(1,1),
                     coef.constraints = c(1,1),
                     contrasts = list(f2 = "contr.sum"))

res.summary <- summary(res, short = FALSE)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate,
  c(-0.90552506790519915469,  1.00420547521745429087, -0.90552506790519915469,  1.00420547521745429087), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate,
  c(0.64793126338493944871, -0.42893128334048874484), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate,
  c(0.73667859723415496376, 0.91829826305234418804, 0.83710611971260706632), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`,
  c(0.18265647646915181279, 0.17467164336594759311, 0.18265647646915181279, 0.17467164336594759311), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`,
  c(0.14153585035881874332, 0.13989943737309898375), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`,
  c(0.17823839430849069965, 0.06885865619995423792, 0.13103608299872890330), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -134.17046176709035876, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 283.39468697504094052, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 303.00349482656417877, tolerance = tolerance))

##################################################################################
res <- mvord::mvord(formula = MMO(Y, i, j) ~  0 + X1 + X2,
                     data = df,
                     link = mvlogit(df = 8L),
                     error.structure = cor_general(~1))
                     #threshold.constraints = c(1,1),
                     #coef.constraints = c(1,1))
res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
# mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-1.6170817306633420, 1.7855897338762188, -1.6170817306633420, 1.7855897338762188), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$coefficients$Estimate, c(1.07242200717115987, 1.07242200717115987, -0.76715925377701732, -0.76715925377701732), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.85317690560688531), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.28997405649500174, 0.27427389826231802, 0.28997405649500174, 0.27427389826231802), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.24111402270822993, 0.24111402270822993, 0.24156664773225886, 0.24156664773225886), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.063316529381183581), tolerance = tolerance))
# mvord:::check(all.equal(logLik(res), -135.41665313840898, tolerance = tolerance))
# mvord:::check(all.equal(AIC(res), 281.35962206629165, tolerance = tolerance))
# mvord:::check(all.equal(BIC(res), 295.07104409780789, tolerance = tolerance))

##################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     #index = c("i", "j"),
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cov_general(~1),
                     threshold.constraints = c(1,1),
                     threshold.values = list(c(-1,NA),
                                             c(-1,NA)),
                     coef.constraints = c(1,1))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-1.000000000000000000000,  1.082667980083374503764, -1.000000000000000000000,  1.082667980083374503764), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.6883494614209317852271, -0.4540651765505133163892), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.8561536378031743277361, 1.0019345838689188710191, 1.0898281993157663549709), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.0000000000000000000000, 0.2492604341815170820862, 0.0000000000000000000000, 0.2492604341815170820862), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1512803674185349001036, 0.1570413444091093568833), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.06182222146981963123435, 0.18518796068392143205905, 0.19759328972385112321852), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -134.6391112878561102661, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 282.0441800225207202857, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 298.6729258905298252103, tolerance = tolerance))
########################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     #index = c("i", "j"),
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_equi(~1),
                     threshold.constraints = c(1,1),
                     coef.constraints = cbind(c(1,1),c(1,2)))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.9627566523584676350112,  1.0338802913596250032668, -0.9627566523584676350112,  1.0338802913596247812222), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.6382046741095841468905, -0.4467656111955374820255, -0.4075015551913261924177), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$error.structure$Estimate, c(1.272900571254629964457), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.1670386574395647250046, 0.1520309160654948199554, 0.1670386574395647250046, 0.1520309160654948199554), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1367049280926301846328, 0.1590087966266460206555, 0.1385097918969951935608), tolerance = tolerance))
# mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.2322161184793045674013), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -134.8432321319771745038, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 282.4524217107628487611, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 299.0811675787719536856, tolerance = tolerance))

res2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                      data = df,
                      #index = c("i", "j"),
                      link = mvprobit(),
                      control = mvord.control(solver = "BFGS"),
                      error.structure = cor_equi(~1),
                      threshold.constraints = c(1,1),
                      coef.constraints = list(X2 = cbind(c(1,1,0,0), c(0,0,1,1)),
                        X1 = matrix(rep(1,4), ncol = 1)))

mvord:::check(all.equal(res$beta, res2$beta, tolerance = tolerance))
mvord:::check(all.equal(res$sebeta, res2$sebeta, tolerance = tolerance))

########################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_ar1(~ 1 + X1),
                     threshold.constraints = c(1,1),
                     coef.constraints = c(1,1))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.9572295115261755249492,  1.0374679396833161870717, -0.9572295115261755249492,  1.0374679396833161870717), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.6530420867428219366957, -0.4227313679354475217664), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(1.298915541971320308789, 0.293508923629326790028), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.1655828816811988069002, 0.1530695261669812456962, 0.1655828816811988069002, 0.1530695261669812456962), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1358401620200463122412, 0.1382021771627104100855), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.2451746600755327076815, 0.2941698864686516090572), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -133.938827675498004055, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 280.6436127978045078635, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 297.272358665813612788, tolerance = tolerance))

res2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                      data = df,
                      #index = c("i", "j"),
                      link = mvprobit(),
                      control = mvord.control(solver = "BFGS"),
                      #se = TRUE,
                      error.structure = cor_ar1(~1 + X1),
                      threshold.constraints = c(1,1),
                      coef.constraints = list(matrix(rep(1,4), ncol = 1), matrix(rep(1,4),
                        ncol = 1)))

mvord:::check(all.equal(res$beta, res2$beta, tolerance = tolerance))
mvord:::check(all.equal(res$sebeta, res2$sebeta, tolerance = tolerance))

########################################################################################
res <- mvord(formula = MMO2(Y1,Y2) ~ 0 + X1 + X2,
                      data = data_toy_example,
                      link = mvprobit(),
                      control = mvord.control(solver = "BFGS"),
                      #se = TRUE,
                      error.structure = cor_general(~1),
                      threshold.constraints = c(1,1),
                      coef.constraints = cbind(c(1,2),c(NA,1)))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.8989294252736504953205,  0.9841336594057152886705, -0.8989294252736503842982,  0.9841336594057153996928), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.68194485928859438494953,  0.46837658310856356003171, -0.05217980198794509860694), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.8903906571446290607597), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.1624607776227008359982, 0.1598165721396140603883, 0.1624607776227008359982, 0.1598165721396140603883), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.17496189292260000414103, 0.15710600669362786985239, 0.08979821541610534529898), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.0631755330240382911855), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -137.6494615446075613363, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 288.064880536023622426, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 304.6936264040327273506, tolerance = tolerance))


res2 <- mvord::mvord(formula = MMO2(Y1,Y2) ~ 0 + X1 + X2,
                       data = data_toy_example,
                       link = mvprobit(),
                       control = mvord.control(solver = "BFGS"),
                       #se = TRUE,
                       error.structure = cor_general(~1),
                       threshold.constraints = c(1,1),
                      coef.constraints = list(X1 = cbind(c(1,1,0,0), c(0,0,1,1)),
                        X2 = matrix(c(rep(0,2),rep(1,2)), ncol = 1)))

mvord:::check(all.equal(res$beta, res2$beta, tolerance = tolerance))
mvord:::check(all.equal(res$sebeta, res2$sebeta, tolerance = tolerance))

########################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + offset(X2),
                      data = df,
                      #index = c("i", "j"),
                      link = mvprobit(),
                      control = mvord.control(solver = "BFGS"),
                      #se = TRUE,
                      error.structure = cor_general(~1),
                      threshold.constraints = c(1,2),
                      coef.constraints = list(X1 = cbind(c(1,1,0,0), c(0,0,1,1))))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.7672652623027305107684,  0.8919048301000485068357, -0.7967534462687086982413,  0.8729408047614937160574), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.5337008195901541407480, 0.3438557606717420056519), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.9408085278204662005308), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.1202334189234334493879, 0.1305894297520299174309, 0.1117234735991832761393, 0.1366223811642423013879), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1363643050543318635537, 0.1063416779401736012023), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.03369720676949919391241), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -194.3258468555947047207, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 403.705457152049632441, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 423.314265003572927526, tolerance = tolerance))

res2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1,
                     offset = list(df$X2[1:100], df$X2[101:200]),
                     data = df,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,2),
                     coef.constraints = list(X1 = cbind(c(1,1,0,0), c(0,0,1,1))))

mvord:::check(all.equal(res$beta, res2$beta, tolerance = tolerance))
mvord:::check(all.equal(res$sebeta, res2$sebeta, tolerance = tolerance))

res3 <- mvord::mvord(formula = MMO2(Y1,Y2) ~ 0 + X1 + offset(X2),
                       data = data_toy_example,
                       link = mvprobit(),
                       control = mvord.control(solver = "BFGS"),
                       #se = TRUE,
                       error.structure = cor_general(~1),
                       threshold.constraints = c(1,2),
                       coef.constraints = list(X1 = cbind(c(1,1,0,0), c(0,0,1,1))))
mvord:::check(all.equal(res$beta, res3$beta, tolerance = tolerance))
mvord:::check(all.equal(res$sebeta, res3$sebeta, tolerance = tolerance))

########################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
              data = df,
              #index = c("i", "j"),
              link = mvlogit(),
              #solver = "newuoa",
              #se = TRUE,
              error.structure = cor_general(~1),
              threshold.constraints = c(1,1),
              coef.constraints = cbind(c(NA, NA), c(1,2)),
              coef.values = cbind(c(1, 1), c(NA,NA)))
#rho <- res$rho


res.summary <- summary(res, short = FALSE)

options(digits = 22)

tolerance2 <- 1e-4
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-1.583288915548024533564,  1.758642501421955106622, -1.583288915548024533564,  1.758642501421955106622), tolerance = tolerance2))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(-0.7730608423417825170176, -0.7243032675736262859800), tolerance = tolerance2))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.8558944465934688050623), tolerance = tolerance2))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.2585108482987530376107, 0.2425159370850978601819, 0.2585108482987530376107, 0.2425159370850978601819), tolerance = tolerance2))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.2673280928737744588375, 0.2440398419906193439033), tolerance = tolerance2))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.06059919526643173237623), tolerance = tolerance2))
mvord:::check(all.equal(logLik(res)[[1]], -135.5210996495467554723, tolerance = tolerance2))
mvord:::check(all.equal(AIC(res), 281.568515088567210114, tolerance = tolerance2))
mvord:::check(all.equal(BIC(res), 295.2799371200834457341, tolerance = tolerance2))

# res1 <- mvord::mvord(formula = MMO(Y) ~ 1 + X1 + X2,
#                       data = df,
#                       index = c("i", "j"),
#                       link = mvprobit(),
#                       solver = "BFGS",
#                       se = TRUE,
#                       error.structure = cor_equi(~1),
#                       threshold.constraints = c(1,1),
#                       coef.constraints = list(Intercept = matrix(rep(1,6), ncol = 1),
#                                               X2 = cbind(c(1,1,1,0,0,0), c(0,0,0,1,1,1)), X1 = matrix(rep(1,6), ncol = 1)))


########################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X3,
                     data = df,
                     #index = c("i", "j"),
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     #se = T,
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,1),
                     coef.constraints = c(1,1))
res.summary <- summary(res, short = FALSE)

options(digits = 22)

mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-1.212706987250618206886,  0.758451539011639419563, -1.212706987250618206886,  0.758451539011639419563), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.58234630409261545214150,  0.06594956465951176682871, -0.63986755805831896370961), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.8579757656611409766256), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.2113140787381464491546, 0.1886079031355630819533, 0.2113140787381464491546, 0.1886079031355630819533), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1338742749509412854891, 0.3942507844634823288565, 0.2395962995094961522913), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.06239182988606906121731), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -136.0526219854373835005, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 284.8712014176832667545, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 301.499947285692371679, tolerance = tolerance))
########################################################################################
res <- mvord::mvord(formula = MMO(Y) ~ 1 + X1 * X2,
                     data = df,
                     #index = c("i", "j"),
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     #se = T,
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,1),
                     threshold.values = list(c(-1,NA),
                                             c(-1,NA)),
                     coef.constraints = c(1,1))
res.summary <- summary(res, short = TRUE)

options(digits = 22)

mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-1.000000000000000000000,  1.002695135787324165477), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(-0.02191358273802538475516,  0.62162104624596714597118, -0.42762435428426853745165, -0.10256725599972128792903), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.8535748052245384354109), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.0000000000000000000000, 0.2278165916815008962271), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1671006248895306467439, 0.1385408652510198845853, 0.1406676052717547309445, 0.1503815732862690546234), tolerance = tolerance))
#mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.063097735777156036), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.06309098922138993426056), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -134.6255603757719541136, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 282.0170781983524079806, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 298.6458240663615129051, tolerance = tolerance))
########################################################################################
df_NA <- df[-c(1,90:110),]


res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df_NA,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,2))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate,
  c(-1.0201263773742383911269,  1.1418036492524186176212, -0.9066282230038859024646,  0.9994511726914282467860), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.8373361664658567349306,  0.4982107530258053085248, -0.4474027944949834356692, -0.3539038996339812781500),
 tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.9106503078209118307029), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.1999664122961394840949, 0.2237909069393826111405, 0.1905965159570501277209, 0.1740306974019930341679), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.1889847005044749950198, 0.1663038426085674814647, 0.1876357516513920686840, 0.1448962459360401411335), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.07055356607861483497768), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -119.4526280916837492896, tolerance = tolerance))
mvord:::check(all.equal(AIC(res),  258.7052561833675099479, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 284.3969426996999345647, tolerance = tolerance))

#weights
df_NA$weights <- 0.5
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                    data = df_NA,
                    link = mvprobit(),
                    weights.name = "weights",
                    error.structure = cor_general(~1),
                    threshold.constraints = c(1,2))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate,
                        c(-1.0201283793713198377873,  1.1418048707357395521456, -0.9066161533747665313143,  0.9994581570190742558779), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.8373332969208180376341,  0.4982136050008418859392, -0.4474066175154422508875, -0.3539061566757252808024),
                        tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.9106469196938988819312), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.3999331442401471981007, 0.4475833313547766811880, 0.3811935020996289336104, 0.3480606771313752845209), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.3779704499331774103510, 0.3326083385001772918521, 0.3752721890043577146479, 0.2897924558296670061175), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.1411107414052604758226), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -59.72631403934757088336, tolerance = tolerance))
mvord:::check(all.equal(AIC(res),  139.2526280786951531354, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 164.944314595027606174, tolerance = tolerance))
##########################################################################################
df_23 <- df

##For the first response (first 100 obs), replace the level 3 by 2.
## Now the first response has 2 levels and the second one 3 levels
df_23$Y[which(df_23$Y[1:100] == 3)] <- 2

## Must be converted to integer. Ordered factor is not OK, as we have diff number of responses.
df_23$Y <- as.integer(df_23$Y)

res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df_23,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_general(~1),
                     threshold.constraints = c(1,2))

res.summary <- summary(res, short = FALSE)


mvord:::check(all.equal(res.summary$thresholds$Estimate,
                        c(-0.9724336904108278334391, -0.9773415102989740921302,  1.0122464288213952610107), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.8296955208405573101160,  0.5333263657626993170524, -0.3947931985275920929723, -0.3862597650564197349077),
                        tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.8529580497839671648919), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.2010157484670624383760, 0.1873095626070425123721, 0.1571342875452053644558), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.2240929429677874784588, 0.1463835117300714416810, 0.1893139616655838686210, 0.1397639987481304746364), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.1022697891863128538681), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -109.147957798495554016, tolerance = tolerance))
mvord:::check(all.equal(AIC(res),  235.6872199448172011671, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 258.3408737360180111864, tolerance = tolerance))

##########################################################################################
## make three responses
suppressWarnings(RNGversion("3.5.0"))
set.seed(100)
Y3 <- sample(1:3, nobs, replace = TRUE)
dat_3 <- cbind(Y3, data_toy_example)
res <- mvord::mvord(formula = MMO2(Y1, Y2, Y3) ~ 0,
                     data = dat_3,
                     link = mvprobit(),
                     error.structure = cor_general(~1),
                     control= mvord.control(solver="BFGS",se=TRUE))
res.summary <- summary(res, short = FALSE)


mvord:::check(all.equal(res.summary$thresholds$Estimate,
                        c(-0.70654808839454852354, 0.84166389562778798350 ,-0.80564017944211274713,
                          0.87834342403465937021, -0.52550677773899212575,  0.43872841951666374793), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate,
                        c(0.908534877536098850470 , -0.068886066520331828977, -0.150586382347683067628), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`,
                        c(0.14413615638981053246, 0.15094934798830694778, 0.14777311388303066009,
                          0.15259796098859162994, 0.13805473124275577379, 0.13625729487947765839), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.039421737990101658744,0.133978933777893494117,0.129319546675785962409), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]],  -565.1857975044365502981, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 1163.199625517649110407, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 1205.960928690734590418, tolerance = tolerance))
##########################################
### fixed coefs fixed thresholds
###########################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df_NA,
                     link = mvprobit(),
                     error.structure = cor_general(~1),
                     threshold.values = list(c(-1,1), c(-1,1)),
                     coef.values  = matrix(c(0.8,-0.5, 0.8,-0.5),ncol=2, byrow=TRUE),
                     control= mvord.control(solver="BFGS",se=TRUE))

res.summary <- summary(res, short = FALSE)


mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.854364908267473466275), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c( 0.064212481459284348473), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -124.6298356200302066554, tolerance = tolerance))
mvord:::check(all.equal(AIC(res),  251.280079403325714793, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 253.9016800682576047166, tolerance = tolerance))


##########################################
### fixed cor_general
###########################################
#polychor
res <- mvord::mvord(formula = MMO(Y) ~ 1,
                     data = df,
                     link = mvprobit(),
                     error.structure = cor_general(~1, value = rep(0.5, 100), fixed = TRUE),
                     threshold.constraints = c(1,1),
                     control= mvord.control(solver="BFGS",se=TRUE))

res.summary <- summary(res, short = FALSE)
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(0.0000000000000000000, 1.75287953185495593011, 0.0000000000000000000, 1.75287953185495593011), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.80138487301048844103, 0.83916182795533666994), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.00000000000000000000, 0.12271683314006323617, 0.00000000000000000000, 0.12271683314006323617), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.13627526573459616821, 0.14543826409239238306), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -167.9806162341531887705 , tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 342.1467994786156623377, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 350.2040268579602866339, tolerance = tolerance))

##########################################
### fixed cor_equi
###########################################
res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     #index = c("i", "j"),
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_equi(~0+X1, fixed = TRUE),
                     threshold.constraints = c(1,1),
                     coef.constraints = cbind(c(1,1),c(1,2)))

res.summary <- summary(res, short = FALSE)

options(digits = 22)

# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.971882533677633775326,  1.049862299350415195676,  -0.971882533677633775326,   1.049862299350415195676), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(0.641833309819990649459, -0.450733949174842496443, -0.406790878681064282940), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.106998427476494920030 ,0.099201905823191174894, 0.106998427476494920030, 0.099201905823191174894), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.085721117636038332566,0.188620706054737241608, 0.157255593575882451907), tolerance = tolerance))
mvord:::check(all.equal(logLik(res)[[1]], -161.6316436136202128182, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 333.7896030167141248057, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 347.5010250482303604258, tolerance = tolerance))
##########################################
### try all fixed
###########################################
res <- mvord::mvord(formula = MMO(Y) ~ 0,
                     data = df,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_equi(~1, value = 0.3, fixed = TRUE),
                     threshold.values = list(c(-1,1), c(-1,1)))


mvord:::check(all.equal(logLik(res)[[1]], -179.9334335736217553858, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 359.8668671472435107717, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 359.8668671472435107717, tolerance = tolerance))

res <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2,
                     data = df,
                     link = mvprobit(),
                     control = mvord.control(solver = "BFGS"),
                     error.structure = cor_equi(~1, value = 0.8, fixed = TRUE),
                     threshold.values = list(c(-1,1), c(-1,1)),
                     coef.values = rbind(c(0.8,-0.5), c(0.8,-0.5)))

mvord:::check(all.equal(logLik(res)[[1]],  -136.4789503795024359079, tolerance = tolerance))
mvord:::check(all.equal(AIC(res), 272.9579007590048718157, tolerance = tolerance))
mvord:::check(all.equal(BIC(res), 272.9579007590048718157, tolerance = tolerance))

## make three responses
suppressWarnings(RNGversion("3.5.0"))
set.seed(100)
Y3 <- factor(sample(1:3, nobs, replace = TRUE), ordered = TRUE)
dat_3 <- cbind(Y3, data_toy_example)
Rfixed <- matrix(c(1,0.9,0.2, 0.9, 1, 0.5, 0.2,0.5,1), nrow = 3)
res <- mvord::mvord(formula = MMO2(Y1, Y2, Y3) ~ 0,
                     data = dat_3,
                     link = mvprobit(),
                     error.structure = cor_general(~1, value = Rfixed[lower.tri(Rfixed)], fixed = TRUE),
                     threshold.values = list(c(-1,1), c(-1,1), c(-1,1)),
                     control= mvord.control(solver="BFGS",se=TRUE))

mvord:::check(all.equal(logLik(res)[[1]], -647.8423017485074524302, tolerance = tolerance))
mvord:::check(all.equal(AIC(res),  1295.68460349701490486, tolerance = tolerance))
mvord:::check(all.equal(BIC(res),  1295.68460349701490486, tolerance = tolerance))

###############################
# tests for predict functions #
###############################
data_toy_example$X3 <- cut(data_toy_example$X2, c(-Inf, -0.2, 0.2, Inf))

dat_train <- data_toy_example[1:90, ]
dat_test  <- data_toy_example[-c(1:90), ]

## mvord2
### example with factor covariates
res_train <- mvord::mvord(formula = MMO2(Y1, Y2) ~ 0 + X1 + X2 + X3,
                     data = dat_train,
                     link = mvprobit(),
                     contrasts = list(X3 = "contr.sum"))
tolerance <- 0.00001
#### predict in sample
phat_in_mvord2 <- predict(res_train)
mvord:::check(all.equal(
  unname(phat_in_mvord2[1:10]), c(0.1464235557313708635530, 0.4321397714843485116099, 0.5968270247741616074677,
                        0.0166994091681152423412, 0.5854390261665686212567, 0.1371096523508007203329,
                        0.1675153858158427988556, 0.5777761575951960715258, 0.5680673832081621910106,
                        0.2403876149104018367098),  tolerance = tolerance))

pjhat_in_mvord2 <- joint_probabilities(res_train, response.cat = c(1,2))
pjhat_in_mvord2_1 <- joint_probabilities(res_train, response.cat = c(1,1))

p1 <- joint_probabilities(res_train, response.cat = c(1,NA))
p2 <- joint_probabilities(res_train, response.cat = c(2,NA))
min(rowSums(cbind(pjhat_in_mvord2/p1, pjhat_in_mvord2_1/p1)))

mvord:::check(all.equal(
  unname(pjhat_in_mvord2[1:10]), c(0.0414494114970503091388565, 0.0007763727128762831775077, 0.0258791356876044154056160,
                                   0.0166994147132800141442033, 0.0468724891285079159342075, 0.0021647187725574168482012,
                                   0.0172638911663273442176347, 0.0094837831334130262561644, 0.0190099916431371585012755,
                                   0.0051418125225030086866695),  tolerance = tolerance))

pmhat_in_mvord2 <- marginal_predict(res_train)
mvord:::check(all.equal(
  unname(c(pmhat_in_mvord2[1:5,1:2])), c(0.16369555938725077748330, 0.51206442964897258551815, 0.70085752018010938346748,
                                         0.07780827569973108870371, 0.68584056337028265204481, 0.23658538526249328626250,
                                         0.62104312126876426436439, 0.65391318126703112945108, 0.62425836768740627924501,
                                         0.64830276740170678095865),  tolerance = tolerance))



#### predict with newdata
phat_new_mvord2 <- predict(res_train, newdata = dat_test)
mvord:::check(all.equal(
  unname(phat_new_mvord2), c(0.50788676686680833682885, 0.53701228897073405299523, 0.41931860775657509021741,
                  0.56837451159716412263379, 0.46158948077333511461617, 0.02998639798438534898040,
                  0.06314648572131303927435, 0.21981855723244725364651, 0.57800150278644024659513,
                  0.57873631730860186639376),  tolerance = tolerance))

pjhat_new_mvord2 <- joint_probabilities(res_train, response.cat = c(1,2), newdata = dat_test)
mvord:::check(all.equal(
  unname(pjhat_new_mvord2), c(0.098893552409077684073324, 0.007494474951204033175145, 0.223630828339767884216371,
                              0.012572935869544998865877, 0.125624085603978530301106, 0.021774380504568757732642,
                              0.063146525118494334360975, 0.000117581041524150720079, 0.017153034070206046868279,
                              0.049461049315582095164956),  tolerance = tolerance))

pmhat_new_mvord2 <- marginal_predict(res_train, newdata = dat_test)
mvord:::check(all.equal(
  unname(c(pmhat_new_mvord2[1:5, 1:2])), c(0.5879858728748059704117, 0.6245339414959365509361, 0.6429603932672183219665,
                                           0.6645168757663386660539, 0.5872139451846044577721, 0.6110788197600012239263,
                                           0.6185525371757517598681, 0.4254630210115252220149, 0.6355172734548747426331,
                                           0.4875151835646630571475),  tolerance = tolerance))


### example with offsets
res_train2 <- mvord::mvord(formula = MMO2(Y1, Y2) ~ 0 + X1 + offset(X2) + X3,
                            data = dat_train,
                            link = mvprobit(),
                            contrasts = list(X3 = "contr.sum"))
#### predict in sample
phat_in_offset_mvord2 <- predict(res_train2)
mvord:::check(all.equal(
  unname(phat_in_offset_mvord2[1:10]), c(0.04651187641465470007374, 0.37252604618539558734014, 0.52685355820976886853657,
                                         0.01765062793667723783919, 0.56828750653284465510495, 0.13097924136563171559899,
                                         0.01252594308629018104995, 0.55895559865178667813268, 0.57808854222406069744977,
                                         0.47885651300164144839044),  tolerance = tolerance))

pjhat_in_mvord2 <- joint_probabilities(res_train, response.cat = c(1,2))
mvord:::check(all.equal(
  unname(pjhat_in_mvord2[1:10]), c(0.0414494114970503091388565, 0.0007763727128762831775077, 0.0258791356876044154056160,
                                  0.0166994147132800141442033, 0.0468724891285079159342075, 0.0021647187725574168482012,
                                  0.0172638911663273442176347, 0.0094837831334130262561644, 0.0190099916431371585012755,
                                  0.0051418125225030086866695),  tolerance = tolerance))

pmhat_in_mvord2 <- marginal_predict(res_train)
mvord:::check(all.equal(
  unname(c(pmhat_in_mvord2[1:5,1:2])), c(0.16369555938725077748330, 0.51206442964897258551815, 0.70085752018010938346748,
                                        0.07780827569973108870371, 0.68584056337028265204481, 0.23658538526249328626250,
                                        0.62104312126876426436439, 0.65391318126703112945108, 0.62425836768740627924501,
                                        0.64830276740170678095865),  tolerance = tolerance))



#### predict with newdata
phat_new_offset_mvord2 <- predict(res_train2, newdata = dat_test)
mvord:::check(all.equal(
  unname(phat_new_offset_mvord2), c(0.53788568792038271570988, 0.47226651862936941395077, 0.59736906923858912321634,
                                    0.51730764069679358030385, 0.18750767549430336078586, 0.01993645212535863353587,
                                    0.05609361865051418205574, 0.15695387935814331115125, 0.57079407463923570453801,
                                    0.46447453401393973271283),  tolerance = tolerance))


pjhat_new_offset_mvord2 <- joint_probabilities(res_train2, response.cat = c(1,2),
                                               newdata = dat_test)
mvord:::check(all.equal(
  unname(pjhat_new_offset_mvord2), c(0.0839907843603347747940546, 0.0053737551172480801930931, 0.1657125298469289131908511,
                                      0.0095039724989203211436006, 0.1110141336202859763115924, 0.0308028085793187883512090,
                                      0.0560936429791563384572584, 0.0000659201973896017534571, 0.0161492922521473680763648,
                                      0.0556615607647627450016437),  tolerance = tolerance))

pmhat_new_offset_mvord2 <- marginal_predict(res_train2, newdata = dat_test)
mvord:::check(all.equal(
  unname(c(pmhat_new_offset_mvord2[1:5, 1:2])), c(0.6410297683731078777214, 0.5530240404054447278526, 0.7630836817892157064591,
                                                  0.6103504462453783752096, 0.2985228989322880055468, 0.6269027239883980806567,
                                                  0.5461645309760024824541, 0.6047542519330054711091, 0.5786243958547179211394,
                                                  0.2011002687374978670221),  tolerance = tolerance))

## mvord
dat_train <- df[df$i %in% 1:90, ]
dat_test  <- df[!(df$i %in% 1:90), ]

res_train <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + X2 + X3,
                           data = dat_train,
                           link = mvprobit(),
                           contrasts = list(X3 = "contr.sum"))
### example with factor covariates
#### predictions in sample
phat_in_mvord <- predict(res_train)
mvord:::check(all.equal(
  unname(phat_in_mvord[1:10]), c(0.1464235557313708635530, 0.4321397714843485116099, 0.5968270247741616074677,
                                 0.0166994091681152423412, 0.5854390261665686212567, 0.1371096523508007203329,
                                 0.1675153858158427988556, 0.5777761575951960715258, 0.5680673832081621910106,
                                 0.2403876149104018367098),  tolerance = tolerance))

pjhat_in_mvord <- joint_probabilities(res_train, response.cat = c(1,2))
mvord:::check(all.equal(
  unname(pjhat_in_mvord[1:10]), c(0.0414494114970503091388565, 0.0007763727128762831775077, 0.0258791356876044154056160,
                                  0.0166994147132800141442033, 0.0468724891285079159342075, 0.0021647187725574168482012,
                                  0.0172638911663273442176347, 0.0094837831334130262561644, 0.0190099916431371585012755,
                                  0.0051418125225030086866695),  tolerance = tolerance))

pmhat_in_mvord <- marginal_predict(res_train)
mvord:::check(all.equal(
  unname(c(pmhat_in_mvord[1:5,1:2])), c(0.16369555938725077748330, 0.51206442964897258551815, 0.70085752018010938346748,
                                     0.07780827569973108870371, 0.68584056337028265204481, 0.23658538526249328626250,
                                     0.62104312126876426436439, 0.65391318126703112945108, 0.62425836768740627924501,
                                     0.64830276740170678095865),  tolerance = tolerance))


#### predict with newdata
phat_new_mvord <- predict(res_train, newdata = dat_test)
mvord:::check(all.equal(
  unname(phat_new_mvord), c(0.50788676686680833682885, 0.53701228897073405299523, 0.41931860775657509021741,
                            0.56837451159716412263379, 0.46158948077333511461617, 0.02998639798438534898040,
                            0.06314648572131303927435, 0.21981855723244725364651, 0.57800150278644024659513,
                            0.57873631730860186639376),  tolerance = tolerance))

pjhat_new_mvord <- joint_probabilities(res_train, response.cat = c(1,2), newdata = dat_test)
mvord:::check(all.equal(
  unname(pjhat_new_mvord), c(0.098893552409077684073324, 0.007494474951204033175145, 0.223630828339767884216371,
                             0.012572935869544998865877, 0.125624085603978530301106, 0.021774380504568757732642,
                             0.063146525118494334360975, 0.000117581041524150720079, 0.017153034070206046868279,
                             0.049461049315582095164956),  tolerance = tolerance))

pmhat_new_mvord <- marginal_predict(res_train, newdata = dat_test)
mvord:::check(all.equal(
  unname(c(pmhat_new_mvord[1:5,1:2])),
  c(0.5879858728748059704117, 0.6245339414959365509361, 0.6429603932672183219665,
    0.6645168757663386660539, 0.5872139451846044577721, 0.6110788197600012239263,
    0.6185525371757517598681, 0.4254630210115252220149, 0.6355172734548747426331,
    0.4875151835646630571475),  tolerance = tolerance))






### example with fixed error structure
res_train2_vec <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                           data = dat_train,
                           link = mvprobit(),
                           cor_ar1(~1, value = rep(0,length(unique(dat_train$i))), fixed = TRUE))

mvord:::check(all.equal(logLik(res_train2_vec)[[1]], -148.0213745871528772113, tolerance = tolerance))


res_train2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                           data = dat_train,
                           link = mvprobit(),
                           cor_ar1(~1, value = 0, fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2_vec)[[1]], -148.0213745871528772113, tolerance = tolerance))

phat_mvord <- predict(res_train2)

phat_mp_mvord <- marginal_predict(res_train2)

pjhat_mvord <- joint_probabilities(res_train2, response.cat = c(1,2))

# res_train2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
#                            data = dat_train,
#                            link = mvprobit(),
#                            cor_ar1(~1))
#
# res_train2_free$rho$optpar
# res_train2_free$rho$varGamma

res_train2_vec <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                               data = dat_train,
                               link = mvprobit(),
                               cor_equi(~1, value = rep(0,length(unique(dat_train$i))), fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2_vec)[[1]], -148.0213745871528772113, tolerance = tolerance))

res_train2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                           data = dat_train,
                           link = mvprobit(),
                           cor_equi(~1, value = 0, fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2)[[1]], -148.0213745871528772113, tolerance = tolerance))

phat_mvord <- predict(res_train2)

phat_mp_mvord <- marginal_predict(res_train2)

pjhat_mvord <- joint_probabilities(res_train2, response.cat = c(1,2))


# cor_general
res_train2_vec <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                               data = dat_train,
                               link = mvprobit(),
                               cor_general(~1, value = matrix(rep(0,length(unique(dat_train$i))), ncol = 1), fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2_vec)[[1]], -148.0213745871528772113, tolerance = tolerance))

res_train2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                           data = dat_train,
                           link = mvprobit(),
                           cor_general(~1, value = 0, fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2)[[1]], -148.0213745871528772113, tolerance = tolerance))

phat_mvord <- predict(res_train2)

phat_mp_mvord <- marginal_predict(res_train2)

pjhat_mvord <- joint_probabilities(res_train2, response.cat = c(1,2))

#cov_general
res_train2_vec <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                               data = dat_train,
                               link = mvprobit(),
                               cov_general(~1, value = matrix(rep(c(1,0,1),length(unique(dat_train$i))), ncol = 3, byrow= TRUE), fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2_vec)[[1]], -187.3253828363251329847, tolerance = tolerance))

res_train2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 +X2,
                           data = dat_train,
                           link = mvprobit(),
                           error.structure = cov_general(~1, value = c(1,0,1), fixed = TRUE))
mvord:::check(all.equal(logLik(res_train2)[[1]], -187.3253828363251329847, tolerance = tolerance))

error_structure(res_train2, type ="sigmas")


phat_mvord <- predict(res_train2)
#warnings()

phat_mp_mvord <- marginal_predict(res_train2)

pjhat_mvord <- joint_probabilities(res_train2, response.cat = c(1,2))

### example with offsets
res_train2 <- mvord::mvord(formula = MMO(Y) ~ 0 + X1 + offset(X2) + X3,
                           data = dat_train,
                           link = mvprobit(),
                           contrasts = list(X3 = "contr.sum"))
#### predict in sample
phat_in_offset_mvord <- predict(res_train2)
mvord:::check(all.equal(
  unname(phat_in_offset_mvord[1:10]), c(0.04651187641465470007374, 0.37252604618539558734014, 0.52685355820976886853657,
                                        0.01765062793667723783919, 0.56828750653284465510495, 0.13097924136563171559899,
                                        0.01252594308629018104995, 0.55895559865178667813268, 0.57808854222406069744977,
                                        0.47885651300164144839044),  tolerance = tolerance))

pjhat_in_offset_mvord <- joint_probabilities(res_train2, response.cat = c(1,2))
mvord:::check(all.equal(
  unname(pjhat_in_offset_mvord[1:10]), c(0.0472419001589241827065990, 0.0005045621079537432329687, 0.0160654878216455188066902,
                                         0.0176506444778379567583926, 0.0349580018370856671072744, 0.0017616969113226077503498,
                                         0.0043762331129306719645911, 0.0079872301725372740754949, 0.0211556317124290127473785,
                                         0.0023057907016680312395351),  tolerance = tolerance))

pmhat_in_offset_mvord <- marginal_predict(res_train2)
mvord:::check(all.equal(
  unname(c(pmhat_in_offset_mvord[1:5, 1:2])), c(0.05687801748367393717132, 0.43314394646591181103901, 0.63410255320208053220199,
                                                0.11659868247577427624595, 0.68287284549406501721336, 0.07832694375364634975512,
                                                0.55065457504773962504885, 0.57624707519535223188001, 0.62708167599782682621878,
                                                0.62089949532970512002805),  tolerance = tolerance))


## predict with newdata
phat_new_offset_mvord <- predict(res_train2, newdata = dat_test)
mvord:::check(all.equal(
  unname(phat_new_offset_mvord), c(0.53788568792038271570988, 0.47226651862936941395077, 0.59736906923858912321634,
                                   0.51730764069679358030385, 0.18750767549430336078586, 0.01993645212535863353587,
                                   0.05609361865051418205574, 0.15695387935814331115125, 0.57079407463923570453801,
                                   0.46447453401393973271283),  tolerance = tolerance))

pjhat_new_offset_mvord <- joint_probabilities(res_train2, response.cat = c(1,2), newdata = dat_test)
mvord:::check(all.equal(
  unname(pjhat_new_offset_mvord), c(0.0839907843603347747940546, 0.0053737551172480801930931, 0.1657125298469289131908511,
                                    0.0095039724989203211436006, 0.1110141336202859763115924, 0.0308028085793187883512090,
                                    0.0560936429791563384572584, 0.0000659201973896017534571, 0.0161492922521473680763648,
                                    0.0556615607647627450016437),  tolerance = tolerance))

pmhat_new_offset_mvord <- marginal_predict(res_train2, newdata = dat_test)
mvord:::check(all.equal(
  unname(c(pmhat_new_offset_mvord[1:5,1:2])), c(
    0.6410297683731078777214, 0.5530240404054447278526, 0.7630836817892157064591,
    0.6103504462453783752096, 0.2985228989322880055468, 0.6269027239883980806567,
    0.5461645309760024824541, 0.6047542519330054711091, 0.5786243958547179211394,
    0.2011002687374978670221),  tolerance = tolerance))


