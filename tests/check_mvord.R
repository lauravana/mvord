library(mvord)

#load data
data("data_mvord")
head(data_mvord)

#-------------
# cor_general
#-------------
## Not run:
# approx 1 min
res_cor <- mvord:::mvord(formula = MMO(rating, firm_id, rater_id) ~ 0 + X1 + X2 + X3 + X4 + X5,
                   #formula ~ 0 ... without intercept
                   #index = c("firm_id", "rater_id"),
                   #not necessary if firmID is first column and rater is second column in data
                   data = data_mvord, #choose data
                   response.levels = list(c("G","F","E", "D", "C", "B", "A"),
                                          c("G","F","E", "D", "C", "B", "A"),
                                          c("O","N","M","L", "K", "J", "I", "H")),
                   #list for each rater;
                   #need to be set if specific levels/labels are desired (not in natural ordering)
                   #response.names = c("rater1", "rater2", "rater3"),
                   # set if not all raters are used and specifies ordering
                   link = mvprobit(), #probit or alogit
                   error.structure = cor_general(~1), #different error structures
                   coef.constraints = cbind(c(1,2,2),
                                            c(1,1,2),
                                            c(NA,1,2),
                                            c(NA,NA,NA),
                                            c(1,1,2)),#either a vector or a matrix
                   coef.values = cbind(c(NA,NA,NA),
                                       c(NA,NA,NA),
                                       c(0,NA,NA),
                                       c(1,1,1),
                                       c(NA,NA,NA)),
                   #matrix (possible if coef.constraints is a matrix)
                   threshold.constraints = c(1,1,2),
                   #solver = "newuoa") #BFGS is faster
)
print(res_cor, call = FALSE)
res_cor.summary <- summary(res_cor, short = FALSE)
# res_cor.summary
# thresholds(res_cor)
# coefficients(res_cor)
# get.error.struct(res_cor)

#check par
options(digits = 22)
# paste(format(res_cor.summary$Estimate), collapse = ",")
mvord:::check(all.equal(c(res_cor.summary$thresholds$Estimate,
  res_cor.summary$coefficients$Estimate, res_cor.summary$error.structure$Estimate),
                c(-4.0875998663039574410050, -2.4604528485010810356926, -1.0242005011600983088726,  1.0341012929138195808321,
                  2.5374262532354427968073,  4.2860166974636388914632, -4.0875998663039574410050, -2.4604528485010814797818,
                  -1.0242005011600980868280,  1.0341012929138198028767, 2.5374262532354423527181,  4.2860166974636388914632,
                  -4.7740156535556401706799, -3.4516236175055778900855, -1.8314457828246328841004, -0.3928419217472290947413,
                  1.2299327843176965924243,  2.6060983851102483832562,  4.5595822763854227943625,
                  0.8664911462254030194430, 0.6363387387318012455850, -0.5439454055658743403256, -1.0715577249351615485295,
                  0.1840960060925358188211, -0.2125229066909701258759, -0.6290716517330824375520, -0.8363720960361500367242,
                  0.8708985232957860977976, 0.8316844611287248500986, 0.6992995912049327911220), tolerance = 0.001))


#check se
#paste(format(res_cor.summary$`Std. Error`), collapse = ",")
mvord:::check(all.equal(c(res_cor.summary$thresholds$`Std. Error`, res_cor.summary$coefficients$`Std. Error`, res_cor.summary$error.structure$`Std. Error`),
                c(0.13888976956270382423497, 0.07005019934680442594832, 0.05088894731787228786768, 0.04835893519039383137148,
                  0.07369578409058712142876, 0.14939028227977182994302, 0.13888976956270382423497, 0.07005019934680442594832,
                  0.05088894731787228786768, 0.04835893519039383137148, 0.07369578409058712142876, 0.14939028227977182994302,
                  0.19462260360636776668208, 0.09993848248356777819179, 0.06526136869387920425023, 0.05351395618005048981924,
                  0.05968587266688644282775, 0.08268507331077515343232, 0.17303260905267009439612,
                  0.03758737492257332307721, 0.03152574715668300664451, 0.03360806296478834054309, 0.03811269880849090657682,
                  0.02565728423846840716704, 0.02682951404421570623660, 0.03371348291873306479705, 0.03804459680325304737902,
                  0.01387001816700043389796, 0.01712987494537120125582, 0.02259124348297502152261), tolerance = 0.001))

#check loglik

mvord:::check(all.equal(logLik(res_cor)[[1]], -5222.586652085715286375, tolerance = 0.001))
# #check AIC

mvord:::check(all.equal(AIC(res_cor), 10542.26761241730491747, tolerance = 0.001))
# #check BIC

mvord:::check(all.equal(BIC(res_cor), 10780.52516434370772913, tolerance = 0.001))

#-------------
# cov_general
#-------------
#approx 4 min
res_cov <- mvord(formula = MMO(rating) ~ 1 + X1 + X2 + X3 + X4 + X5,
                   #formula ~ 0 ... without intercept
                   #index = c("firm_id", "rater_id"),
                   #not necessary if firmID is first column and rater is second column in data
                   data = data_mvord, #choose data
                   # response.levels = list(c("G","F","E", "D", "C", "B", "A"),
                   #                        c("G","F","E", "D", "C", "B", "A"),
                   #                        c("O","N","M","L", "K", "J", "I", "H")),
                   #list for each rater;
                   #need to be set if specific levels/labels are desired
                   #response.names = c("rater1", "rater2", "rater3"),
                   # set if not all raters are used and specifies ordering
                   link = mvprobit(), #probit or alogit
                   error.structure = cov_general(~1), #different error structures
                   threshold.constraints = NULL, #vector
                   threshold.values = list(c(-4,NA,NA,NA,NA,4.5),
                                           c(-4,NA,NA,NA,NA,4),
                                           c(-5,NA,NA,NA,NA,NA,4.5))
                   #list for each rater
                 #  solver = "newuoa"
                 ) #does not converge with BFGS
print(res_cov)
res_cov.summary <- summary(res_cov, short = FALSE)
# res_cov.summary
# thresholds(res_cov)
# coefficients(res_cov)
# get.error.struct(res_cov)

#check par
options(digits = 22)
# paste(format(res_cov.summary$Estimate), collapse = ",")
mvord:::check(all.equal(c(res_cov.summary$thresholds$Estimate, res_cov.summary$coefficients$Estimate, res_cov.summary$error.structure$Estimate),
                        c(-3.99999999999999955591079, -1.85982365447351649656582, -0.43470909142659736046355,  1.58641774049757722231391,
                          2.96037994068721710405612,  4.50000000000000000000000, -4.00000000000000000000000, -2.61392958294945998432013,
                          -1.10140861829631964141640, 0.88229858712515785157393, 2.30290944979030598460668, 4.00000000000000000000000,
                          -5.00000000000000000000000,-2.98889276400575454317732, -1.59289755742616012668122,  0.05477441508605492004325,
                          1.51490119327636074686438,  3.15733614802702655666167, 4.50000000000000000000000,
                          0.57475565987818355573324, -0.10035177605432644976080, -0.34544103523249458653765, -0.84865793302510805773409,
                          -0.62532200123506509470417, -0.64036762687110748704100, 0.54908289697688772434958,  0.51281144422778190961054,
                          1.09088937657879370135561,  0.06274731838029297403825, -0.14435692152398768572930,  0.24903885203080558530253,
                          -0.97825733605997255981634, -0.95005949054885707738549, -1.00191424874327172922506,  0.60733724097710339862033,
                          0.61457931743825822135108,  0.84840756313303988811469,
                          0.8771199640568675404140, 0.8333523304563162925618, 0.7053205703262825920774, 0.9636634829420597236904,
                          0.9864635447453774519389, 1.0272486989106017762907), tolerance = 0.001))

#check se
#paste(format(res_cov.summary$`Std. Error`), collapse = ",")
mvord:::check(all.equal(c(res_cov.summary$thresholds$`Std. Error`, res_cov.summary$coefficients$`Std. Error`, res_cov.summary$error.structure$`Std. Error`),
                        c(0.0000000000000000000000, 0.1707198082928001736658, 0.1415648215532093745495, 0.1196828594539791484896,
                          0.1197226677694545049491, 0.0000000000000000000000, 0.0000000000000000000000, 0.1419637699778572847986,
                          0.1261649786447358867481, 0.1302280055940158998151, 0.1482844005184661972940, 0.0000000000000000000000,
                          0.0000000000000000000000, 0.1462049787534562272917, 0.1345254787833224108251, 0.1344466539934549287327,
                          0.1458571506341596601963, 0.1663486287296551957571, 0.0000000000000000000000,
                          0.11938247752028616210929, 0.11664096558459673136365, 0.12433968158786336266619, 0.04129312264023522055512,
                          0.03855178362525488872103, 0.03707871312395525503769, 0.03514809893473093416194, 0.03563359073101304147491,
                          0.04064811566529103609158, 0.03700380043276757513482, 0.03560741648478967652514, 0.03751080307137479030732,
                          0.04177285509039119215657, 0.04263046247676769567869, 0.04240261466557804920230, 0.03628413680128383250745,
                          0.03718012557457621436452, 0.03967242720303063713283,
                          0.01368858253874749575374, 0.01749697857211561166646, 0.02272151989271832886463, 0.03524061684745433825627,
                          0.03615424177173217484826, 0.03572289119325788159243), tolerance = 0.001))

#check loglik
mvord:::check(all.equal(logLik(res_cov)[[1]], -5202.614799284418950265, tolerance = 0.001))
#check AIC
mvord:::check(all.equal(AIC(res_cov), 70214653818038641, tolerance = 0.001))
#check BIC
mvord:::check(all.equal(BIC(res_cov), 10904.2201015204518626, tolerance = 0.001))


#-------------
# cor_ar1
#-------------
#approx 4min
data(data_mvord_panel)
head(data_mvord_panel)
mult.obs = 5
res_AR1 <- mvord(formula = MMO(rating) ~ 0 + X1 + X2 + X3 + X4 + X5,
                   #formula ~ 0 ... without intercept
                  # index = c("firm_id", "year"),
                   #not necessary if firmID is first column and rater is second column in data
                   data = data_mvord_panel[data_mvord_panel$year %in% c("year3", "year4", "year5", "year6", "year7"),], #choose data
                   #response.levels = rep(list(c("G","F","E", "D", "C", "B", "A")), mult.obs),
                   #list for each rater;
                   #need to be set if specific levels/labels are desired (not in natural ordering)
                   #response.names = c("year3", "year4", "year5", "year6", "year7"),
                   # set if not all raters are used and specifies ordering
                   link = mvprobit(), #probit or aalogit
                   error.structure = cor_ar1(~1), #different error structures
                   threshold.constraints = c(1,1,1,2,2),
                   coef.constraints = c(1,1,1,2,2)
                 #  solver = "BFGS"
                 )
print(res_AR1)
res_AR1.summary <- summary(res_AR1, short = FALSE)
# res_AR1.summary
# thresholds(res_AR1)
# coefficients(res_AR1)
# get.error.struct(res_AR1)
# get.error.struct(res_AR1, type = "corr")

#coef.constraints
all.equal(as.vector(unlist(res_AR1$rho$constraints)), as.vector(rep(cbind(c(rep(1,18), rep(0,12)), c(rep(0,18), rep(1,12))),5)))


#check par
options(digits = 22)
# paste(format(res_AR1.summary$Estimate), collapse = ",")
mvord:::check(all.equal(c(res_AR1.summary$thresholds$Estimate, res_AR1.summary$coefficients$Estimate, res_AR1.summary$error.structure$Estimate),
                        c(-4.3809404612238251885969, -2.6014487978132363465988, -1.0378979572353754790015,  0.9504962937178964565987,
                          2.4210536148501784481368,  3.8370726903362868398517, -4.3809404612238251885969, -2.6014487978132359025096,
                          -1.0378979572353754790015,  0.9504962937178963455764,  2.4210536148501788922260,  3.8370726903362872839409,
                          -4.3809404612238251885969, -2.6014487978132359025096, -1.0378979572353754790015,  0.9504962937178963455764,
                          2.4210536148501788922260,  3.8370726903362872839409, -4.4597432947972990291419, -2.4256366480242990135707,
                          -1.1413565451747142986960,  0.5221627867500505670861,  1.8351142770044337471091,  4.5244731435131306795938,
                          -4.4597432947972990291419, -2.4256366480242990135707, -1.1413565451747142986960,  0.5221627867500505670861,
                          1.8351142770044337471091,  4.5244731435131306795938,
                          -0.775910372581866236707526, -0.608304608756956910475822,  0.496222943081394807229856,  0.954819462523996520531000,
                          0.002047146347684709270093,  0.185836762421712348158920, -0.986354543211104717315152, -0.964865124883243052700266,
                          0.587702469993099940737125, 0.791062503299648667187682), tolerance = 0.001))

#check se
#paste(format(res_AR1.summary$`Std. Error`), collapse = ",")
mvord:::check(all.equal(c(res_AR1.summary$thresholds$`Std. Error`, res_AR1.summary$coefficients$`Std. Error`, res_AR1.summary$error.structure$`Std. Error`),
                        c(0.13836231110364854979267, 0.07295803082496739311313, 0.04211686864255907714050, 0.04260253632122401601379,
                          0.06549954279731404205300, 0.10989716258897458400767, 0.13836231110364854979267, 0.07295803082496739311313,
                          0.04211686864255907714050, 0.04260253632122401601379, 0.06549954279731404205300, 0.10989716258897458400767,
                          0.13836231110364854979267, 0.07295803082496739311313, 0.04211686864255907714050, 0.04260253632122401601379,
                          0.06549954279731404205300, 0.10989716258897458400767, 0.15236538869407589835703, 0.07582416010103604220305,
                          0.05103025585296900451526, 0.04106759083997663095644, 0.05830682416009386886957, 0.13082366188528546380176,
                          0.15236538869407589835703, 0.07582416010103604220305, 0.05103025585296900451526, 0.04106759083997663095644,
                          0.05830682416009386886957, 0.13082366188528546380176,
                          0.02175120297630910026765, 0.02232476848318381765224, 0.01797340879042605074623, 0.02754107431822683596523,
                          0.01358302836775518267209, 0.01790948518849638793071, 0.02512315989817485994973, 0.02807609101903460468996,
                          0.01911438661776309819174, 0.02599194462760747934005,
                          0.03212181813757332415893), tolerance = 0.001))

#check loglik
mvord:::check(all.equal(logLik(res_AR1)[[1]], -17305.52075838709788513, tolerance = 0.001))
#check AIC
mvord:::check(all.equal(AIC(res_AR1), 34785.50813512730383081, tolerance = 0.001))
#check BIC
mvord:::check(all.equal(BIC(res_AR1), 35213.62786874162702588, tolerance = 0.001))





library(mvord)

#load data
data("data_mvord2")
head(data_mvord2)

#-------------
# cor_general
#-------------
## Not run:
res_cor <- mvord(formula = MMO2(rater1, rater2, rater3) ~ 0 + X1 + X2 + X3 + X4 + X5,
                        #formula ~ 0 ... without intercept
                        data = data_mvord2, #choose data
                        link = mvprobit(), #probit or alogit
                        error.structure = cor_general(~1), #different error structures
                    coef.constraints = cbind(c(1,2,2),
                                             c(1,1,2),
                                             c(NA,1,2),
                                             c(NA,NA,NA),
                                             c(1,1,2)),#either a vector or a matrix
                    coef.values = cbind(c(NA,NA,NA),
                                        c(NA,NA,NA),
                                        c(0,NA,NA),
                                        c(1,1,1),
                                        c(NA,NA,NA)),
                    #matrix (possible if coef.constraints is a matrix)
                    threshold.constraints = c(1,1,2)
                  #  solver = "BFGS"
                  )
print(res_cor)
res_cor.summary <- summary(res_cor, short = FALSE)
# res_cor.summary
# thresholds(res_cor)
# coefficients(res_cor)
# get.error.struct(res_cor)

#check par
options(digits = 22)
# paste(format(res_cor.summary$Estimate), collapse = ",")
mvord:::check(all.equal(c(res_cor.summary$thresholds$Estimate,
                          res_cor.summary$coefficients$Estimate, res_cor.summary$error.structure$Estimate),
                        c(-4.0875998663039574410050, -2.4604528485010810356926, -1.0242005011600983088726,  1.0341012929138195808321,
                          2.5374262532354427968073,  4.2860166974636388914632, -4.0875998663039574410050, -2.4604528485010814797818,
                          -1.0242005011600980868280,  1.0341012929138198028767, 2.5374262532354423527181,  4.2860166974636388914632,
                          -4.7740156535556401706799, -3.4516236175055778900855, -1.8314457828246328841004, -0.3928419217472290947413,
                          1.2299327843176965924243,  2.6060983851102483832562,  4.5595822763854227943625,
                          0.8664911462254030194430, 0.6363387387318012455850, -0.5439454055658743403256, -1.0715577249351615485295,
                          0.1840960060925358188211, -0.2125229066909701258759, -0.6290716517330824375520, -0.8363720960361500367242,
                          0.8708985232957860977976, 0.8316844611287248500986, 0.6992995912049327911220), tolerance = 0.001))


#check se
#paste(format(res_cor.summary$`Std. Error`), collapse = ",")
mvord:::check(all.equal(c(res_cor.summary$thresholds$`Std. Error`, res_cor.summary$coefficients$`Std. Error`, res_cor.summary$error.structure$`Std. Error`),
                        c(0.13888976956270382423497, 0.07005019934680442594832, 0.05088894731787228786768, 0.04835893519039383137148,
                          0.07369578409058712142876, 0.14939028227977182994302, 0.13888976956270382423497, 0.07005019934680442594832,
                          0.05088894731787228786768, 0.04835893519039383137148, 0.07369578409058712142876, 0.14939028227977182994302,
                          0.19462260360636776668208, 0.09993848248356777819179, 0.06526136869387920425023, 0.05351395618005048981924,
                          0.05968587266688644282775, 0.08268507331077515343232, 0.17303260905267009439612,
                          0.03758737492257332307721, 0.03152574715668300664451, 0.03360806296478834054309, 0.03811269880849090657682,
                          0.02565728423846840716704, 0.02682951404421570623660, 0.03371348291873306479705, 0.03804459680325304737902,
                          0.01387001816700043389796, 0.01712987494537120125582, 0.02259124348297502152261), tolerance = 0.001))

#check loglik

mvord:::check(all.equal(logLik(res_cor)[[1]], -5222.586652085715286375, tolerance = 0.001))
# #check AIC

mvord:::check(all.equal(AIC(res_cor), 10542.26761241730491747, tolerance = 0.001))
# #check BIC

mvord:::check(all.equal(BIC(res_cor), 10780.52516434370772913, tolerance = 0.001))

###############################################################################################
res_cor.logit <- mvord(formula = MMO2(rater1, rater2, rater3) ~ 0 + X1 + X2 + X3 + X4 + X5,
                  #formula ~ 0 ... without intercept
                  data = data_mvord2, #choose data
                  link = mvlogit(), #probit or alogit
                  error.structure = cor_general(~1), #different error structures
                  coef.constraints = cbind(c(1,2,2),
                                           c(1,1,2),
                                           c(NA,1,2),
                                           c(NA,NA,NA),
                                           c(1,1,2)),#either a vector or a matrix
                  coef.values = cbind(c(NA,NA,NA),
                                      c(NA,NA,NA),
                                      c(0,NA,NA),
                                      c(1,1,1),
                                      c(NA,NA,NA)),
                  #matrix (possible if coef.constraints is a matrix)
                  threshold.constraints = c(1,1,2)
                  #solver = "newuoa"
                  )
print(res_cor.logit)
res_cor.logit.summary <- summary(res_cor.logit, short = FALSE)
#TODO
