############
### INLA ###
############
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
require(INLA)
require(dplyr)

#Create functions for outputs
inla_emarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mean of variance"))
  results
}

inla_mmarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.mmarginal(inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mode of variance"))
  results
}

#Load Data
RoostAll.final <- read.csv("Roostfinal_HMM.csv")
Stationary.final <- read.csv("Stationaryfinal_HMM.csv")
Mobile.final <- read.csv("Mobilefinal_HMM.csv")
AllPoints.final <- read.csv("AllPoints.csv")

#Number of points for each track
numpoints <- AllPoints.final %>% group_by(ID) %>% summarize(Total = n()/11) %>% arrange(Total)
numpoints
mean(numpoints$Total)

#Set mean and precision priors for slope coefficients
mean.beta <- 0
prec.beta <- 1e-4  


################################################################################### 
### 1st Scale - Roost Points
sink("Results/RoostResults_HMM.csv")
cat("Roost ~ Basal Area + Wind Chill + BA*WC + BA|IND")
cat('\n')
formula.random <- Use ~  -1 + BA.Z + WC_prev + BA.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,BA.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.BA.WCprev.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)),
                              control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.BA.WCprev.roost$summary.fixed)
write.csv(model.BA.WCprev.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.BA.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.BA.WCprev.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.BA.WCprev.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.BA.WCprev.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.BA.WCprev.roost$mlik[2,1])
cat("\n")
cat("\n")


cat("Roost ~ Basal Area + SnowDepth + BA*SD + BA|IND")
cat('\n')
formula.random <- Use ~ -1 + BA.Z + SD + BA.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,BA.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.BA.SD.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)),
                          control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.BA.SD.roost$summary.fixed)
write.csv(model.BA.SD.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.BA.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.BA.SD.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.BA.SD.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.BA.SD.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.BA.SD.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Distance to Edge + Wind Chill + DtE*WC + DtE|IND')
cat('\n')
formula.random <- Use ~  -1 + DtFE.Z + WC_prev + DtFE.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,DtFE.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.WCprev.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)),
                                control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.DtFE.WCprev.roost$summary.fixed)
write.csv(model.DtFE.WCprev.roost$summary.hyperpar)
cat("Posterior Mean of Variance: ")
cat('\n')
cat(inla_emarginal(model.DtFE.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.WCprev.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.DtFE.WCprev.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.DtFE.WCprev.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.DtFE.WCprev.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Distance to Edge + SnowDepth + DtE*SD + DtE|IND')
cat('\n')
formula.random <- Use ~ -1 + DtFE.Z + SD + DtFE.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,DtFE.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.SD.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                            control.fixed = list(
                              mean = mean.beta,
                              prec = list(default = prec.beta)),
                            control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.DtFE.SD.roost$summary.fixed)
write.csv(model.DtFE.SD.roost$summary.hyperpar)
cat('Posterior Meand of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.SD.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.DtFE.SD.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.DtFE.SD.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.DtFE.SD.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Percent Softwood + Wind Chill + SW*WC + SW|IND')
cat('\n')
formula.random <- Use ~  -1 + SW.Z + WC_prev + SW.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,SW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.SW.WCprev.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)),
                              control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.SW.WCprev.roost$summary.fixed)
write.csv(model.SW.WCprev.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.SW.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.SW.WCprev.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.SW.WCprev.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.SW.WCprev.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.SW.WCprev.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Percent Softwood + SnowDepth + SW*SD + SW|IND')
cat('\n')
formula.random <- Use ~ -1 + SW.Z + SD + SW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,SW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.SW.SD.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)),
                          control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.SW.SD.roost$summary.fixed)
write.csv(model.SW.SD.roost$summary.hyperpar)
cat("Posterior Mean of Variance")
cat('\n')
cat(inla_emarginal(model.SW.SD.roost))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.SW.SD.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.SW.SD.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.SW.SD.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.SW.SD.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Mean Height + Wind Chill + Ht*WC + Ht|IND')
cat('\n')
formula.random <- Use ~  -1 + Ht.Z + WC_prev + Ht.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Ht.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.Ht.WCprev.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)),
                              control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.Ht.WCprev.roost$summary.fixed)
write.csv(model.Ht.WCprev.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.WCprev.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.Ht.WCprev.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.Ht.WCprev.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.Ht.WCprev.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Mean Height + SnowDepth + Ht*SD + Ht|IND')
cat('\n')
formula.random <- Use ~ -1 + Ht.Z + SD + Ht.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Ht.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.Ht.SD.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)),
                          control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.Ht.SD.roost$summary.fixed)
write.csv(model.Ht.SD.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.SD.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.Ht.SD.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.Ht.SD.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.Ht.SD.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Wind Exposure + Wind Chill + WindExp*WC + WindExp|IND')
cat('\n')
formula.random <- Use ~  -1 + Wind.Exp.Z + WC_prev + Wind.Exp.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Wind.Exp.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.WCprev.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.WindExp.WCprev.roost$summary.fixed)
write.csv(model.WindExp.WCprev.roost$summary.hyperpar)
cat('Posterio Mean of Variance')
cat('\n')
cat(inla_emarginal(model.WindExp.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.WindExp.WCprev.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.WindExp.WCprev.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.WindExp.WCprev.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.WindExp.WCprev.roost$mlik[2,1])
cat("\n")
cat("\n")


cat('Roost ~ Wind Exposure + SnowDepth + WindExp*SD + WindExp|IND')
cat('\n')
formula.random <- Use ~ -1 + Wind.Exp.Z + SD + Wind.Exp.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Wind.Exp.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.SD.roost <- inla(formula.random, family ="Poisson", data=RoostAll.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)),
                               control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.WindExp.SD.roost$summary.fixed)
write.csv(model.WindExp.SD.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.SD.roost))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.WindExp.SD.roost$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.WindExp.SD.roost$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.WindExp.SD.roost$mlik[2,1])
cat("\n")
cat("\n")

sink()


################################################################################### 
### 2nd Scale - Stationary
sink('Results/StationaryResults_HMM.csv')
cat('Stationary ~ Basal Area + Wind Chill + BA*WC + BA|IND')
cat('\n')
formula.random <- Use ~  -1 + BA.Z + WC + BA.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,BA.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.BA.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.BA.WCprev.Stationary$summary.fixed)
write.csv(model.BA.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.BA.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.BA.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.BA.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.BA.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.BA.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Basal Area + SnowDepth + BA*SD + BA|IND')
cat('\n')
formula.random <- Use ~ -1 + BA.Z + SD + BA.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,BA.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.BA.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)),
                               control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.BA.SD.Stationary$summary.fixed)
write.csv(model.BA.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.BA.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.BA.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.BA.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.BA.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.BA.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Distance to Edge + Wind Chill + DtE*WC + DtE|IND')
cat('\n')
formula.random <- Use ~  -1 + DtFE.Z + WC + DtFE.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,DtFE.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                     control.fixed = list(
                                       mean = mean.beta,
                                       prec = list(default = prec.beta)),
                                     control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.DtFE.WCprev.Stationary$summary.fixed)
write.csv(model.DtFE.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.DtFE.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.DtFE.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.DtFE.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Distance to Edge + SnowDepth + DtE*SD + DtE|IND')
cat('\n')
formula.random <- Use ~ -1 + DtFE.Z + SD + DtFE.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,DtFE.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                 control.fixed = list(
                                   mean = mean.beta,
                                   prec = list(default = prec.beta)),
                                 control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.DtFE.SD.Stationary$summary.fixed)
write.csv(model.DtFE.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.DtFE.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.DtFE.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.DtFE.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


# cat('Stationary ~ Percent Softwood + Wind Chill + SW*WC + SW|IND')
# cat('\n')
# formula.random <- Use ~  -1 + SW.Z + WC + SW.Z*WC + #Fixed effects
#   StepLength.Z + #step length
#   f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
#   f(ID,SW.Z,values=1:26,model="iid", #Random Slope
#     hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
# 
# model.SW.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
#                                    control.fixed = list(
#                                      mean = mean.beta,
#                                      prec = list(default = prec.beta)),
# control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
# write.csv(model.SW.WCprev.Stationary$summary.fixed)
# write.csv(model.SW.WCprev.Stationary$summary.hyperpar)
# cat('Posterior Mean of Variance: ')
# cat('\n')
# cat(inla_emarginal(model.SW.WCprev.Stationary))
# cat('\n')
# cat('Posterior Mode of Variance: ')
# cat('\n')
# cat(inla_mmarginal(model.SW.WCprev.Stationary))
# cat("\n")
# cat('WAIC: ')
# cat('\n')
# cat(model.SW.WCprev.Stationary$waic$waic)
# cat("\n")
# cat('DIC: ')
# cat('\n')
# cat(model.SW.WCprev.Stationary$dic$dic)
# cat("\n")
# cat('Marginal Likelihood: ')
# cat('\n')
# cat(model.SW.WCprev.Stationary$mlik[2,1])
# cat("\n")
# cat("\n")
# 
# 
# cat('Stationary ~ Percent Softwood + SnowDepth + SW*SD + SW|IND')
# cat('\n')
# formula.random <- Use ~ -1 + SW.Z + SD + SW.Z*SD +
#   StepLength.Z + 
#   f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
#   f(ID,SW.Z,values=1:26,model="iid",
#     hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
# 
# model.SW.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
#                                control.fixed = list(
#                                  mean = mean.beta,
#                                  prec = list(default = prec.beta)),
# control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
# write.csv(model.SW.SD.Stationary$summary.fixed)
# write.csv(model.SW.SD.Stationary$summary.hyperpar)
# cat('Posterior Mean of Variance')
# cat('\n')
# cat(inla_emarginal(model.SW.SD.Stationary))
# cat('\n')
# cat('Posterior Mode of Variance')
# cat('\n')
# cat(inla_mmarginal(model.SW.SD.Stationary))
# cat("\n")
# cat('WAIC: ')
# cat('\n')
# cat(model.SW.SD.Stationary$waic$waic)
# cat("\n")
# cat('DIC: ')
# cat('\n')
# cat(model.SW.SD.Stationary$dic$dic)
# cat("\n")
# cat('Marginal Likelihood: ')
# cat('\n')
# cat(model.SW.SD.Stationary$mlik[2,1])
# cat("\n")
# cat("\n")


cat('Stationary ~ Mean Height + Wind Chill + Ht*WC + Ht|IND')
cat('\n')
formula.random <- Use ~  -1 + Ht.Z + WC + Ht.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Ht.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.Ht.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.Ht.WCprev.Stationary$summary.fixed)
write.csv(model.Ht.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.Ht.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.Ht.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.Ht.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Mean Height + SnowDepth + Ht*SD + Ht|IND')
cat('\n')
formula.random <- Use ~ -1 + Ht.Z + SD + Ht.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Ht.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.Ht.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)),
                               control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.Ht.SD.Stationary$summary.fixed)
write.csv(model.Ht.SD.Stationary$summary.hyperpar)
cat('Posteruir Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.Ht.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.Ht.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.Ht.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Wind Exposure + Wind Chill + WindExp*WC + WindExp|IND')
cat('\n')
formula.random <- Use ~  -1 + Wind.Exp.Z + WC + Wind.Exp.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Wind.Exp.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                        control.fixed = list(
                                          mean = mean.beta,
                                          prec = list(default = prec.beta)),
                                        control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.WindExp.WCprev.Stationary$summary.fixed)
write.csv(model.WindExp.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.WindExp.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.WindExp.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.WindExp.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Wind Exposure + SnowDepth + WindExp*SD + WindExp|IND')
cat('\n')
formula.random <- Use ~ -1 + Wind.Exp.Z + SD + Wind.Exp.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Wind.Exp.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)),
                                    control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.WindExp.SD.Stationary$summary.fixed)
write.csv(model.WindExp.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.WindExp.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.WindExp.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.WindExp.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion Ag + Wind Chill + PropAg*WC + PropAg|IND')
cat('\n')
formula.random <- Use ~  -1 + PropAg.Z + WC + PropAg.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropAg.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                       control.fixed = list(
                                         mean = mean.beta,
                                         prec = list(default = prec.beta)),
                                       control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropAg.WCprev.Stationary$summary.fixed)
write.csv(model.PropAg.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropAg.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropAg.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropAg.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion Ag + SnowDepth + PropAg*SD + PropAg|IND')
cat('\n')
formula.random <- Use ~ -1 + PropAg.Z + SD + PropAg.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropAg.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropAg.SD.Stationary$summary.fixed)
write.csv(model.PropAg.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropAg.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropAg.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropAg.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion Dev + Wind Chill + PropDev*WC + PropDev|IND')
cat('\n')
formula.random <- Use ~  -1 + PropDev.Z + WC + PropDev.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropDev.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                        control.fixed = list(
                                          mean = mean.beta,
                                          prec = list(default = prec.beta)),
                                        control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropDev.WCprev.Stationary$summary.fixed)
write.csv(model.PropDev.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropDev.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropDev.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropDev.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropDev.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion Dev + SnowDepth + PropDev*SD + PropDev|IND')
cat('\n')
formula.random <- Use ~ -1 + PropDev.Z + SD + PropDev.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropDev.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)),
                                    control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropDev.SD.Stationary$summary.fixed)
write.csv(model.PropDev.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropDev.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropDev.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropDev.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropDev.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion SW + Wind Chill + PropSW*WC + PropSW|IND')
cat('\n')
formula.random <- Use ~  -1 + PropSW.Z + WC + PropSW.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropSW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                       control.fixed = list(
                                         mean = mean.beta,
                                         prec = list(default = prec.beta)),
                                       control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropSW.WCprev.Stationary$summary.fixed)
write.csv(model.PropSW.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.PropSW.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropSW.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropSW.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropSW.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropSW.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion SW + SnowDepth + PropSW*SD + PropSW|IND')
cat('\n')
formula.random <- Use ~ -1 + PropSW.Z + SD + PropSW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropSW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropSW.SD.Stationary$summary.fixed)
write.csv(model.PropSW.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropSW.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropSW.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropSW.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropSW.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropSW.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion FoodSub + Wind Chill + PropFoodSub*WC + PropFoodSub|IND')
cat('\n')
formula.random <- Use ~  -1 + PropFoodSub.Z + WC + PropFoodSub.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropFoodSub.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.WCprev.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                            control.fixed = list(
                                              mean = mean.beta,
                                              prec = list(default = prec.beta)),
                                            control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropFoodSub.WCprev.Stationary$summary.fixed)
write.csv(model.PropFoodSub.WCprev.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.WCprev.Stationary))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.WCprev.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropFoodSub.WCprev.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropFoodSub.WCprev.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropFoodSub.WCprev.Stationary$mlik[2,1])
cat("\n")
cat("\n")


cat('Stationary ~ Proportion FoodSub + SnowDepth + PropFoodSub*SD + PropFoodSub|IND')
cat('\n')
formula.random <- Use ~ -1 + PropFoodSub.Z + SD + PropFoodSub.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropFoodSub.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.SD.Stationary <- inla(formula.random, family ="Poisson", data=Stationary.final, 
                                        control.fixed = list(
                                          mean = mean.beta,
                                          prec = list(default = prec.beta)),
                                        control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropFoodSub.SD.Stationary$summary.fixed)
write.csv(model.PropFoodSub.SD.Stationary$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.SD.Stationary))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.SD.Stationary))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropFoodSub.SD.Stationary$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropFoodSub.SD.Stationary$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropFoodSub.SD.Stationary$mlik[2,1])
cat("\n")
cat("\n")



sink()


################################################################################### 
### 3rd Scale - Day Use, All Points
sink('Results/MobileResults_HMM.csv')


cat('Mobile ~ Distance to Edge + Wind Chill + DtE*WC + DtE|IND')
cat('\n')
formula.random <- Use ~  -1 + DtFE.Z + WC + DtFE.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,DtFE.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                 control.fixed = list(
                                   mean = mean.beta,
                                   prec = list(default = prec.beta)),
                                 control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.DtFE.WCprev.Mobile$summary.fixed)
write.csv(model.DtFE.WCprev.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.WCprev.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.WCprev.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.DtFE.WCprev.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.DtFE.WCprev.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.DtFE.WCprev.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Distance to Edge + SnowDepth + DtE*SD + DtE|IND')
cat('\n')
formula.random <- Use ~ -1 + DtFE.Z + SD + DtFE.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,DtFE.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                             control.fixed = list(
                               mean = mean.beta,
                               prec = list(default = prec.beta)),
                             control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.DtFE.SD.Mobile$summary.fixed)
write.csv(model.DtFE.SD.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.SD.Mobile))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.DtFE.SD.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.DtFE.SD.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.DtFE.SD.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.DtFE.SD.Mobile$mlik[2,1])
cat("\n")
cat("\n")


# cat('Mobile ~ Percent Softwood + Wind Chill + SW*WC + SW|IND')
# cat('\n')
# formula.random <- Use ~  -1 + SW.Z + WC + SW.Z*WC + #Fixed effects
#   StepLength.Z + #step length
#   f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
#   f(ID,SW.Z,values=1:26,model="iid", #Random Slope
#     hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
# 
# model.SW.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
#                                control.fixed = list(
#                                  mean = mean.beta,
#                                  prec = list(default = prec.beta)),
# control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
# write.csv(model.SW.WCprev.Mobile$summary.fixed)
# write.csv(model.SW.WCprev.Mobile$summary.hyperpar)
# cat('Posterior Mean of Variance: ')
# cat('\n')
# cat(inla_emarginal(model.SW.WCprev.Mobile))
# cat('\n')
# cat('Posterior Mode of Variance')
# cat('\n')
# cat(inla_mmarginal(model.SW.WCprev.Mobile))
# cat("\n")
# cat('WAIC: ')
# cat('\n')
# cat(model.SW.WCprev.Mobile$waic$waic)
# cat("\n")
# cat('DIC: ')
# cat('\n')
# cat(model.SW.WCprev.Mobile$dic$dic)
# cat("\n")
# cat('Marginal Likelihood: ')
# cat('\n')
# cat(model.SW.WCprev.Mobile$mlik[2,1])
# cat("\n")
# cat("\n")
# 
# 
# cat('Mobile ~ Percent Softwood + SnowDepth + SW*SD + SW|IND')
# cat('\n')
# formula.random <- Use ~ -1 + SW.Z + SD + SW.Z*SD +
#   StepLength.Z + 
#   f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
#   f(ID,SW.Z,values=1:26,model="iid",
#     hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
# 
# model.SW.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
#                            control.fixed = list(
#                              mean = mean.beta,
#                              prec = list(default = prec.beta)),
# control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
# write.csv(model.SW.SD.Mobile$summary.fixed)
# write.csv(model.SW.SD.Mobile$summary.hyperpar)
# cat('Posterior Mean of Variance: ')
# cat('\n')
# cat(inla_emarginal(model.SW.SD.Mobile))
# cat('\n')
# cat('Posterior Mode of Variance: ')
# cat('\n')
# cat(inla_mmarginal(model.SW.SD.Mobile))
# cat("\n")
# cat('WAIC: ')
# cat('\n')
# cat(model.SW.SD.Mobile$waic$waic)
# cat("\n")
# cat('DIC: ')
# cat('\n')
# cat(model.SW.SD.Mobile$dic$dic)
# cat("\n")
# cat('Marginal Likelihood: ')
# cat('\n')
# cat(model.SW.SD.Mobile$mlik[2,1])
# cat("\n")
# cat("\n")


cat('Mobile ~ Wind Exposure + Wind Chill + WindExp*WC + WindExp|IND')
cat('\n')
formula.random <- Use ~  -1 + Wind.Exp.Z + WC + Wind.Exp.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Wind.Exp.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)),
                                    control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.WindExp.WCprev.Mobile$summary.fixed)
write.csv(model.WindExp.WCprev.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.WCprev.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.WCprev.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.WindExp.WCprev.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.WindExp.WCprev.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.WindExp.WCprev.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Wind Exposure + SnowDepth + WindExp*SD + WindExp|IND')
cat('\n')
formula.random <- Use ~ -1 + Wind.Exp.Z + SD + Wind.Exp.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Wind.Exp.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)),
                                control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.WindExp.SD.Mobile$summary.fixed)
write.csv(model.WindExp.SD.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.SD.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.SD.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.WindExp.SD.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.WindExp.SD.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.WindExp.SD.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion Ag + Wind Chill + PropAg*WC + PropAg|IND')
cat('\n')
formula.random <- Use ~  -1 + PropAg.Z + WC + PropAg.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropAg.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropAg.WCprev.Mobile$summary.fixed)
write.csv(model.PropAg.WCprev.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.WCprev.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.WCprev.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropAg.WCprev.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropAg.WCprev.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropAg.WCprev.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion Ag + SnowDepth + PropAg*SD + PropAg|IND')
cat('\n')
formula.random <- Use ~ -1 + PropAg.Z + SD + PropAg.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropAg.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)),
                               control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropAg.SD.Mobile$summary.fixed)
write.csv(model.PropAg.SD.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.SD.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.SD.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropAg.SD.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropAg.SD.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropAg.SD.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion Dev + Wind Chill + PropDev*WC + PropDev|IND')
cat('\n')
formula.random <- Use ~  -1 + PropDev.Z + WC + PropDev.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropDev.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)),
                                    control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropDev.WCprev.Mobile$summary.fixed)
write.csv(model.PropDev.WCprev.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.PropDev.WCprev.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.WCprev.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropDev.WCprev.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropDev.WCprev.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropDev.WCprev.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion Dev + SnowDepth + PropDev*SD + PropDev|IND')
cat('\n')
formula.random <- Use ~ -1 + PropDev.Z + SD + PropDev.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropDev.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)),
                                control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropDev.SD.Mobile$summary.fixed)
write.csv(model.PropDev.SD.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropDev.SD.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.SD.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropDev.SD.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropDev.SD.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropDev.SD.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion SW + Wind Chill + PropSW*WC + PropSW|IND')
cat('\n')
formula.random <- Use ~  -1 + PropSW.Z + WC + PropSW.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropSW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)),
                                   control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropSW.WCprev.Mobile$summary.fixed)
write.csv(model.PropSW.WCprev.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropSW.WCprev.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropSW.WCprev.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropSW.WCprev.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropSW.WCprev.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropSW.WCprev.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion SW + SnowDepth + PropSW*SD + PropSW|IND')
cat('\n')
formula.random <- Use ~ -1 + PropSW.Z + SD + PropSW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropSW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)),
                               control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropSW.SD.Mobile$summary.fixed)
write.csv(model.PropSW.SD.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropSW.SD.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropSW.SD.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropSW.SD.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropSW.SD.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropSW.SD.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion FoodSubsidy + Wind Chill + PropFS*WC + PropFS|IND')
cat('\n')
formula.random <- Use ~  -1 + PropFoodSub.Z + WC + PropFoodSub.Z*WC + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropFoodSub.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.WCprev.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                        control.fixed = list(
                                          mean = mean.beta,
                                          prec = list(default = prec.beta)),
                                        control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropFoodSub.WCprev.Mobile$summary.fixed)
write.csv(model.PropFoodSub.WCprev.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.WCprev.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.WCprev.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropFoodSub.WCprev.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropFoodSub.WCprev.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropFoodSub.WCprev.Mobile$mlik[2,1])
cat("\n")
cat("\n")


cat('Mobile ~ Proportion FoodSub + SnowDepth + PropFoodSub*SD + PropFoodSub|IND')
cat('\n')
formula.random <- Use ~ -1 + PropFoodSub.Z + SD + PropFoodSub.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropFoodSub.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.SD.Mobile <- inla(formula.random, family ="Poisson", data=Mobile.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)),
                                    control.compute = list(mlik = TRUE, dic = TRUE, waic = TRUE))
write.csv(model.PropFoodSub.SD.Mobile$summary.fixed)
write.csv(model.PropFoodSub.SD.Mobile$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.SD.Mobile))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.SD.Mobile))
cat("\n")
cat('WAIC: ')
cat('\n')
cat(model.PropFoodSub.SD.Mobile$waic$waic)
cat("\n")
cat('DIC: ')
cat('\n')
cat(model.PropFoodSub.SD.Mobile$dic$dic)
cat("\n")
cat('Marginal Likelihood: ')
cat('\n')
cat(model.PropFoodSub.SD.Mobile$mlik[2,1])
cat("\n")
cat("\n")


sink()
