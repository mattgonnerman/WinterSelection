############
### INLA ###
############

setwd("E:/Maine Drive/Analysis/Kaj Thesis")

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
RoostAll.final <- read.csv("RoostAll.csv")
DayForest.final <- read.csv("DayForest.csv")
DayAll.final <- read.csv("DayAll.csv")

#Set mean and precision priors for slope coefficients
mean.beta <- 0
prec.beta <- 1e-4  


################################################################################### 
### 1st Scale - Roost Points
sink("RoostResults.csv")
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
                            prec = list(default = prec.beta)))
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
                            prec = list(default = prec.beta)))
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
                            prec = list(default = prec.beta)))
write.csv(model.DtFE.WCprev.roost$summary.fixed)
write.csv(model.DtFE.WCprev.roost$summary.hyperpar)
cat("Posterior Mean of Variance: ")
cat('\n')
cat(inla_emarginal(model.DtFE.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.WCprev.roost))
cat('\n')
cat('\n')


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
                            prec = list(default = prec.beta)))
write.csv(model.DtFE.SD.roost$summary.fixed)
write.csv(model.DtFE.SD.roost$summary.hyperpar)
cat('Posterior Meand of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.SD.roost))
cat('\n')
cat('\n')


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
                              prec = list(default = prec.beta)))
write.csv(model.SW.WCprev.roost$summary.fixed)
write.csv(model.SW.WCprev.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.SW.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.SW.WCprev.roost))
cat('\n')
cat('\n')


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
                            prec = list(default = prec.beta)))
write.csv(model.SW.SD.roost$summary.fixed)
write.csv(model.SW.SD.roost$summary.hyperpar)
cat("Posterior Mean of Variance")
cat('\n')
cat(inla_emarginal(model.SW.SD.roost))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.SW.SD.roost))
cat('\n')
cat('\n')


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
                                prec = list(default = prec.beta)))
write.csv(model.Ht.WCprev.roost$summary.fixed)
write.csv(model.Ht.WCprev.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.WCprev.roost))
cat('\n')
cat('\n')


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
                            prec = list(default = prec.beta)))
write.csv(model.Ht.SD.roost$summary.fixed)
write.csv(model.Ht.SD.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.SD.roost))
cat('\n')
cat('\n')


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
                                prec = list(default = prec.beta)))
write.csv(model.WindExp.WCprev.roost$summary.fixed)
write.csv(model.WindExp.WCprev.roost$summary.hyperpar)
cat('Posterio Mean of Variance')
cat('\n')
cat(inla_emarginal(model.WindExp.WCprev.roost))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.WindExp.WCprev.roost))
cat('\n')
cat('\n')


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
                            prec = list(default = prec.beta)))
write.csv(model.WindExp.SD.roost$summary.fixed)
write.csv(model.WindExp.SD.roost$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.SD.roost))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.SD.roost))


sink()


################################################################################### 
### 2nd Scale - Day Use, Forested Points only
sink('DayForestResults.csv')
cat('DayTree ~ Basal Area + Wind Chill + BA*WC + BA|IND')
cat('\n')
formula.random <- Use ~  -1 + BA.Z + WC_prev + BA.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,BA.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.BA.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)))
write.csv(model.BA.WCprev.daytree$summary.fixed)
write.csv(model.BA.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.BA.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.BA.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Basal Area + SnowDepth + BA*SD + BA|IND')
cat('\n')
formula.random <- Use ~ -1 + BA.Z + SD + BA.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,BA.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.BA.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)))
write.csv(model.BA.SD.daytree$summary.fixed)
write.csv(model.BA.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.BA.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.BA.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Distance to Edge + Wind Chill + DtE*WC + DtE|IND')
cat('\n')
formula.random <- Use ~  -1 + DtFE.Z + WC_prev + DtFE.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,DtFE.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)))
write.csv(model.DtFE.WCprev.daytree$summary.fixed)
write.csv(model.DtFE.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Distance to Edge + SnowDepth + DtE*SD + DtE|IND')
cat('\n')
formula.random <- Use ~ -1 + DtFE.Z + SD + DtFE.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,DtFE.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                            control.fixed = list(
                              mean = mean.beta,
                              prec = list(default = prec.beta)))
write.csv(model.DtFE.SD.daytree$summary.fixed)
write.csv(model.DtFE.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Percent Softwood + Wind Chill + SW*WC + SW|IND')
cat('\n')
formula.random <- Use ~  -1 + SW.Z + WC_prev + SW.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,SW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.SW.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)))
write.csv(model.SW.WCprev.daytree$summary.fixed)
write.csv(model.SW.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.SW.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.SW.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Percent Softwood + SnowDepth + SW*SD + SW|IND')
cat('\n')
formula.random <- Use ~ -1 + SW.Z + SD + SW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,SW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.SW.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)))
write.csv(model.SW.SD.daytree$summary.fixed)
write.csv(model.SW.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.SW.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.SW.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Mean Height + Wind Chill + Ht*WC + Ht|IND')
cat('\n')
formula.random <- Use ~  -1 + Ht.Z + WC_prev + Ht.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Ht.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.Ht.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)))
write.csv(model.Ht.WCprev.daytree$summary.fixed)
write.csv(model.Ht.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Mean Height + SnowDepth + Ht*SD + Ht|IND')
cat('\n')
formula.random <- Use ~ -1 + Ht.Z + SD + Ht.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Ht.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.Ht.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)))
write.csv(model.Ht.SD.daytree$summary.fixed)
write.csv(model.Ht.SD.daytree$summary.hyperpar)
cat('Posteruir Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.Ht.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.Ht.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Wind Exposure + Wind Chill + WindExp*WC + WindExp|IND')
cat('\n')
formula.random <- Use ~  -1 + Wind.Exp.Z + WC_prev + Wind.Exp.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Wind.Exp.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)))
write.csv(model.WindExp.WCprev.daytree$summary.fixed)
write.csv(model.WindExp.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Wind Exposure + SnowDepth + WindExp*SD + WindExp|IND')
cat('\n')
formula.random <- Use ~ -1 + Wind.Exp.Z + SD + Wind.Exp.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Wind.Exp.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)))
write.csv(model.WindExp.SD.daytree$summary.fixed)
write.csv(model.WindExp.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion Ag + Wind Chill + PropAg*WC + PropAg|IND')
cat('\n')
formula.random <- Use ~  -1 + PropAg.Z + WC_prev + PropAg.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropAg.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                     control.fixed = list(
                                       mean = mean.beta,
                                       prec = list(default = prec.beta)))
write.csv(model.PropAg.WCprev.daytree$summary.fixed)
write.csv(model.PropAg.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion Ag + SnowDepth + PropAg*SD + PropAg|IND')
cat('\n')
formula.random <- Use ~ -1 + PropAg.Z + SD + PropAg.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropAg.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                 control.fixed = list(
                                   mean = mean.beta,
                                   prec = list(default = prec.beta)))
write.csv(model.PropAg.SD.daytree$summary.fixed)
write.csv(model.PropAg.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion Dev + Wind Chill + PropDev*WC + PropDev|IND')
cat('\n')
formula.random <- Use ~  -1 + PropDev.Z + WC_prev + PropDev.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropDev.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)))
write.csv(model.PropDev.WCprev.daytree$summary.fixed)
write.csv(model.PropDev.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropDev.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion Dev + SnowDepth + PropDev*SD + PropDev|IND')
cat('\n')
formula.random <- Use ~ -1 + PropDev.Z + SD + PropDev.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropDev.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)))
write.csv(model.PropDev.SD.daytree$summary.fixed)
write.csv(model.PropDev.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropDev.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion SW + Wind Chill + PropSW*WC + PropSW|IND')
cat('\n')
formula.random <- Use ~  -1 + PropSW.Z + WC_prev + PropSW.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropSW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                     control.fixed = list(
                                       mean = mean.beta,
                                       prec = list(default = prec.beta)))
write.csv(model.PropSW.WCprev.daytree$summary.fixed)
write.csv(model.PropSW.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.PropSW.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropSW.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion SW + SnowDepth + PropSW*SD + PropSW|IND')
cat('\n')
formula.random <- Use ~ -1 + PropSW.Z + SD + PropSW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropSW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                 control.fixed = list(
                                   mean = mean.beta,
                                   prec = list(default = prec.beta)))
write.csv(model.PropSW.SD.daytree$summary.fixed)
write.csv(model.PropSW.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropSW.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropSW.SD.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion FoodSub + Wind Chill + PropFoodSub*WC + PropFoodSub|IND')
cat('\n')
formula.random <- Use ~  -1 + PropFoodSub.Z + WC_prev + PropFoodSub.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropFoodSub.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.WCprev.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)))
write.csv(model.PropFoodSub.WCprev.daytree$summary.fixed)
write.csv(model.PropFoodSub.WCprev.daytree$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.WCprev.daytree))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.WCprev.daytree))
cat('\n')
cat('\n')


cat('DayTree ~ Proportion FoodSub + SnowDepth + PropFoodSub*SD + PropFoodSub|IND')
cat('\n')
formula.random <- Use ~ -1 + PropFoodSub.Z + SD + PropFoodSub.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropFoodSub.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.SD.daytree <- inla(formula.random, family ="Poisson", data=DayForest.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)))
write.csv(model.PropFoodSub.SD.daytree$summary.fixed)
write.csv(model.PropFoodSub.SD.daytree$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.SD.daytree))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.SD.daytree))
cat('\n')
cat('\n')



sink()


################################################################################### 
### 3rd Scale - Day Use, All Points
sink('DayAllResults.csv')


cat('AllDay ~ Distance to Edge + Wind Chill + DtE*WC + DtE|IND')
cat('\n')
formula.random <- Use ~  -1 + DtFE.Z + WC_prev + DtFE.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,DtFE.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                  control.fixed = list(
                                    mean = mean.beta,
                                    prec = list(default = prec.beta)))
write.csv(model.DtFE.WCprev.allday$summary.fixed)
write.csv(model.DtFE.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.DtFE.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Distance to Edge + SnowDepth + DtE*SD + DtE|IND')
cat('\n')
formula.random <- Use ~ -1 + DtFE.Z + SD + DtFE.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,DtFE.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.DtFE.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                              control.fixed = list(
                                mean = mean.beta,
                                prec = list(default = prec.beta)))
write.csv(model.DtFE.SD.allday$summary.fixed)
write.csv(model.DtFE.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.DtFE.SD.allday))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.DtFE.SD.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Percent Softwood + Wind Chill + SW*WC + SW|IND')
cat('\n')
formula.random <- Use ~  -1 + SW.Z + WC_prev + SW.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,SW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.SW.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)))
write.csv(model.SW.WCprev.allday$summary.fixed)
write.csv(model.SW.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.SW.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance')
cat('\n')
cat(inla_mmarginal(model.SW.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Percent Softwood + SnowDepth + SW*SD + SW|IND')
cat('\n')
formula.random <- Use ~ -1 + SW.Z + SD + SW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,SW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.SW.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                            control.fixed = list(
                              mean = mean.beta,
                              prec = list(default = prec.beta)))
write.csv(model.SW.SD.allday$summary.fixed)
write.csv(model.SW.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.SW.SD.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.SW.SD.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Wind Exposure + Wind Chill + WindExp*WC + WindExp|IND')
cat('\n')
formula.random <- Use ~  -1 + Wind.Exp.Z + WC_prev + Wind.Exp.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,Wind.Exp.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                     control.fixed = list(
                                       mean = mean.beta,
                                       prec = list(default = prec.beta)))
write.csv(model.WindExp.WCprev.allday$summary.fixed)
write.csv(model.WindExp.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Wind Exposure + SnowDepth + WindExp*SD + WindExp|IND')
cat('\n')
formula.random <- Use ~ -1 + Wind.Exp.Z + SD + Wind.Exp.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,Wind.Exp.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.WindExp.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                 control.fixed = list(
                                   mean = mean.beta,
                                   prec = list(default = prec.beta)))
write.csv(model.WindExp.SD.allday$summary.fixed)
write.csv(model.WindExp.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.WindExp.SD.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.WindExp.SD.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion Ag + Wind Chill + PropAg*WC + PropAg|IND')
cat('\n')
formula.random <- Use ~  -1 + PropAg.Z + WC_prev + PropAg.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropAg.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)))
write.csv(model.PropAg.WCprev.allday$summary.fixed)
write.csv(model.PropAg.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion Ag + SnowDepth + PropAg*SD + PropAg|IND')
cat('\n')
formula.random <- Use ~ -1 + PropAg.Z + SD + PropAg.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropAg.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropAg.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)))
write.csv(model.PropAg.SD.allday$summary.fixed)
write.csv(model.PropAg.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropAg.SD.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropAg.SD.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion Dev + Wind Chill + PropDev*WC + PropDev|IND')
cat('\n')
formula.random <- Use ~  -1 + PropDev.Z + WC_prev + PropDev.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropDev.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                     control.fixed = list(
                                       mean = mean.beta,
                                       prec = list(default = prec.beta)))
write.csv(model.PropDev.WCprev.allday$summary.fixed)
write.csv(model.PropDev.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance')
cat('\n')
cat(inla_emarginal(model.PropDev.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion Dev + SnowDepth + PropDev*SD + PropDev|IND')
cat('\n')
formula.random <- Use ~ -1 + PropDev.Z + SD + PropDev.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropDev.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropDev.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                 control.fixed = list(
                                   mean = mean.beta,
                                   prec = list(default = prec.beta)))
write.csv(model.PropDev.SD.allday$summary.fixed)
write.csv(model.PropDev.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropDev.SD.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropDev.SD.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion SW + Wind Chill + PropSW*WC + PropSW|IND')
cat('\n')
formula.random <- Use ~  -1 + PropSW.Z + WC_prev + PropSW.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropSW.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                    control.fixed = list(
                                      mean = mean.beta,
                                      prec = list(default = prec.beta)))
write.csv(model.PropSW.WCprev.allday$summary.fixed)
write.csv(model.PropSW.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropSW.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropSW.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion SW + SnowDepth + PropSW*SD + PropSW|IND')
cat('\n')
formula.random <- Use ~ -1 + PropSW.Z + SD + PropSW.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropSW.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropSW.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                control.fixed = list(
                                  mean = mean.beta,
                                  prec = list(default = prec.beta)))
write.csv(model.PropSW.SD.allday$summary.fixed)
write.csv(model.PropSW.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropSW.SD.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropSW.SD.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion FoodSubsidy + Wind Chill + PropFS*WC + PropFS|IND')
cat('\n')
formula.random <- Use ~  -1 + PropFoodSub.Z + WC_prev + PropFoodSub.Z*WC_prev + #Fixed effects
  StepLength.Z + #step length
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) + #Conditional
  f(ID,PropFoodSub.Z,values=1:26,model="iid", #Random Slope
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.WCprev.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                                   control.fixed = list(
                                     mean = mean.beta,
                                     prec = list(default = prec.beta)))
write.csv(model.PropFoodSub.WCprev.allday$summary.fixed)
write.csv(model.PropFoodSub.WCprev.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.WCprev.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.WCprev.allday))
cat('\n')
cat('\n')


cat('AllDay ~ Proportion FoodSub + SnowDepth + PropFoodSub*SD + PropFoodSub|IND')
cat('\n')
formula.random <- Use ~ -1 + PropFoodSub.Z + SD + PropFoodSub.Z*SD +
  StepLength.Z + 
  f(StratID, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ID,PropFoodSub.Z,values=1:26,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

model.PropFoodSub.SD.allday <- inla(formula.random, family ="Poisson", data=DayAll.final, 
                               control.fixed = list(
                                 mean = mean.beta,
                                 prec = list(default = prec.beta)))
write.csv(model.PropFoodSub.SD.allday$summary.fixed)
write.csv(model.PropFoodSub.SD.allday$summary.hyperpar)
cat('Posterior Mean of Variance: ')
cat('\n')
cat(inla_emarginal(model.PropFoodSub.SD.allday))
cat('\n')
cat('Posterior Mode of Variance: ')
cat('\n')
cat(inla_mmarginal(model.PropFoodSub.SD.allday))
cat('\n')
cat('\n')


sink()
