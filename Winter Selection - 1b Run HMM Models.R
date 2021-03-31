### HMM Options Common Between Models ###
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 3 # Number of states
stateNames <- c("roost", "loafing", "foraging") # label states
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))} #prevents working parameters from straying along boundary


##############################################################
### MODEL 0 - Nothing Specified
# initial parameters
Par0_m0.zm <- list(step=c(2, 50,120, #mean
                          2, 50,100, #sd
                          .99, .005, .005), #zero mass
                   angle = c(.01, 0.2 ,0.3))

# fit model
turk_m0.zm <- fitHMM(data = turkeyData.zm,
                     retryFits = 1,
                     nbStates = nbStates,
                     dist = dist,
                     Par0 = Par0_m0.zm,
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames,
                     ncores = (detectCores()/2))

##############################################################
### MODEL 1
stepDM<-matrix(c(
  1,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
                c("mean_123:(Intercept)", "mean_2","mean_3",
                  paste0("sd_",1:nbStates,":(Intercept)"),
                  paste0("zero_",1:nbStates,":(Intercept)"))))
#define the directions of the differences
stepworkBounds <- matrix(c(-Inf, 0,0,-Inf, 0,0,rep(-Inf,3),rep(Inf,9)),nrow = 3*nbStates,
                         dimnames=list(colnames(stepDM),c("lower","upper")))

stepBounds <- matrix(c(0,5,
                       0,Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       .98,1,
                       0,.02,
                       0,.02), nrow = 3*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters:
# Concentration -> searching < dispersal
angleDM<-matrix(c(1,0,0,
                  1,1,0,
                  1,1,1),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_1:(Intercept)","concentration_2","concentration_3")))

#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
                          nrow = ncol(angleDM),
                          dimnames=list(colnames(angleDM),c("lower","upper")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,0.94,
                        0,0.94,
                        0,0.94),nrow = nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper")))

#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,
                 angle=angleworkBounds)
userBounds <- list(step = stepBounds,
                   angle = angleBounds)


Par_m1.zm <- getPar0(model = turk_m0.zm,
                     formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
                     DM = DM,
                     workBounds = workBounds,
                     userBounds = userBounds)


# fit model
turk_m1.zm <- fitHMM(data = turkeyData.zm,
                     retryFits = 1,
                     nbStates = nbStates,
                     dist = dist,
                     Par0 = Par_m1.zm$Par,
                     beta0 = Par_m1.zm$beta,
                     DM = DM,
                     workBounds = workBounds,
                     userBounds = userBounds,
                     estAngleMean = list(angle=FALSE),
                     prior = prior,
                     stateNames = stateNames,
                     ncores = (detectCores()/2),
                     formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID
)


##############################################################
### MODEL 2 - Restriction on Stationary Mean Step length
#Step Length Quantiles - Non Roost movements
quantile(turkeyData.zm[turkeyData.zm$step > 0,]$step, seq(0, .9, .1), na.rm = T)

stepDM<-matrix(c(
  1,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
                c("mean_123:(Intercept)", "mean_2","mean_3",
                  paste0("sd_",1:nbStates,":(Intercept)"),
                  paste0("zero_",1:nbStates,":(Intercept)"))))
#define the directions of the differences
stepworkBounds <- matrix(c(-Inf, 0,0,-Inf, 0,0,rep(-Inf,3),rep(Inf,9)),nrow = 3*nbStates,
                         dimnames=list(colnames(stepDM),c("lower","upper")))

stepBounds <- matrix(c(0,5,
                       0,20,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       .98,1,
                       0,.02,
                       0,.02), nrow = 3*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters:
# Concentration -> searching < dispersal
angleDM<-matrix(c(1,0,0,
                  1,1,0,
                  1,1,1),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_1:(Intercept)","concentration_2","concentration_3")))

#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
                          nrow = ncol(angleDM),
                          dimnames=list(colnames(angleDM),c("lower","upper")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,0.94,
                        0,0.94,
                        0,0.94),nrow = nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper")))

#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,
                 angle=angleworkBounds)
userBounds <- list(step = stepBounds,
                   angle = angleBounds)


Par_m2.zm <- getPar0(model = turk_m1.zm,
                     formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
                     DM = DM,
                     workBounds = workBounds,
                     userBounds = userBounds)


# fit model
turk_m2.zm <- fitHMM(data = turkeyData.zm,
                     retryFits = retryFits,
                     nbStates = nbStates,
                     dist = dist,
                     Par0 = Par_m2.zm$Par,
                     beta0 = Par_m2.zm$beta,
                     DM = DM,
                     workBounds = workBounds,
                     userBounds = userBounds,
                     estAngleMean = list(angle=FALSE),
                     prior = prior,
                     stateNames = stateNames,
                     ncores = (detectCores()/2),
                     formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID
)




# ##############################################################
# ### MODEL 1 - Daily Cycle in Transition Probability
# # initial parameters
# Par0_m1.zm <- getPar0(model=turk_m0.zm,
#                       formula=~cosinor(hour, period = 24))
# 
# # fit model
# turk_m1.zm <- fitHMM(data = turkeyData.zm,
#                      retryFits = retryFits,
#                      nbStates = nbStates,
#                      dist = dist,
#                      Par0 = Par0_m1.zm$Par,
#                      beta0 = Par0_m1.zm$beta,
#                      estAngleMean = list(angle=FALSE),
#                      stateNames = stateNames,
#                      ncores = 5,
#                      formula = ~cosinor(hour, period = 24)
#                      )
# 
# 
# ##############################################################
# ### MODEL 2 - Daily Cycle and Weather Effects on Transition Probability
# # initial parameters
# Par0_m2.zm <- getPar0(model=turk_m1.zm,
#                       formula= ~ WC.Z + SD.Z + cosinor(hour, period = 24),
#                       covNames=c("SD.Z", "WC.Z"))
# 
# # fit model
# turk_m2.zm <- fitHMM(data = turkeyData.zm,
#                      retryFits = retryFits,
#                      nbStates = nbStates,
#                      dist = dist,
#                      Par0 = Par0_m2.zm$Par,
#                      beta0 = Par0_m2.zm$beta,
#                      estAngleMean = list(angle=FALSE),
#                      stateNames = stateNames,
#                      ncores = 5,
#                      formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24)
#                      )
# 
# 
# ##############################################################
# ### MODEL 3 - Full Transition Probability Formula, ID and Daily Cycle affect Step Length
# # formulas for parameters of state-dependent observation distributions
# DM <- list(step = list(mean = ~ ID + cosinor(hour, period = 24),
#                        sd = ~ ID + cosinor(hour, period = 24),
#                        zeromass = ~ 1))
# 
# # initial parameters
# Par0_m3.zm <- getPar0(model = turk_m2.zm,
#                        formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
#                        DM = DM,
#                        covNames=c("SD.Z", "WC.Z"))
# 
# # fit model
# turk_m3.zm <- fitHMM(data = turkeyData.zm,
#                       retryFits = retryFits,
#                       nbStates = nbStates,
#                       dist = dist,
#                       Par0 = Par0_m3.zm$Par,
#                       beta0 = Par0_m3.zm$beta,
#                       estAngleMean = list(angle=FALSE),
#                       stateNames = stateNames,
#                       DM = DM,
#                       ncores = 5,
#                       formula = ~WC.Z + SD.Z + cosinor(hour, period = 24) + ID
# )

# ##############################################################
# ### MODEL 4 
# # DOES NOT WORK WELL
# DM <- list(step = list(mean = ~ SD.Z + WC.Z + ID + cosinor(hour, period = 24),
#                        sd = ~ SD.Z + WC.Z + ID + cosinor(hour, period = 24),
#                        zeromass = ~ 1))
# 
# # initial parameters
# Par0_m4.zm <- getPar0(model = turk_m0.zm,
#                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
#                       DM = DM,
#                       covNames=c("SD.Z", "WC.Z"))
# 
# # fit model
# turk_m4.zm <- fitHMM(data = turkeyData.zm,
#                      retryFits = retryFits,
#                      nbStates = nbStates,
#                      dist = dist,
#                      Par0 = Par0_m4.zm$Par,
#                      beta0 = Par0_m4.zm$beta,
#                      estAngleMean = list(angle=FALSE),
#                      stateNames = stateNames,
#                      DM = DM,
#                      ncores = 5,
#                      formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID
# )




# Par_m5b.zm <- getPar0(model = turk_m0.zm,
#                      formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24),
#                      DM = DM,
#                      workBounds = workBounds,
#                      userBounds = userBounds)
# 
# 
# # fit model
# turk_m5b.zm <- fitHMM(data = turkeyData.zm,
#                      # retryFits = retryFits,
#                      nbStates = nbStates,
#                      dist = dist,
#                      Par0 = Par_m5b.zm$Par,
#                      beta0 = Par_m5b.zm$beta,
#                      DM = DM,
#                      workBounds = workBounds,
#                      userBounds = userBounds,
#                      estAngleMean = list(angle=FALSE),
#                      prior = prior,
#                      stateNames = stateNames,
#                      ncores = 5,
#                      formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24)
# )

# ##############################################################
# ### MODEL 
# # formulas for parameters of state-dependent observation distributions
# stepBounds <- matrix(c(0,5,
#                        0,50,
#                        50, Inf,
#                        0, Inf,
#                        0, Inf,
#                        0, Inf,
#                        .98,1,
#                        0,.02,
#                        0,.02), nrow = 3*nbStates, byrow = T,
#                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# 
# ## constrain turning angle concentration parameters:
# # Concentration -> searching < dispersal
# angleDM<-matrix(c(1,0,0,
#                   1,1,0,
#                   1,1,1),nrow = nbStates,byrow=TRUE,
#                 dimnames=list(paste0("concentration_",1:nbStates),
#                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# 
# #define direction of differences
# angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
#                           nrow = ncol(angleDM),
#                           dimnames=list(colnames(angleDM),c("lower","upper")))
# 
# #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# #Userbound contraint on angle
# angleBounds <- matrix(c(0,0.94,
#                         0,0.94,
#                         0,0.94),nrow = nbStates,
#                       byrow=TRUE, dimnames=list(rownames(angleDM),
#                                                 c("lower","upper")))
# 
# #Bundle individual parameter DM and workbounds
# DM<-list(step=list(mean = ~ cosinor(hour, period = 24),
#                     sd = ~ 1,
#                     zeromass = ~ 1),
#          angle=angleDM)
# workBounds<-list(angle=angleworkBounds)
# userBounds <- list(step = stepBounds,
#                    angle = angleBounds)
# # initial parameters
# Par0_m4.zm <- getPar0(model = turk_m2.zm,
#                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
#                       DM = DM,
#                       workBounds = workBounds,
#                       userBounds = userBounds,
#                       covNames=c("SD.Z", "WC.Z"))
# 
# # fit model
# turk_m4.zm <- fitHMM(data = turkeyData.zm,
#                      retryFits = retryFits,
#                      nbStates = nbStates,
#                      dist = dist,
#                      Par0 = Par0_m4.zm$Par,
#                      beta0 = Par0_m4.zm$beta,
#                      estAngleMean = list(angle=FALSE),
#                      stateNames = stateNames,
#                      DM = DM,
#                      workBounds = workBounds,
#                      userBounds = userBounds,
#                      ncores = 5,
#                      formula = ~WC.Z + SD.Z + cosinor(hour, period = 24) + ID
# ) 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# {
# 
# ### CANDIDATE MODELS ###
# ### CANDIDATE MODELS ###
# ## A = Parameter Restrictions
# # Step DM - Step Mean and SD, Roost < Loafing < Foraging
# # Step Bound - Step Mean for Roost < 5; Zero Mass for roost >.98
# # Angle DM - Roost < Loafing < Foraging
# # Angle Bounds - All 0< & <.94 
# ## B = Step Length
# # ID in all models
# ## C = Transition Probabilities
# 
# #Step Length Quantiles - Non Roost movements
# quantile(turkeyData.zm[turkeyData.zm$step > 0,]$step, seq(0.5, .9, .1), na.rm = T)
# 
# ########################################################################
# ################################## A ###################################
# ########################################################################
# 
# ### Same for all A models
# stepDM<-matrix(c(
#   1,0,0,0,0,0,0,0,0,
#   1,1,0,0,0,0,0,0,0,
#   1,1,1,0,0,0,0,0,0,
#   0,0,0,1,0,0,0,0,0,
#   0,0,0,1,1,0,0,0,0,
#   0,0,0,1,1,1,0,0,0,
#   0,0,0,0,0,0,1,0,0,
#   0,0,0,0,0,0,0,1,0,
#   0,0,0,0,0,0,0,0,1),
#   nrow = 3*nbStates,byrow=TRUE,
#   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
#                 c("mean_123:(Intercept)", "mean_2","mean_3",
#                   paste0("sd_",1:nbStates,":(Intercept)"),
#                   paste0("zero_",1:nbStates,":(Intercept)"))))
# #define the directions of the differences
# stepworkBounds <- matrix(c(-Inf, 0,0,-Inf, 0,0,rep(-Inf,3),rep(Inf,9)),nrow = 3*nbStates,
#                          dimnames=list(colnames(stepDM),c("lower","upper")))
# 
# ## constrain turning angle concentration parameters:
# # Concentration -> searching < dispersal
# angleDM<-matrix(c(1,0,0,
#                   1,1,0,
#                   1,1,1),nrow = nbStates,byrow=TRUE,
#                 dimnames=list(paste0("concentration_",1:nbStates),
#                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# 
# #define direction of differences
# angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
#                           nrow = ncol(angleDM),
#                           dimnames=list(colnames(angleDM),c("lower","upper")))
# 
# #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# #Userbound contraint on angle
# angleBounds <- matrix(c(0,0.94,
#                         0,0.94,
#                         0,0.94),nrow = nbStates,
#                       byrow=TRUE, dimnames=list(rownames(angleDM),
#                                                 c("lower","upper")))
# 
# ##############################################################
# ### MODEL A3 - 70% quantile for loafing step length
# #Userbound constraint on step
# stepBounds <- matrix(c(0,5,
#                        0,109.69861,
#                        0, Inf,
#                        0, Inf,
#                        0, Inf,
#                        0, Inf,
#                        .98,1,
#                        0,.02,
#                        0,.02), nrow = 3*nbStates, byrow = T,
#                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# 
# #Bundle individual parameter DM and workbounds
# DM<-list(step=stepDM,angle=angleDM)
# workBounds<-list(step=stepworkBounds,
#                  angle=angleworkBounds)
# userBounds <- list(step = stepBounds,
#                    angle = angleBounds)
# 
# # initial parameters
# Par <- list(step=c(2, 75,150, #mean
#                    2, 10,25, #sd
#                    .99, .005, .005), #zero mass
#             angle = c(.01, 0.2 ,0.3))
# 
# Par0_mA3.zm <- getParDM(data = turkeyData.zm,
#                         nbStates = nbStates,
#                         zeroInflation = list(step = T,
#                                              angle = T),
#                         dist = dist,
#                         Par = Par,
#                         DM = DM,
#                         workBounds = workBounds,
#                         userBounds = userBounds,
#                         estAngleMean = list(angle = FALSE))
# 
# # fit model
# turk_mA3.zm <- fitHMM(data = turkeyData.zm,
#                       retryFits = retryFits,
#                       nbStates = nbStates,
#                       dist = dist,
#                       Par0 = Par0_mA3.zm,
#                       DM = DM,
#                       workBounds = workBounds,
#                       userBounds = userBounds,
#                       estAngleMean = list(angle=FALSE),
#                       prior = prior,
#                       stateNames = stateNames,
#                       ncores = 5
# )
# {
#   # ##############################################################
#   # ### MODEL A0 - No limit for loafing step length
#   # #Userbound constraint on step
#   # stepBounds <- matrix(c(0, 5,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        .98, 1,
#   #                        0, .02,
#   #                        0, .02), nrow = 3*nbStates, byrow = T,
#   #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
#   # 
#   # #Bundle individual parameter DM and workbounds
#   # DM<-list(step=stepDM,angle=angleDM)
#   # workBounds<-list(step=stepworkBounds,
#   #                  angle=angleworkBounds)
#   # userBounds <- list(step = stepBounds,
#   #                    angle = angleBounds)
#   # 
#   # # initial parameters
#   # Par <- list(step=c(2, 45,200, #mean
#   #                    2, 10,25, #sd
#   #                    .99, .005, .005), #zero mass
#   #             angle = c(.01, 0.2 ,0.3))
#   # 
#   # Par0_mA0.zm <- getParDM(data = turkeyData.zm,
#   #                         nbStates = nbStates,
#   #                         zeroInflation = list(step = T,
#   #                                              angle = T),
#   #                         dist = dist,
#   #                         Par = Par,
#   #                         DM = DM,
#   #                         workBounds = workBounds,
#   #                         userBounds = userBounds,
#   #                         estAngleMean = list(angle = FALSE))
#   # 
#   # # fit model
#   # turk_mA0.zm <- fitHMM(data = turkeyData.zm,
#   #                       retryFits = retryFits,
#   #                       nbStates = nbStates,
#   #                       dist = dist,
#   #                       Par0 = Par0_mA0.zm,
#   #                       DM = DM,
#   #                       workBounds = workBounds,
#   #                       userBounds = userBounds,
#   #                       estAngleMean = list(angle=FALSE),
#   #                       prior = prior,
#   #                       stateNames = stateNames,
#   #                       ncores = 5
#   # )
#   # 
#   # ##############################################################
#   # ### MODEL A1 - 50% quantile for loafing step length
#   # #Userbound constraint on step
#   # stepBounds <- matrix(c(0,5,
#   #                        0,51.04313,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        .98,1,
#   #                        0,.02,
#   #                        0,.02), nrow = 3*nbStates, byrow = T,
#   #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
#   # 
#   # #Bundle individual parameter DM and workbounds
#   # DM<-list(step=stepDM,angle=angleDM)
#   # workBounds<-list(step=stepworkBounds,
#   #                  angle=angleworkBounds)
#   # userBounds <- list(step = stepBounds,
#   #                    angle = angleBounds)
#   # 
#   # # initial parameters
#   # Par <- list(step=c(2, 25,150, #mean
#   #                    2, 10,25, #sd
#   #                    .99, .005, .005), #zero mass
#   #             angle = c(.01, 0.2 ,0.3))
#   # 
#   # Par0_mA1.zm <- getParDM(data = turkeyData.zm,
#   #                         nbStates = nbStates,
#   #                         zeroInflation = list(step = T,
#   #                                              angle = T),
#   #                         dist = dist,
#   #                         Par = Par,
#   #                         DM = DM,
#   #                         workBounds = workBounds,
#   #                         userBounds = userBounds,
#   #                         estAngleMean = list(angle = FALSE))
#   # 
#   # # fit model
#   # turk_mA1.zm <- fitHMM(data = turkeyData.zm,
#   #                       retryFits = retryFits,
#   #                       nbStates = nbStates,
#   #                       dist = dist,
#   #                       Par0 = Par0_mA1.zm,
#   #                       DM = DM,
#   #                       workBounds = workBounds,
#   #                       userBounds = userBounds,
#   #                       estAngleMean = list(angle=FALSE),
#   #                       prior = prior,
#   #                       stateNames = stateNames,
#   #                       ncores = 5
#   # )
#   # 
#   # ##############################################################
#   # ### MODEL A2 - 60% quantile for loafing step length
#   # #Userbound constraint on step
#   # stepBounds <- matrix(c(0,5,
#   #                        0,76.35757,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        0, Inf,
#   #                        .98,1,
#   #                        0,.02,
#   #                        0,.02), nrow = 3*nbStates, byrow = T,
#   #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
#   # 
#   # #Bundle individual parameter DM and workbounds
#   # DM<-list(step=stepDM,angle=angleDM)
#   # workBounds<-list(step=stepworkBounds,
#   #                  angle=angleworkBounds)
#   # userBounds <- list(step = stepBounds,
#   #                    angle = angleBounds)
#   # 
#   # # initial parameters
#   # Par <- list(step=c(2, 40,150, #mean
#   #                    2, 10,25, #sd
#   #                    .99, .005, .005), #zero mass
#   #             angle = c(.01, 0.2 ,0.3))
#   # 
#   # Par0_mA2.zm <- getParDM(data = turkeyData.zm,
#   #                         nbStates = nbStates,
#   #                         zeroInflation = list(step = T,
#   #                                              angle = T),
#   #                         dist = dist,
#   #                         Par = Par,
#   #                         DM = DM,
#   #                         workBounds = workBounds,
#   #                         userBounds = userBounds,
#   #                         estAngleMean = list(angle = FALSE))
#   # 
#   # # fit model
#   # turk_mA2.zm <- fitHMM(data = turkeyData.zm,
#   #                       retryFits = retryFits,
#   #                       nbStates = nbStates,
#   #                       dist = dist,
#   #                       Par0 = Par0_mA2.zm,
#   #                       DM = DM,
#   #                       workBounds = workBounds,
#   #                       userBounds = userBounds,
#   #                       estAngleMean = list(angle=FALSE),
#   #                       prior = prior,
#   #                       stateNames = stateNames,
#   #                       ncores = 5
#   # )
# } #A0-2
# {
# # ##############################################################
# # ### MODEL A4 - 80% quantile for loafing step length
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,5,
# #                        0,160.15313,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.02,
# #                        0,.02), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 100,200, #mean
# #                    2, 10,25, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_mA4.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         zeroInflation = list(step = T,
# #                                              angle = T),
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE))
# # 
# # # fit model
# # turk_mA4.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_mA4.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5
# # )
# # 
# # ##############################################################
# # ### MODEL A5 - 90% quantile for loafing step length
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,5,
# #                        0,254.72528,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.02,
# #                        0,.02), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 150,300, #mean
# #                    2, 10,25, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_mA5.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         zeroInflation = list(step = T,
# #                                              angle = T),
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE))
# # 
# # # fit model
# # turk_mA5.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_mA5.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5
# # )
# } #A4-5
# 
# ### Top Model = MODEL A3 ###
# 
# 
# 
# 
# ########################################################################
# ################################## B ###################################
# ########################################################################
# stepDM<-matrix(c(
#   "WC.Z + SD.Z",0,0,0,0,0,0,0,0,
#   "WC.Z + SD.Z",1,0,0,0,0,0,0,0,
#   "WC.Z + SD.Z",1,1,0,0,0,0,0,0,
#   0,0,0,"WC.Z + SD.Z",0,0,0,0,0,
#   0,0,0,"WC.Z + SD.Z",1,0,0,0,0,
#   0,0,0,"WC.Z + SD.Z",1,1,0,0,0,
#   0,0,0,0,0,0,1,0,0,
#   0,0,0,0,0,0,0,1,0,
#   0,0,0,0,0,0,0,0,1),
#   nrow = 3*nbStates,byrow=TRUE,
#   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
#                 c("mean_123:(Intercept)", "mean_2","mean_3",
#                   paste0("sd_",1:nbStates,":(Intercept)"),
#                   paste0("zero_",1:nbStates,":(Intercept)"))))
# #define the directions of the differences
# stepworkBounds <- matrix(c(-Inf, 0,0,-Inf, 0,0,rep(-Inf,3),rep(Inf,9)),nrow = 3*nbStates,
#                          dimnames=list(colnames(stepDM),c("lower","upper")))
# 
# stepBounds <- matrix(c(0,5,
#                        0,109.69861,
#                        0, Inf,
#                        0, Inf,
#                        0, Inf,
#                        0, Inf,
#                        .98,1,
#                        0,.02,
#                        0,.02), nrow = 3*nbStates, byrow = T,
#                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# 
# ## constrain turning angle concentration parameters:
# # Concentration -> searching < dispersal
# angleDM<-matrix(c(1,0,0,
#                   1,1,0,
#                   1,1,1),nrow = nbStates,byrow=TRUE,
#                 dimnames=list(paste0("concentration_",1:nbStates),
#                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# 
# #define direction of differences
# angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
#                           nrow = ncol(angleDM),
#                           dimnames=list(colnames(angleDM),c("lower","upper")))
# 
# #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# #Userbound contraint on angle
# angleBounds <- matrix(c(0,0.94,
#                         0,0.94,
#                         0,0.94),nrow = nbStates,
#                       byrow=TRUE, dimnames=list(rownames(angleDM),
#                                                 c("lower","upper")))
# 
# #Bundle individual parameter DM and workbounds
# DM<-list(step=stepDM,angle=angleDM)
# workBounds<-list(step=stepworkBounds,
#                  angle=angleworkBounds)
# userBounds <- list(step = stepBounds,
#                    angle = angleBounds)
# 
# 
# Par_mA3.zm <- getPar0(model = turk_mA3.zm,
#                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID, 
#                       DM = DM,
#                       workBounds = workBounds,
#                       userBounds = userBounds)
# 
# 
# # fit model
# turk_mB0.zm <- fitHMM(data = turkeyData.zm,
#                       retryFits = retryFits,
#                       nbStates = nbStates,
#                       dist = dist,
#                       Par0 = Par0_mB0.zm$Par,
#                       beta0 = Par0_mB0.zm$beta,
#                       DM = DM,
#                       workBounds = workBounds,
#                       userBounds = userBounds,
#                       estAngleMean = list(angle=FALSE),
#                       prior = prior,
#                       stateNames = stateNames,
#                       ncores = 5,
#                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID
# )
# 
# 
# }
# {
# 
# ### ##Test if including parameter bounds improves model
# # # 13 - DM and Parameter Bounds, No other specifications
# # # 14 - DM and Parameter Bounds, Hour of Day effect on transition probability
# # # 15 - DM and Parameter Bounds, Wind Chill effect on transition probability
# # # 16 - DM and Parameter Bounds, Snow Depth effect on transition probability
# # # 17 - DM and Parameter Bounds, Hour of Day*Wind Chill effect on transition probability
# # # 18 - DM and Parameter Bounds, Hour of Day*Snow Depth effect on transition probability
# # # 19 - DM and Parameter Bounds, Hour of Day+Snow Depth+Wind Chill effect on transition probability
# # 
# # 
# # ### HMM Options Common Between Models ###
# # retryFits <- 5 # number attempt to re-fit based on random perturbation
# # nbStates <- 3 # Number of states
# # stateNames <- c("roost", "loafing", "foraging") # label states
# # dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes
# # prior <- function(par){sum(dnorm(par,0,10,log=TRUE))} #prevents working parameters from straying along boundary
# # 
# # 
# # ### NEED TO FIND THE BEST BOUNDS/DM FOR BASE MODEL
# # # ##############################################################
# # # ### MODEL 0A 
# # # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # # ## Step Bound - Zero Mass for roost >.98
# # # ## Angle DM - Roost < Loafing < Foraging
# # # ## Angle Bounds - All 0< & <.94
# # # stepDM<-matrix(c(
# # #   1,0,0,0,0,0,0,0,0,
# # #   1,1,0,0,0,0,0,0,0,
# # #   1,1,1,0,0,0,0,0,0,
# # #   0,0,0,1,0,0,0,0,0,
# # #   0,0,0,0,1,0,0,0,0,
# # #   0,0,0,0,0,1,0,0,0,
# # #   0,0,0,0,0,0,1,0,0,
# # #   0,0,0,0,0,0,0,1,0,
# # #   0,0,0,0,0,0,0,0,1),
# # #   nrow = 3*nbStates,byrow=TRUE,
# # #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# # #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# # #                   paste0("sd_",1:nbStates,":(Intercept)"),
# # #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # # 
# # # #define the directions of the differences
# # # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# # #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # # #Userbound constraint on step
# # # stepBounds <- matrix(c(0,Inf,
# # #                        0,Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        .98,1,
# # #                        0,.01,
# # #                        0,.01), nrow = 3*nbStates, byrow = T,
# # #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # # 
# # # ## constrain turning angle concentration parameters:
# # # # Concentration -> searching < dispersal
# # # angleDM<-matrix(c(1,0,0,
# # #                   1,1,0,
# # #                   1,1,1),nrow = nbStates,byrow=TRUE,
# # #                 dimnames=list(paste0("concentration_",1:nbStates),
# # #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # # 
# # # #define direction of differences
# # # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# # #                           nrow = ncol(angleDM),
# # #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # # 
# # # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # # #Userbound contraint on angle
# # # angleBounds <- matrix(c(0,0.94,
# # #                         0,0.94,
# # #                         0,0.94),nrow = nbStates,
# # #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# # #                                                 c("lower","upper")))
# # # 
# # # #Bundle individual parameter DM and workbounds
# # # DM<-list(step=stepDM,angle=angleDM)
# # # workBounds<-list(step=stepworkBounds,
# # #                  angle=angleworkBounds)
# # # userBounds <- list(step = stepBounds,
# # #                    angle = angleBounds)
# # # 
# # # # initial parameters
# # # Par <- list(step=c(2, 50,120, #mean
# # #                    2, 15,15, #sd
# # #                    .99, .005, .005), #zero mass
# # #             angle = c(.01, 0.2 ,0.3))
# # # 
# # # Par0_m0A.zm <- getParDM(data = turkeyData.zm,
# # #                         nbStates = nbStates,
# # #                         dist = dist,
# # #                         Par = Par,
# # #                         DM = DM,
# # #                         workBounds = workBounds,
# # #                         userBounds = userBounds,
# # #                         estAngleMean = list(angle = FALSE))
# # # 
# # # # fit model
# # # turk_m0A.zm <- fitHMM(data = turkeyData.zm,
# # #                       retryFits = retryFits,
# # #                       nbStates = nbStates,
# # #                       dist = dist,
# # #                       Par0 = Par0_m0A.zm,
# # #                       DM = DM,
# # #                       workBounds = workBounds,
# # #                       userBounds = userBounds,
# # #                       estAngleMean = list(angle=FALSE),
# # #                       prior = prior,
# # #                       stateNames = stateNames,
# # #                       ncores = 5
# # # )
# # # 
# # # ##############################################################
# # # ### MODEL 0B 
# # # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # # ## Angle DM - Roost < Loafing < Foraging
# # # ## Angle Bounds - All 0< & <.94
# # # stepDM<-matrix(c(
# # #   1,0,0,0,0,0,0,0,0,
# # #   1,1,0,0,0,0,0,0,0,
# # #   1,1,1,0,0,0,0,0,0,
# # #   0,0,0,1,0,0,0,0,0,
# # #   0,0,0,0,1,0,0,0,0,
# # #   0,0,0,0,0,1,0,0,0,
# # #   0,0,0,0,0,0,1,0,0,
# # #   0,0,0,0,0,0,0,1,0,
# # #   0,0,0,0,0,0,0,0,1),
# # #   nrow = 3*nbStates,byrow=TRUE,
# # #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# # #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# # #                   paste0("sd_",1:nbStates,":(Intercept)"),
# # #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # # 
# # # #define the directions of the differences
# # # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# # #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # # #Userbound constraint on step
# # # stepBounds <- matrix(c(0,10,
# # #                        0,Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        .98,1,
# # #                        0,.02,
# # #                        0,.02), nrow = 3*nbStates, byrow = T,
# # #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # # 
# # # ## constrain turning angle concentration parameters:
# # # # Concentration -> searching < dispersal
# # # angleDM<-matrix(c(1,0,0,
# # #                   1,1,0,
# # #                   1,1,1),nrow = nbStates,byrow=TRUE,
# # #                 dimnames=list(paste0("concentration_",1:nbStates),
# # #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # # 
# # # #define direction of differences
# # # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# # #                           nrow = ncol(angleDM),
# # #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # # 
# # # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # # #Userbound contraint on angle
# # # angleBounds <- matrix(c(0,0.94,
# # #                         0,0.94,
# # #                         0,0.94),nrow = nbStates,
# # #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# # #                                                 c("lower","upper")))
# # # 
# # # #Bundle individual parameter DM and workbounds
# # # DM<-list(step=stepDM,angle=angleDM)
# # # workBounds<-list(step=stepworkBounds,
# # #                  angle=angleworkBounds)
# # # userBounds <- list(step = stepBounds,
# # #                    angle = angleBounds)
# # # 
# # # # initial parameters
# # # Par <- list(step=c(2, 50,120, #mean
# # #                    2, 15,15, #sd
# # #                    .99, .005, .005), #zero mass
# # #             angle = c(.01, 0.2 ,0.3))
# # # 
# # # Par0_m0B.zm <- getParDM(data = turkeyData.zm,
# # #                         nbStates = nbStates,
# # #                         dist = dist,
# # #                         Par = Par,
# # #                         DM = DM,
# # #                         workBounds = workBounds,
# # #                         userBounds = userBounds,
# # #                         estAngleMean = list(angle = FALSE))
# # # 
# # # # fit model
# # # turk_m0B.zm <- fitHMM(data = turkeyData.zm,
# # #                       retryFits = retryFits,
# # #                       nbStates = nbStates,
# # #                       dist = dist,
# # #                       Par0 = Par0_m0B.zm,
# # #                       DM = DM,
# # #                       workBounds = workBounds,
# # #                       userBounds = userBounds,
# # #                       estAngleMean = list(angle=FALSE),
# # #                       prior = prior,
# # #                       stateNames = stateNames,
# # #                       ncores = 5
# # # )
# # # 
# # # ##############################################################
# # # ### MODEL 0C - Different Initial parameters 
# # # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # # ## Angle DM - Roost < Loafing < Foraging
# # # ## Angle Bounds - All 0< & <.94
# # # stepDM<-matrix(c(
# # #   1,0,0,0,0,0,0,0,0,
# # #   1,1,0,0,0,0,0,0,0,
# # #   1,1,1,0,0,0,0,0,0,
# # #   0,0,0,1,0,0,0,0,0,
# # #   0,0,0,0,1,0,0,0,0,
# # #   0,0,0,0,0,1,0,0,0,
# # #   0,0,0,0,0,0,1,0,0,
# # #   0,0,0,0,0,0,0,1,0,
# # #   0,0,0,0,0,0,0,0,1),
# # #   nrow = 3*nbStates,byrow=TRUE,
# # #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# # #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# # #                   paste0("sd_",1:nbStates,":(Intercept)"),
# # #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # # 
# # # #define the directions of the differences
# # # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# # #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # # #Userbound constraint on step
# # # stepBounds <- matrix(c(0,10,
# # #                        0,Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        .98,1,
# # #                        0,.02,
# # #                        0,.02), nrow = 3*nbStates, byrow = T,
# # #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # # 
# # # ## constrain turning angle concentration parameters:
# # # # Concentration -> searching < dispersal
# # # angleDM<-matrix(c(1,0,0,
# # #                   1,1,0,
# # #                   1,1,1),nrow = nbStates,byrow=TRUE,
# # #                 dimnames=list(paste0("concentration_",1:nbStates),
# # #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # # 
# # # #define direction of differences
# # # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# # #                           nrow = ncol(angleDM),
# # #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # # 
# # # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # # #Userbound contraint on angle
# # # angleBounds <- matrix(c(0,0.94,
# # #                         0,0.94,
# # #                         0,0.94),nrow = nbStates,
# # #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# # #                                                 c("lower","upper")))
# # # 
# # # #Bundle individual parameter DM and workbounds
# # # DM<-list(step=stepDM,angle=angleDM)
# # # workBounds<-list(step=stepworkBounds,
# # #                  angle=angleworkBounds)
# # # userBounds <- list(step = stepBounds,
# # #                    angle = angleBounds)
# # # 
# # # # initial parameters
# # # Par <- list(step=c(6, 50,120, #mean
# # #                    6, 15,15, #sd
# # #                    .99, .005, .005), #zero mass
# # #             angle = c(.01, 0.2 ,0.3))
# # # 
# # # Par0_m0C.zm <- getParDM(data = turkeyData.zm,
# # #                         nbStates = nbStates,
# # #                         dist = dist,
# # #                         Par = Par,
# # #                         DM = DM,
# # #                         workBounds = workBounds,
# # #                         userBounds = userBounds,
# # #                         estAngleMean = list(angle = FALSE))
# # # 
# # # # fit model
# # # turk_m0C.zm <- fitHMM(data = turkeyData.zm,
# # #                       retryFits = retryFits,
# # #                       nbStates = nbStates,
# # #                       dist = dist,
# # #                       Par0 = Par0_m0C.zm,
# # #                       DM = DM,
# # #                       workBounds = workBounds,
# # #                       userBounds = userBounds,
# # #                       estAngleMean = list(angle=FALSE),
# # #                       prior = prior,
# # #                       stateNames = stateNames,
# # #                       ncores = 5
# # # )
# # # 
# # # ##############################################################
# # # ### MODEL 0D 
# # # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # # ## Angle DM - Roost < Loafing < Foraging
# # # ## Angle Bounds - All 0< & <.94
# # # stepDM<-matrix(c(
# # #   1,0,0,0,0,0,0,0,0,
# # #   1,1,0,0,0,0,0,0,0,
# # #   1,1,1,0,0,0,0,0,0,
# # #   0,0,0,1,0,0,0,0,0,
# # #   0,0,0,0,1,0,0,0,0,
# # #   0,0,0,0,0,1,0,0,0,
# # #   0,0,0,0,0,0,1,0,0,
# # #   0,0,0,0,0,0,0,1,0,
# # #   0,0,0,0,0,0,0,0,1),
# # #   nrow = 3*nbStates,byrow=TRUE,
# # #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# # #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# # #                   paste0("sd_",1:nbStates,":(Intercept)"),
# # #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # # 
# # # #define the directions of the differences
# # # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# # #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # # #Userbound constraint on step
# # # stepBounds <- matrix(c(0,15,
# # #                        0,Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        .98,1,
# # #                        0,.02,
# # #                        0,.02), nrow = 3*nbStates, byrow = T,
# # #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # # 
# # # ## constrain turning angle concentration parameters:
# # # # Concentration -> searching < dispersal
# # # angleDM<-matrix(c(1,0,0,
# # #                   1,1,0,
# # #                   1,1,1),nrow = nbStates,byrow=TRUE,
# # #                 dimnames=list(paste0("concentration_",1:nbStates),
# # #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # # 
# # # #define direction of differences
# # # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# # #                           nrow = ncol(angleDM),
# # #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # # 
# # # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # # #Userbound contraint on angle
# # # angleBounds <- matrix(c(0,0.94,
# # #                         0,0.94,
# # #                         0,0.94),nrow = nbStates,
# # #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# # #                                                 c("lower","upper")))
# # # 
# # # #Bundle individual parameter DM and workbounds
# # # DM<-list(step=stepDM,angle=angleDM)
# # # workBounds<-list(step=stepworkBounds,
# # #                  angle=angleworkBounds)
# # # userBounds <- list(step = stepBounds,
# # #                    angle = angleBounds)
# # # 
# # # # initial parameters
# # # Par <- list(step=c(2, 50,120, #mean
# # #                    2, 15,15, #sd
# # #                    .99, .005, .005), #zero mass
# # #             angle = c(.01, 0.2 ,0.3))
# # # 
# # # Par0_m0D.zm <- getParDM(data = turkeyData.zm,
# # #                         nbStates = nbStates,
# # #                         dist = dist,
# # #                         Par = Par,
# # #                         DM = DM,
# # #                         workBounds = workBounds,
# # #                         userBounds = userBounds,
# # #                         estAngleMean = list(angle = FALSE))
# # # 
# # # # fit model
# # # turk_m0D.zm <- fitHMM(data = turkeyData.zm,
# # #                       retryFits = retryFits,
# # #                       nbStates = nbStates,
# # #                       dist = dist,
# # #                       Par0 = Par0_m0D.zm,
# # #                       DM = DM,
# # #                       workBounds = workBounds,
# # #                       userBounds = userBounds,
# # #                       estAngleMean = list(angle=FALSE),
# # #                       prior = prior,
# # #                       stateNames = stateNames,
# # #                       ncores = 5
# # # )
# # # 
# # # ##############################################################
# # # ### MODEL 0E - Different Initial parameters 
# # # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # # ## Angle DM - Roost < Loafing < Foraging
# # # ## Angle Bounds - All 0< & <.94
# # # stepDM<-matrix(c(
# # #   1,0,0,0,0,0,0,0,0,
# # #   1,1,0,0,0,0,0,0,0,
# # #   1,1,1,0,0,0,0,0,0,
# # #   0,0,0,1,0,0,0,0,0,
# # #   0,0,0,0,1,0,0,0,0,
# # #   0,0,0,0,0,1,0,0,0,
# # #   0,0,0,0,0,0,1,0,0,
# # #   0,0,0,0,0,0,0,1,0,
# # #   0,0,0,0,0,0,0,0,1),
# # #   nrow = 3*nbStates,byrow=TRUE,
# # #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# # #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# # #                   paste0("sd_",1:nbStates,":(Intercept)"),
# # #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # # 
# # # #define the directions of the differences
# # # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# # #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # # #Userbound constraint on step
# # # stepBounds <- matrix(c(0,15,
# # #                        0,Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        .98,1,
# # #                        0,.02,
# # #                        0,.02), nrow = 3*nbStates, byrow = T,
# # #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # # 
# # # ## constrain turning angle concentration parameters:
# # # # Concentration -> searching < dispersal
# # # angleDM<-matrix(c(1,0,0,
# # #                   1,1,0,
# # #                   1,1,1),nrow = nbStates,byrow=TRUE,
# # #                 dimnames=list(paste0("concentration_",1:nbStates),
# # #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # # 
# # # #define direction of differences
# # # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# # #                           nrow = ncol(angleDM),
# # #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # # 
# # # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # # #Userbound contraint on angle
# # # angleBounds <- matrix(c(0,0.94,
# # #                         0,0.94,
# # #                         0,0.94),nrow = nbStates,
# # #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# # #                                                 c("lower","upper")))
# # # 
# # # #Bundle individual parameter DM and workbounds
# # # DM<-list(step=stepDM,angle=angleDM)
# # # workBounds<-list(step=stepworkBounds,
# # #                  angle=angleworkBounds)
# # # userBounds <- list(step = stepBounds,
# # #                    angle = angleBounds)
# # # 
# # # # initial parameters
# # # Par <- list(step=c(6, 50,120, #mean
# # #                    6, 15,15, #sd
# # #                    .99, .005, .005), #zero mass
# # #             angle = c(.01, 0.2 ,0.3))
# # # 
# # # Par0_m0E.zm <- getParDM(data = turkeyData.zm,
# # #                         nbStates = nbStates,
# # #                         dist = dist,
# # #                         Par = Par,
# # #                         DM = DM,
# # #                         workBounds = workBounds,
# # #                         userBounds = userBounds,
# # #                         estAngleMean = list(angle = FALSE))
# # # 
# # # # fit model
# # # turk_m0E.zm <- fitHMM(data = turkeyData.zm,
# # #                       retryFits = retryFits,
# # #                       nbStates = nbStates,
# # #                       dist = dist,
# # #                       Par0 = Par0_m0E.zm,
# # #                       DM = DM,
# # #                       workBounds = workBounds,
# # #                       userBounds = userBounds,
# # #                       estAngleMean = list(angle=FALSE),
# # #                       prior = prior,
# # #                       stateNames = stateNames,
# # #                       ncores = 5
# # # )
# # # 
# # # ##############################################################
# # # ### MODEL 0F 
# # # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # # ## Angle DM - Roost < Loafing < Foraging
# # # ## Angle Bounds - All 0< & <.94
# # # stepDM<-matrix(c(
# # #   1,0,0,0,0,0,0,0,0,
# # #   1,1,0,0,0,0,0,0,0,
# # #   1,1,1,0,0,0,0,0,0,
# # #   0,0,0,1,0,0,0,0,0,
# # #   0,0,0,0,1,0,0,0,0,
# # #   0,0,0,0,0,1,0,0,0,
# # #   0,0,0,0,0,0,1,0,0,
# # #   0,0,0,0,0,0,0,1,0,
# # #   0,0,0,0,0,0,0,0,1),
# # #   nrow = 3*nbStates,byrow=TRUE,
# # #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# # #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# # #                   paste0("sd_",1:nbStates,":(Intercept)"),
# # #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # # 
# # # #define the directions of the differences
# # # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# # #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # # #Userbound constraint on step
# # # stepBounds <- matrix(c(0,5,
# # #                        0,Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        0, Inf,
# # #                        .98,1,
# # #                        0,.02,
# # #                        0,.02), nrow = 3*nbStates, byrow = T,
# # #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # # 
# # # ## constrain turning angle concentration parameters:
# # # # Concentration -> searching < dispersal
# # # angleDM<-matrix(c(1,0,0,
# # #                   1,1,0,
# # #                   1,1,1),nrow = nbStates,byrow=TRUE,
# # #                 dimnames=list(paste0("concentration_",1:nbStates),
# # #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # # 
# # # #define direction of differences
# # # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# # #                           nrow = ncol(angleDM),
# # #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # # 
# # # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # # #Userbound contraint on angle
# # # angleBounds <- matrix(c(0,0.94,
# # #                         0,0.94,
# # #                         0,0.94),nrow = nbStates,
# # #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# # #                                                 c("lower","upper")))
# # # 
# # # #Bundle individual parameter DM and workbounds
# # # DM<-list(step=stepDM,angle=angleDM)
# # # workBounds<-list(step=stepworkBounds,
# # #                  angle=angleworkBounds)
# # # userBounds <- list(step = stepBounds,
# # #                    angle = angleBounds)
# # # 
# # # # initial parameters
# # # Par <- list(step=c(2, 50,120, #mean
# # #                    2, 15,15, #sd
# # #                    .99, .005, .005), #zero mass
# # #             angle = c(.01, 0.2 ,0.3))
# # # 
# # # Par0_m0F.zm <- getParDM(data = turkeyData.zm,
# # #                         nbStates = nbStates,
# # #                         dist = dist,
# # #                         Par = Par,
# # #                         DM = DM,
# # #                         workBounds = workBounds,
# # #                         userBounds = userBounds,
# # #                         estAngleMean = list(angle = FALSE))
# # # 
# # # # fit model
# # # turk_m0F.zm <- fitHMM(data = turkeyData.zm,
# # #                       retryFits = retryFits,
# # #                       nbStates = nbStates,
# # #                       dist = dist,
# # #                       Par0 = Par0_m0F.zm,
# # #                       DM = DM,
# # #                       workBounds = workBounds,
# # #                       userBounds = userBounds,
# # #                       estAngleMean = list(angle=FALSE),
# # #                       prior = prior,
# # #                       stateNames = stateNames,
# # #                       ncores = 5
# # # )
# # 
# # ##############################################################
# # ### MODEL 0G - Different Initial parameters 
# # ## Step DM - only Step Mean, Roost < Loafing < Foraging
# # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # ## Angle DM - Roost < Loafing < Foraging
# # ## Angle Bounds - All 0< & <.94
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,0,1,0,0,0,0,
# #   0,0,0,0,0,1,0,0,0,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0,
# #   0,0,0,0,0,0,0,0,1),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,5,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.02,
# #                        0,.02), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(4, 50,120, #mean
# #                    6, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m0G.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE))
# # 
# # # fit model
# # turk_m0G.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m0G.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5
# # )
# # 
# # ##############################################################
# # ### MODEL 1A - Different Initial parameters 
# # ## Step DM - Step Mean and SD, Roost < Loafing < Foraging
# # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # ## Angle DM - Roost < Loafing < Foraging
# # ## Angle Bounds - All 0< & <.94
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0,
# #   0,0,0,0,0,0,0,0,1),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,-Inf, 0,0,rep(-Inf,3),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,5,
# #                        6,85,
# #                        86, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.02,
# #                        0,.02), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 45,200, #mean
# #                    2, 10,25, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m1A.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         zeroInflation = list(step = T,
# #                                              angle = T),
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE))
# # 
# # # fit model
# # turk_m1A.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m1A.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5
# # )
# ###NO PAREMETER RESTRICTIONS SPECIFIED###
# ## All models will have 3 states = Roosting, Loafing, Foraging
# # 0 - Null Model = no specifications
# # 1 - No Parameter Bounds, Hour of day effect on transition probability
# # 2 - No Parameter Bounds, Hour of day effect on transition probability and step length
# # 3 - No Parameter Bounds, Wind Chill effect on transition probability
# # 4 - No Parameter Bounds, Wind Chill effect on transition probability and step length
# # 5 - No Parameter Bounds, Snow Depth effect on transition probability
# # 6 - No Parameter Bounds, Snow Depth effect on transition probability and step length
# # 7 - No Parameter Bounds, Hour of day*Wind Chill effect on transition probability
# # 8 - No Parameter Bounds, Hour of day*Wind Chill effect on transition probability and step length
# # 9 - No Parameter Bounds, Hour of day*Snow Depth effect on transition probability
# # 10 - No Parameter Bounds, Hour of day*Snow Depth effect on transition probability and step length
# # 11 - No Parameter Bounds, Hour of day + Snow Depth + Wind Chill effect on transition probability
# # 12 - No Parameter Bounds, Hour of day + Snow Depth + Wind Chill effect on transition probability and step length
# 
# # ##############################################################
# # ### MODEL 0 - Nothing Specified
# # # initial parameters
# # Par0_m0.zm <- list(step=c(2, 50,120, #mean
# #                           2, 50,100, #sd
# #                           .99, .005, .005), #zero mass
# #                    angle = c(.01, 0.2 ,0.3))
# # 
# # # fit model
# # turk_m0.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m0.zm,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5)
# # 
# # 
# # ##############################################################
# # ### MODEL 1 - No Parameter Bounds, Hour of day effect on transition probability
# # # initial parameters
# # Par0_m1.zm <- getPar0(model=turk_m0.zm, 
# #                       formula=~cosinor(hour, period = 24))
# # 
# # # fit model
# # turk_m1.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m1.zm$Par,
# #                      beta0 = Par0_m1.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5,
# #                      formula = ~cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 2 - No Parameter Bounds, Hour of day effect on transition probability and step length
# # # formulas for parameters of state-dependent observation distributions
# # DM <- list(step = list(mean = ~ cosinor(hour, period = 24),
# #                        sd = ~ cosinor(hour, period = 24),
# #                        zeromass = ~ 1))
# # 
# # # initial parameters
# # Par0_m2.zm <- getPar0(model=turk_m1.zm,
# #                       formula=~cosinor(hour, period = 24),
# #                       DM = DM)
# # 
# # # fit model
# # turk_m2.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m2.zm$Par,
# #                      beta0 = Par0_m2.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      DM = DM,
# #                      ncores = 5,
# #                      formula = ~cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 3 - No Parameter Bounds, Wind Chill effect on transition probability
# # # initial parameters
# # Par0_m3.zm <- getPar0(model=turk_m0.zm, 
# #                       formula= ~ WC.Z,
# #                       covNames=c("WC.Z"))
# # 
# # # fit model
# # turk_m3.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m3.zm$Par,
# #                      beta0 = Par0_m3.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5,
# #                      formula = ~ WC.Z
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 4 - No Parameter Bounds, Wind Chill effect on transition probability and step length
# # # formulas for parameters of state-dependent observation distributions
# # DM <- list(step = list(mean = ~ WC.Z,
# #                        sd = ~ WC.Z,
# #                        zeromass = ~ 1))
# # 
# # # initial parameters
# # Par0_m4.zm <- getPar0(model = turk_m3.zm,
# #                       formula = ~ WC.Z,
# #                       DM = DM,
# #                       covNames=c("WC.Z"))
# # 
# # # fit model
# # turk_m4.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m4.zm$Par,
# #                      beta0 = Par0_m4.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      DM = DM,
# #                      ncores = 5,
# #                      formula = ~WC.Z
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 5 - No Parameter Bounds, Snow Depth effect on transition probability
# # # initial parameters
# # Par0_m5.zm <- getPar0(model=turk_m0.zm, 
# #                       formula= ~ SD.Z,
# #                       covNames=c("SD.Z"))
# # 
# # # fit model
# # turk_m5.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m5.zm$Par,
# #                      beta0 = Par0_m5.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5,
# #                      formula = ~ SD.Z
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 6 - No Parameter Bounds, Snow Depth effect on transition probability and step length
# # # formulas for parameters of state-dependent observation distributions
# # DM <- list(step = list(mean = ~ SD.Z,
# #                        sd = ~ SD.Z,
# #                        zeromass = ~ 1))
# # 
# # # initial parameters
# # Par0_m6.zm <- getPar0(model = turk_m5.zm,
# #                       formula = ~ SD.Z,
# #                       DM = DM,
# #                       covNames=c("SD.Z"))
# # 
# # # fit model
# # turk_m6.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m6.zm$Par,
# #                      beta0 = Par0_m6.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      DM = DM,
# #                      ncores = 5,
# #                      formula = ~SD.Z
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 7 - No Parameter Bounds, Hour of Day * Wind Chill effect on transition probability
# # # initial parameters
# # Par0_m7.zm <- getPar0(model=turk_m3.zm, 
# #                       formula= ~ WC.Z*cosinor(hour, period = 24),
# #                       covNames=c("WC.Z"))
# # 
# # # fit model
# # turk_m7.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m7.zm$Par,
# #                      beta0 = Par0_m7.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5,
# #                      formula = ~ WC.Z*cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 8 - No Parameter Bounds, Hour of Day * Wind Chill effect on transition probability and step length
# # # formulas for parameters of state-dependent observation distributions
# # DM <- list(step = list(mean = ~ WC.Z*cosinor(hour, period = 24),
# #                        sd = ~ WC.Z*cosinor(hour, period = 24),
# #                        zeromass = ~ 1))
# # 
# # # initial parameters
# # Par0_m8.zm <- getPar0(model = turk_m7.zm,
# #                       formula = ~ WC.Z*cosinor(hour, period = 24),
# #                       DM = DM,
# #                       covNames=c("WC.Z"))
# # 
# # # fit model
# # turk_m8.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m8.zm$Par,
# #                      beta0 = Par0_m8.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      DM = DM,
# #                      ncores = 5,
# #                      formula = ~WC.Z*cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 9 - No Parameter Bounds, Hour of Day * Snow Depth effect on transition probability
# # # initial parameters
# # Par0_m9.zm <- getPar0(model=turk_m1.zm, 
# #                       formula= ~ SD.Z*cosinor(hour, period = 24),
# #                       covNames=c("SD.Z"))
# # 
# # # fit model
# # turk_m9.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m9.zm$Par,
# #                      beta0 = Par0_m9.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5,
# #                      formula = ~ SD.Z*cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 10 - No Parameter Bounds, Hour of day * Snow Depth effect on transition probability and step length
# # # formulas for parameters of state-dependent observation distributions
# # DM <- list(step = list(mean = ~ SD.Z*cosinor(hour, period = 24),
# #                        sd = ~ SD.Z*cosinor(hour, period = 24),
# #                        zeromass = ~ 1))
# # 
# # # initial parameters
# # Par0_m10.zm <- getPar0(model = turk_m9.zm,
# #                       formula = ~ SD.Z*cosinor(hour, period = 24),
# #                       DM = DM,
# #                       covNames=c("SD.Z"))
# # 
# # # fit model
# # turk_m10.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m10.zm$Par,
# #                      beta0 = Par0_m10.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      DM = DM,
# #                      ncores = 5,
# #                      formula = ~SD.Z*cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# # ### MODEL 11 - No Parameter Bounds, Hour of Day + Snow Depth + Wind Chill effect on transition probability
# # # initial parameters
# # Par0_m11.zm <- getPar0(model=turk_m1.zm, 
# #                       formula= ~ WC.Z + SD.Z + cosinor(hour, period = 24),
# #                       covNames=c("SD.Z", "WC.Z"))
# # 
# # # fit model
# # turk_m11.zm <- fitHMM(data = turkeyData.zm, 
# #                      retryFits = retryFits,
# #                      nbStates = nbStates, 
# #                      dist = dist, 
# #                      Par0 = Par0_m11.zm$Par,
# #                      beta0 = Par0_m11.zm$beta,
# #                      estAngleMean = list(angle=FALSE),
# #                      stateNames = stateNames,
# #                      ncores = 5,
# #                      formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24)
# #                      )
# # 
# # 
# # ##############################################################
# ### MODEL 12 - No Parameter Bounds, Hour of day + Snow Depth + Wind Chill effect on transition probability and step length
# # formulas for parameters of state-dependent observation distributions
# DM <- list(step = list(mean = ~ ID + cosinor(hour, period = 24),
#                        sd = ~ ID + cosinor(hour, period = 24),
#                        zeromass = ~ 1))
# 
# # initial parameters
# Par0_m12.zm <- getPar0(model = turk_m11.zm,
#                        formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24),
#                        DM = DM,
#                        covNames=c("SD.Z", "WC.Z"))
# 
# # fit model
# turk_m12.zm <- fitHMM(data = turkeyData.zm,
#                       retryFits = retryFits,
#                       nbStates = nbStates,
#                       dist = dist,
#                       Par0 = Par0_m12.zm$Par,
#                       beta0 = Par0_m12.zm$beta,
#                       estAngleMean = list(angle=FALSE),
#                       stateNames = stateNames,
#                       DM = DM,
#                       ncores = 5,
#                       formula = ~WC.Z + SD.Z + cosinor(hour, period = 24) + ID
#                       )
# 
# # ### Includes Weather and Time Covariates
# # ##############################################################
# # ### MODEL 13 - DM and Workbounds added; No covariates
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m13.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE))
# # 
# # # fit model
# # turk_m13.zm <- fitHMM(data = turkeyData.zm,
# #                      retryFits = retryFits,
# #                      nbStates = nbStates,
# #                      dist = dist,
# #                      Par0 = Par0_m13.zm,
# #                      DM = DM,
# #                      workBounds = workBounds,
# #                      userBounds = userBounds,
# #                      estAngleMean = list(angle=FALSE),
# #                      prior = prior,
# #                      stateNames = stateNames,
# #                      ncores = 5
# # )
# # 
# # 
# # ##############################################################
# # ### MODEL 14 - DM and Workbounds added; Hour of Day Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m14.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~cosinor(hour, period = 24)
# #                         )
# # 
# # # fit model
# # turk_m14.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m14.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~cosinor(hour, period = 24)
# # )
# # 
# # 
# # ##############################################################
# # ### MODEL 15 - DM and Workbounds added; Wind Chill Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m15.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ WC.Z,
# #                         covNames = c("WC.Z")
# # )
# # 
# # # fit model
# # turk_m15.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m15.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ WC.Z
# # )
# # 
# # 
# # ##############################################################
# # ### MODEL 16 - DM and Workbounds added; Snow Depth Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m16.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ SD.Z,
# #                         covNames = c("SD.Z")
# # )
# # 
# # # fit model
# # turk_m16.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m16.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ SD.Z
# # )
# # 
# # 
# # ##############################################################
# # ### MODEL 17 - DM and Workbounds added; Wind Chill*Hour Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m17.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ WC.Z*cosinor(hour, period = 24),
# #                         covNames = c("WC.Z")
# # )
# # 
# # # fit model
# # turk_m17.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m17.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ WC.Z*cosinor(hour, period = 24)
# # )
# # 
# # 
# # ##############################################################
# # ### MODEL 18 - DM and Workbounds added; SnowDepth*Hour Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m18.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ SD.Z*cosinor(hour, period = 24),
# #                         covNames = c("SD.Z")
# # )
# # 
# # # fit model
# # turk_m18.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m18.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ SD.Z*cosinor(hour, period = 24)
# # )
# # 
# # 
# # ##############################################################
# # ### MODEL 19 - DM and Workbounds added; Wind Chill + SnowDepth + Hour Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,1,1,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,10,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.01,
# #                        0,.01), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 50,120, #mean
# #                    2, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m19.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24),
# #                         covNames = c("SD.Z", "WC.Z")
# # )
# # 
# # # fit model
# # turk_m19.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m19.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24)
# # )
# # 
# # ##############################################################
# # ### MODEL 20 - DM and Workbounds added; ID + Wind Chill + SnowDepth + Hour Effect on Transition Probability
# # ## constrain step length parameters:
# # # Mean/SD -> foraging>loafing>roosting
# # # Zero Mass -> roost>loafing=foraging
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0,
# #   0,0,0,0,0,0,0,0,1),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,5,
# #                        0,Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.02,
# #                        0,.02), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(4, 50,120, #mean
# #                    6, 15,15, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m20.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
# #                         covNames = c("SD.Z", "WC.Z")
# # )
# # 
# # # fit model
# # turk_m20.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m20.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID
# # )
# # 
# # ##############################################################
# # ### MODEL 21 - Different Initial parameters 
# # ## Step DM - Step Mean and SD, Roost < Loafing < Foraging
# # ## Step Bound - Step Mean for Roost < 10; Zero Mass for roost >.98
# # ## Angle DM - Roost < Loafing < Foraging
# # ## Angle Bounds - All 0< & <.94
# # stepDM<-matrix(c(
# #   1,0,0,0,0,0,0,0,0,
# #   1,1,0,0,0,0,0,0,0,
# #   1,1,1,0,0,0,0,0,0,
# #   0,0,0,1,0,0,0,0,0,
# #   0,0,0,1,1,0,0,0,0,
# #   0,0,0,1,1,1,0,0,0,
# #   0,0,0,0,0,0,1,0,0,
# #   0,0,0,0,0,0,0,1,0,
# #   0,0,0,0,0,0,0,0,1),
# #   nrow = 3*nbStates,byrow=TRUE,
# #   dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
# #                 c("mean_123:(Intercept)", "mean_2","mean_3",
# #                   paste0("sd_",1:nbStates,":(Intercept)"),
# #                   paste0("zero_",1:nbStates,":(Intercept)"))))
# # 
# # #define the directions of the differences
# # stepworkBounds <- matrix(c(-Inf, 0,0,-Inf, 0,0,rep(-Inf,3),rep(Inf,9)),nrow = 3*nbStates,
# #                          dimnames=list(colnames(stepDM),c("lower","upper")))
# # #Userbound constraint on step
# # stepBounds <- matrix(c(0,5,
# #                        6,85,
# #                        86, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        0, Inf,
# #                        .98,1,
# #                        0,.02,
# #                        0,.02), nrow = 3*nbStates, byrow = T,
# #                      dimnames = list(rownames(stepDM), c("lower", "upper")))
# # 
# # ## constrain turning angle concentration parameters:
# # # Concentration -> searching < dispersal
# # angleDM<-matrix(c(1,0,0,
# #                   1,1,0,
# #                   1,1,1),nrow = nbStates,byrow=TRUE,
# #                 dimnames=list(paste0("concentration_",1:nbStates),
# #                               c("concentration_1:(Intercept)","concentration_2","concentration_3")))
# # 
# # #define direction of differences
# # angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
# #                           nrow = ncol(angleDM),
# #                           dimnames=list(colnames(angleDM),c("lower","upper")))
# # 
# # #Restrict angle concentration such that dispersal > .75 while dispersal>localized
# # #Userbound contraint on angle
# # angleBounds <- matrix(c(0,0.94,
# #                         0,0.94,
# #                         0,0.94),nrow = nbStates,
# #                       byrow=TRUE, dimnames=list(rownames(angleDM),
# #                                                 c("lower","upper")))
# # 
# # #Bundle individual parameter DM and workbounds
# # DM<-list(step=stepDM,angle=angleDM)
# # workBounds<-list(step=stepworkBounds,
# #                  angle=angleworkBounds)
# # userBounds <- list(step = stepBounds,
# #                    angle = angleBounds)
# # 
# # # initial parameters
# # Par <- list(step=c(2, 45,200, #mean
# #                    2, 10,25, #sd
# #                    .99, .005, .005), #zero mass
# #             angle = c(.01, 0.2 ,0.3))
# # 
# # Par0_m21.zm <- getParDM(data = turkeyData.zm,
# #                         nbStates = nbStates,
# #                         zeroInflation = list(step = T,
# #                                              angle = T),
# #                         dist = dist,
# #                         Par = Par,
# #                         DM = DM,
# #                         workBounds = workBounds,
# #                         userBounds = userBounds,
# #                         estAngleMean = list(angle = FALSE),
# #                         formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID,
# #                         covNames = c("SD.Z", "WC.Z"))
# # 
# # # fit model
# # turk_m21.zm <- fitHMM(data = turkeyData.zm,
# #                       retryFits = retryFits,
# #                       nbStates = nbStates,
# #                       dist = dist,
# #                       Par0 = Par0_m21.zm,
# #                       DM = DM,
# #                       workBounds = workBounds,
# #                       userBounds = userBounds,
# #                       estAngleMean = list(angle=FALSE),
# #                       prior = prior,
# #                       stateNames = stateNames,
# #                       ncores = 5,
# #                       formula = ~ WC.Z + SD.Z + cosinor(hour, period = 24) + ID
# # )
#   
# } #Previous Model Attempts