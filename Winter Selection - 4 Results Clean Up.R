require(dplyr)

#####################
### ROOST RESULTS ###
#####################
roostlist <- list()
roostlist[[1]] <- model.BA.WCprev.roost
roostlist[[2]] <- model.BA.SD.roost
roostlist[[3]] <- model.DtFE.WCprev.roost
roostlist[[4]] <- model.DtFE.SD.roost
roostlist[[5]] <- model.SW.WCprev.roost
roostlist[[6]] <- model.SW.SD.roost
roostlist[[7]] <- model.Ht.WCprev.roost
roostlist[[8]] <- model.Ht.SD.roost
roostlist[[9]] <- model.WindExp.WCprev.roost
roostlist[[10]] <- model.WindExp.SD.roost
  
rm(roostresults)
rm(roostmodelselection)
rm(info.cow)
for(i in 1:length(roostlist)){
  results.matrix <- summary(roostlist[[i]])$fixed
  indresult <- as.data.frame(results.matrix) %>%
    dplyr::select(mean, sd, "0.025quant", "0.975quant") %>%
    mutate(WeatherCov = rownames(results.matrix)[2]) %>%
    mutate(HabitatCov = rownames(results.matrix)[1]) %>%
    mutate(CovName = rownames(results.matrix))
  
  info.cow.ind <- data.frame(Beh_State = "Roosting",
                                LC_Cov = indresult[1,6],
                                LC_Coef = indresult[1,1],
                                Weath_Cov = indresult[1,5],
                                Weath_Coef = indresult[2,1],
                                Interaction = indresult[4,1])
  
  ind.ms <- data.frame(Beh_State = "Roosting",
                       LC_Cov = indresult[1,6],
                       Weath_Cov = indresult[1,5],
                       WAIC = roostlist[[i]]$waic$waic,
                       DIC = roostlist[[i]]$dic$dic,
                       MLik = roostlist[[i]]$mlik[2,1]
                       )
  
  if(exists("roostresults")){
    roostresults <- rbind(roostresults, indresult)
    info.cow <- rbind(info.cow, info.cow.ind)
    roostmodelselection <-  rbind(roostmodelselection, ind.ms)
  }else{
    roostresults <- indresult
    info.cow <- info.cow.ind
    roostmodelselection <- ind.ms
  }
}

write.csv(roostresults, "Results/roostresults.full.csv", row.names = F)
write.csv(roostmodelselection %>% arrange(Weath_Cov, DIC), "Results/roostmodelselection.csv", row.names = F)
roostinteractions <- roostresults %>%
  filter(grepl(":", CovName, fixed = T)) %>%
  mutate(Analysis = "Roosting")

#######################
### LOAFING RESULTS ###
#######################
loaflist <- list()
loaflist[[1]] <- model.BA.WCprev.Stationary
loaflist[[2]] <- model.BA.SD.Stationary
loaflist[[3]] <- model.DtFE.WCprev.Stationary
loaflist[[4]] <- model.DtFE.SD.Stationary
loaflist[[5]] <- model.Ht.WCprev.Stationary
loaflist[[6]] <- model.Ht.SD.Stationary
loaflist[[7]] <- model.WindExp.WCprev.Stationary
loaflist[[8]] <- model.WindExp.SD.Stationary
loaflist[[9]] <- model.PropAg.WCprev.Stationary
loaflist[[10]] <- model.PropAg.SD.Stationary
loaflist[[11]] <- model.PropDev.WCprev.Stationary
loaflist[[12]] <- model.PropDev.SD.Stationary
loaflist[[13]] <- model.PropSW.WCprev.Stationary
loaflist[[14]] <- model.PropSW.SD.Stationary
loaflist[[15]] <- model.PropFoodSub.WCprev.Stationary
loaflist[[16]] <- model.PropFoodSub.SD.Stationary

rm(loafresults)
rm(loafmodelselection)
for(i in 1:length(loaflist)){
  results.matrix <- summary(loaflist[[i]])$fixed
  indresult <- as.data.frame(results.matrix) %>%
    dplyr::select(mean, sd, "0.025quant", "0.975quant") %>%
    mutate(WeatherCov = rownames(results.matrix)[2]) %>%
    mutate(HabitatCov = rownames(results.matrix)[1]) %>%
    mutate(CovName = rownames(results.matrix))
  
  info.cow.ind <- data.frame(Beh_State = "Loafing",
                             LC_Cov = indresult[1,6],
                             LC_Coef = indresult[1,1],
                             Weath_Cov = indresult[1,5],
                             Weath_Coef = indresult[2,1],
                             Interaction = indresult[4,1])
  info.cow <- rbind(info.cow, info.cow.ind)
  
  ind.ms <- data.frame(Beh_State = "Loafing",
                       LC_Cov = indresult[1,6],
                       Weath_Cov = indresult[1,5],
                       WAIC = loaflist[[i]]$waic$waic,
                       DIC = loaflist[[i]]$dic$dic,
                       MLik = loaflist[[i]]$mlik[2,1]
  )
  
  if(exists("loafresults")){
    loafresults <- rbind(loafresults, indresult)
    loafmodelselection <-  rbind(loafmodelselection, ind.ms)
  }else{
    loafresults <- indresult
    loafmodelselection <- ind.ms
  }
}

write.csv(loafresults, "Results/loafresults.full.csv", row.names = F)
write.csv(loafmodelselection %>% arrange(Weath_Cov, DIC), "Results/loafmodelselection.csv", row.names = F)
loafinteractions <- loafresults %>%
  filter(grepl(":", CovName, fixed = T)) %>%
  mutate(Analysis = "Loafing")

fullinteractions <- rbind(roostinteractions, loafinteractions)




########################
### FORAGING RESULTS ###
########################
foragelist <- list()
foragelist[[1]] <- model.DtFE.WCprev.Mobile
foragelist[[2]] <- model.DtFE.SD.Mobile
foragelist[[3]] <- model.WindExp.WCprev.Mobile
foragelist[[4]] <- model.WindExp.SD.Mobile
foragelist[[5]] <- model.PropAg.WCprev.Mobile
foragelist[[6]] <- model.PropAg.SD.Mobile
foragelist[[7]] <- model.PropDev.WCprev.Mobile
foragelist[[8]] <- model.PropDev.SD.Mobile
foragelist[[9]] <- model.PropSW.WCprev.Mobile
foragelist[[10]] <- model.PropSW.SD.Mobile
foragelist[[11]] <- model.PropFoodSub.WCprev.Mobile
foragelist[[12]] <- model.PropFoodSub.SD.Mobile

rm(forageresults)
rm(foragemodelselection)
for(i in 1:length(foragelist)){
  results.matrix <- summary(foragelist[[i]])$fixed
  indresult <- as.data.frame(results.matrix) %>%
    dplyr::select(mean, sd, "0.025quant", "0.975quant") %>%
    mutate(WeatherCov = rownames(results.matrix)[2]) %>%
    mutate(HabitatCov = rownames(results.matrix)[1]) %>%
    mutate(CovName = rownames(results.matrix))
  
  info.cow.ind <- data.frame(Beh_State = "Foraging",
                             LC_Cov = indresult[1,6],
                             LC_Coef = indresult[1,1],
                             Weath_Cov = indresult[1,5],
                             Weath_Coef = indresult[2,1],
                             Interaction = indresult[4,1])
  info.cow <- rbind(info.cow, info.cow.ind)
  
  ind.ms <- data.frame(Beh_State = "Foraging",
                       LC_Cov = indresult[1,6],
                       Weath_Cov = indresult[1,5],
                       WAIC = foragelist[[i]]$waic$waic,
                       DIC = foragelist[[i]]$dic$dic,
                       MLik = foragelist[[i]]$mlik[2,1]
  )
  
  if(exists("forageresults")){
    forageresults <- rbind(forageresults, indresult)
    foragemodelselection <-  rbind(foragemodelselection, ind.ms)
  }else{
    forageresults <- indresult
    foragemodelselection <- ind.ms
  }
}

write.csv(forageresults, "Results/forageresults.full.csv", row.names = F)
forageinteractions <- forageresults %>%
  filter(grepl(":", CovName, fixed = T)) %>%
  mutate(Analysis = "Foraging")


fullinteractions <- rbind(fullinteractions, forageinteractions)
write.csv(fullinteractions, 'Results/InteractionResults.csv', row.names = F)
write.csv(foragemodelselection %>% arrange(Weath_Cov, DIC), "Results/foragemodelselection.csv", row.names = F)
write.csv(info.cow, 'Results/CowplotData.csv', row.names = F)
