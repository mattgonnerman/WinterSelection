require(dplyr)
LCnames <- c("Dist. to Forest Edge", "Agriculture", "Developed", "Food Subsidy","Percent Softwood",
              "Wind Exposure", "Basal Area", "Mean Tree Height", "Percent Softwood")

Roost.SSF.results <- read.csv("Results/roostresults.full.csv") %>% 
  mutate(CovName = ifelse(CovName == "DtFE.Z", "Dist. to Forest Edge",
                          ifelse(CovName == "PropAg.Z", "Agriculture",
                                 ifelse(CovName == "PropDev.Z", "Developed",
                                        ifelse(CovName == "PropFoodSub.Z", "Food Subsidy",
                                               ifelse(CovName == "PropSW.Z", "Percent Softwood",
                                                      ifelse(CovName == "Wind.Exp.Z", "Wind Exposure",
                                                             ifelse(CovName == "BA.Z", "Basal Area",
                                                                    ifelse(CovName == "Ht.Z", "Mean Tree Height",
                                                                           ifelse(CovName == "SW.Z", "Percent Softwood",
                                                                                  CovName)))))))))) %>%
  mutate(HabitatCov = ifelse(HabitatCov == "DtFE.Z", "Dist. to Forest Edge",
                   ifelse(HabitatCov == "PropAg.Z", "Agriculture",
                                 ifelse(HabitatCov == "PropDev.Z", "Developed",
                                        ifelse(HabitatCov == "PropFoodSub.Z", "Food Subsidy",
                                               ifelse(HabitatCov == "PropSW.Z", "Percent Softwood",
                                                      ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                                                             ifelse(HabitatCov == "BA.Z", "Basal Area",
                                                                    ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                                                                           ifelse(HabitatCov == "SW.Z", "Percent Softwood",
                                                                                  HabitatCov)))))))))) %>%
  mutate(WeatherCov = ifelse(WeatherCov == "SD", "Snow Depth", "Wind Chill")) %>%
  filter(CovName %in% LCnames) %>% 
  mutate('Behavioral State' = "Roosting")

Stationary.SSF.results <- read.csv("Results/stationaryresults.full.csv") %>% 
  mutate(CovName = ifelse(CovName == "DtFE.Z", "Dist. to Forest Edge",
                          ifelse(CovName == "PropAg.Z", "Agriculture",
                                 ifelse(CovName == "PropDev.Z", "Developed",
                                        ifelse(CovName == "PropFoodSub.Z", "Food Subsidy",
                                               ifelse(CovName == "PropSW.Z", "Percent Softwood",
                                                      ifelse(CovName == "Wind.Exp.Z", "Wind Exposure",
                                                             ifelse(CovName == "BA.Z", "Basal Area",
                                                                    ifelse(CovName == "Ht.Z", "Mean Tree Height",
                                                                           ifelse(CovName == "SW.Z", "Percent Softwood",
                                                                                  CovName)))))))))) %>%
  mutate(HabitatCov = ifelse(HabitatCov == "DtFE.Z", "Dist. to Forest Edge",
                             ifelse(HabitatCov == "PropAg.Z", "Agriculture",
                                    ifelse(HabitatCov == "PropDev.Z", "Developed",
                                           ifelse(HabitatCov == "PropFoodSub.Z", "Food Subsidy",
                                                  ifelse(HabitatCov == "PropSW.Z", "Percent Softwood",
                                                         ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                                                                ifelse(HabitatCov == "BA.Z", "Basal Area",
                                                                       ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                                                                              ifelse(HabitatCov == "SW.Z", "Percent Softwood",
                                                                                     HabitatCov)))))))))) %>%
  mutate(WeatherCov = ifelse(WeatherCov == "SD", "Snow Depth", "Wind Chill")) %>%
  filter(CovName %in% LCnames) %>% 
  mutate('Behavioral State' = "Stationary")

Mobile.SSF.results <- read.csv("Results/mobileresults.full.csv") %>% 
  mutate(CovName = ifelse(CovName == "DtFE.Z", "Dist. to Forest Edge",
                          ifelse(CovName == "PropAg.Z", "Agriculture",
                                 ifelse(CovName == "PropDev.Z", "Developed",
                                        ifelse(CovName == "PropFoodSub.Z", "Food Subsidy",
                                               ifelse(CovName == "PropSW.Z", "Percent Softwood",
                                                      ifelse(CovName == "Wind.Exp.Z", "Wind Exposure",
                                                             ifelse(CovName == "BA.Z", "Basal Area",
                                                                    ifelse(CovName == "Ht.Z", "Mean Tree Height",
                                                                           ifelse(CovName == "SW.Z", "Percent Softwood",
                                                                                  CovName)))))))))) %>%
  mutate(HabitatCov = ifelse(HabitatCov == "DtFE.Z", "Dist. to Forest Edge",
                             ifelse(HabitatCov == "PropAg.Z", "Agriculture",
                                    ifelse(HabitatCov == "PropDev.Z", "Developed",
                                           ifelse(HabitatCov == "PropFoodSub.Z", "Food Subsidy",
                                                  ifelse(HabitatCov == "PropSW.Z", "Percent Softwood",
                                                         ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                                                                ifelse(HabitatCov == "BA.Z", "Basal Area",
                                                                       ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                                                                              ifelse(HabitatCov == "SW.Z", "Percent Softwood",
                                                                                     HabitatCov)))))))))) %>%
  mutate(WeatherCov = ifelse(WeatherCov == "SD", "Snow Depth", "Wind Chill")) %>%
  filter(CovName %in% LCnames) %>% 
  mutate('Behavioral State' = "Mobile")

allLCresults <- rbind(Roost.SSF.results, Stationary.SSF.results, Mobile.SSF.results) %>%
  dplyr::select(-CovName)
write.csv(allLCresults, "Results/LC_Coefficients.csv", row.names = F)
