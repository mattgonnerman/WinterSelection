require(dplyr)
require(ggplot2)
require(forcats)
require(tidyr)

### Bring in the original data to unscale
Roost.raw <- read.csv("Roostfinal_HMM.csv")
Roost.final <- Roost.raw %>%
  summarize(Mean.DtFE = mean(DtFE),
            SD.DtFE = sd(DtFE),
            Mean.BA = mean(BA),
            SD.BA = sd(BA),
            Mean.SW = mean(SW),
            SD.SW = sd(SW),
            Mean.WE = mean(Wind_Exp),
            SD.WE = sd(Wind_Exp),
            Mean.AG = mean(PropAg),
            SD.AG = sd(PropAg),
            Mean.DEV = mean(PropDev),
            SD.DEV = sd(PropDev),
  )
Loafing.raw <- read.csv("Stationaryfinal_HMM.csv")
Loafing.final <- Roost.raw %>%
  summarize(Mean.DtFE = mean(DtFE),
            SD.DtFE = sd(DtFE),
            Mean.BA = mean(BA),
            SD.BA = sd(BA),
            Mean.SW = mean(SW),
            SD.SW = sd(SW),
            Mean.WE = mean(Wind_Exp),
            SD.WE = sd(Wind_Exp),
            Mean.AG = mean(PropAg),
            SD.AG = sd(PropAg),
            Mean.DEV = mean(PropDev),
            SD.DEV = sd(PropDev),
  )
Foraging.raw <- read.csv("Stationaryfinal_HMM.csv")
Foraging.final <- Foraging.raw %>%
  summarize(Mean.DtFE = mean(DtFE),
            SD.DtFE = sd(DtFE),
            Mean.BA = mean(BA),
            SD.BA = sd(BA),
            Mean.SW = mean(SW),
            SD.SW = sd(SW),
            Mean.WE = mean(Wind_Exp),
            SD.WE = sd(Wind_Exp),
            Mean.AG = mean(PropAg),
            SD.AG = sd(PropAg),
            Mean.DEV = mean(PropDev),
            SD.DEV = sd(PropDev),
  )


################################################################################################
### Plot Specific Resource Selection Relationships For Each Land Cover Covariate
Roost.SSF.results <- read.csv("Results/roostresults.full.csv") %>%
  mutate(Beh_State = "Roosting")
Stationary.SSF.results <- read.csv("Results/stationaryresults.full.csv") %>%
  mutate(Beh_State = "Stationary")
Mobile.SSF.results <- read.csv("Results/mobileresults.full.csv") %>%
  mutate(Beh_State = "Mobile")

interactions.raw <- rbind(Roost.SSF.results, Stationary.SSF.results, Mobile.SSF.results) %>%
  mutate(CovName = ifelse(CovName == HabitatCov, "LC_Coef", 
                          ifelse(CovName == WeatherCov, "W_Coef",
                                 ifelse(CovName == "StepLength.Z", "SL_Coef",
                                        "Int_Coef"))))  %>%
  mutate(WeatherCov = ifelse(WeatherCov == "WC_prev", "WC", WeatherCov)) %>%
  dplyr::select(-sd) %>%
  mutate(X0.025quant = ifelse(CovName == "LC_Coef", X0.025quant, NA)) %>%
  mutate(X0.975quant = ifelse(CovName == "LC_Coef", X0.975quant, NA)) %>%
  group_by(Beh_State, HabitatCov, WeatherCov) %>%
  pivot_wider(names_from = CovName, values_from = c(mean, X0.025quant, X0.975quant))%>% 
  ungroup() %>%
  dplyr::select(where(~!all(is.na(.x)))) %>%
  arrange(Beh_State, HabitatCov) %>% 
  mutate(HabitatCov = ifelse(HabitatCov == "DtFE.Z", "Dist. to Forest Edge",
                  ifelse(HabitatCov == "PropAg.Z", "Agriculture",
                  ifelse(HabitatCov == "PropDev.Z", "Developed",
                  ifelse(HabitatCov == "PropFoodSub.Z", "Food Subsidy",
                  ifelse(HabitatCov == "PropSW.Z", "Proportion Conifer",
                  ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                  ifelse(HabitatCov == "BA.Z", "Basal Area",
                  ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                  ifelse(HabitatCov == "SW.Z", "Percent Conifer",
                         HabitatCov)))))))))) %>%
  mutate(WeatherCov = ifelse(WeatherCov == "WC", "Wind Chill", "Snow Depth")) %>%
  dplyr::select(Beh_State, Weath_Cov = WeatherCov, LC_Cov = HabitatCov, LC_Coef = mean_LC_Coef, 
         Weath_Coef = mean_W_Coef, Int_Coef = mean_Int_Coef, LCL = X0.025quant_LC_Coef, UCL = X0.975quant_LC_Coef)


#Wind Exposure
WE.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Wind Exposure") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


WE.SD.df <- data.frame(Behavior = c(rep("Roosting", 21),
                                    rep("Stationary", 21),
                                    rep("Mobile", 21)),
                      LC = WE.SD.raw$LC_Cov[1],
                      LC.Coef = c(rep(WE.SD.raw$LC_Coef[1], 21),
                                  rep(WE.SD.raw$LC_Coef[2], 21),
                                  rep(WE.SD.raw$LC_Coef[3], 21)),
                      LC.LCL = c(rep(WE.SD.raw$LCL[1], 21),
                                  rep(WE.SD.raw$LCL[2], 21),
                                  rep(WE.SD.raw$LCL[3], 21)),
                      LC.UCL = c(rep(WE.SD.raw$UCL[1], 21),
                                  rep(WE.SD.raw$UCL[2], 21),
                                  rep(WE.SD.raw$UCL[3], 21)),
                      W.Coef = c(rep(WE.SD.raw$Weath_Coef[1], 21),
                                 rep(WE.SD.raw$Weath_Coef[2], 21),
                                 rep(WE.SD.raw$Weath_Coef[3], 21)),
                      Int.Coef = c(rep(WE.SD.raw$Int_Coef[1], 21),
                                   rep(WE.SD.raw$Int_Coef[2], 21),
                                   rep(WE.SD.raw$Int_Coef[3], 21)),
                      LC.Val = c(seq(min(Roost.raw$Wind.Exp.Z), max(Roost.raw$Wind.Exp.Z),length.out = 21),
                                 seq(min(Loafing.raw$Wind.Exp.Z), max(Loafing.raw$Wind.Exp.Z),length.out = 21),
                                 seq(min(Foraging.raw$Wind.Exp.Z), max(Foraging.raw$Wind.Exp.Z),length.out = 21)),
                      W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) %>%
  mutate(LC.Val = ifelse(Behavior == "Roosting", (LC.Val*Roost.final$SD.WE[1]) + Roost.final$Mean.WE[1],
                         ifelse(Behavior == "Stationary", (LC.Val*Loafing.final$SD.WE[1]) + Loafing.final$Mean.WE[1],
                                (LC.Val*Foraging.final$SD.WE[1]) + Foraging.final$Mean.WE[1])))

WE.SD.input <- WE.SD.df %>%
  # mutate(Check1 = (LC.Coef*LC.Val)) %>%
  # mutate(Check2 = (Int.Coef*LC.Val*W.Val)) %>%
  # mutate(Check3 = Check1 + Check2) %>%
  # mutate(Check4 = (LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)) %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))

WE.SD.graph <- ggplot(data = WE.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .2) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(WE.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                   # labels = c("Roosting", "Stationary", "Mobile"),
                   values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("longdash", "solid", "twodash")) +
  scale_y_continuous(trans = "log10")

#Basal Area
BA.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Basal Area") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


BA.SD.df <- data.frame(Behavior = c(rep("Roosting", 21),
                                    rep("Stationary", 21),
                                    rep("Mobile", 21)),
                       LC = BA.SD.raw$LC_Cov[1],
                       LC.Coef = c(rep(BA.SD.raw$LC_Coef[1], 21),
                                   rep(BA.SD.raw$LC_Coef[2], 21),
                                   rep(BA.SD.raw$LC_Coef[3], 21)),
                       LC.LCL = c(rep(BA.SD.raw$LCL[1], 21),
                                  rep(BA.SD.raw$LCL[2], 21),
                                  rep(BA.SD.raw$LCL[3], 21)),
                       LC.UCL = c(rep(BA.SD.raw$UCL[1], 21),
                                  rep(BA.SD.raw$UCL[2], 21),
                                  rep(BA.SD.raw$UCL[3], 21)),
                       W.Coef = c(rep(BA.SD.raw$Weath_Coef[1], 21),
                                  rep(BA.SD.raw$Weath_Coef[2], 21),
                                  rep(BA.SD.raw$Weath_Coef[3], 21)),
                       Int.Coef = c(rep(BA.SD.raw$Int_Coef[1], 21),
                                    rep(BA.SD.raw$Int_Coef[2], 21),
                                    rep(BA.SD.raw$Int_Coef[3], 21)),
                       LC.Val = c(seq(min(Roost.raw$BA.Z), max(Roost.raw$BA.Z),length.out = 21),
                                  seq(min(Loafing.raw$BA.Z), max(Loafing.raw$BA.Z),length.out = 21),
                                  seq(min(Foraging.raw$BA.Z), max(Foraging.raw$BA.Z),length.out = 21)),
                       W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) %>%
  mutate(LC.Val = ifelse(Behavior == "Roosting", (LC.Val*Roost.final$SD.BA[1]) + Roost.final$Mean.BA[1],
                         ifelse(Behavior == "Stationary", (LC.Val*Loafing.final$SD.BA[1]) + Loafing.final$Mean.BA[1],
                                (LC.Val*Foraging.final$SD.BA[1]) + Foraging.final$Mean.BA[1])))


BA.SD.input <- BA.SD.df %>%
  # mutate(Check1 = (LC.Coef*LC.Val)) %>%
  # mutate(Check2 = (Int.Coef*LC.Val*W.Val)) %>%
  # mutate(Check3 = Check1 + Check2) %>%
  # mutate(Check4 = (LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)) %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
BA.SD.graph <- ggplot(data = BA.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(BA.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("longdash", "solid", "twodash"))

#Distance to Forest Edge
DtFE.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Dist. to Forest Edge") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


DtFE.SD.df <- data.frame(Behavior = c(rep("Roosting", 21),
                                    rep("Stationary", 21),
                                    rep("Mobile", 21)),
                       LC = DtFE.SD.raw$LC_Cov[1],
                       LC.Coef = c(rep(DtFE.SD.raw$LC_Coef[1], 21),
                                   rep(DtFE.SD.raw$LC_Coef[2], 21),
                                   rep(DtFE.SD.raw$LC_Coef[3], 21)),
                       LC.LCL = c(rep(DtFE.SD.raw$LCL[1], 21),
                                  rep(DtFE.SD.raw$LCL[2], 21),
                                  rep(DtFE.SD.raw$LCL[3], 21)),
                       LC.UCL = c(rep(DtFE.SD.raw$UCL[1], 21),
                                  rep(DtFE.SD.raw$UCL[2], 21),
                                  rep(DtFE.SD.raw$UCL[3], 21)),
                       W.Coef = c(rep(DtFE.SD.raw$Weath_Coef[1], 21),
                                  rep(DtFE.SD.raw$Weath_Coef[2], 21),
                                  rep(DtFE.SD.raw$Weath_Coef[3], 21)),
                       Int.Coef = c(rep(DtFE.SD.raw$Int_Coef[1], 21),
                                    rep(DtFE.SD.raw$Int_Coef[2], 21),
                                    rep(DtFE.SD.raw$Int_Coef[3], 21)),
                       LC.Val = c(seq(min(Roost.raw$DtFE.Z), max(Roost.raw$DtFE.Z),length.out = 21),
                                  seq(min(Loafing.raw$DtFE.Z), max(Loafing.raw$DtFE.Z),length.out = 21),
                                  seq(min(Foraging.raw$DtFE.Z), max(Foraging.raw$DtFE.Z),length.out = 21)),
                       W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) %>%
  mutate(LC.Val = ifelse(Behavior == "Roosting", (LC.Val*Roost.final$SD.DtFE[1]) + Roost.final$Mean.DtFE[1],
                         ifelse(Behavior == "Stationary", (LC.Val*Loafing.final$SD.DtFE[1]) + Loafing.final$Mean.DtFE[1],
                                (LC.Val*Foraging.final$SD.DtFE[1]) + Foraging.final$Mean.DtFE[1])))

DtFE.SD.input <- DtFE.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val))) %>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
DtFE.SD.graph <- ggplot(data = DtFE.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(DtFE.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("longdash", "solid", "twodash")) +
  scale_y_continuous(trans = "log10")

#Percent Softwood
PropSW.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Percent Conifer") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


PropSW.SD.df <- data.frame(Behavior = c(rep("Roosting", 21),
                                      rep("Stationary", 21),
                                      rep("Mobile", 21)),
                         LC = PropSW.SD.raw$LC_Cov[1],
                         LC.Coef = c(rep(PropSW.SD.raw$LC_Coef[1], 21),
                                     rep(PropSW.SD.raw$LC_Coef[2], 21),
                                     rep(PropSW.SD.raw$LC_Coef[3], 21)),
                         LC.LCL = c(rep(PropSW.SD.raw$LCL[1], 21),
                                    rep(PropSW.SD.raw$LCL[2], 21),
                                    rep(PropSW.SD.raw$LCL[3], 21)),
                         LC.UCL = c(rep(PropSW.SD.raw$UCL[1], 21),
                                    rep(PropSW.SD.raw$UCL[2], 21),
                                    rep(PropSW.SD.raw$UCL[3], 21)),
                         W.Coef = c(rep(PropSW.SD.raw$Weath_Coef[1], 21),
                                    rep(PropSW.SD.raw$Weath_Coef[2], 21),
                                    rep(PropSW.SD.raw$Weath_Coef[3], 21)),
                         Int.Coef = c(rep(PropSW.SD.raw$Int_Coef[1], 21),
                                      rep(PropSW.SD.raw$Int_Coef[2], 21),
                                      rep(PropSW.SD.raw$Int_Coef[3], 21)),
                         LC.Val = c(seq(min(Roost.raw$SW.Z), max(Roost.raw$SW.Z),length.out = 21),
                                    seq(min(Loafing.raw$SW.Z), max(Loafing.raw$SW.Z),length.out = 21),
                                    seq(min(Foraging.raw$SW.Z), max(Foraging.raw$SW.Z),length.out = 21)),
                         W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) %>%
  mutate(LC.Val = ifelse(Behavior == "Roosting", (LC.Val*Roost.final$SD.SW[1]) + Roost.final$Mean.SW[1],
                         ifelse(Behavior == "Stationary", (LC.Val*Loafing.final$SD.SW[1]) + Loafing.final$Mean.SW[1],
                                (LC.Val*Foraging.final$SD.SW[1]) + Foraging.final$Mean.SW[1])))
PropSW.SD.input <- PropSW.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val))) %>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
PropSW.SD.graph <- ggplot(data = PropSW.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(PropSW.SD.df$LC[1]) + ylab("Relative Probability of Selection") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("longdash", "solid", "twodash")) +
  scale_y_continuous(trans = "log10")

#Agriculture
FS.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Agriculture") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


FS.SD.df <- data.frame(Behavior = c(rep("Stationary", 21),
                                      rep("Mobile", 21)),
                         LC = FS.SD.raw$LC_Cov[1],
                         LC.Coef = c(rep(FS.SD.raw$LC_Coef[1], 21),
                                     rep(FS.SD.raw$LC_Coef[2], 21)),
                         LC.LCL = c(rep(FS.SD.raw$LCL[1], 21),
                                    rep(FS.SD.raw$LCL[2], 21)),
                         LC.UCL = c(rep(FS.SD.raw$UCL[1], 21),
                                    rep(FS.SD.raw$UCL[2], 21)),
                         W.Coef = c(rep(FS.SD.raw$Weath_Coef[1], 21),
                                    rep(FS.SD.raw$Weath_Coef[2], 21)),
                         Int.Coef = c(rep(FS.SD.raw$Int_Coef[1], 21),
                                      rep(FS.SD.raw$Int_Coef[2], 21)),
                       LC.Val = c(seq(min(Loafing.raw$PropAg.Z), max(Loafing.raw$PropAg.Z),length.out = 21),
                                  seq(min(Foraging.raw$PropAg.Z), max(Foraging.raw$PropAg.Z),length.out = 21)),
                         W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile")))  %>%
  mutate(LC.Val = ifelse(Behavior == "Roosting", (LC.Val*Roost.final$SD.AG[1]) + Roost.final$Mean.AG[1],
                         ifelse(Behavior == "Stationary", (LC.Val*Loafing.final$SD.AG[1]) + Loafing.final$Mean.AG[1],
                                (LC.Val*Foraging.final$SD.AG[1]) + Foraging.final$Mean.AG[1])))

FS.SD.input <- FS.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val))) %>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
FS.SD.graph <- ggplot(data = FS.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(FS.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#f1a806", "#46e300")) +
  scale_linetype_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("solid", "twodash"))


FS.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Agriculture") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)

#Developds
DEV.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Developed") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


DEV.SD.df <- data.frame(Behavior = c(rep("Stationary", 21),
                                    rep("Mobile", 21)),
                       LC = DEV.SD.raw$LC_Cov[1],
                       LC.Coef = c(rep(DEV.SD.raw$LC_Coef[1], 21),
                                   rep(DEV.SD.raw$LC_Coef[2], 21)),
                       LC.LCL = c(rep(DEV.SD.raw$LCL[1], 21),
                                  rep(DEV.SD.raw$LCL[2], 21)),
                       LC.UCL = c(rep(DEV.SD.raw$UCL[1], 21),
                                  rep(DEV.SD.raw$UCL[2], 21)),
                       W.Coef = c(rep(DEV.SD.raw$Weath_Coef[1], 21),
                                  rep(DEV.SD.raw$Weath_Coef[2], 21)),
                       Int.Coef = c(rep(DEV.SD.raw$Int_Coef[1], 21),
                                    rep(DEV.SD.raw$Int_Coef[2], 21)),
                       LC.Val = c(seq(min(Loafing.raw$PropDev.Z), max(Loafing.raw$PropDev.Z),length.out = 21),
                                  seq(min(Foraging.raw$PropDev.Z), max(Foraging.raw$PropDev.Z),length.out = 21)),
                       W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile")))   %>%
  mutate(LC.Val = ifelse(Behavior == "Roosting", (LC.Val*Roost.final$SD.DEV[1]) + Roost.final$Mean.DEV[1],
                         ifelse(Behavior == "Stationary", (LC.Val*Loafing.final$SD.DEV[1]) + Loafing.final$Mean.DEV[1],
                                (LC.Val*Foraging.final$SD.DEV[1]) + Foraging.final$Mean.DEV[1])))
DEV.SD.input <- DEV.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val))) %>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
DEV.SD.graph <- ggplot(data = DEV.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(DEV.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#f1a806", "#46e300")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("solid", "twodash"))


DEV.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Agriculture") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)

require(cowplot)
legend <- get_legend(WE.SD.graph + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(3,"inch")))
WE.SD.graph <- WE.SD.graph + theme(legend.position = "none")
BA.SD.graph <- BA.SD.graph + theme(legend.position = "none")
DtFE.SD.graph <- DtFE.SD.graph + theme(legend.position = "none")
PropSW.SD.graph <- PropSW.SD.graph + theme(legend.position = "none")
FS.SD.graph <- FS.SD.graph + theme(legend.position = "none")
DEV.SD.graph <- DEV.SD.graph + theme(legend.position = "none")

LC_grid <- plot_grid(plotlist = list(DtFE.SD.graph, BA.SD.graph,
                          PropSW.SD.graph, WE.SD.graph,
                          FS.SD.graph, DEV.SD.graph),
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)



jpeg('Results/Land Cover Grid.jpg', height = 35, width= 31, units = "in", res = 300)
plot_grid(plotlist = list(LC_grid, legend),
          nrow = 2, 
          rel_heights = c(1,.1)
)
dev.off()
