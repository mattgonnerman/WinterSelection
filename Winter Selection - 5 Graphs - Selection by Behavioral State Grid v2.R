require(dplyr)
require(ggplot2)
require(forcats)
require(tidyr)


################################################################################################
### Plot Specific Resource Selection Relationships For Each Land Cover Covariate
Roost.SSF.results <- read.csv("Results/roostresults.full.csv") %>%
  mutate(Beh_State = "Roosting") %>%
  mutate(WeatherCov = ifelse(WeatherCov == WC_prev, "WC", WeatherCov))
Stationary.SSF.results <- read.csv("Results/stationaryresults.full.csv") %>%
  mutate(Beh_State = "Stationary") %>%
  mutate(WeatherCov = ifelse(WeatherCov == WC_prev, "WC", WeatherCov))
Mobile.SSF.results <- read.csv("Results/mobileresults.full.csv") %>%
  mutate(Beh_State = "Mobile") %>%
  mutate(WeatherCov = ifelse(WeatherCov == WC_prev, "WC", WeatherCov))

interactions.raw <- rbind(Roost.SSF.results, Stationary.SSF.results, Mobile.SSF.results) %>%
  mutate(CovName = ifelse(CovName == HabitatCov, "LC_Coef", 
                          ifelse(CovName == WeatherCov, "W_Coef",
                                 ifelse(CovName == "StepLength.Z", "SL_Coef",
                                        "Int_Coef")))) %>%
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
                  ifelse(HabitatCov == "PropSW.Z", "Softwoods",
                  ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                  ifelse(HabitatCov == "BA.Z", "Basal Area",
                  ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                  ifelse(HabitatCov == "SW.Z", "Percent Softwood",
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
                      LC.Val = rep(seq(-2, 2,.2),3),
                      W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
WE.SD.input <- WE.SD.df %>%
  # mutate(Check1 = (LC.Coef*LC.Val)) %>%
  # mutate(Check2 = (Int.Coef*LC.Val*W.Val)) %>%
  # mutate(Check3 = Check1 + Check2) %>%
  # mutate(Check4 = (LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)) %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
WE.SD.graph <- ggplot(data = WE.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
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
                       LC.Val = rep(seq(-2, 2,.2),3),
                       W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
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
                        values = c("longdash", "solid", "twodash"))

#Proportion Softwood
PropSW.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  filter(LC_Cov == "Softwoods" | LC_Cov == "Percent Softwood") %>%
  mutate(LC_Cov = "Percent Softwood") %>%
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
                         LC.Val = rep(seq(-2, 2,.2),3),
                         W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
PropSW.SD.input <- PropSW.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val))) %>%
  mutate(Est.LCL = exp((LC.LCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))%>%
  mutate(Est.UCL = exp((LC.UCL*LC.Val) + (Int.Coef*LC.Val*W.Val)))
PropSW.SD.graph <- ggplot(data = PropSW.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.LCL, ymax = Est.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 55) +
  xlab(PropSW.SD.df$LC[1]) + ylab("") +
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
                         LC.Val = rep(seq(-2, 2,.2),2),
                         W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
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

require(cowplot)
legend <- get_legend(WE.SD.graph + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(3,"inch")))
WE.SD.graph <- WE.SD.graph + theme(legend.position = "none")
DtFE.SD.graph <- DtFE.SD.graph + theme(legend.position = "none")
PropSW.SD.graph <- PropSW.SD.graph + theme(legend.position = "none")
FS.SD.graph <- FS.SD.graph + theme(legend.position = "none")

LC_grid <- plot_grid(plotlist = list(WE.SD.graph,DtFE.SD.graph,
                          PropSW.SD.graph, FS.SD.graph),
          nrow = 2,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)



jpeg('Results/Land Cover Grid.jpg', width = 2600, height = 1900)
plot_grid(plotlist = list(LC_grid, legend),
          nrow = 2, 
          rel_heights = c(1,.1)
)
dev.off()
