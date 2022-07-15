require(dplyr)
require(ggplot2)
require(forcats)
require(tidyr)

'%notin%' <- Negate('%in%')

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
            Mean.HT = mean(Ht, na.rm = T),
            SD.HT = sd(Ht, na.rm = T)
  )
Loafing.raw <- read.csv("Stationaryfinal_HMM.csv")
Loafing.final <- Loafing.raw %>%
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
            Mean.HT = mean(Ht, na.rm = T),
            SD.HT = sd(Ht, na.rm = T)
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
            Mean.HT = mean(Ht),
            SD.HT = sd(Ht)
  )

################################################################################################
### Plot Specific Resource Selection Relationships For Each Land Cover Covariate
Stationary.SSF.results <- read.csv("Results/stationaryresults.full.csv") %>%
  mutate(Beh_State = "Stationary")
Mobile.SSF.results <- read.csv("Results/mobileresults.full.csv") %>%
  mutate(Beh_State = "Mobile")


interactions.DtFE <- rbind(Stationary.SSF.results, Mobile.SSF.results) %>%
  filter(grepl("DtFE", HabitatCov, fixed = T)) %>%
  # mutate(Int.Switch = ifelse(CovName == "DtFE.Z", 0, 1)) %>%
  mutate(WeatherCov = ifelse(WeatherCov == "WC_prev", "WC", WeatherCov)) %>%
  dplyr::select(-sd) %>%
  mutate(CovName = gsub("WC", "W", CovName)) %>%
  mutate(CovName = gsub("SD", "W", CovName)) %>%
  filter(CovName %notin% c("W:InFor", "W:InFor", "DtFE.Z:W:InFor")) %>%
  group_by(Beh_State, WeatherCov) %>%
  rename(Est = mean, LCL = X0.025quant, UCL = X0.975quant) %>%
  pivot_wider(names_from = c(CovName), values_from = c(Est, LCL, UCL)) %>% 
  ungroup() %>%
  dplyr::select(where(~!all(is.na(.x)))) %>%
  arrange(Beh_State, HabitatCov) %>% 
  mutate(HabitatCov = "Dist. to Forest Edge") %>%
  mutate(WeatherCov = ifelse(WeatherCov == "WC", "Wind Chill", "Snow Depth")) %>%
  dplyr::select(Beh_State, Weath_Cov = WeatherCov, LC_Cov = HabitatCov, 
                LC_Coef = Est_DtFE.Z, Weath_Coef = Est_W, Int_Coef = "Est_DtFE.Z:W",
                LCL = LCL_DtFE.Z, UCL = UCL_DtFE.Z,
                Est_InFor, Int_InFor = "Est_DtFE.Z:InFor", LCL_InFor = "LCL_InFor", UCL_InFor = "UCL_InFor")


#Separate Distance to Forest Edge Figure
DtFE.SD.raw <- interactions.DtFE %>%
  filter(LC_Cov == "Dist. to Forest Edge") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


DtFE.SD.df <- data.frame(Beh_State = c(rep("Stationary", 21),
                                       rep("Mobile", 21)),
                         LC.Val = c(seq(min(Loafing.raw$DtFE.Z), max(Loafing.raw$DtFE.Z),length.out = 21),
                                    seq(min(Foraging.raw$DtFE.Z), max(Foraging.raw$DtFE.Z),length.out = 21)),
                         W.Val = 0) 
DtFE.df <-  merge(DtFE.SD.df, DtFE.SD.raw, by = "Beh_State") 



DtFE.SD.input <- DtFE.df %>%
  filter(Weath_Cov == "Snow Depth") %>% 
  mutate(Est.out = exp(LC_Coef*LC.Val),
         Est.out.LCL = exp(LCL*LC.Val),
         Est.out.UCL = exp(UCL*LC.Val)) %>%
  mutate(Est.in = exp((LC_Coef*LC.Val) + Est_InFor + (LC.Val*Int_InFor)),
         Est.in.LCL = exp((LCL*LC.Val) + Est_InFor + (LC.Val*Int_InFor)),
         Est.in.UCL = exp((UCL*LC.Val) + Est_InFor + (LC.Val*Int_InFor)))%>%
  mutate(LC.Val = ifelse(Beh_State == "Roosting", (LC.Val*Roost.final$SD.DtFE[1]) + Roost.final$Mean.DtFE[1],
                         ifelse(Beh_State == "Stationary", (LC.Val*Loafing.final$SD.DtFE[1]) + Loafing.final$Mean.DtFE[1],
                                (LC.Val*Foraging.final$SD.DtFE[1]) + Foraging.final$Mean.DtFE[1]))) %>%
  rename(Behavior = Beh_State) 
DtFE.SD.out.graph <- ggplot(data = DtFE.SD.input, aes(x = LC.Val, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.out.LCL, ymax = Est.out.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(y = Est.out, color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 30) +
  xlab("Dist. to Forest Edge (m)") + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#46e300", "#f1a806")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#46e300", "#f1a806")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("twodash","solid")) +
  coord_cartesian(ylim = c(0,10)) +
  labs(subtitle = "Outside Forest Stand (Wind Chill)")
# DtFE.SD.out.graph

DtFE.SD.in.graph <- ggplot(data = DtFE.SD.input, aes(x = LC.Val, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.in.LCL, ymax = Est.in.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(y = Est.in, color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 30) +
  xlab("Dist. to Forest Edge (m)") + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#46e300", "#f1a806")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#46e300", "#f1a806")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("twodash","solid")) +
  coord_cartesian(ylim = c(0,30), xlim = c(0,200))+
  labs(subtitle = "Inside Forest Stand (Wind Chill)")
# DtFE.SD.in.graph

DtFE.WC.input <- DtFE.df %>%
  filter(Weath_Cov == "Wind Chill") %>% 
  mutate(Est.out = exp(LC_Coef*LC.Val),
         Est.out.LCL = exp(LCL*LC.Val),
         Est.out.UCL = exp(UCL*LC.Val)) %>%
  mutate(Est.in = exp((LC_Coef*LC.Val) + Est_InFor + (LC.Val*Int_InFor)),
         Est.in.LCL = exp((LCL*LC.Val) + Est_InFor + (LC.Val*Int_InFor)),
         Est.in.UCL = exp((UCL*LC.Val) + Est_InFor + (LC.Val*Int_InFor)))%>%
  mutate(LC.Val = ifelse(Beh_State == "Roosting", (LC.Val*Roost.final$SD.DtFE[1]) + Roost.final$Mean.DtFE[1],
                         ifelse(Beh_State == "Stationary", (LC.Val*Loafing.final$SD.DtFE[1]) + Loafing.final$Mean.DtFE[1],
                                (LC.Val*Foraging.final$SD.DtFE[1]) + Foraging.final$Mean.DtFE[1]))) %>%
  rename(Behavior = Beh_State) 
DtFE.WC.out.graph <- ggplot(data = DtFE.SD.input, aes(x = LC.Val, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.out.LCL, ymax = Est.out.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(y = Est.out, color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 30) +
  xlab("Dist. to Forest Edge (m)") + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#46e300", "#f1a806")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#46e300", "#f1a806")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("twodash","solid")) +
  coord_cartesian(ylim = c(0,10)) +
  labs(subtitle = "Outside Forest Stand (Wind Chill)")
# DtFE.WC.out.graph

DtFE.WC.in.graph <- ggplot(data = DtFE.SD.input, aes(x = LC.Val, group = Behavior)) +
  geom_ribbon(aes(ymin = Est.in.LCL, ymax = Est.in.UCL, fill = Behavior), alpha = .4) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(y = Est.in, color = Behavior, linetype = Behavior), size = 4) +
  theme_classic(base_size = 30) +
  xlab("Dist. to Forest Edge (m)") + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#46e300", "#f1a806")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roosting", "Stationary", "Mobile"),
                    values = c("#46e300", "#f1a806")) +
  scale_linetype_manual(name = "Behavioral\nState",
                        # labels = c("Roosting", "Stationary", "Mobile"),
                        values = c("twodash","solid")) +
  coord_cartesian(ylim = c(0,30), xlim = c(0,200)) +
  labs(subtitle = "Inside Forest Stand (Wind Chill)")
# DtFE.WC.in.graph


require(cowplot)
legend <- get_legend(DtFE.WC.in.graph + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(3,"inch")))
DTFE1.graph <- DtFE.SD.out.graph + theme(legend.position = "none")
DTFE2.graph <- DtFE.SD.in.graph + theme(legend.position = "none")
DTFE3.graph <- DtFE.WC.out.graph + theme(legend.position = "none")
DTFE4.graph <- DtFE.WC.in.graph + theme(legend.position = "none")


LC_grid <- plot_grid(plotlist = list(DTFE1.graph, DTFE2.graph,
                                     DTFE3.graph, DTFE4.graph),
                     nrow = 2,
                     # labels = "auto",
                     # label_size = 35,
                     align = "hv",
                     axis = "lb"
)



jpeg('Results/Land Cover Grid - DTFE.jpg', height = 20, width= 20, units = "in", res = 300)
plot_grid(plotlist = list(LC_grid, legend),
          nrow = 2, 
          rel_heights = c(1,.1)
)
dev.off()
