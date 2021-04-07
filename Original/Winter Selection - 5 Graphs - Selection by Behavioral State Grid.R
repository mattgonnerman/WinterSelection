require(dplyr)
require(ggplot2)
require(forcats)

################################################################################################
### Plot Specific Resource Selection Relationships For Each Land Cover Covariate
interactions.raw <- read.csv("Results/CowplotData.csv") %>%
  arrange(Beh_State, LC_Cov) %>% 
  mutate(LC_Cov = ifelse(LC_Cov == "DtFE.Z", "Dist. to Forest Edge",
                  ifelse(LC_Cov == "PropAg.Z", "Agriculture",
                  ifelse(LC_Cov == "PropDev.Z", "Developed",
                  ifelse(LC_Cov == "PropFoodSub.Z", "Food Subsidy",
                  ifelse(LC_Cov == "PropSW.Z", "Softwoods",
                  ifelse(LC_Cov == "Wind.Exp.Z", "Wind Exposure",
                  ifelse(LC_Cov == "BA.Z", "Basal Area",
                  ifelse(LC_Cov == "Ht.Z", "Mean Tree Height",
                  ifelse(LC_Cov == "SW.Z", "Percent Softwood",
                  LC_Cov))))))))))


#Wind Exposure
WE.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "SD") %>% 
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
                      W.Coef = c(rep(WE.SD.raw$Weath_Coef[1], 21),
                                 rep(WE.SD.raw$Weath_Coef[2], 21),
                                 rep(WE.SD.raw$Weath_Coef[3], 21)),
                      Int.Coef = c(rep(WE.SD.raw$Interaction[1], 21),
                                   rep(WE.SD.raw$Interaction[2], 21),
                                   rep(WE.SD.raw$Interaction[3], 21)),
                      LC.Val = rep(seq(-2, 2,.2),3),
                      W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
WE.SD.input <- WE.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
WE.SD.graph <- ggplot(data = WE.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior), size = 3) +
  theme_classic(base_size = 55) +
  xlab(WE.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300"))

#Distance to Forest Edge
DtFE.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "SD") %>% 
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
                       W.Coef = c(rep(DtFE.SD.raw$Weath_Coef[1], 21),
                                  rep(DtFE.SD.raw$Weath_Coef[2], 21),
                                  rep(DtFE.SD.raw$Weath_Coef[3], 21)),
                       Int.Coef = c(rep(DtFE.SD.raw$Interaction[1], 21),
                                    rep(DtFE.SD.raw$Interaction[2], 21),
                                    rep(DtFE.SD.raw$Interaction[3], 21)),
                       LC.Val = rep(seq(-2, 2,.2),3),
                       W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
DtFE.SD.input <- DtFE.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
DtFE.SD.graph <- ggplot(data = DtFE.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior), size = 3) +
  theme_classic(base_size = 55) +
  xlab(DtFE.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300"))

#Softwood
PropSW.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "SD") %>% 
  filter(LC_Cov == "Softwoods") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)

PercSWSD.raw <- interactions.raw %>%
  filter(Weath_Cov == "SD") %>% 
  filter(LC_Cov == "Percent Softwood") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) 

PropSW.SD.raw <- rbind(PropSW.SD.raw, PercSWSD.raw)%>%
  arrange(Beh_State)


PropSW.SD.df <- data.frame(Behavior = c(rep("Roosting", 21),
                                        rep("Stationary", 21),
                                      rep("Mobile", 21)),
                         LC = PropSW.SD.raw$LC_Cov[1],
                         LC.Coef = c(rep(PropSW.SD.raw$LC_Coef[1], 21),
                                     rep(PropSW.SD.raw$LC_Coef[2], 21),
                                     rep(PropSW.SD.raw$LC_Coef[3], 21)),
                         W.Coef = c(rep(PropSW.SD.raw$Weath_Coef[1], 21),
                                    rep(PropSW.SD.raw$Weath_Coef[2], 21),
                                    rep(PropSW.SD.raw$Weath_Coef[3], 21)),
                         Int.Coef = c(rep(PropSW.SD.raw$Interaction[1], 21),
                                      rep(PropSW.SD.raw$Interaction[2], 21),
                                      rep(PropSW.SD.raw$Interaction[3], 21)),
                         LC.Val = rep(seq(-2, 2,.2),3),
                         W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) 
PropSW.SD.input <- PropSW.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
PropSW.SD.graph <- ggplot(data = PropSW.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior), size = 3) +
  theme_classic(base_size = 55) +
  xlab(PropSW.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#1cade4", "#f1a806", "#46e300"))

#Food Subsidy
FS.SD.raw <- interactions.raw %>%
  filter(Weath_Cov == "SD") %>% 
  filter(LC_Cov == "Food Subsidy") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Beh_State)


FS.SD.df <- data.frame(Behavior = c(rep("Stationary", 21),
                                    rep("Mobile", 21)),
                       LC = FS.SD.raw$LC_Cov[1],
                       LC.Coef = c(rep(FS.SD.raw$LC_Coef[1], 21),
                                   rep(FS.SD.raw$LC_Coef[2], 21)),
                       W.Coef = c(rep(FS.SD.raw$Weath_Coef[1], 21),
                                  rep(FS.SD.raw$Weath_Coef[2], 21)),
                       Int.Coef = c(rep(FS.SD.raw$Interaction[1], 21),
                                    rep(FS.SD.raw$Interaction[2], 21)),
                       LC.Val = rep(seq(-2, 2,.2),2),
                       W.Val = 4) %>%
  mutate(Behavior = factor(Behavior, levels = c("Roosting","Stationary","Mobile"))) %>%
  arrange(Behavior)
FS.SD.input <- FS.SD.df %>%
  mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
FS.SD.graph <- ggplot(data = FS.SD.input, aes(x = LC.Val, y = Est, group = Behavior)) +
  geom_hline(yintercept = 1, color = "black", linetype = 2, size = 1.5) +
  geom_line(aes(color = Behavior), size = 3) +
  theme_classic(base_size = 55) +
  xlab(FS.SD.df$LC[1]) + ylab("") +
  theme(legend.key.width=unit(.6,"inch")) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roosting", "Stationary", "Mobile"),
                     values = c("#f1a806", "#46e300"))

require(cowplot)
legend <- get_legend(WE.SD.graph + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(1,"inch")))
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





