require(dplyr)
require(ggplot2)
require(forcats)

### Comparison of Magnitude of Interactions
interactions.raw <- read.csv("Results/InteractionResults.csv") %>% 
  mutate(HabitatCov = ifelse(HabitatCov == "DtFE.Z", "Dist. to Forest Edge",
                  ifelse(HabitatCov == "PropAg.Z", "Agriculture",
                  ifelse(HabitatCov == "PropDev.Z", "Developed",
                  ifelse(HabitatCov == "PropFoodSub.Z", "Food Subsidy",
                  ifelse(HabitatCov == "PropSW.Z", "Percent Softwood",
                  ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                  ifelse(HabitatCov == "BA.Z", "Basal Area",
                  ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                  ifelse(HabitatCov == "SW.Z", "Percent Softwood",
                         HabitatCov))))))))))%>%
  mutate(Beh_State = factor(Analysis, levels = c("Roosting", "Stationary", "Mobile"))) %>%
  mutate(LC_Cov = as.factor(HabitatCov)) %>%
  rename(Interaction = mean, SD = sd)
# interactions.raw$Beh_State <- factor(interactions.raw$Beh_State,
#                                      levels = c("Roost", "Loafing", "Foraging"))
# interactions.raw$LC_Cov <- factor(interactions.raw$LC_Cov,
#                                   levels = c("Wind Exposure", "Distance to Edge", "Proportion Ag",
#                                              "Proportion Dev", "Proportion SW", "% Softwood",
#                                              "Mean Tree Height", "Basal Area"))

#Graph showing interaction terms for snow depth
int.Snow <- interactions.raw %>% filter(WeatherCov == "SD")
ggplot(data = int.Snow, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2, size = 1.5) +
    geom_point(size = 3,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant),
                width = 0, size = 1.2,
                position = position_dodge(width = .4)) +
theme_bw(base_size = 20) + 
  xlab("Coefficient Estimate") +
  ylab("") +
  ggtitle("Snow Depth") +
  labs(color = "Behavioral\nState") + 
  theme(legend.title.align=0.5) + 
  scale_colour_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#1cade4", "#f1a806", "#46e300")) +   
  scale_shape_manual(name = "Behavioral\nState",
                     # labels = c("Roost", "Stationary", "Mobile"),
                     values = c(15, 19, 17))
ggsave("Results/SnowDepth_InteractionComp.jpeg", width = 10, height = 12, units = "in")


#Graph showing interaction terms for Previous days wind chill
int.Wind <- interactions.raw %>% filter(WeatherCov == "WC_prev")
ggplot(data = int.Wind, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2, size = 1.5) +
  geom_point(size = 3,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant),
                width = 0, size = 1.2,
                position = position_dodge(width = .4)) +
  theme_bw(base_size = 20) + 
  xlab("Coefficient Estimate") +
  ylab("") +
  ggtitle("Wind Chill") +
  labs(color = "Behavioral\nState") + 
  theme(legend.title.align=0.5) + 
  scale_colour_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#1cade4", "#f1a806", "#46e300")) +   
  scale_shape_manual(name = "Behavioral\nState",
                     # labels = c("Roost", "Stationary", "Mobile"),
                     values = c(15, 19, 17))
ggsave("Results/WindChill_InteractionComp.jpeg", width = 10, height = 12, units = "in")

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




################################################################################################
### Plot Matrix showing selection at Poor, Average, and Good Weather
require(cowplot)

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
  

int.Snow.Roost <- interactions.raw %>% filter(Weath_Cov == "SD") %>% filter(Beh_State == "Roosting")
int.Wind.Roost <- interactions.raw %>% filter(Weath_Cov == "WC_prev") %>% filter(Beh_State == "Roosting")
int.Snow.Loaf <- interactions.raw %>% filter(Weath_Cov == "SD") %>% filter(Beh_State == "Stationary") %>% distinct()
int.Wind.Loaf <- interactions.raw %>% filter(Weath_Cov == "WC_prev") %>% filter(Beh_State == "Stationary") %>% distinct()
int.Snow.Forage <- interactions.raw %>% filter(Weath_Cov == "SD") %>% filter(Beh_State == "Mobile") %>% distinct()
int.Wind.Forage <- interactions.raw %>% filter(Weath_Cov == "WC_prev") %>% filter(Beh_State == "Mobile") %>% distinct()
int.Wind <- interactions.raw %>% filter(Weath_Cov == "WC_prev")
int.Snow <- interactions.raw %>% filter(Weath_Cov == "SD")

# Condition Thresholds (Used summary on raw data and chose near 1st/3rd Quantile and Mean)
###Roosting 
##Snow Depth 
snow.list <- list()
snow.plots <- list()
for(i in 1:length(int.Snow.Roost$LC_Cov)){
  snow.df <- data.frame(Behavior = int.Snow.Roost$Beh_State[i],
                        LC = int.Snow.Roost$LC_Cov[i],
                        LC.Coef = int.Snow.Roost$LC_Coef[i],
                        W.Coef = int.Snow.Roost$Weath_Coef[i],
                        Int.Coef = int.Snow.Roost$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(-0.8962,0,6.2074), each = 21),
                        W.Condition = rep(c("Favorable","Average","Poor"), each = 21)) %>%
    mutate(W.Condition = factor(W.Condition, levels = c("Favorable","Average","Poor"))) %>%
    arrange(W.Condition)
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.6) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 35) +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))

jpeg('Results/Roost_Snow.jpg', width = 1500, height = 1500)
plot_grid(plotlist = snow.plots,
          legend,
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
          )
dev.off()

##Wind Chill
wind.list <- list()
wind.plots <- list()
for(i in 1:length(int.Wind.Roost$LC_Cov)){
  wind.df <- data.frame(Behavior = int.Wind.Roost$Beh_State[i],
                        LC = int.Wind.Roost$LC_Cov[i],
                        LC.Coef = int.Wind.Roost$LC_Coef[i],
                        W.Coef = int.Wind.Roost$Weath_Coef[i],
                        Int.Coef = int.Wind.Roost$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(1.83460,0,-2.15336), each = 21),
                        W.Condition = rep(c("Favorable","Average","Poor"), each = 21)) %>%
    mutate(W.Condition = factor(W.Condition, levels = c("Favorable","Average","Poor"))) %>%
    arrange(W.Condition)
  wind.list[[i]] <- wind.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
  wind.plot <- ggplot(data = wind.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.6) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 35) +
    xlab(wind.df$LC[i]) + ylab("")
  
  wind.plots[[i]] <- wind.plot + theme(legend.position="none")
}

legend <- get_legend(wind.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))


jpeg('Results/Roost_Wind.jpg', width = 1500, height = 1500)
plot_grid(plotlist = wind.plots,
          legend,
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)
dev.off()

###Loafing 
##Snow Depth 
snow.list <- list()
snow.plots <- list()
for(i in 1:length(int.Snow.Loaf$LC_Cov)){
  snow.df <- data.frame(Behavior = int.Snow.Loaf$Beh_State[i],
                        LC = int.Snow.Loaf$LC_Cov[i],
                        LC.Coef = int.Snow.Loaf$LC_Coef[i],
                        W.Coef = int.Snow.Loaf$Weath_Coef[i],
                        Int.Coef = int.Snow.Loaf$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(-1.13211,0,4.04389), each = 21),
                        W.Condition = rep(c("Favorable","Average","Poor"), each = 21)) %>%
    mutate(W.Condition = factor(W.Condition, levels = c("Favorable","Average","Poor")))
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.6) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 35) +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))

jpeg('Results/Loafing_Snow.jpg', width = 1800, height = 1500)
plot_grid(plotlist = snow.plots,
          legend,
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)
dev.off()

##Wind Chill
wind.list <- list()
wind.plots <- list()
for(i in 1:length(int.Wind.Loaf$LC_Cov)){
  wind.df <- data.frame(Behavior = int.Wind.Loaf$Beh_State[i],
                        LC = int.Wind.Loaf$LC_Cov[i],
                        LC.Coef = int.Wind.Loaf$LC_Coef[i],
                        W.Coef = int.Wind.Loaf$Weath_Coef[i],
                        Int.Coef = int.Wind.Loaf$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(2.04987,0,-1.97040), each = 21),
                        W.Condition = rep(c("Favorable","Average","Poor"), each = 21)) %>%
    mutate(W.Condition = factor(W.Condition, levels = c("Favorable","Average","Poor")))
  wind.list[[i]] <- wind.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
  wind.plot <- ggplot(data = wind.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.6) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 35) +
    xlab(wind.df$LC[i]) + ylab("")
  
  wind.plots[[i]] <- wind.plot + theme(legend.position="none")
}

legend <- get_legend(wind.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))

jpeg('Results/Loafing_Wind.jpg', width = 1800, height = 1500)
plot_grid(plotlist = wind.plots,
          legend,
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)
dev.off()

###Foraging
##Snow Depth 
snow.list <- list()
snow.plots <- list()
for(i in 1:length(int.Snow.Forage$LC_Cov)){
  snow.df <- data.frame(Behavior = int.Snow.Forage$Beh_State[i],
                        LC = int.Snow.Forage$LC_Cov[i],
                        LC.Coef = int.Snow.Forage$LC_Coef[i],
                        W.Coef = int.Snow.Forage$Weath_Coef[i],
                        Int.Coef = int.Snow.Forage$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(-.8782,0,7.1013), each = 21),
                        W.Condition = rep(c("Favorable","Average","Poor"), each = 21)) %>%
    mutate(W.Condition = factor(W.Condition, levels = c("Favorable","Average","Poor")))
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.6) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 35) +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))


jpeg('Results/Foraging_Snow.jpg', width = 1800, height = 1500)
plot_grid(plotlist = snow.plots,
          legend,
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)
dev.off()

##Wind Chill
wind.list <- list()
wind.plots <- list()
for(i in 1:length(int.Wind.Forage$LC_Cov)){
  wind.df <- data.frame(Behavior = int.Wind.Forage$Beh_State[i],
                        LC = int.Wind.Forage$LC_Cov[i],
                        LC.Coef = int.Wind.Forage$LC_Coef[i],
                        W.Coef = int.Wind.Forage$Weath_Coef[i],
                        Int.Coef = int.Wind.Forage$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(1.7557,0,-2.2648), each = 21),
                        W.Condition = rep(c("Favorable","Average","Poor"), each = 21)) %>%
    mutate(W.Condition = factor(W.Condition, levels = c("Favorable","Average","Poor")))
  wind.list[[i]] <- wind.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (Int.Coef*LC.Val*W.Val)))
  wind.plot <- ggplot(data = wind.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.6) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 35) +
    xlab(wind.df$LC[i]) + ylab("")
  
  wind.plots[[i]] <- wind.plot + theme(legend.position="none")
}

legend <- get_legend(wind.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))


jpeg('Results/Foraging_Wind.jpg', width = 1800, height = 1500)
plot_grid(plotlist = wind.plots,
          legend,
          nrow = 3,
          # labels = "auto",
          # label_size = 35,
          align = "hv",
          axis = "lb"
)
dev.off()

