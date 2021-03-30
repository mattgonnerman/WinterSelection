require(dplyr)
require(ggplot2)
require(forcats)

### Comparison of Magnitude of Interactions
interactions.raw <- read.csv("Results/InteractionResults.csv") %>%
  mutate(Beh_State = as.factor(Analysis)) %>%
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
  geom_point(size = 1.5,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant),
                width = .2,
                position = position_dodge(width = .4)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2) +
  theme_bw() + 
  xlab("Coefficient Estimate") +
  ylab("") +
  ggtitle("Snow Depth") +
  labs(color = "Behavioral\nState") + 
  theme(legend.title.align=0.5) + 
  scale_colour_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#46e300", "#f1a806", "#1cade4")) +   
  scale_shape_manual(name = "Behavioral\nState",
                     # labels = c("Roost", "Stationary", "Mobile"),
                     values = c(15, 19, 17))
ggsave("Results/SnowDepth_InteractionComp.jpeg", width = 8, height = 7, units = "in")


#Graph showing interaction terms for Previous days wind chill
int.Wind <- interactions.raw %>% filter(WeatherCov == "WC_prev")
ggplot(data = int.Wind, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_point(size = 1.5,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant),
                width = .2,
                position = position_dodge(width = .4)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2) +
  theme_bw() + 
  xlab("Coefficient Estimate") +
  ylab("") +
  ggtitle("Wind Chill") +
  labs(color = "Behavioral\nState") + 
  theme(legend.title.align=0.5) + 
  scale_colour_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#46e300", "#f1a806", "#1cade4")) +   
  scale_shape_manual(name = "Behavioral\nState",
                     # labels = c("Roost", "Stationary", "Mobile"),
                     values = c(15, 19, 17))
ggsave("Results/WindChill_InteractionComp.jpeg", width = 8, height = 7, units = "in")


#make big points
#remove endcaps on error bars
#thicker error lines

################################################################################################
### Plot Matrix showing selection at Poor, Average, and Good Weather
require(cowplot)

interactions.raw <- read.csv("Results/CowplotData.csv") %>%
  # mutate(Beh_State = factor(Beh_State, levels = c("Roost", "Stationary", "Mobile"))) %>%
  # mutate(LC_Cov = factor(LC_Cov,
  #                        levels = c("Distance to Edge", "Wind Exposure", "Proportion Ag",
  #                                   "Proportion Dev", "Proportion SW",
  #                                   "Mean Tree Height", "Basal Area", "% Softwood"))) %>%
  arrange(Beh_State, LC_Cov)

int.Snow.Roost <- interactions.raw %>% filter(Weath_Cov == "SD") %>% filter(Beh_State == "Roosting")
int.Wind.Roost <- interactions.raw %>% filter(Weath_Cov == "WC_prev") %>% filter(Beh_State == "Roosting")
int.Snow.Loaf <- interactions.raw %>% filter(Weath_Cov == "SD") %>% filter(Beh_State == "Loafing") %>% distinct()
int.Wind.Loaf <- interactions.raw %>% filter(Weath_Cov == "WC_prev") %>% filter(Beh_State == "Loafing") %>% distinct()
int.Snow.Forage <- interactions.raw %>% filter(Weath_Cov == "SD") %>% filter(Beh_State == "Foraging") %>% distinct()
int.Wind.Forage <- interactions.raw %>% filter(Weath_Cov == "WC_prev") %>% filter(Beh_State == "Foraging") %>% distinct()
int.Wind <- interactions.raw %>% filter(Weath_Cov == "WC_prev")
int.Snow <- interactions.raw %>% filter(Weath_Cov == "SD")
# Condition Thresholds (Used summary on raw data and chose near 1st/3rd Quantile and Mean)
# Wind Chill/Roost = 4, 15, 27
# Snow Depth/Roost = 0, 4, 8

# #This is all plots as one.
# i = 1
# snow.list <- list()
# snow.plots <- list()
# for(i in 1:length(int.Snow$LC_Cov)){
#   snow.df <- data.frame(Behavior = int.Snow$Beh_State[i],
#                         LC = int.Snow$LC_Cov[i],
#                         LC.Coef = int.Snow$LC_Coef[i],
#                         W.Coef = int.Snow$Weath_Coef[i],
#                         Int.Coef = int.Snow$Interaction[i],
#                         LC.Val = rep(seq(-2, 2,.2),3),
#                         W.Val = rep(c(0,4,8), each = 21),
#                         W.Condition = rep(c("Good","Average","Poor"), each = 21))
#   snow.list[[i]] <- snow.df %>%
#     mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
#   snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
#     geom_line(aes(linetype = W.Condition), size = 1.4) +
#     scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
#     theme_classic() +
#     xlab(snow.df$LC[i]) + ylab("")
#   
#   snow.plots[[i]] <- snow.plot + theme(legend.position="none")
# }
# 
# legend <- get_legend(snow.plot + theme(legend.position = "bottom"))
# plot_grid(plotlist = snow.plots,
#           legend,
#           labels = "auto",
#           nrow = 3,
#           align = "hv",
#           axis = "lb")

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
                        W.Val = rep(c(0,4,8), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.4) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 25) +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))

jpeg('Results/Roost_Snow.jpg', width = 1300, height = 1500)
plot_grid(plotlist = snow.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
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
                        W.Val = rep(c(27,15,4), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  wind.list[[i]] <- wind.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  wind.plot <- ggplot(data = wind.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.4) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 25) +
    xlab(wind.df$LC[i]) + ylab("")
  
  wind.plots[[i]] <- wind.plot + theme(legend.position="none")
}

legend <- get_legend(wind.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))


jpeg('Results/Roost_Wind.jpg', width = 1300, height = 1500)
plot_grid(plotlist = wind.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
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
                        W.Val = rep(c(0,4,8), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.4) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 25) +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))

jpeg('Results/Loafing_Snow.jpg', width = 1800, height = 1500)
plot_grid(plotlist = snow.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
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
                        W.Val = rep(c(27,15,4), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  wind.list[[i]] <- wind.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  wind.plot <- ggplot(data = wind.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.4) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 25) +
    xlab(wind.df$LC[i]) + ylab("")
  
  wind.plots[[i]] <- wind.plot + theme(legend.position="none")
}

legend <- get_legend(wind.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))

jpeg('Results/Loafing_Wind.jpg', width = 1800, height = 1500)
plot_grid(plotlist = wind.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
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
                        W.Val = rep(c(0,4,8), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.4) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 25) +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))


jpeg('Results/Foraging_Snow.jpg', width = 1800, height = 1500)
plot_grid(plotlist = snow.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
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
                        W.Val = rep(c(27,15,4), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  wind.list[[i]] <- wind.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  wind.plot <- ggplot(data = wind.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition), size = 1.4) +
    scale_linetype_manual(values = c("solid", 'dotdash', "dotted")) +
    theme_classic(base_size = 25) +
    xlab(wind.df$LC[i]) + ylab("")
  
  wind.plots[[i]] <- wind.plot + theme(legend.position="none")
}

legend <- get_legend(wind.plot + theme(legend.position = "right", legend.key.width=unit(1,"inch")))


jpeg('Results/Foraging_Wind.jpg', width = 1800, height = 1500)
plot_grid(plotlist = wind.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
dev.off()

#### Momentuhmm Transition Probability graphs

