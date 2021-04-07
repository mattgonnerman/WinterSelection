require(dplyr)
require(ggplot2)
require(forcats)

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
