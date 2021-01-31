require(dplyr)
require(ggplot2)
require(forcats)

### Comparison of Magnitude of Interactions
interactions.raw <- read.csv("InteractionResults.csv")
interactions.raw$Beh_State <- factor(interactions.raw$Beh_State,
                                     levels = c("Roost", "Stationary", "Mobile"))
interactions.raw$LC_Cov <- factor(interactions.raw$LC_Cov,
                                  levels = c("Wind Exposure", "Distance to Edge", "Proportion Ag",
                                             "Proportion Dev", "Proportion SW", "% Softwood",
                                             "Mean Tree Height", "Basal Area"))

int.Snow <- interactions.raw %>% filter(Weath_Cov == "Snow Depth")
ggplot(data = int.Snow, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_point(size = 1.5,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = Interaction - (1.96*SD), xmax = Interaction + (1.96*SD)),
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
                      labels = c("Roost", "Stationary", "Mobile"),
                      values = c("yellow4", "violetred4", "royalblue4")) +   
  scale_shape_manual(name = "Behavioral\nState",
                     labels = c("Roost", "Stationary", "Mobile"),
                     values = c(15, 19, 17))
ggsave("SnowDepth_InteractionComp.jpeg", width = 8, height = 7, units = "in")

int.Wind <- interactions.raw %>% filter(Weath_Cov == "Wind Chill")
ggplot(data = int.Wind, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_point(size = 1.5,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = Interaction - (1.96*SD), xmax = Interaction + (1.96*SD)),
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
                      labels = c("Roost", "Stationary", "Mobile"),
                      values = c("yellow4", "violetred4", "royalblue4")) +   
  scale_shape_manual(name = "Behavioral\nState",
                     labels = c("Roost", "Stationary", "Mobile"),
                     values = c(15, 19, 17))
ggsave("WindChill_InteractionComp.jpeg", width = 8, height = 7, units = "in")


#make big points
#remove endcaps on error bars
#thicker error lines

################################################################################################
### Plot Matrix showing selection at Poor, Average, and Good Weather
require(cowplot)

interactions.raw <- read.csv("InteractionResults.csv") %>%
  mutate(Beh_State = factor(Beh_State, levels = c("Roost", "Stationary", "Mobile"))) %>%
  mutate(LC_Cov = factor(LC_Cov,
                         levels = c("Distance to Edge", "Wind Exposure", "Proportion Ag",
                                    "Proportion Dev", "Proportion SW",
                                    "Mean Tree Height", "Basal Area", "% Softwood"))) %>%
  arrange(Beh_State, LC_Cov)

int.Snow <- interactions.raw %>% filter(Weath_Cov == "Snow Depth")
int.Wind <- interactions.raw %>% filter(Weath_Cov == "Wind Chill")
# Condition Thresholds (Used summary on raw data and chose near 1st/3rd Quantile and Mean)
# Wind Chill/Roost = 4, 15, 27
# Snow Depth/Roost = 0, 4, 8


i = 1
snow.list <- list()
snow.plots <- list()
for(i in 1:length(int.Snow$LC_Cov)){
  snow.df <- data.frame(Behavior = int.Snow$Beh_State[i],
                        LC = int.Snow$LC_Cov[i],
                        LC.Coef = int.Snow$LC_Coef[i],
                        W.Coef = int.Snow$Weath_Coef[i],
                        Int.Coef = int.Snow$Interaction[i],
                        LC.Val = rep(seq(-2, 2,.2),3),
                        W.Val = rep(c(0,4,8), each = 21),
                        W.Condition = rep(c("Good","Average","Poor"), each = 21))
  snow.list[[i]] <- snow.df %>%
    mutate(Est = exp((LC.Coef*LC.Val) + (W.Coef*W.Val) + (Int.Coef*LC.Val*W.Val)))
  snow.plot <- ggplot(data = snow.list[[i]], aes(x = LC.Val, y = Est, group = W.Condition)) +
    geom_line(aes(linetype = W.Condition)) +
    theme_classic() +
    xlab(snow.df$LC[i]) + ylab("")
  
  snow.plots[[i]] <- snow.plot + theme(legend.position="none")
}

legend <- get_legend(snow.plot + theme(legend.position = "bottom"))
plot_grid(plotlist = snow.plots,
          legend,
          labels = "auto",
          nrow = 3,
          align = "hv",
          axis = "lb")
