require(dplyr)
require(ggplot2)
require(forcats)

LCnames <- c("Dist. to Forest Edge", 
             "Basal Area", 
             "Mean Tree Height", 
             "Percent Softwood",
             "Wind Exposure", 
             "Agriculture", 
             "Developed", 
             "Food Subsidy", 
             "Proportion Softwood")

### Comparison of Magnitude of Interactions
interactions.raw <- read.csv("Results/InteractionResults.csv") %>% 
  mutate(HabitatCov = ifelse(HabitatCov == "DtFE.Z", "Dist. to Forest Edge",
                             ifelse(HabitatCov == "PropAg.Z", "Agriculture",
                                    ifelse(HabitatCov == "PropDev.Z", "Developed",
                                           ifelse(HabitatCov == "PropFoodSub.Z", "Food Subsidy",
                                                  ifelse(HabitatCov == "PropSW.Z", "Proportion Softwood",
                                                         ifelse(HabitatCov == "Wind.Exp.Z", "Wind Exposure",
                                                                ifelse(HabitatCov == "BA.Z", "Basal Area",
                                                                       ifelse(HabitatCov == "Ht.Z", "Mean Tree Height",
                                                                              ifelse(HabitatCov == "SW.Z", "Percent Softwood",
                                                                                     HabitatCov))))))))))%>%
  mutate(Beh_State = factor(Analysis, levels = c("Roosting", "Stationary", "Mobile"))) %>%
  mutate(LC_Cov = factor(HabitatCov, levels = LCnames)) %>%
  arrange(Beh_State, LC_Cov, WeatherCov) %>%
  rename(Interaction = mean, SD = sd)
# interactions.raw$Beh_State <- factor(interactions.raw$Beh_State,
#                                      levels = c("Roost", "Loafing", "Foraging"))
# interactions.raw$LC_Cov <- factor(interactions.raw$LC_Cov,
#                                   levels = c("Wind Exposure", "Distance to Edge", "Proportion Ag",
#                                              "Proportion Dev", "Proportion SW", "% Softwood",
#                                              "Mean Tree Height", "Basal Area"))

#Graph showing interaction terms for snow depth
int.Snow <- interactions.raw %>% filter(WeatherCov == "SD") %>%
  mutate(LC_Cov = fct_reorder(LC_Cov, rev(LCnames)))
int.snow.graph <- ggplot(data = int.Snow, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2, size = 1.5) +
  geom_point(size = 8,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant),
                width = 0, size = 2,
                position = position_dodge(width = .4)) +
  theme_classic(base_size = 50) + 
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
# int.snow.graph
# ggsave("Results/SnowDepth_InteractionComp.jpeg", width = 10, height = 12, units = "in")


#Graph showing interaction terms for Previous days wind chill
int.Wind <- interactions.raw %>% 
  mutate(WeatherCov = ifelse(WeatherCov == "WC_prev", "WC", WeatherCov)) %>%
  filter(WeatherCov == "WC") %>%
  mutate(LC_Cov = fct_reorder(LC_Cov, rev(LCnames)))
int.wind.graph <- ggplot(data = int.Wind, aes(y = LC_Cov, x = Interaction, shape = Beh_State, color = Beh_State)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2, size = 1.5) +
  geom_point(size = 8,
             position = position_dodge(width = .4)) +
  geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant),
                width = 0, size = 2,
                position = position_dodge(width = .4)) +
  theme_classic(base_size = 50) + 
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
                     values = c(15, 19, 17)) +
  theme(axis.text.y = element_blank())
# int.wind.graph
# ggsave("Results/WindChill_InteractionComp.jpeg", width = 10, height = 12, units = "in")

#Save as a combined graph
require(cowplot)
legend <- get_legend(int.wind.graph + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(3,"inch")))
int.snow.grid <- int.snow.graph + theme(legend.position = "none")
int.wind.grid <- int.wind.graph + theme(legend.position = "none")


# INT_grid <- plot_grid(plotlist = list(int.snow.grid, int.wind.grid),
#                      nrow = 1,
#                      # labels = "auto",
#                      # label_size = 35,
#                      align = "hv",
#                      axis = "lb"
# )

require(patchwork)
INT_grid <- int.snow.grid + int.wind.grid

jpeg('Results/SSF Interactions Grid.png', width = 2600, height = 1900)
plot_grid(plotlist = list(INT_grid, legend),
          nrow = 2, 
          rel_heights = c(1,.1)
)
dev.off()
