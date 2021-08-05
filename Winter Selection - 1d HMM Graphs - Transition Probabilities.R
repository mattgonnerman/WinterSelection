require(ggplot2)
require(dplyr)
require(tidyverse)
require(momentuHMM)
load("hmmtopmodel.RData")
m <- hmm.top.model


##################################################
## Plot the t.p. as functions of the covariates ##
##################################################
#Load MLE of Transition Porbability betas
TPbetas <- read.csv('Results/HMM - MLE of betas.csv')
row.names(TPbetas) <- TPbetas$X
TPbetas <- TPbetas %>% select(-X)
TPbetas <- as.data.frame(t(TPbetas)) %>%
  dplyr::select(Intercept.coef = '(Intercept)', WC.coef = WC.Z, 
                SD.coef = SD.Z, Cos.coef = 'cosinorCos(hour, period = 24)', Sin.coef = 'cosinorSin(hour, period = 24)') %>%
  mutate(Start_State = c(1,1,2,2,3,3) ) %>%
  mutate(End_State = c(2,3,1,3,1,2) )
TPLCL <- as.data.frame(t(hmm.top.model$CIbeta$beta$lower)) %>%
  dplyr::select(Intercept.coef = '(Intercept)', WC.coef = WC.Z, 
                SD.coef = SD.Z, Cos.coef = 'cosinorCos(hour, period = 24)', Sin.coef = 'cosinorSin(hour, period = 24)') %>%
  mutate(Start_State = c(1,1,2,2,3,3) ) %>%
  mutate(End_State = c(2,3,1,3,1,2) )
TPUCL <- as.data.frame(t(hmm.top.model$CIbeta$beta$upper)) %>%
  dplyr::select(Intercept.coef = '(Intercept)', WC.coef = WC.Z, 
                SD.coef = SD.Z, Cos.coef = 'cosinorCos(hour, period = 24)', Sin.coef = 'cosinorSin(hour, period = 24)') %>%
  mutate(Start_State = c(1,1,2,2,3,3) ) %>%
  mutate(End_State = c(2,3,1,3,1,2) )

TPdata <- expand.grid(
  WC.cov = c(0,seq(min(m$data$WC.Z), max(m$data$WC.Z), .1)),
  SD.cov = c(0,seq(min(m$data$SD.Z), max(m$data$SD.Z), .1)),
  Hour.cov = c(5, 12, 15, 19),
  Start_State = 1:3,
  End_State = 1:3) %>%
  filter(Start_State != End_State)
TP.full <- merge(TPdata, TPbetas, by = c("Start_State", "End_State"), all.x = T) %>%
  mutate(TP_reg = Intercept.coef + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*Hour.cov/24) + Sin.coef*sin(2*pi*Hour.cov/24)) %>%
  mutate(TP_exp = exp(TP_reg))
TP.full.LCL <- merge(TPdata, TPLCL, by = c("Start_State", "End_State"), all.x = T) %>%
  mutate(TP_reg = Intercept.coef + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*Hour.cov/24) + Sin.coef*sin(2*pi*Hour.cov/24)) %>%
  mutate(TP_exp = exp(TP_reg))
TP.full.UCL <- merge(TPdata, TPUCL, by = c("Start_State", "End_State"), all.x = T) %>%
  mutate(TP_reg = Intercept.coef + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*Hour.cov/24) + Sin.coef*sin(2*pi*Hour.cov/24)) %>%
  mutate(TP_exp = exp(TP_reg))

TP1 <- TP.full %>%
  filter(Start_State == 1) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp12 = '1_2', exp13 = '1_3') %>%
  mutate(TP1_2 = exp12/(1+exp12+exp13))%>%
  mutate(TP1_3 = exp13/(1+exp12+exp13))%>%
  mutate(TP1_1 = 1 - TP1_2 - TP1_3)

TP2 <- TP.full %>%
  filter(Start_State == 2) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp21 = '2_1', exp23 = '2_3') %>%
  mutate(TP2_1 = exp21/(1+exp21+exp23))%>%
  mutate(TP2_3 = exp23/(1+exp21+exp23))%>%
  mutate(TP2_2 = 1 - TP2_1 - TP2_3)

TP3 <- TP.full %>%
  filter(Start_State == 3) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp31 = '3_1', exp32 = '3_2') %>%
  mutate(TP3_1 = exp31/(1+exp31+exp32))%>%
  mutate(TP3_2 = exp32/(1+exp31+exp32))%>%
  mutate(TP3_3 = 1 - TP3_1 - TP3_2)

TP1.lcl <- TP.full.LCL %>%
  filter(Start_State == 1) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp12 = '1_2', exp13 = '1_3') %>%
  mutate(TP1_2 = exp12/(1+exp12+exp13))%>%
  mutate(TP1_3 = exp13/(1+exp12+exp13))%>%
  mutate(TP1_1 = 1 - TP1_2 - TP1_3)

TP2.lcl <- TP.full.LCL %>%
  filter(Start_State == 2) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp21 = '2_1', exp23 = '2_3') %>%
  mutate(TP2_1 = exp21/(1+exp21+exp23))%>%
  mutate(TP2_3 = exp23/(1+exp21+exp23))%>%
  mutate(TP2_2 = 1 - TP2_1 - TP2_3)

TP3.lcl <- TP.full.LCL %>%
  filter(Start_State == 3) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp31 = '3_1', exp32 = '3_2') %>%
  mutate(TP3_1 = exp31/(1+exp31+exp32))%>%
  mutate(TP3_2 = exp32/(1+exp31+exp32))%>%
  mutate(TP3_3 = 1 - TP3_1 - TP3_2)

TP1.ucl <- TP.full.UCL %>%
  filter(Start_State == 1) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp12 = '1_2', exp13 = '1_3') %>%
  mutate(TP1_2 = exp12/(1+exp12+exp13))%>%
  mutate(TP1_3 = exp13/(1+exp12+exp13))%>%
  mutate(TP1_1 = 1 - TP1_2 - TP1_3)

TP2.ucl <- TP.full.UCL %>%
  filter(Start_State == 2) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp21 = '2_1', exp23 = '2_3') %>%
  mutate(TP2_1 = exp21/(1+exp21+exp23))%>%
  mutate(TP2_3 = exp23/(1+exp21+exp23))%>%
  mutate(TP2_2 = 1 - TP2_1 - TP2_3)

TP3.ucl <- TP.full.UCL %>%
  filter(Start_State == 3) %>%
  select(Start_State, End_State, WC.cov, SD.cov, Hour.cov, Start_State, TP_exp) %>%
  pivot_wider(names_from = c(Start_State, End_State), values_from = TP_exp) %>%
  rename(exp31 = '3_1', exp32 = '3_2') %>%
  mutate(TP3_1 = exp31/(1+exp31+exp32))%>%
  mutate(TP3_2 = exp32/(1+exp31+exp32))%>%
  mutate(TP3_3 = 1 - TP3_1 - TP3_2)

TP1.CL <- merge(TP1.lcl, TP1.ucl, by = c("WC.cov", "SD.cov", "Hour.cov"))
TP2.CL <- merge(TP2.lcl, TP2.ucl, by = c("WC.cov", "SD.cov", "Hour.cov"))
TP3.CL <- merge(TP3.lcl, TP3.ucl, by = c("WC.cov", "SD.cov", "Hour.cov"))

### Make Plots and then format as grid using cowplot
## 2 --> ? | Wind Chill
TP2.WC <- TP2 %>%
  filter(SD.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
TP2.WC.CL <- TP2.CL %>%
  filter(SD.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
WC2 <- ggplot(TP2.WC, aes(x = WC.cov)) +
  geom_line(aes(y = TP2_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP2_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP2_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP2.WC.CL, aes(ymin = TP2_1.x, ymax = TP2_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP2.WC.CL, aes(ymin = TP2_2.x, ymax = TP2_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP2.WC.CL, aes(ymin = TP2_3.x, ymax = TP2_3.y),alpha = .3, fill = "#46e300") +
  theme_classic(base_size = 50) +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend")  +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") +
  xlab(element_blank()) +
  ylab(paste("P(Stationary ", sprintf("\u2192"), " X)"))

## 3 --> ? | Wind Chill
TP3.WC <- TP3 %>%
  filter(SD.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
TP3.WC.CL <- TP3.CL %>%
  filter(SD.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
WC3 <- ggplot(TP3.WC, aes(x = WC.cov)) +
  geom_line(aes(y = TP3_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP3_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP3_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP3.WC.CL, aes(ymin = TP3_1.x, ymax = TP3_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP3.WC.CL, aes(ymin = TP3_2.x, ymax = TP3_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP3.WC.CL, aes(ymin = TP3_3.x, ymax = TP3_3.y),alpha = .3, fill = "#46e300") +
  theme_classic(base_size = 50) +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend")  +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") +
  xlab("Wind Chill") +
  ylab(paste("P(Mobile ", sprintf("\u2192"), " X)"))

## 2 --> ? | Snow Depth
TP2.SD <- TP2 %>%
  filter(WC.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
TP2.SD.CL <- TP2.CL %>%
  filter(WC.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
SD2 <- ggplot(TP2.SD, aes(x = SD.cov)) +
  geom_line(aes(y = TP2_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP2_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP2_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP2.SD.CL, aes(ymin = TP2_1.x, ymax = TP2_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP2.SD.CL, aes(ymin = TP2_2.x, ymax = TP2_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP2.SD.CL, aes(ymin = TP2_3.x, ymax = TP2_3.y),alpha = .3, fill = "#46e300") +
  theme_classic(base_size = 50) +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend")  +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") +
  xlab(element_blank()) +
  xlim(-1, 2.5) +
  ylab(element_blank())

## 3 --> ? | Snow Depth
TP3.SD <- TP3 %>%
  filter(WC.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
TP3.SD.CL <- TP3.CL %>%
  filter(WC.cov == 0) %>%
  filter(Hour.cov == 12) %>%
  mutate(Hour.cov = as.factor(Hour.cov))
SD3 <- ggplot(TP3.SD, aes(x = SD.cov)) +
  geom_line(aes(y = TP3_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP3_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP3_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP3.SD.CL, aes(ymin = TP3_1.x, ymax = TP3_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP3.SD.CL, aes(ymin = TP3_2.x, ymax = TP3_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP3.SD.CL, aes(ymin = TP3_3.x, ymax = TP3_3.y),alpha = .3, fill = "#46e300") +
  theme_classic(base_size = 50) +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend")  +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") +
  xlab("Snow Depth") +
  xlim(-1, 2.5) +
  ylab(element_blank())



### Time Budget Component
tb.data <- read.csv("HMMBehavioralStates_output.csv") %>%
  dplyr::select(State, WC.Z, SD.Z) %>%
  mutate(State = factor(State, levels = c(1,2,3), labels = c("Roosting", "Stationary", "Mobile")))

tb.WC <- ggplot(data = tb.data, aes(x = State, y = WC.Z)) +
  geom_boxplot(aes(fill = State), alpha = .7, lwd = 1.5) +
  theme_classic(base_size = 50) +
  ylim(-2.5,2.5) +
  scale_fill_manual(values = c("#1cade4", "#f1a806", "#46e300")) + 
  theme(legend.position = "none") +
  labs(y = "Wind Chill", y = element_blank())

tb.SD <- ggplot(data = tb.data, aes(x = State, y = SD.Z)) +
  geom_boxplot(aes(fill = State), alpha = .7, outlier.shape = NA, lwd = 1.5) +
  theme_classic(base_size = 50) +
  ylim(-2.5,2.5) +
  scale_fill_manual(values = c("#1cade4", "#f1a806", "#46e300")) + 
  theme(legend.position = "none") +
  labs(y = "Snow Depth", x = element_blank())

tb.WC1 <- ggplot(data = tb.data, aes(x = WC.Z, group = State)) +
  geom_density(aes(fill = State, color = State), alpha = .3, lwd = 1.5) +
  theme_classic(base_size = 50) +
  scale_fill_manual(values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_color_manual(values = c("#1cade4", "#f1a806", "#46e300")) +
  theme(legend.position = "none") +
  labs(y = "Density", x = element_blank())

tb.SD1 <- ggplot(data = tb.data, aes(x = SD.Z, group = State)) +
  geom_density(aes(fill = State, color = State), alpha = .3, lwd = 1.5) +
  theme_classic(base_size = 50) +
  scale_fill_manual(values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_color_manual(values = c("#1cade4", "#f1a806", "#46e300")) +
  theme(legend.position = "none") +
  xlim(-1, 2.5) +
  labs(y = "Density", x = element_blank())

#Grid Plot
require(cowplot)
legend <- get_legend(SD3 + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(3,"inch")))
WC2.graph <- WC2 + theme(legend.position = "none")
WC3.graph <- WC3 + theme(legend.position = "none")
SD2.graph <- SD2 + theme(legend.position = "none")
SD3.graph <- SD3 + theme(legend.position = "none")

TP_grid <- plot_grid(plotlist = list(tb.WC1, tb.SD1,
                                     WC2.graph, SD2.graph,
                                     WC3.graph, SD3.graph),
                     nrow = 3,
                     labels = "AUTO",
                     label_size = 45,
                     align = "hv",
                     axis = "lb"
)



jpeg('Results/Transition Probability Grid.png', width = 2600, height = 2500)
plot_grid(plotlist = list(TP_grid, legend),
          nrow = 2, 
          rel_heights = c(1,.05)
)
dev.off()
