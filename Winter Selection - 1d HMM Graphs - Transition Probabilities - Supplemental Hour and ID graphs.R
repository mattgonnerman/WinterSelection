require(ggplot2)
require(dplyr)
require(tidyverse)
require(momentuHMM)
load("hmmtopmodel.RData")
hmm.top.model

###############################################
## Plot the t.p. as functions of Hour of Day ##
###############################################
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
  WC.cov = 0,
  SD.cov = 0,
  Hour.cov = seq(0,24, .1),
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
## Stationary
hour.stationary <- ggplot(TP2, aes(x = Hour.cov)) +
  geom_line(aes(y = TP2_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP2_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP2_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP2.CL, aes(ymin = TP2_1.x, ymax = TP2_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP2.CL, aes(ymin = TP2_2.x, ymax = TP2_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP2.CL, aes(ymin = TP2_3.x, ymax = TP2_3.y),alpha = .3, fill = "#46e300") +
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
  ylab(paste("P(Stationary ", sprintf("\u2192"), " X)")) +
  scale_x_continuous(breaks = (seq(0,24,4)),
                     labels = c(20, seq(0,20,4)))

## Mobile
hour.mobile <- ggplot(TP3, aes(x = Hour.cov)) +
  geom_line(aes(y = TP3_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP3_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP3_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP3.CL, aes(ymin = TP3_1.x, ymax = TP3_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP3.CL, aes(ymin = TP3_2.x, ymax = TP3_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP3.CL, aes(ymin = TP3_3.x, ymax = TP3_3.y),alpha = .3, fill = "#46e300") +
  theme_classic(base_size = 50) +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend")  +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") +
  xlab("Hour") +
  ylab(paste("P(Mobile ", sprintf("\u2192"), " X)")) +
  scale_x_continuous(breaks = (seq(0,24,4)),
                     labels = c(20, seq(0,20,4)))

##Roosting
hour.roosting <- ggplot(TP1, aes(x = Hour.cov)) +
  geom_line(aes(y = TP1_1, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = TP1_2, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = TP1_3, color = "#46e300", linetype = "twodash"), size = 3) +
  geom_ribbon(data = TP1.CL, aes(ymin = TP1_1.x, ymax = TP1_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = TP1.CL, aes(ymin = TP1_2.x, ymax = TP1_2.y),alpha = .3, fill = "#f1a806") +
  geom_ribbon(data = TP1.CL, aes(ymin = TP1_3.x, ymax = TP1_3.y),alpha = .3, fill = "#46e300") +
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
  ylab(paste("P(Roosting ", sprintf("\u2192"), " X)"))  +
  scale_x_continuous(breaks = (seq(0,24,4)),
                     labels = c(20, seq(0,20,4)))


#Grid Plot
require(cowplot)
legend <- get_legend(hour.roosting + theme(legend.title = element_blank() ,legend.position = "bottom", legend.key.width=unit(2.4,"inch")))
mobile.graph <- hour.mobile + theme(legend.position = "none")
stationary.graph <- hour.stationary + theme(legend.position = "none")
roost.graph <- hour.roosting + theme(legend.position = "none")

TP_grid <- plot_grid(plotlist = list(roost.graph, 
                                     stationary.graph, 
                                     mobile.graph),
                     nrow = 3,
                     # labels = "auto",
                     # label_size = 35,
                     align = "hv",
                     axis = "lb"
)



jpeg('Results/Transition Probability Grid - Hour.png', 
     units = "in", width = 17, height = 30, res = 300)
plot_grid(plotlist = list(TP_grid, legend),
          nrow = 2, 
          rel_heights = c(1,.1)
)
dev.off()


###############################################
## Plot the t.p. as functions of Hour of Day ##
###############################################
ind.covs.est <- read.csv('Results/HMM - MLE of betas.csv') %>% 
  filter(substr(X, 1,3) == "IDX") %>%
  pivot_longer(names_to = "TransStates", cols = 2:7) %>%
  mutate(SS = substr(TransStates, 2,2),
         ES = substr(TransStates, 7,7)) %>%
  rename(TP = value, ID = X) %>% dplyr::select(-TransStates) %>%
  mutate(ID = as.factor(ID))

ind.covs.lcl <- as.data.frame(hmm.top.model$CIbeta$beta$lower) %>% 
  mutate(X = row.names(.)) %>%
  filter(substr(X, 1,3) == "IDX") %>%
  pivot_longer(names_to = "TransStates", cols = 1:6) %>%
  mutate(SS = substr(TransStates, 1,1),
         ES = substr(TransStates, 6,6)) %>%
  rename(LCL = value, ID = X) %>% dplyr::select(-TransStates) %>%
  mutate(ID = as.factor(ID))

ind.covs.ucl <- as.data.frame(hmm.top.model$CIbeta$beta$upper) %>% 
  mutate(X = row.names(.)) %>%
  filter(substr(X, 1,3) == "IDX") %>%
  pivot_longer(names_to = "TransStates", cols = 1:6) %>%
  mutate(SS = substr(TransStates, 1,1),
         ES = substr(TransStates, 6,6)) %>%
  rename(UCL = value, ID = X) %>% dplyr::select(-TransStates) %>%
  mutate(ID = as.factor(ID))

ind.covs <- merge(ind.covs.est,ind.covs.lcl, by = c("ID", "SS", "ES")) %>%
  merge(., ind.covs.ucl, by = c("ID", "SS", "ES"))

ind.12 <- ggplot(ind.covs %>% filter(SS == 1 & ES == 2), aes(x = ID)) +
  geom_point(aes(y = TP), size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), size = 1, width = .1) +
  theme_classic(base_size = 50) +
  ylab(paste("Beta(Roosting ", sprintf("\u2192"), " Stationary)"))  +
  theme(axis.text.x = element_blank())

ind.13 <- ggplot(ind.covs %>% filter(SS == 1 & ES == 3), aes(x = ID)) +
  geom_point(aes(y = TP), size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), size = 1, width = .1) +
  theme_classic(base_size = 50) +
  ylab(paste("Beta(Roosting ", sprintf("\u2192"), " Mobile)"))  +
  theme(axis.text.x = element_blank())

ind.21 <- ggplot(ind.covs %>% filter(SS == 2 & ES == 1), aes(x = ID)) +
  geom_point(aes(y = TP), size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), size = 1, width = .1) +
  theme_classic(base_size = 50) +
  ylab(paste("Beta(Stationary ", sprintf("\u2192"), " Roosting)"))  +
  theme(axis.text.x = element_blank())

ind.23 <- ggplot(ind.covs %>% filter(SS == 2 & ES == 3), aes(x = ID)) +
  geom_point(aes(y = TP), size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), size = 1, width = .1) +
  theme_classic(base_size = 50) +
  ylab(paste("Beta(Stationary ", sprintf("\u2192"), " Mobile)"))  +
  theme(axis.text.x = element_blank())

ind.32 <- ggplot(ind.covs %>% filter(SS == 3 & ES == 2), aes(x = ID)) +
  geom_point(aes(y = TP), size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), size = 1, width = .1) +
  theme_classic(base_size = 50) +
  ylab(paste("Beta(Mobile ", sprintf("\u2192"), " Stationary)"))  +
  theme(axis.text.x = element_blank())

ind.31 <- ggplot(ind.covs %>% filter(SS == 3 & ES == 1), aes(x = ID)) +
  geom_point(aes(y = TP), size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), size = 1, width = .1) +
  theme_classic(base_size = 50) +
  ylab(paste("Beta(Mobile ", sprintf("\u2192"), " Roosting)"))  +
  theme(axis.text.x = element_blank())

ind_grid <- plot_grid(plotlist = list(ind.12,ind.13,ind.21,ind.23,ind.31,ind.32),
                     nrow = 3,
                     # labels = "auto",
                     # label_size = 35,
                     align = "hv",
                     axis = "lb"
)

jpeg('Results/Transition Probability Grid - Ind.png', 
     units = "in", width = 25, height = 32, res = 600)
ind_grid
dev.off()
