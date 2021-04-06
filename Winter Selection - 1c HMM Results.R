#Final Model
hmm.top.model <- turk_m2.zm

require(momentuHMM)
load("hmmtopmodel.RData")
hmm.top.model

#Associate behavioral states from best model with original data 
turk_top_states <- viterbi(hmm.top.model)
turkey_states <- turkeyData.zm
turkey_states$State <- turk_top_states

#Diagnostic Plots
plotPR(hmm.top.model)
hist(turk_top_states)
plot(hmm.top.model,legend.pos="right")

save(hmm.top.model, file = "hmmtopmodel.RData")

#Save HMM Outputs as CSV
write.csv(hmm.top.model$rawCovs, 'Results/HMM - RawCovs.csv')
write.csv(hmm.top.model$mle$beta, 'Results/HMM - MLE of betas.csv')
write.csv(hmm.top.model$mle$angle, 'Results/HMM - MLE of angle.csv')
write.csv(hmm.top.model$mle$step, 'Results/HMM - MLE of step.csv')

#Transition probabilities with Confidence Intervals
require(dplyr)
Beta_TP_beta_lower <- as.data.frame(hmm.top.model$CIbeta$beta$lower) %>%
  mutate(CI = "lower") %>%
  mutate(Type = "beta")
Beta_TP_beta_lower$ID <- row.names(Beta_TP_beta_lower)
row.names(Beta_TP_beta_lower) <- c()
Beta_TP_beta_upper <- as.data.frame(hmm.top.model$CIbeta$beta$upper)  %>%
  mutate(CI = "upper") %>%
  mutate(Type = "beta")
Beta_TP_beta_upper$ID <- row.names(Beta_TP_beta_upper)
row.names(Beta_TP_beta_upper) <- c()
Beta_TP_est_beta <- as.data.frame(hmm.top.model$mle$beta)  %>%
  mutate(CI = "estimate") %>%
  mutate(Type = "beta")
Beta_TP_est_beta$ID <- row.names(Beta_TP_est_beta)
row.names(Beta_TP_est_beta) <- c()

require(tidyr)
TP_beta_w_CI <- rbind(Beta_TP_beta_lower, Beta_TP_beta_upper, Beta_TP_est_beta)
colnames(TP_beta_w_CI) <- c("TP1_2", "TP1_3", "TP2_1", "TP2_3", "TP3_1", "TP3_2", "CI", "Type", "ID")
TP_beta_w_CI_wide <- TP_beta_w_CI %>%
  select(-Type) %>%
  mutate(CI = factor(CI, levels = c("estimate", "lower", "upper"))) %>%
  group_by(ID) %>%
  arrange(CI) %>%
  pivot_wider(names_from = CI, values_from = c(TP1_2, TP1_3, TP2_1, TP2_3, TP3_1, TP3_2))%>% 
  ungroup()
write.csv(TP_beta_w_CI_wide, "Results/HMM - TP_with_CI.csv", row.names = F)

# 
# Beta_TP_real_lower <- as.data.frame(hmm.top.model$CIbeta$real$lower) %>%
#   mutate(CI = "lower") %>%
#   mutate(Type = "real")
# Beta_TP_real_upper <- as.data.frame(hmm.top.model$CIbeta$real$upper)  %>%
#   mutate(CI = "upper") %>%
#   mutate(Type = "real")



#Write output to csv to bring in, simplify, and combine with SSF data
write.csv(turkey_states, "Results/HMMBehavioralStates_output.csv")

#Checkout Leaflet

###########################################################################
########################
### CHECKING OUTPUTS ###
########################
# #Where are the cases localized temporally?
# state_times <- turkey_states %>%
#   mutate(State = as.factor(State)) %>%
#   group_by(State, hour) %>%
#   summarize(Total = n())
# 
# require(tidyr)
# state_days <- turkey_states %>%
#   mutate(State = as.factor(State)) %>%
#   group_by(ID, YDay, State) %>%
#   summarize(Total = n()) %>%
#   group_by(YDay, State) %>%
#   summarize(Avg = mean(Total)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = "State", values_from = "Avg") %>%
#   rename(Roosting = `1`, Loafing = `2`, Foraging = `3`) %>%
#   mutate(Ratio = Loafing/Foraging)%>%
#   mutate(L_propday = Loafing/(24-Roosting))%>%
#   mutate(F_propday = Foraging/(24-Roosting))
# 
# 
# ggplot(data = state_times, aes(x = hour, y = Total, group = as.factor(State))) +
#   geom_point(aes(shape = State, color = as.factor(State), size = 3))
# 
# ggplot(data = state_days, aes(x = YDay)) +
#   geom_point(aes(y = L_propday, shape = '21', size = 3)) +
#   geom_point(aes(y = F_propday, shape = '25', size = 3)) +
#   xlim(0,100)


# turkey_states1 <- turkey_states %>%
#   mutate(Date = as.Date(Timestamp)) %>%
#   #filter(hour(Timestamp) %in% ) #Maybe consider filtering out nightime hours but will need to use suncalc
#   group_by(ID, Date) %>%
#   summarize(State_Day_Avg = mean(State))
# 
# DailyStates <- ggplot(data=turkey_states1, aes(x = Date, y=State_Day_Avg)) +
#   geom_point(size = 1.5) +
#   xlab("Date") +
#   ylab("Daily Average State") +
#   theme_grey(base_size = 18) +
#   geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1) + #best fit line
#   theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
# DailyStates
# 
# id_list = unique(turkey_states$ID)
# id = id_list[3]
# states_onebird <- turkey_states1 %>% filter( ID == id)
# graph_onebird <- ggplot(data = states_onebird, aes(x = Date, y = State_Day_Avg)) +
#   geom_point(size = 1.5) +
#   labs(x = "Date", y = "State", title = id) +
#   theme_grey(base_size = 18) +
#   geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1)  #best fit line
# graph_onebird  


#Need to figure out a way to get the State values onto the points
#Need to figure out a way to incorproate individual variation into transition/mean etc.
#induce some sort of variation in the crawl wrap