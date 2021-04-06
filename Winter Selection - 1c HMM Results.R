require(momentuHMM)
load("hmmtopmodel.RData")
hmm.top.model

#Associate behavioral states from best model with original data 
turk_top_states <- viterbi(hmm.top.model)
turkey_states <- turkeyData.zm
turkey_states$State <- turk_top_states

#Write output to csv to bring in, simplify, and combine with SSF data
write.csv(turkey_states, "Results/HMMBehavioralStates_output.csv")

#Diagnostic Plots
plotPR(hmm.top.model)
hist(turk_top_states)
plot(hmm.top.model,legend.pos="right")

#Save Raw HMM Outputs as CSV
write.csv(hmm.top.model$rawCovs, 'Results/HMM - RawCovs.csv')
write.csv(hmm.top.model$mle$beta, 'Results/HMM - MLE of betas.csv')
write.csv(hmm.top.model$mle$angle, 'Results/HMM - MLE of angle.csv')
write.csv(hmm.top.model$mle$step, 'Results/HMM - MLE of step.csv')

#Save Betas for Transition probabilities with Confidence Intervals
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





#Step Length Estimates and CI
step.est <- as.data.frame(hmm.top.model$CIreal$step$est) %>%
  mutate(CI = "Estimate")
step.est$Type <- row.names(step.est)
step.lower <- as.data.frame(hmm.top.model$CIreal$step$lower) %>%
  mutate(CI = "Lower")
step.lower$Type <- row.names(step.est)
step.upper <- as.data.frame(hmm.top.model$CIreal$step$upper) %>%
  mutate(CI = "Upper")
step.upper$Type <- row.names(step.est)

step.all <- rbind(step.est, step.lower, step.upper) %>%
  rename(Roosting = roost, Stationary = loafing, Mobile = foraging) %>%
  group_by(Type) %>%
  pivot_wider(names_from = CI, values_from = c(Roosting, Stationary, Mobile))
write.csv(step.all, "Results/HMM - Step Length with CI.csv", row.names = F)

#Angle Concentration Estimates and CI
angle.est <- as.data.frame(hmm.top.model$CIreal$angle$est) %>%
  mutate(CI = "Estimate")
angle.est$Type <- row.names(angle.est)
angle.lower <- as.data.frame(hmm.top.model$CIreal$angle$lower) %>%
  mutate(CI = "Lower")
angle.lower$Type <- row.names(angle.est)
angle.upper <- as.data.frame(hmm.top.model$CIreal$angle$upper) %>%
  mutate(CI = "Upper")
angle.upper$Type <- row.names(angle.est)

angle.all <- rbind(angle.est, angle.lower, angle.upper) %>%
  rename(Roosting = roost, Stationary = loafing, Mobile = foraging)# %>%
  # group_by(Type) %>%
  # pivot_wider(names_from = CI, values_from = c(Roosting, Stationary, Mobile))
write.csv(angle.all, "Results/HMM - Angle Concentration with CI.csv", row.names = F)




#Time Spent in each State
require(lubridate)
turkey_states <-  read.csv("Results/HMMBehavioralStates_output.csv")
time_spent <- turkey_states %>%
  select(ID, step, angle, Timestamp,WC.Z, SD.Z, State) %>%
  mutate(Timestamp_UTC = as.POSIXct(Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  mutate(Timestamp = with_tz(Timestamp_UTC, tzone = "America/New_York")) %>%
  mutate(hour = hour(Timestamp)) %>%
  mutate(OrdinalDay = yday(Timestamp)) %>%
  mutate(Date = as.Date(Timestamp)) %>%
  mutate(Year = year(Timestamp)) %>%
  mutate(State1 = ifelse(State == 1, 1, 0))  %>%
  mutate(State2 = ifelse(State == 2, 1, 0)) %>%
  mutate(State3 = ifelse(State == 3, 1, 0))

#By hour of day
time_spent_hour <- time_spent %>%
  group_by(hour) %>%
  summarize(PropState1 = mean(State1), 
            PropState2 = mean(State2),
            PropState3 = mean(State3)) %>%
  pivot_longer(2:4, names_to = "State")

ggplot(time_spent_hour, aes(x = hour, y = value, group = State)) +
  geom_line(aes(color = State))

#By ordinal Day
time_spent_ordinalday <- time_spent %>%
  group_by(Year, OrdinalDay) %>%
  summarize(PropState1 = mean(State1), 
            PropState2 = mean(State2),
            PropState3 = mean(State3)) %>%
  pivot_longer(3:5, names_to = "State") %>%
  filter(OrdinalDay < 300)

ggplot(time_spent_ordinalday[time_spent_ordinalday$Year == 2019,], aes(x = OrdinalDay, y = value, group = State)) +
  geom_line(aes(color = State)) +
  geom_line(data = weather.raw, aes(y = WC, group = State)) +
  scale_y_continuous("Probability in State", sec.axis = sec_axis(trans=~.*coeff, name = "Wind Chill"))

#By Weather covariates
weather.raw <- read.csv("WeatherVariables.csv") %>%
  mutate(Date = as.Date(DATE, format = "%m/%d/%Y"))%>%
  mutate(OrdinalDay = yday(Date)) %>%
  filter(OrdinalDay < 300) %>%
  filter(Date %in% time_spent$Date) %>%
  mutate(State = "Weather") %>%
  mutate(AWND = ifelse(is.na(AWND), 0, AWND)) %>% #If AWND is NA, that means no wind
  dplyr::select(Date, OrdinalDay, State, TMIN, AWND, WDF2, SD = SNWD)%>%
  mutate(WC = 13.12 + (.6215*TMIN)-(11.37*(AWND^0.16))+(.3965*TMIN*(AWND^0.16))) #calculate windchill 


time_spent_weather <- merge(time_spent, weather.raw, by = "Date", all.x = T)










