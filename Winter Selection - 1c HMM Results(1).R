# ### Model Selection
# # Finding best parameter restrictions
# hmm.aic.A <- AIC(turk_mA0.zm,
#                  turk_mA1.zm,
#                  turk_mA2.zm,
#                  turk_mA3.zm,
#                  turk_mA4.zm,
#                  turk_mA5.zm
#                )
# hmm.aic.A
# write.csv(hmm.aic.A, "HMM_AIC_A.csv", row.names = F)
# hmm.top.model <- get(hmm.aic.A[1,1])
# 
# #Best Model With Covariates
# hmm.aic.B <- AIC(get(hmm.aic.A[1,1]),
#                  
# )
# hmm.aic.B
# write.csv(hmm.aic.B, "HMM_AIC_B.csv", row.names = F)
# 
# 
# 
# hmm.aicC[1,1]
# hmm.top.model <- get(hmm.aicC[1,1])

hmm.top.model <- turk_m5.zm

write.csv(turk_m5.zm$rawCovs, 'HMM - RawCovs.csv')
write.csv(turk_m5.zm$mle$beta, 'HMM - MLE of betas.csv')
write.csv(turk_m5.zm$mle$angle, 'HMM - MLE of angle.csv')
write.csv(turk_m5.zm$mle$step, 'HMM - MLE of step.csv')

#Associate behavioral states from best model with original data 
turk_top_states <- viterbi(hmm.top.model)
turkey_states <- turkeyData.zm
turkey_states$State <- turk_top_states
#Diagnostic Plots
# plotPR(hmm.top.model)
# hist(turk_top_states)
# plot(hmm.top.model,legend.pos="right")

#Write output to csv to bring in, simplify, and combine with SSF data
write.csv(turkey_states, "HMMBehavioralStates_output.csv")


#Checkout Leaflet


###########################################################################
########################
### CHECKING OUTPUTS ###
########################
# Look at output maps/graphs
# plot(turk_m1.zm,legend.pos="right")
# hist(turk_m1_states)

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