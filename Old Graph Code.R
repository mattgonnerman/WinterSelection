{  
  ###############################################################################################
  # # 1->?
  # TP1_2 <- data.frame(Intercept = TPbetas[1,1],
  #                     WC.coef = TPbetas[1,2],
  #                     SD.coef = TPbetas[1,3],
  #                     Cos.coef = TPbetas[1,4],
  #                     Sin.coef = TPbetas[1,5],
  #                     hour = 9,
  #                     SD.cov = 0,
  #                     WC.cov = seq(min(hmm.top.model$data$WC.Z),max(hmm.top.model$data$WC.Z), .1)
  #                     ) %>%
  #   mutate(TP1_2_reg = Intercept + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*hour/24) + Sin.coef*sin(2*pi*hour/24)) %>%
  #   mutate(TP1_2_exp = exp(TP1_2_reg)) %>%
  #   select(WC.cov, TP1_2_reg, TP1_2_exp)
  # TP1_3 <- data.frame(Intercept = TPbetas[2,1],
  #                     WC.coef = TPbetas[2,2],
  #                     SD.coef = TPbetas[2,3],
  #                     Cos.coef = TPbetas[2,4],
  #                     Sin.coef = TPbetas[2,5],
  #                     hour = 9,
  #                     SD.cov = 0,
  #                     WC.cov = seq(min(hmm.top.model$data$WC.Z),max(hmm.top.model$data$WC.Z), .1)
  # ) %>%
  #   mutate(TP1_3_reg = Intercept + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*hour/24) + Sin.coef*sin(2*pi*hour/24)) %>%
  #   mutate(TP1_3_exp = exp(TP1_3_reg)) %>%
  #   select(TP1_3_reg, TP1_3_exp)
  # TP1 <- cbind(TP1_2, TP1_3) %>%
  #   mutate(TP1_2 = TP1_2_exp/(1+TP1_2_exp+TP1_3_exp)) %>%
  #   mutate(TP1_3 = TP1_3_exp/(1+TP1_2_exp+TP1_3_exp)) %>%
  #   mutate(TP1_1 = 1-TP1_3-TP1_2)
  # 
  # ggplot(TP1) +
  #   geom_line(aes(x = WC.cov, y = TP1_1)) +
  #   geom_line(aes(x = WC.cov, y = TP1_2)) +
  #   geom_line(aes(x = WC.cov, y = TP1_3))
  # 
  # # 2->?
  # TP2_1 <- data.frame(Intercept = TPbetas[3,1],
  #                     WC.coef = TPbetas[3,2],
  #                     SD.coef = TPbetas[3,3],
  #                     Cos.coef = TPbetas[3,4],
  #                     Sin.coef = TPbetas[3,5],
  #                     hour = 11.54,
  #                     SD.cov = 0,
  #                     WC.cov = seq(min(hmm.top.model$data$WC.Z),max(hmm.top.model$data$WC.Z), .1)
  # ) %>%
  #   mutate(TP2_1_reg = Intercept + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*hour/24) + Sin.coef*sin(2*pi*hour/24)) %>%
  #   mutate(TP2_1_exp = exp(TP2_1_reg)) %>%
  #   select(WC.cov, TP2_1_reg, TP2_1_exp)
  # TP2_3 <- data.frame(Intercept = TPbetas[4,1],
  #                     WC.coef = TPbetas[4,2],
  #                     SD.coef = TPbetas[4,3],
  #                     Cos.coef = TPbetas[4,4],
  #                     Sin.coef = TPbetas[4,5],
  #                     hour = 11.54,
  #                     SD.cov = 0,
  #                     WC.cov = seq(min(hmm.top.model$data$WC.Z),max(hmm.top.model$data$WC.Z), .1)
  # ) %>%
  #   mutate(TP2_3_reg = Intercept + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*hour/24) + Sin.coef*sin(2*pi*hour/24)) %>%
  #   mutate(TP2_3_exp = exp(TP2_3_reg)) %>%
  #   select(TP2_3_reg, TP2_3_exp)
  # TP2 <- cbind(TP2_1, TP2_3) %>%
  #   mutate(TP2_1 = TP2_1_exp/(1+TP2_1_exp+TP2_3_exp)) %>%
  #   mutate(TP2_3 = TP2_3_exp/(1+TP2_1_exp+TP2_3_exp)) %>%
  #   mutate(TP2_2 = 1-TP2_3-TP2_1)
  # 
  # ggplot(TP2) +
  #   geom_line(aes(x = WC.cov, y = TP2_2)) +
  #   geom_line(aes(x = WC.cov, y = TP2_1)) +
  #   geom_line(aes(x = WC.cov, y = TP2_3))
  # 
  # # 3->?
  # TP3_1 <- data.frame(Intercept = TPbetas[5,1],
  #                     WC.coef = TPbetas[5,2],
  #                     SD.coef = TPbetas[5,3],
  #                     Cos.coef = TPbetas[5,4],
  #                     Sin.coef = TPbetas[5,5],
  #                     hour = 11.54,
  #                     SD.cov = 0,
  #                     WC.cov = seq(min(hmm.top.model$data$WC.Z),max(hmm.top.model$data$WC.Z), .1)
  # ) %>%
  #   mutate(TP3_1_reg = Intercept + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*hour/24) + Sin.coef*sin(2*pi*hour/24)) %>%
  #   mutate(TP3_1_exp = exp(TP3_1_reg)) %>%
  #   select(WC.cov, TP3_1_reg, TP3_1_exp)
  # TP3_2 <- data.frame(Intercept = TPbetas[6,1],
  #                     WC.coef = TPbetas[6,2],
  #                     SD.coef = TPbetas[6,3],
  #                     Cos.coef = TPbetas[6,4],
  #                     Sin.coef = TPbetas[6,5],
  #                     hour = 11.54,
  #                     SD.cov = 0,
  #                     WC.cov = seq(min(hmm.top.model$data$WC.Z),max(hmm.top.model$data$WC.Z), .1)
  # ) %>%
  #   mutate(TP3_2_reg = Intercept + WC.coef*WC.cov + SD.coef*SD.cov + Cos.coef*cos(2*pi*hour/24) + Sin.coef*sin(2*pi*hour/24)) %>%
  #   mutate(TP3_2_exp = exp(TP3_2_reg)) %>%
  #   select(TP3_2_reg, TP3_2_exp)
  # TP3 <- cbind(TP3_1, TP3_2) %>%
  #   mutate(TP3_1 = TP3_1_exp/(1+TP3_1_exp+TP3_2_exp)) %>%
  #   mutate(TP3_2 = TP3_2_exp/(1+TP3_1_exp+TP3_2_exp)) %>%
  #   mutate(TP3_3 = 1-TP3_2-TP3_1)
  # 
  # ggplot(TP3) +
  #   geom_line(aes(x = WC.cov, y = TP3_3)) +
  #   geom_line(aes(x = WC.cov, y = TP3_1)) +
  #   geom_line(aes(x = WC.cov, y = TP3_2))
  # 
  # 
  # 
  # #####################################################
  # ###Trying to think of cool ways to visualize data ###
  # #####################################################
  # #Weather Data
  # TPdata <- read.csv("Results/HMMBehavioralStates_output.csv") %>%
  #   dplyr::select(ID, Timestamp, step, angle, WC.Z, SD.Z, hour, State) %>%
  #   mutate(State = as.factor(State))
  # 
  # Sdata <- TPdata %>% filter(State == 2)
  # Mdata <- TPdata %>% filter(State == 3)
  # 
  # ggplot(Sdata, aes(x = step, y = angle)) +
  #   geom_hex(bins = 50) +
  #   scale_fill_viridis_c(option = "magma") +
  #   theme_classic()
  # 
  # ggplot(Mdata, aes(x = step, y = angle)) +
  #   geom_hex(bins = 50) +
  #   scale_fill_viridis_c(option = "magma") +
  #   theme_classic() + 
  #   xlim(0, 1000)
  # 
  # 
  # ggplot(TPdata, aes(x = log(step), group = State)) +
  #   geom_density(aes(color = State, fill = State), alpha = .4) +
  #   theme_classic() 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # #Checkout Leaflet
  # ###########################################################################
  # ########################
  # ### CHECKING OUTPUTS ###
  # ########################
  # # #Where are the cases localized temporally?
  # # state_times <- turkey_states %>%
  # #   mutate(State = as.factor(State)) %>%
  # #   group_by(State, hour) %>%
  # #   summarize(Total = n())
  # # 
  # # require(tidyr)
  # # state_days <- turkey_states %>%
  # #   mutate(State = as.factor(State)) %>%
  # #   group_by(ID, YDay, State) %>%
  # #   summarize(Total = n()) %>%
  # #   group_by(YDay, State) %>%
  # #   summarize(Avg = mean(Total)) %>%
  # #   ungroup() %>%
  # #   pivot_wider(names_from = "State", values_from = "Avg") %>%
  # #   rename(Roosting = `1`, Loafing = `2`, Foraging = `3`) %>%
  # #   mutate(Ratio = Loafing/Foraging)%>%
  # #   mutate(L_propday = Loafing/(24-Roosting))%>%
  # #   mutate(F_propday = Foraging/(24-Roosting))
  # # 
  # # 
  # # ggplot(data = state_times, aes(x = hour, y = Total, group = as.factor(State))) +
  # #   geom_point(aes(shape = State, color = as.factor(State), size = 3))
  # # 
  # # ggplot(data = state_days, aes(x = YDay)) +
  # #   geom_point(aes(y = L_propday, shape = '21', size = 3)) +
  # #   geom_point(aes(y = F_propday, shape = '25', size = 3)) +
  # #   xlim(0,100)
  # 
  # 
  # # turkey_states1 <- turkey_states %>%
  # #   mutate(Date = as.Date(Timestamp)) %>%
  # #   #filter(hour(Timestamp) %in% ) #Maybe consider filtering out nightime hours but will need to use suncalc
  # #   group_by(ID, Date) %>%
  # #   summarize(State_Day_Avg = mean(State))
  # # 
  # # DailyStates <- ggplot(data=turkey_states1, aes(x = Date, y=State_Day_Avg)) +
  # #   geom_point(size = 1.5) +
  # #   xlab("Date") +
  # #   ylab("Daily Average State") +
  # #   theme_grey(base_size = 18) +
  # #   geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1) + #best fit line
  # #   theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  # #   theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
  # # DailyStates
  # # 
  # # id_list = unique(turkey_states$ID)
  # # id = id_list[3]
  # # states_onebird <- turkey_states1 %>% filter( ID == id)
  # # graph_onebird <- ggplot(data = states_onebird, aes(x = Date, y = State_Day_Avg)) +
  # #   geom_point(size = 1.5) +
  # #   labs(x = "Date", y = "State", title = id) +
  # #   theme_grey(base_size = 18) +
  # #   geom_smooth(method = "loess", span = .25, se=F, col = "red", lwd = 1)  #best fit line
  # # graph_onebird  
  # 
  # 
  # #Need to figure out a way to get the State values onto the points
  # #Need to figure out a way to incorproate individual variation into transition/mean etc.
  # #induce some sort of variation in the crawl wrap
}  #Old Code
