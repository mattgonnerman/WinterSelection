require(momentuHMM)
require(dplyr)
require(ggplot2)


load("hmmtopmodel.RData")
hmm.top.model

turkey_states <- read.csv("Results/HMMBehavioralStates_output.csv")

AllPoints.final <- read.csv("AllPoints.csv")

hist(AllPoints.final$State)

#Proportion In each state for all locations
nrow(AllPoints.final[which(AllPoints.final$State == 1),])/nrow(AllPoints.final)
nrow(AllPoints.final[which(AllPoints.final$State == 2),])/nrow(AllPoints.final)
nrow(AllPoints.final[which(AllPoints.final$State == 3),])/nrow(AllPoints.final)

#Time Spent in each State
require(lubridate)
AllPoints.final <-   read.csv("AllPoints.csv")
time_spent <- AllPoints.final %>%
  filter(case_ == T) %>%
  select(ID, t2_, WC.Z , SD.Z, State) %>%
  mutate(Timestamp = as.POSIXct(t2_, format = "%Y-%m-%d %H:%M:%S", tz = "EST")) %>%
  mutate(hour = hour(Timestamp)) %>%
  mutate(OrdinalDay = yday(Timestamp)) %>%
  mutate(Date = as.Date(Timestamp)) %>%
  mutate(Year = year(Timestamp)) %>%
  mutate(State1 = ifelse(State == 1, 1, 0))  %>%
  mutate(State2 = ifelse(State == 2, 1, 0)) %>%
  mutate(State3 = ifelse(State == 3, 1, 0))

weather.raw <- read.csv("WeatherVariables.csv") %>%
  mutate(DATE = as.Date(DATE, format = "%m/%d/%Y")) %>%
  mutate(AWND = ifelse(is.na(AWND), 0, AWND)) %>% #If AWND is NA, that means no wind
  dplyr::select(Date = DATE, TMIN, AWND, WDF2, SD = SNWD) %>%
  mutate(WC = 13.12 + (.6215*TMIN)-(11.37*(AWND^0.16))+(.3965*TMIN*(AWND^0.16))) %>% #calculate windchill
  mutate(WC_prev = lag(WC, 1)) #want column with previous days wind chill

time_spent.weather <- merge(time_spent, weather.raw, by = "Date", all.x = T)

weathermeans <- time_spent.weather %>%
  group_by(State) %>%
  summarize(Mean_WC = mean(WC),
            Mean_SD = mean(SD))

mean(time_spent.weather$WC)
mean(time_spent.weather$SD)

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










