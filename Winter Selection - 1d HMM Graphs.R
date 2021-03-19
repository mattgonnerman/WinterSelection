require(ggplot2)
require(dplyr)
require(tidyverse)
require(hexbin)

#Load MLE of Transition Porbability betas
TPbetas <- read.csv('Results/HMM - MLE of betas.csv')

#Weather Data
TPdata <- read.csv("Results/HMMBehavioralStates_output.csv") %>%
  dplyr::select(ID, Timestamp, step, angle, WC.Z, SD.Z, hour, State) %>%
  mutate(State = as.factor(State))

Sdata <- TPdata %>% filter(State == 2)
Mdata <- TPdata %>% filter(State == 3)
  
ggplot(Sdata, aes(x = step, y = angle)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(option = "magma") +
  theme_classic()

ggplot(Mdata, aes(x = step, y = angle)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(option = "magma") +
  theme_classic() + 
  xlim(0, 1000)


ggplot(TPdata, aes(x = log(step), group = State)) +
  geom_density(aes(color = State, fill = State), alpha = .4) +
  theme_classic() 

