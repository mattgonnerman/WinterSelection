require(ggplot2)
require(dplyr)
require(tidyverse)
require(hexbin)

######################
## Prepare the data ##
######################
# Raw data
data <- turkeyData.zm
step <- data$step
angle <- data$angle

# MomentuHMM model
m <- hmm.top.model

# Animal Index
nbAnimals <- length(unique(m$data$ID))
animalsInd <- 1:nbAnimals
ID <- unique(m$data$ID)[animalsInd]

# States decoding with Viterbi 
states <- viterbi(m)

#######################################
## Plot the Step Length Density/Hist ##
#######################################
turkey_states <- turkey_states %>%
  filter(!is.na(location_lat)) %>%
  mutate(State = ifelse(State == 1, "Roosting", 
                 ifelse(State == 2, "Stationary",
                 ifelse(State == 3, "Mobile", NA)))) %>%
  mutate(State = factor(State, levels = c("Roosting", "Stationary", "Mobile"))) 


ggplot(turkey_states, aes(x = step, group = State)) +
  geom_density(aes(color = State, fill = State), alpha = .7, size = 1.2) +
  theme_classic(base_size = 25) +
  scale_y_continuous(limits=c(0,.04), oob = scales::squish) + 
  scale_x_continuous(limits=c(0,500)) +
  scale_color_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#1cade4", "#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                      # labels = c("Roost", "Stationary", "Mobile"),
                      values = c("#1cade4", "#f1a806", "#46e300")) + 
  theme(legend.position = c(0.77, 0.77)) +
  xlab("Step Length") +
  ylab("Density")
ggsave("Results/StepLengthDensity.jpeg", width = 7, height = 7, units = "in")

ggplot(turkey_states, aes(x = angle, group = State)) +
  geom_density(aes(color = State), alpha = .7, size = 1.2) +
  theme_classic(base_size = 25) +
  # scale_y_continuous(limits=c(0,.04), oob = scales::squish) + 
  scale_x_continuous(limits=c(-pi,pi)) +
  scale_color_manual(name = "Behavioral\nState",
                     # labels = c("Roost", "Stationary", "Mobile"),
                     values = c("#f1a806", "#46e300")) +
  scale_fill_manual(name = "Behavioral\nState",
                    # labels = c("Roost", "Stationary", "Mobile"),
                    values = c("#f1a806", "#46e300")) + 
  theme(legend.position = c(0.5, 0.3)) +
  xlab("Turning Angle") +
  ylab("Density")
ggsave("Results/TurningAngleDensity.jpeg", width = 7, height = 7, units = "in")


##################################################
## Plot the t.p. as functions of the covariates ##
##################################################
#Load MLE of Transition Porbability betas
TPbetas <- read.csv('Results/HMM - MLE of betas.csv')
row.names(TPbetas) <- TPbetas$X
TPbetas <- TPbetas %>% select(-X)
TPbetas <- as.data.frame(t(TPbetas)) %>%
  dplyr::select('(Intercept)', WC.Z, SD.Z, 'cosinorCos(hour, period = 24)', 'cosinorSin(hour, period = 24)')

X.wc <- seq(-2,.5,.1)
X.sd <- seq(-.5,6,.1)
X.coshour <- cos(2*pi*seq(0,24, .25)/24)
X.sinhour <- sin(2*pi*seq(0,24, .25)/24)

x <-6

eq <- TPbetas[x,1] + TPbetas[x,2]*-1 + TPbetas[x,3]*X.sd + TPbetas[x,4]*cos(2*pi*11.54/24) + TPbetas[x,5]*sin(2*pi*11.54/24)
plot((exp(eq)/(1+exp(eq)))~X.sd)







#####################################################
###Trying to think of cool ways to visualize data ###
#####################################################
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