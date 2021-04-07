require(ggplot2)
require(dplyr)
require(tidyverse)


#######################################
## Plot the Step Length Density/Hist ##
#######################################
m <- hmm.top.model
# Estimated step length parameters
stepMean <- m$mle$step["mean",]
stepSD <- m$mle$step["sd",]

# Estimated turning angle parameters
angleMean <- m$mle$angle["mean",]
angleCon <- m$mle$angle["concentration",]

#Derive Shape and Rate Parameters of Gamma
stepShape <- stepMean^2/stepSD^2
stepRate <- stepMean/stepSD^2

gamma.data <- data.frame(
  Roosting = dgamma(seq(1, 1000, .1), shape = stepShape[1], scale = 1/stepRate[1]),
  Stationary = dgamma(seq(1, 1000, .1), shape = stepShape[2], scale = 1/stepRate[2]),
  Mobile = dgamma(seq(1, 1000, .1), shape = stepShape[3], scale = 1/stepRate[3]),
  X = seq(1, 1000, .1)
)

ggplot(gamma.data, aes(x = X)) +
  geom_line(aes(y = Roosting, color = "#1cade4"), size = 2) +
  geom_line(aes(y = Stationary, color = "#f1a806"), size = 2) +
  geom_line(aes(y = Mobile, color = "#46e300"), size = 2) +
  geom_vline(aes(xintercept = stepMean[1]), color = "#1cade4", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = stepMean[2]), color = "#f1a806", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = stepMean[3]), color = "#46e300", size = 1, linetype = "dashed") +
  scale_y_continuous(limits=c(0,.05)) + 
  scale_x_continuous(limits=c(0,200)) +
  theme_classic(base_size = 25) +
  xlab("Step Length") +
  ylab("Density") +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend") + 
  theme(legend.position = c(0.85, 0.8)) 
ggsave("Results/StepLengthDensity.png", width = 9, height = 7, units = "in")

require(circular)
wrpcauchy.data <- data.frame(
  Roosting = dwrappedcauchy(seq(-pi, pi, .01), mu = angleMean[1], rho = angleCon[1]),
  Stationary = dwrappedcauchy(seq(-pi, pi, .01), mu = angleMean[2], rho = angleCon[2]),
  Mobile = dwrappedcauchy(seq(-pi, pi, .01), mu = angleMean[3], rho = angleCon[3]),
  X = seq(-pi, pi, .01)
)

ggplot(wrpcauchy.data, aes(x = X)) +
  geom_line(aes(y = Roosting, color = "#1cade4"), size = 2) +
  geom_line(aes(y = Stationary, color = "#f1a806"), size = 2) +
  geom_line(aes(y = Mobile, color = "#46e300"), size = 2) +
  scale_y_continuous(limits=c(0,.3)) + 
  scale_x_continuous(limits=c(-pi,pi)) +
  theme_classic(base_size = 25) +
  xlab("Step Length") +
  ylab("Density") +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend") + 
  theme(legend.position = c(0.85, 0.8)) 

ggsave("Results/TurningAngleDensity.png", width = 7, height = 7, units = "in")


