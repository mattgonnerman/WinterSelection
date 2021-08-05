require(ggplot2)
require(dplyr)
require(tidyverse)
require(momentuHMM)
load("hmmtopmodel.RData")
hmm.top.model

#######################################
## Plot the Step Length Density/Hist ##
#######################################
m <- hmm.top.model
nbStates <- 3
# Estimated step length parameters
stepMean <- m$mle$step["mean",]
stepSD <- m$mle$step["sd",]

# Estimated turning angle parameters
angleMean <- m$mle$angle["mean",]
angleCon <- m$mle$angle["concentration",]

#Derive Shape and Rate Parameters of Gamma
stepShape <- stepMean^2/stepSD^2
stepRate <- stepMean/stepSD^2

#Will weight densities according to points
# states <- viterbi(m)
hmm_data.raw <- read.csv("Results/HMMBehavioralStates_output.csv") 
states <- hmm_data.raw$State
states1 <- states[which(m$data$step > 0)]
w <- c()
for(i in 1:3){
  w[i] <- length(which(states1==i))/length(states1)
}

gamma.data <- data.frame(
  Roosting = dgamma(seq(0, 1000, .1), shape = stepShape[1], scale = 1/stepRate[1]),
  Stationary = dgamma(seq(0, 1000, .1), shape = stepShape[2], scale = 1/stepRate[2]),
  Mobile = dgamma(seq(0, 1000, .1), shape = stepShape[3], scale = 1/stepRate[3]),
  X = seq(0, 1000, .1)
)

for(i in 1:nbStates){
  gamma.data[,i] <- gamma.data[,i]*w[i]
}

step.graph <- ggplot(gamma.data, aes(x = X)) +
  geom_line(aes(y = Roosting, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = Stationary, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = Mobile, color = "#46e300", linetype = "twodash"), size = 3) +
  # geom_vline(aes(xintercept = stepMean[1]), color = "#1cade4", size = 1, linetype = "dashed") +
  # geom_vline(aes(xintercept = stepMean[2]), color = "#f1a806", size = 1, linetype = "dashed") +
  # geom_vline(aes(xintercept = stepMean[3]), color = "#46e300", size = 1, linetype = "dashed") +
  scale_y_continuous(limits=c(0,max(gamma.data[,1:nbStates]))) + 
  scale_x_continuous(limits=c(0,200)) +
  theme_classic(base_size = 45) +
  xlab("Step Length") +
  ylab("Density") +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend") + 
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") + 
  theme(legend.position = "none") 
# step.graph
# ggsave("Results/StepLengthDensity.png", width = 9, height = 7, units = "in")

require(circular)
wrpcauchy.data <- data.frame(
  Roosting = dwrappedcauchy(seq(-pi, pi, .01), mu = circular(angleMean[1]), rho = angleCon[1]),
  Stationary = dwrappedcauchy(seq(-pi, pi, .01), mu = circular(angleMean[2]), rho = angleCon[2]),
  Mobile = dwrappedcauchy(seq(-pi, pi, .01), mu = circular(angleMean[3]), rho = angleCon[3]),
  X = seq(-pi, pi, .01)
)

for(i in 1:nbStates){
  wrpcauchy.data[,i] <- wrpcauchy.data[,i]*w[i]
}

anglecon.graph <- ggplot(wrpcauchy.data, aes(x = X)) +
  geom_line(aes(y = Roosting, color = "#1cade4", linetype = "longdash"), size = 3) +
  geom_line(aes(y = Stationary, color = "#f1a806", linetype = "solid"), size = 3) +
  geom_line(aes(y = Mobile, color = "#46e300", linetype = "twodash"), size = 3) +
  scale_y_continuous(limits=c(0,max(wrpcauchy.data[,1:nbStates]))) + 
  scale_x_continuous(limits=c(-pi,pi), 
                     breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = c(expression(paste("-",pi, sep = "")),
                                expression(paste("-",pi,"/2", sep = "")),
                                "0",
                                expression(paste(pi,"/2", sep = "")),
                                expression(pi))) +
  theme_classic(base_size = 45) +
  xlab("Turning Angle") +
  ylab("Density") +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806", "#46e300"),
                       labels = c("Roosting", "Stationary", "Mobile"),
                       guide = "legend") +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid", "twodash"),
                          labels = c("Roosting", "Stationary", "Mobile"),
                          guide = "legend") + 
  theme(legend.position = "none") 
# anglecon.graph
# ggsave("Results/TurningAngleDensity.png", width = 7, height = 7, units = "in")

###################################
### Example Movement Tracks Map ###
###################################
require(dplyr)
require(sf)

setwd("E:/GitHub/WinterSelection")

behpoints <- st_read("./Overview Figure/BehaviorPoints_412_2019.shp") %>%
  mutate(State = factor(State, levels = c(3,2,1)))

n <- nrow(behpoints) - 1
linestrings <- lapply(X = 1:n, FUN = function(x) {
  
  pair <- st_combine(behpoints[x:(x+1),])
  line <- st_as_sf(st_cast(pair, "LINESTRING"))
  return(line)
  
})

# One MULTILINESTRING object with all the LINESTRINGS
multilinestring <- do.call("rbind", linestrings)
multilinestring$State <- behpoints$State[2:(n+1)]

## get the geom coordinates as data.frame
geomdf <- st_coordinates(multilinestring)

## reverse Y coords
geomdf[,"Y"] <- geomdf[,"Y"]*-1

simplemap <- ggplot() +
  geom_sf(data = multilinestring, aes(color = State), size = 1.5) +
  geom_sf(data = behpoints, aes(color = State), size = 5) + 
  scale_color_manual(values = c("#46e300", "#f1a806", "#1cade4")) +
  theme_void() +
  theme(legend.position = "none") 
  coord_flip()

#Combine into grid
require(cowplot)
require(patchwork)
legend <- get_legend(step.graph + theme(legend.title = element_blank() ,
                                        legend.position = "bottom", 
                                        legend.key.width=unit(3,"inch")))
Density_grid <- (step.graph / anglecon.graph) | simplemap
Density_grid <- Density_grid + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 45))
png('Results/Density Grid 2.png', width = 1500, height = 1450)
plot_grid(plotlist = list(Density_grid, legend),
          nrow = 2,
          rel_heights = c(20,1)
)
dev.off()
