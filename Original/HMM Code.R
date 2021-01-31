### Set Working Directory
setwd("E:/Maine Drive/Analysis/Kaj Thesis")

#Load Packages
require(move)
require(moveud)
require(dplyr)
require(momentuHMM)
require(ggplot2)
require(lubridate)

################################################
### Load Movement Data Used in INLA Analysis ###
################################################
#Login to download data from Movebank directly
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

#Load list with BirdID and associated time frame
birdlist <- read.csv("BirdList.csv")

### Loop to download points, create dBBMM home range, and create df of all locations
### Only Looking at data between January 1 through March 15
'%notin%' <- Negate('%in%') #Need the opposite of %in%

rm(full_all)#Clear objects so you don't duplicate data (for if not starting new session)
for(i in 1:nrow(birdlist)){
  BirdID <- as.character(birdlist[i,1])
  Year <- birdlist[i,2]
  StartDate <- paste(Year, "0101000000000", sep = "")
  EndDate <- paste(Year, "0316000000000", sep = "")
  
  #Get all pointsa from Jan1 through March15
  winterpoints_ind <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login, 
                                      animalName = BirdID, timestamp_start = StartDate, timestamp_end = EndDate)
  
  winterpoints_ind@data$Time <- hour(winterpoints_ind@data$timestamp)
  winterpoints_ind@data$BirdID <- paste(BirdID, Year, sep = "_")
  winterpoints_ind@data$Year <- Year
  winterpoints_ind@data$TrackID <- paste(BirdID, Year, sep = "_")
  rownames(winterpoints_ind@idData) <- paste(BirdID, Year, sep = "_")

  if(exists("full_all")){ #rbind ind bird data to create one large df
    full_all <- moveStack(c(full_all, winterpoints_ind))
  }else{
    full_all <- winterpoints_ind
  }
}

raw.move.df <- as.data.frame(coordinates(full_all)) %>%
  mutate(Timestamp = round_date(timestamps(full_all), unit = "hour")) %>%
  mutate(ID = as.character(full_all@trackId)) %>%
  rename(location_long = coords.x1, location_lat = coords.x2)

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(raw.move.df[,1:2], proj4string=CRS(projection(full_all)))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=19 ellps=WGS84"))

# add UTM locations to data frame
raw.move.df$x <- attr(utmcoord,"coords")[,1]
raw.move.df$y <- attr(utmcoord,"coords")[,2]

# Duplicate roost locations for all night time hours
# includes random error based on gps location uncertainty (~12.7m)
roosts <- raw.move.df %>%
  rename(lon = location_long, lat = location_lat) %>%
  mutate(date = as.Date(Timestamp)) %>%
  filter(hour(Timestamp) < 6 & hour(Timestamp) > 3)

require(suncalc)
suntimes <- getSunlightTimes(data = roosts, tz = "GMT", keep = c("sunrise", "sunset"))
roosts$sunrise <- suntimes$sunrise
roosts$sunset <- suntimes$sunset

for(j in 1:nrow(roosts)){
  roosttime <- roosts$Timestamp[j]
  hr_before_sun <- floor(as.numeric(roosts$sunrise[j]-roosts$Timestamp[j]))-1
  hr_after_sun <- floor(24 - as.numeric(roosts$sunset[j]-roosts$Timestamp[j]))-1
  
  working.df <- roosts[rep(j, hr_before_sun + hr_after_sun),] %>%
    dplyr::select(location_long = lon, location_lat = lat, Timestamp, ID, x, y)
  
  for(i in 1:hr_after_sun){
    working.df$Timestamp[i] <- working.df$Timestamp[i] - (i*60*60) #change time
    #working.df$x[i] <- working.df$x[i] + rnorm(1,0,6.5) #change x UTM
    #working.df$y[i] <- working.df$y[i] + rnorm(1,0,6.5) #change y UTM
  }
  for(i in 1:hr_before_sun){
    working.df$Timestamp[i+hr_after_sun] <- working.df$Timestamp[i+hr_after_sun] + (i*60*60) #change time
    #working.df$x[i+hr_after_sun] <- working.df$x[i+hr_after_sun] + rnorm(1,0,6.5) #change x UTM
    #working.df$y[i+hr_after_sun] <- working.df$y[i+hr_after_sun] + rnorm(1,0,6.5) #change y UTM
  }
  raw.move.df <- rbind(raw.move.df, working.df) #add to original database
}

raw.move.df <- raw.move.df %>% arrange(ID, Timestamp)

#add a covariate for Ordinal Day
raw.move.df$YDay <- yday(raw.move.df$Timestamp)

#limits optimization from straying outside parameter bounds
ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)

#Use crawl function to fill in missing locations
Turkey.crawl.zm <- crawlWrap(raw.move.df, Time.name = "Timestamp", timeStep = "hour",
                             attempts=20, fixPar = c(NA, NA), prior = ln.prior, retryFits=10)
#plot(Turkey.crawl.zm)

# create momentuHMMData object from crwData object
turkeyData.zm <- prepData(data=Turkey.crawl.zm, covNames = "YDay")
#View(turkeyData.zm)

# add cosinor covariate based on hour of day
turkeyData.zm$hour <- as.numeric(strftime(turkeyData.zm$Timestamp, format = "%H", tz="GMT"))

#determine proportion of zero step lengths
whichzero <- which(turkeyData.zm$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(turkeyData.zm)


###########################################################
####### Stationary, Searching, and Dispersing Models ######
###########################################################
### MODEL 1
### momentuHMM Function Objects
nSims <- 100 # number of imputatons
retryFits <- 10 # number attempt to re-fit based on random perturbation
nbStates <- 3 # Number of states
stateNames <- c("roost", "loafing", "foraging") # label states
dist = list(step = "gamma", angle = "wrpcauchy") # distributions for observation processes

## constrain step length parameters: 
# Mean -> dispersal>searching>roosting
# SD -> no relationship
stepDM<-matrix(c(
  1,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,1),
  nrow = 3*nbStates,byrow=TRUE,
  dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates), paste0("zero_",1:nbStates)),
                c("mean_123:(Intercept)", "mean_2","mean_3",
                  paste0("sd_",1:nbStates,":(Intercept)"),
                  paste0("zero_",1:nbStates,":(Intercept)"))))

#define the directions of the differences
stepworkBounds <- matrix(c(-Inf, 0,0,rep(-Inf,6),rep(Inf,9)),nrow = 3*nbStates,
                         dimnames=list(colnames(stepDM),c("lower","upper")))
#Userbound constraint on step
stepBounds <- matrix(c(0,10,
                       0,Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       0, Inf,
                       .98,1,
                       0,.01,
                       0,.01), nrow = 3*nbStates, byrow = T,
                     dimnames = list(rownames(stepDM), c("lower", "upper")))

## constrain turning angle concentration parameters: 
# Concentration -> searching < dispersal
angleDM<-matrix(c(1,0,0,
                  1,1,0,
                  1,1,1),nrow = nbStates,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_1:(Intercept)","concentration_2","concentration_3")))

#Restrict angle concentration such that dispersal > .75 while dispersal>localized
#Userbound contraint on angle
angleBounds <- matrix(c(0,0.94, 
                        0,0.94, 
                        0,0.94),nrow = nbStates,
                      byrow=TRUE, dimnames=list(rownames(angleDM),
                                                c("lower","upper"))) 
#define direction of differences
angleworkBounds <- matrix(c(-Inf,0,0,rep(Inf,3)),
                          nrow = ncol(angleDM), 
                          dimnames=list(colnames(angleDM),c("lower","upper")))


#Bundle individual parameter DM and workbounds
DM<-list(step=stepDM,angle=angleDM)
workBounds<-list(step=stepworkBounds,
                 angle=angleworkBounds)
userBounds <- list(step = stepBounds, 
                   angle = angleBounds)

#prevents working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

#Fixes "roost" for all, fixPar = fixPar in fitHMM
#fixPar<-list(step=c(rep(NA,nbStates*2),NA, rep(stats::qlogis(1.e-100), 3)))

# initial parameters
Par <- list(step=c(6, 50,120, #mean
                   6, 15,15, #sd
                   .99, .005, .005), #zero mass
            angle = c(.01, 0.2 ,0.3))

Par0_m1.zm <- getParDM(data = turkeyData.zm,
                       nbStates = nbStates, 
                       dist = dist,
                       Par = Par, 
                       DM = DM, 
                       workBounds = workBounds, 
                       userBounds = userBounds,
                       estAngleMean = list(angle = FALSE))

# fit model 1
turk_m1.zm <- fitHMM(data = turkeyData.zm, 
                     nSims = nSims,
                     nbStates = nbStates, 
                     dist = dist, 
                     Par0 = Par0_m1.zm, 
                     DM = DM, 
                     workBounds = workBounds, 
                     userBounds = userBounds,
                     estAngleMean = list(angle=FALSE), 
                     prior = prior,
                     stateNames = stateNames,
                     formula = ~ID,
                     )
turk_m1_states <- viterbi(turk_m1.zm)
turkey_states <- turkeyData.zm
turkey_states$State <- turk_m1_states

#Write output to csv to bring in, simplify, and combine with SSF data
write.csv(turkey_states, "HMMBehavioralStates_output.csv")





###########################################################################
########################
### CHECKING OUTPUTS ###
########################
# Look at output maps/graphs
plot(turk_m1.zm,legend.pos="right")
hist(turk_m1_states)

#Where are the cases localized temporally?
state_times <- turkey_states %>%
  mutate(State = as.factor(State)) %>%
  group_by(State, hour) %>%
  summarize(Total = n())

require(tidyr)
state_days <- turkey_states %>%
  mutate(State = as.factor(State)) %>%
  group_by(ID, YDay, State) %>%
  summarize(Total = n()) %>%
    group_by(YDay, State) %>%
  summarize(Avg = mean(Total)) %>%
  ungroup() %>%
  pivot_wider(names_from = "State", values_from = "Avg") %>%
  rename(Roosting = `1`, Loafing = `2`, Foraging = `3`) %>%
  mutate(Ratio = Loafing/Foraging)%>%
  mutate(L_propday = Loafing/(24-Roosting))%>%
  mutate(F_propday = Foraging/(24-Roosting))


ggplot(data = state_times, aes(x = hour, y = Total, group = as.factor(State))) +
  geom_point(aes(shape = State, color = as.factor(State), size = 3))

ggplot(data = state_days, aes(x = YDay)) +
  geom_point(aes(y = L_propday, shape = '21', size = 3)) +
  geom_point(aes(y = F_propday, shape = '25', size = 3)) +
  xlim(0,100)


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