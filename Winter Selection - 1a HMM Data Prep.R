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
Turkey.crawl.zm1 <- crawlWrap(raw.move.df, Time.name = "Timestamp", timeStep = "hour",
                             attempts=20, fixPar = c(NA, NA), prior = ln.prior, retryFits=10, ncores = 5)
#plot(Turkey.crawl.zm)
# create momentuHMMData object from crwData object
turkeyData.zm1 <- prepData(data=Turkey.crawl.zm1)
#View(turkeyData.zm)

#Add covariate for weather (Wind Chill and Snow Depth)
#Daily weather summaries downloaded for Bangor Airport weather station
weather.raw <- read.csv("WeatherVariables.csv") %>%
  mutate(Date = as.Date(DATE, format = "%m/%d/%Y"))%>%
  mutate(AWND = ifelse(is.na(AWND), 0, AWND)) %>% #If AWND is NA, that means no wind
  dplyr::select(Date, TMIN, AWND, WDF2, SD = SNWD)%>%
  mutate(WC = 13.12 + (.6215*TMIN)-(11.37*(AWND^0.16))+(.3965*TMIN*(AWND^0.16))) #calculate windchill
turkeyData.zm1$Date <- as.Date(turkeyData.zm1$Timestamp)
#Add windchill and snow depth by date
turkeyData.weather <- merge(turkeyData.zm1, weather.raw, by = "Date", all.x = T) %>%
  mutate(WC.Z =  as.numeric(scale(WC, center = T, scale = T))) %>%
  mutate(SD.Z =  as.numeric(scale(SD, center = T, scale = T))) %>%
  dplyr::select(Timestamp, ID, SD.Z, WC.Z) 
Turkey.crawl.zm <- crawlMerge(Turkey.crawl.zm1, turkeyData.weather, Time.name = "Timestamp")

# create momentuHMMData object from crwData object
turkeyData.zm <- prepData(data=Turkey.crawl.zm, covNames = "SD.Z","WC.Z")
#View(turkeyData.zm)

# add cosinor covariate based on hour of day
turkeyData.zm$hour <- as.numeric(strftime(turkeyData.zm$Timestamp, format = "%H", tz="GMT"))

#determine proportion of zero step lengths
whichzero <- which(turkeyData.zm$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(turkeyData.zm)

#Add covariate for time since midday
turkeyData.zm$T2MidDay <- abs(12-turkeyData.zm$hour)
