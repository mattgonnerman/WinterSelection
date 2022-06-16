#Load Packages
require(move)
require(moveud)
require(dplyr)
require(momentuHMM)
require(ggplot2)
require(lubridate)
require(parallel)

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

require(lubridate)
raw.move.df <- raw.move.df %>% arrange(ID, Timestamp) %>%
  mutate(Timestamp = with_tz(Timestamp, tzone = "America/New_York"))

#add a covariate for Ordinal Day
raw.move.df$YDay <- yday(raw.move.df$Timestamp)

#limits optimization from straying outside parameter bounds
ln.prior <- function(theta) dnorm(theta[2],-4,2,log=TRUE)

#Use crawl function to fill in missing locations
Turkey.crawl.zm1 <- crawlWrap(raw.move.df, Time.name = "Timestamp", timeStep = "hour",
                             attempts=10, fixPar = c(NA, NA), prior = ln.prior, theta = c(1,1), retryFits=10, ncores = (detectCores()/2))

#plot(Turkey.crawl.zm)
# create momentuHMMData object from crwData object
turkeyData.zm1 <- prepData(data=Turkey.crawl.zm1)
#View(turkeyData.zm)
save(turkeyData.zm1, file = "./Data/turkeycrawl.RData")
load("./Data/turkeycrawl.RData")

### Loop through points to get locaiton specific weather covariates
source("./Weather Functions.R")
options(prism.path = "E:/GitHub/WinterSelection/Data/prismtmp")
lapply(c("sf", "raster", "prism", "rWind", "gdistance"), require, character.only = T)

dates <- unique(sort(as.Date(turkeyData.zm1$Timestamp, tz = "America/New_York")))
rm("pointsall")
for(i in 1:length(dates)){
  date <- dates[i]
  
  #Used Locations
  points <- turkeyData.zm1 %>% filter(as.Date(Timestamp, tz = "America/New_York") == date) %>%
    dplyr::select(ID, long = x, lat = y, Timestamp) %>%
    mutate(Date = as.Date(Timestamp, tz = "America/New_York")) %>%
    filter(!is.na(long)) %>%
    st_as_sf(., coords = c("long", "lat"), crs =32619)  %>%
    st_transform(4326)
  
  #Snow Raster
  # #Only need commented code if first time running 
  # download_SNODAS(date, out_dir = "C:/Users/mgonn/Downloads/Raw/")
  # unpack_SNODAS()
  # rasterize_SNODAS(data_dir = "C:/Users/mgonn/Downloads/Unpack/", 
  #                  out_dir = "E:/GitHub/WinterSelection/Data/Rasters/", 
  #                  format = "GTiff")
  snow.raster <- raster(list.files("E:/GitHub/WinterSelection/Data/Rasters/", full.names = T)[i])
  
  #Wind Raster
  wind.df <- wind.dl(yyyy = year(date), mm = month(date), dd = day(date), 10,
                     lon1 = floor(min(st_coordinates(points)[,1])-1), lon2 = ceiling(max(st_coordinates(points)[,1])+1),
                     lat1 = floor(min(st_coordinates(points)[,2])-1), lat2 = ceiling(max(st_coordinates(points)[,2])+1))
  wind.raster <- wind2raster(wind.df)
  
  #Temp Raster
  # #Only need commented code if first time running 
  # get_prism_dailys(type = "tmin", minDate = date, maxDate = date, keepZip = F)
  prismdata.list <- ls_prism_data()
  prismlist.pos <- which(grepl(gsub("-", x = date, ""), ls_prism_data()[,1]))
  temp.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos,2])
  
  #Extract Values
  points$SnowDepth <- raster::extract(snow.raster, points)
  points$Temp <- raster::extract(temp.raster, points %>% st_transform(crs(temp.raster)))
  points$WindDir <- raster::extract(wind.raster[[1]], points)
  points$WindSpd <- raster::extract(wind.raster[[2]], points)
  
  #Cleanup Files
  if(exists("pointsall")){
    pointsall <- rbind(pointsall, points)
  }else{
    pointsall <- points
  }
}

st_write(pointsall, dsn = "./Data", layer = "WeatherIndSpec", driver = "ESRI Shapefile", delete_layer = T)

weather.covs <- as.data.frame(pointsall %>%
                                st_drop_geometry()) %>%
  select(ID, Timestamp, SD = SnowDepth, TMIN = Temp, WindDir, WindSpd) %>%
  mutate(AWND = WindSpd*3.6) %>%
  mutate(WC = 13.12 + (.6215*TMIN)-(11.37*(AWND^0.16))+(.3965*TMIN*(AWND^0.16))) %>% #calculate windchill
  select(ID, Timestamp, SD, WC, WindDir)
write.csv(weather.covs, file = "./Data/WeatherVariables.csv")

#Add windchill and snow depth by date
turkeyData.weather <- merge(turkeyData.zm1, weather.covs, by = c("ID", "Timestamp"), all.x = T) %>%
  mutate(WC.Z =  as.numeric(scale(WC, center = T, scale = T))) %>%
  mutate(SD.Z =  as.numeric(scale(SD, center = T, scale = T))) %>%
  dplyr::select(Timestamp, ID, SD.Z, WC.Z) %>%
  mutate(SD.Z = ifelse(is.na(SD.Z), 0, SD.Z),
         WC.Z = ifelse(is.na(WC.Z), 0, WC.Z))
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
