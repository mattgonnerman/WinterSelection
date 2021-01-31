# https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html#5_ssf_with_multiple_animals


pcks <- list("sp", "sf", "dplyr", "raster", "rgdal", "lubridate", "amt")
sapply(pcks, require, char = TRUE)

#Set Working Directory
setwd("E:/Maine Drive/Analysis/Kaj Thesis")


#Load used points
roostlocations.obs.raw <- read.csv("RoostPoints.csv")
daylocations.obs.raw <- read.csv("DaytimeWinterPoints.csv")
winter95UDs.raw <- readOGR(".", "AllWinterUDs_95_polygon")
winter95UDs <- spTransform(winter95UDs.raw, CRS("+proj=longlat +datum=WGS84"))

#Clean up dataframes
roostlocations.obs <- roostlocations.obs.raw %>%
  dplyr::select(BirdID, timestamp, Lat = location_lat, Long = location_long) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(Year = year(timestamp)) %>% mutate(Used = 1) %>% mutate(Roost = 1) %>%
  dplyr::select(-timestamp)

daylocations.obs <- daylocations.obs.raw %>%
  dplyr::select(BirdID, timestamp, Lat = location_lat, Long = location_long) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(Year = year(timestamp)) %>% mutate(Used = 1) %>% mutate(Roost = 0) %>%
  dplyr::select(-timestamp)

### Roost and Daytime Available Locations for RSF ###
rm(roostlocations.random)
rm(daylocations.random)

for(i in 1:length(winter95UDs)){
  Bird <- winter95UDs@data$Bird[i]
  Yearobs <- winter95UDs@data$Year[i]
  
  #How many used roosts/day locations are there?
  observed.roost <- nrow(roostlocations.obs %>% filter(BirdID == Bird) %>% filter(Year == Yearobs)) 
  observed.day <- nrow(daylocations.obs %>% filter(BirdID == Bird) %>% filter(Year == Yearobs))
  
  roost.avail <- as.data.frame(spsample(winter95UDs[i,], observed.roost*10, "random")@coords) %>%
    mutate(BirdID = Bird) %>% mutate(Year = Yearobs) %>% mutate(Used = 0) %>% mutate(Roost = 1) %>%
    rename(Long = x, Lat = y)
    
  day.avail <- as.data.frame(spsample(winter95UDs[i,], observed.day*10, "random")@coords) %>%
    mutate(BirdID = Bird) %>% mutate(Year = Yearobs) %>% mutate(Used = 0) %>% mutate(Roost = 0) %>%
    rename(Long = x, Lat = y)
  
  if(exists("roostlocations.random")){
    roostlocations.random <- rbind(roostlocations.random,roost.avail)
    daylocations.random <- rbind(daylocations.random,day.avail)
  }else{
    roostlocations.random <- roost.avail
    daylocations.random <- day.avail
  }
}

rsf.allpoints <- do.call("rbind", list(roostlocations.obs, daylocations.obs, roostlocations.random, daylocations.random)) %>%
  arrange(Year, BirdID, Used, Roost)
write.csv(rsf.allpoints, "RSF_Used_Avail_points.csv")

### Roost Location SSF Randoms ###
roostlocations.obs <- roostlocations.obs.raw %>%
  dplyr::select(BirdID, timestamp, Lat = location_lat, Long = location_long) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(Year = year(timestamp)) %>%
  mutate(ID = paste(BirdID, Year, sep = "_"))

nested.roosts <- roostlocations.obs %>%
  nest(-ID) %>%
  mutate(track = map(data, ~ mk_track(., Long, Lat, timestamp, crs= CRS("+proj=longlat +datum=WGS84"))))

roost.true.steps <- map(nested.roosts$track, ~ steps)
roost.random.steps <- map(roost.true.steps, random_steps, n=10)
roost.steps <- bind_rows(roost.random.steps, .id="ID")
coordinates(roost.steps) <- ~x2_+y2_
proj4string(roost.steps) <- CRS("+proj=longlat +datum=WGS84")

### Daytime Location SSF Randoms ###

