require(dplyr)
require(sf)
require(move)
require(lubridate)
require(R)
require(maptools)

setwd("E:/GitHub/WinterSelection")

#Login to download data from Movebank directly
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

#Get Movement Track
winterpoints_ind <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login, 
                animalName = "412", timestamp_start = "20190101000000000", timestamp_end = "20190316000000000")
m.dat <- winterpoints_ind@data %>%
  distinct() %>%
  arrange(timestamp)

for(i in 2:nrow(m.dat)){
  m.dat$sincelast[i] <- m.dat$timestamp[i] %--% m.dat$timestamp[i-1]
  m.dat$tillnext[i-1] <- m.dat$timestamp[i] %--% m.dat$timestamp[i-1]
}

m.dat$sincelast[1] <-3600
m.dat$tillnext[nrow(m.dat)] <-3600

for(i in 1:nrow(m.dat)){
  if(abs(m.dat$sincelast[i]) > 120 & abs(m.dat$tillnext[i]) > 120){
    m.dat$keep[i] <- 1
  }else{
    if(minute(m.dat$timestamp[i]) == 0 & second(m.dat$timestamp[i]) == 0){
      m.dat$keep[i] <- 0
      
    }else{
      m.dat$keep[i] <- 1
    }
  }
}

m.dat1 <- m.dat %>% filter(keep == 1) %>%
  dplyr::select(-keep, -tillnext, -sincelast) 
  
  
  
#Get Behavioral States
Bird1points <- read.csv("AllUsedPoints.csv") %>%
  filter(ID == "X412_2019") %>%
  dplyr::select(timestamp, State) %>%
  mutate(timestamp = as.POSIXct(timestamp))

#Merge Behavioral State by nearest date/time
require(data.table)
dt1 <- data.table(m.dat1)
dt2 <- data.table(Bird1points)
setkey(dt1, timestamp)
dt2[,mergerdate:=timestamp]
setkey(dt2, timestamp)
dt3 <- dt2[dt1, roll = 'nearest']
m.dat2 <- setDF(dt3) %>%
  dplyr::select(timestamps = timestamp, location_lat, location_long, State)

m.move <- move(x = m.dat2$location_long,
               y = m.dat2$location_lat,
               time = as.POSIXct(m.dat2$timestamps),
               data = m.dat2,
               proj = projection(winterpoints_ind))
m.move <- spTransform(m.move,  CRSobj="+proj=aeqd", center = T)
m.points <- as(m.move, "SpatialPointsDataFrame")
writeOGR(m.points, 
         dsn = './Overview Figure', 
         layer = "BehaviorPoints_412_2019", 
         driver = "ESRI Shapefile",
         overwrite_layer = T)

m.range <- brownian.bridge.dyn(m.move, 
                    raster = 30, 
                    location.error = 17, 
                    margin = 5, 
                    window.size = 15, 
                    ext = 2.5)
m.range <- raster2contour(m.range, level=c(.95))
writeOGR(m.range, 
         dsn = './Overview Figure', 
         layer = "FullHomeRange_412_2019", 
         driver = "ESRI Shapefile",
         overwrite_layer = T)

m.range3 <- brownian.bridge.dyn(m.move[m.move$State == 3], 
                               raster = 30, 
                               location.error = 17, 
                               margin = 5, 
                               window.size = 15, 
                               ext = 2.5)
m.range3 <- raster2contour(m.range3, level=c(.95))
writeOGR(m.range3, 
         dsn = './Overview Figure', 
         layer = "State3HomeRange_412_2019", 
         driver = "ESRI Shapefile",
         overwrite_layer = T)
