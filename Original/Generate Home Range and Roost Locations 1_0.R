#Load Packages
require(move)
require(ggmap)
require(mapproj)
require(sp)
require(maptools)
require(moveud)
require(dplyr)
require(lubridate)

#Set Working Directory
setwd("E:/Maine Drive/Analysis/Kaj Thesis")


#Login to download data from Movebank directly
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

#Load Bird List with associated Years
birdlist <- read.csv("BirdList.csv")

#Loop to download points, create dBBMM home range, and extract list of roost locations
#Only Looking at data between January 1 through March 15
rm(roostpoints_all)
rm(diurnal_all)
rm(winter.contours)
'%notin%' <- Negate('%in%')
for(i in 1:nrow(birdlist)){
  BirdID <- as.character(birdlist[i,1])
  Year <- birdlist[i,2]
  StartDate <- paste(Year, "0101000000000", sep = "")
  EndDate <- paste(Year, "0316000000000", sep = "")
  
  #Get all points for SSF/RSF
  winterpoints_ind <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login, 
                               animalName = BirdID, timestamp_start = StartDate, timestamp_end = EndDate)
  roostpoints_ind <- winterpoints_ind@data %>% 
    mutate(Time = hour(timestamp)) %>%
    filter(Time %in% 4:7) %>%
    mutate(BirdID = BirdID)
  diurnal_ind <- winterpoints_ind@data %>%  
    mutate(Time = hour(timestamp)) %>%
    filter(Time %notin% 4:7) %>%
    mutate(BirdID = BirdID)
  if(i==1){
    roostpoints_all <- roostpoints_ind
    diurnal_all <- diurnal_ind
  }else{
    roostpoints_all <- rbind(roostpoints_all, roostpoints_ind)
    diurnal_all <- rbind(diurnal_all, diurnal_ind)
  }
  
  #Create UDs for RSF availability
  t.winterpoints_ind <- spTransform(winterpoints_ind, center = T)
  winter.dbbmm <- brownian.bridge.dyn(t.winterpoints_ind, raster = 30, location.error = 17, margin = 5, window.size = 15, ext = 1.7)
  #Convert the dbbmm object to a raster
  #Change level to desired % for countor
  winter.UD <-raster2contour(winter.dbbmm, level=c(.95))
  winter.UD@data$ID <- paste(BirdID, Year, sep="_")
  winter.UD@data$Bird <- BirdID
  winter.UD@data$Season <- "Winter"
  winter.UD@data$Year <- Year
  
  #This if statement binds the individual layers into a new SpatialPolygonDataFrame
  if(exists("winter.contours")){
    winter.UD.t <- spTransform(winter.UD, CRS(proj4string(winter.contours)))
    winter.contours <- rbind(winter.contours, winter.UD.t)
  }else{
    winter.contours <- winter.UD
  }
    
}

write.csv(roostpoints_all, "RoostPoints.csv", row.names = F)
write.csv(diurnal_all, "DaytimeWinterPoints.csv", row.names = F) #Just day points
writeOGR(winter.contours, dsn = 'E:/Maine Drive/Analysis/Kaj Thesis', 
         layer = "AllWinterUDs_95_line", driver = "ESRI Shapefile")

#Convert SpatialLinesDF to SpatialPolygonsDF

##########################################################################
####You will need to Add both the Line and Polygon shapefiles and then####
####    "Join Field" in arcgis using FID to get the final product     ####
##########################################################################

winter.contours.polyset <- SpatialLines2PolySet(winter.contours)
winter.polygons <- PolySet2SpatialPolygons(winter.contours.polyset, close_polys=TRUE)
proj4string(winter.polygons) <- CRS(proj4string(winter.contours))
winter.polygons.df <- as(winter.polygons, "SpatialPolygonsDataFrame")

writeOGR(winter.polygons.df, dsn = 'E:/Maine Drive/Analysis/Kaj Thesis', 
         layer = "AllWinterUDs_95_polygon", driver = "ESRI Shapefile")
