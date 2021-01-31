#### This is to consolidate all the code for the SSF/RSF Prep
#### Need Used and available points for a roost specific and all points analysis
#### For roost specific, get two sets of available points for 2 scales of analysis
#### 1)Available based on roost to roost, 2) Available based on point previous to roost
#### UPDATE: Could not get option 1 to work, and since Erik wanted 2 anyways, didnt push through that issue yet

### Load Packages
require(move)
require(ggmap)
require(mapproj)
require(sp)
require(maptools)
require(moveud)
require(dplyr)
require(lubridate)
require(spdplyr)
require(raster)

### Set Working Directory
setwd("E:/Maine Drive/Analysis/Kaj Thesis")


########################################
### Generate Used Points for SSF/RSF ###
########################################
#Login to download data from Movebank directly
login <- movebankLogin(username = "matthew.gonnerman", password="26qPDLY9YN")

#Load list with BirdID and associated time frame
birdlist <- read.csv("BirdList.csv")

### Loop to download points, create dBBMM home range, and extract df of roost, daytime, and all locations
### Only Looking at data between January 1 through March 15
#Clear objects so you don't duplicate data
rm(roostpoints_all)
rm(diurnal_all)
rm(full_all)
rm(winter.contours)

'%notin%' <- Negate('%in%') #Need the opposite of %in%

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
  
  full_ind <- winterpoints_ind@data %>%  
    mutate(Time = hour(timestamp)) %>%
    mutate(BirdID = BirdID)
  
  if(i==1){
    roostpoints_all <- roostpoints_ind
    diurnal_all <- diurnal_ind
    full_all <- full_ind
  }else{
    roostpoints_all <- rbind(roostpoints_all, roostpoints_ind)
    diurnal_all <- rbind(diurnal_all, diurnal_ind)
    full_all <- rbind(full_all, full_ind)
  }
  
  # #Create UDs for RSF availability
  # t.winterpoints_ind <- spTransform(winterpoints_ind, center = T)
  # winter.dbbmm <- brownian.bridge.dyn(t.winterpoints_ind, raster = 30, location.error = 17, margin = 5, window.size = 15, ext = 1.7)
  # #Convert the dbbmm object to a raster
  # #Change level to desired % for countor
  # winter.UD <-raster2contour(winter.dbbmm, level=c(.95))
  # winter.UD@data$ID <- paste(BirdID, Year, sep="_")
  # winter.UD@data$Bird <- BirdID
  # winter.UD@data$Season <- "Winter"
  # winter.UD@data$Year <- Year
  # 
  # #This if statement binds the individual layers into a new SpatialPolygonDataFrame
  # if(exists("winter.contours")){
  #   winter.UD.t <- spTransform(winter.UD, CRS(proj4string(winter.contours)))
  #   winter.contours <- rbind(winter.contours, winter.UD.t)
  # }else{
  #   winter.contours <- winter.UD
  # }
  
}

write.csv(roostpoints_all, "RoostPoints.csv", row.names = F) #Only roost points
write.csv(diurnal_all, "DaytimeWinterPoints.csv", row.names = F) #Just day points
write.csv(full_all, "AllUsedPoints.csv", row.names = F) #All points 
# writeOGR(winter.contours, dsn = 'E:/Maine Drive/Analysis/Kaj Thesis', #Save ind home ranges as shapefile
#          layer = "AllWinterUDs_95_line", driver = "ESRI Shapefile")

#Convert SpatialLinesDF to SpatialPolygonsDF
####You will need to Add both the Line and Polygon shapefiles and then####
####    "Join Field" in arcgis using FID to get the final product     ####
# winter.contours.polyset <- SpatialLines2PolySet(winter.contours)
# winter.polygons <- PolySet2SpatialPolygons(winter.contours.polyset, close_polys=TRUE)
# proj4string(winter.polygons) <- CRS(proj4string(winter.contours))
# winter.polygons.df <- as(winter.polygons, "SpatialPolygonsDataFrame")
# 
# writeOGR(winter.polygons.df, dsn = 'E:/Maine Drive/Analysis/Kaj Thesis', 
#          layer = "AllWinterUDs_95_polygon", driver = "ESRI Shapefile")


#############################################
### Generate Available Points for SSF/RSF ###
#############################################
# https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html#5_ssf_with_multiple_animals
pcks <- list("sp", "sf", "dplyr", "raster", "rgdal", "lubridate", "amt")
sapply(pcks, require, char = TRUE)

#Load used points
roostlocations.obs.raw <- read.csv("RoostPoints.csv")
daylocations.obs.raw <- read.csv("DaytimeWinterPoints.csv")
fulllocations.obs.raw <- read.csv("AllUsedPoints.csv")
# winter95UDs.raw <- readOGR(".", "AllWinterUDs_95_polygon")
# winter95UDs <- spTransform(winter95UDs.raw, CRS("+proj=longlat +datum=WGS84"))

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

fulllocations.obs <- fulllocations.obs.raw %>%
  dplyr::select(BirdID, timestamp, Lat = location_lat, Long = location_long) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(Year = year(timestamp)) %>% mutate(Used = 1) %>% mutate(Roost = 0) %>%
  dplyr::select(-timestamp)

### Available Locations for RSF ###
### Will end up with a DF that has both used and available points
### Will be able to filter by roost column depending on analysis
# rm(roostlocations.random)
# rm(daylocations.random)
# 
# for(i in 1:length(winter95UDs)){
#   Bird <- winter95UDs@data$Bird[i]
#   Yearobs <- winter95UDs@data$Year[i]
#   
#   #How many used roosts/day locations are there?
#   observed.roost <- nrow(roostlocations.obs %>% filter(BirdID == Bird) %>% filter(Year == Yearobs)) 
#   observed.day <- nrow(daylocations.obs %>% filter(BirdID == Bird) %>% filter(Year == Yearobs))
#   
#   roost.avail <- as.data.frame(spsample(winter95UDs[i,], observed.roost*200, "random")@coords) %>%
#     mutate(BirdID = Bird) %>% mutate(Year = Yearobs) %>% mutate(Used = 0) %>% mutate(Roost = 1) %>%
#     rename(Long = x, Lat = y)
#   
#   day.avail <- as.data.frame(spsample(winter95UDs[i,], observed.day*10, "random")@coords) %>%
#     mutate(BirdID = Bird) %>% mutate(Year = Yearobs) %>% mutate(Used = 0) %>% mutate(Roost = 0) %>%
#     rename(Long = x, Lat = y)
#   
#   if(exists("roostlocations.random")){
#     roostlocations.random <- rbind(roostlocations.random,roost.avail)
#     daylocations.random <- rbind(daylocations.random,day.avail)
#   }else{
#     roostlocations.random <- roost.avail
#     daylocations.random <- day.avail
#   }
# }
# 
# rsf.allpoints <- do.call("rbind", list(roostlocations.obs, daylocations.obs, roostlocations.random, daylocations.random)) %>%
#   arrange(Year, BirdID, Used, Roost)
# write.csv(rsf.allpoints, "RSF_Used_Avail_points.csv", row.names = F)
# rsf.allpoints.sp <- rsf.allpoints
# coordinates(rsf.allpoints.sp) <- ~Long+Lat
# proj4string(rsf.allpoints.sp) <- CRS("+proj=longlat +datum=WGS84")
# rsf.allpoints.sp <- spTransform(rsf.allpoints.sp, CRS = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")


# ### Roost Location SSF Available Points - Based on Previous Roost Location ###
# roostlocations.obs <- roostlocations.obs.raw %>%
#   dplyr::select(BirdID, timestamp, Lat = location_lat, Long = location_long) %>%
#   mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")) %>%
#   mutate(Year = year(timestamp)) %>%
#   mutate(ID = paste(BirdID, Year, sep = "_"))
# 
# nested.roosts <- roostlocations.obs %>%
#   nest(-ID) %>%
#   mutate(track = map(data, ~ mk_track(., Long, Lat, timestamp, crs= CRS("+proj=longlat +datum=WGS84"))))
# 
# roost.true.steps <- map(nested.roosts$track, steps)
# roost.random.steps <- map(roost.true.steps, random_steps, n=10)
# roost.PR.steps <- bind_rows(roost.random.steps, .id="ID")
# coordinates(roost.PR.steps) <- ~x2_+y2_
# proj4string(roost.PR.steps) <- CRS("+proj=longlat +datum=WGS84")


### All Locations SSF Available Points ###
fulllocations.obs <- fulllocations.obs.raw %>%
  dplyr::select(BirdID, timestamp, Lat = location_lat, Long = location_long) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS")) %>%
  mutate(Year = year(timestamp)) %>%
  mutate(ID = paste(BirdID, Year, sep = "_"))

nested.full <- fulllocations.obs %>%
  nest(-ID) %>%
  mutate(track = map(data, ~ mk_track(., Long, Lat, timestamp, crs= CRS("+proj=longlat +datum=WGS84"))))

full.true.steps <- map(nested.full$track, steps)
full.random.steps <- map(full.true.steps, random_steps, n=200, include_observed = T)
full.steps <- bind_rows(full.random.steps, .id="ID")
AllPoints.ssf <- full.steps
coordinates(AllPoints.ssf) <- ~x2_+y2_
proj4string(AllPoints.ssf) <- CRS("+proj=longlat +datum=WGS84")
AllPoints.ssf <- spTransform(AllPoints.ssf, CRS = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

###################################
### Load and Prepare Covariates ###
###################################
require(raster)

#Had to do some Raster editing
# BA.rast.raw <- raster("./Covariate Data/BA_19.tif")
# BA.rast <- projectRaster(BA.rast.raw, crs = crs(AllPoints.ssf))
# DtE.rast.raw <- raster("./Covariate Data/DisttoEdge.tif")
# DtE.rast <- projectRaster(DtE.rast.raw, crs = crs(AllPoints.ssf))
# Ht.rast.raw <- raster("./Covariate Data/HT_Mean_19.tif")
# Ht.rast <- projectRaster(Ht.rast.raw, crs = crs(AllPoints.ssf))
# SW.rast.raw <- raster("./Covariate Data/SW_Percent_19.tif")
# SW.rast <- projectRaster(SW.rast.raw, crs = crs(AllPoints.ssf))
# Ag.rast.raw <- raster("./Covariate Data/Agriculture.tif")
# Ag.rast <- projectRaster(Ag.rast.raw, crs = crs(AllPoints.ssf))
# Ag.rast.crop <- crop(Ag.rast, BA.rast, snap = "near")
# Dev.rast.raw <- raster("./Covariate Data/Developed.tif")
# Dev.rast <- projectRaster(Dev.rast.raw, crs = crs(AllPoints.ssf))
# Dev.rast.crop <- crop(Dev.rast, BA.rast, snap = "near")
# 
# writeRaster(Ag.rast.crop, "Agriculture_crschange.tif", format = "GTiff")
# writeRaster(Dev.rast.crop, "Developed_crschange.tif", format = "GTiff")
# writeRaster(SW.rast, "SW_Percent_19_crschange.tif", format = "GTiff")
# writeRaster(Ht.rast, "HT_Mean_19_crschange.tif", format = "GTiff")
# writeRaster(DtE.rast, "DtE_19_crschange.tif", format = "GTiff")
# writeRaster(BA.rast, "BA_19_crschange.tif", format = "GTiff")

#Load Rasters
BA.rast <- raster("./Covariate Data/BA_final_0.tif")
DtE.rast <- raster("./Covariate Data/DtE_final.tif")
Ht.rast <- raster("./Covariate Data/Ht_final.tif")
SW.rast <- raster("./Covariate Data/SW_final_0.tif")
Ag.rast <- raster("./Covariate Data/AG_final.tif")
Dev.rast <- raster("./Covariate Data/DEV_final.tif")

#Extract Spatial Covariates
#SSF
BA.full.step.cov <- raster::extract(BA.rast, AllPoints.ssf)
BA.full.step.cov[is.na(BA.full.step.cov)] <- 0
DtE.full.step.cov <- raster::extract(DtE.rast, AllPoints.ssf)
DtE.full.step.cov[is.na(DtE.full.step.cov)] <- 0
Ht.full.step.cov <- raster::extract(Ht.rast, AllPoints.ssf)
Ht.full.step.cov[is.na(Ht.full.step.cov)] <- 0
SW.full.step.cov <- raster::extract(SW.rast, AllPoints.ssf)
SW.full.step.cov[is.na(SW.full.step.cov)] <- 0
AG.full.step.cov <- raster::extract(Ag.rast, AllPoints.ssf)
AG.full.step.cov[is.na(AG.full.step.cov)] <- 0
DEV.full.step.cov <- raster::extract(Dev.rast, AllPoints.ssf)
DEV.full.step.cov[is.na(DEV.full.step.cov)] <- 0
# #RSF
# BA.full.rsf.cov <- raster::extract(BA.rast, rsf.allpoints.sp)
# BA.full.rsf.cov[is.na(BA.full.rsf.cov)] <- 0
# DtE.full.rsf.cov <- raster::extract(DtE.rast, rsf.allpoints.sp)
# DtE.full.rsf.cov[is.na(DtE.full.rsf.cov)] <- 0
# Ht.full.rsf.cov <- raster::extract(Ht.rast, rsf.allpoints.sp)
# Ht.full.rsf.cov[is.na(Ht.full.rsf.cov)] <- 0
# SW.full.rsf.cov <- raster::extract(SW.rast, rsf.allpoints.sp)
# SW.full.rsf.cov[is.na(SW.full.rsf.cov)] <- 0
# AG.full.rsf.cov <- raster::extract(Ag.rast, rsf.allpoints.sp)
# AG.full.rsf.cov[is.na(AG.full.rsf.cov)] <- 0
# DEV.full.rsf.cov <- raster::extract(Dev.rast, rsf.allpoints.sp)
# DEV.full.rsf.cov[is.na(DEV.full.rsf.cov)] <- 0


#Add SP Objects
#SSF
AllPoints.ssf$BA <- BA.full.step.cov
AllPoints.ssf$DtE <- DtE.full.step.cov
AllPoints.ssf$Ht <- Ht.full.step.cov
AllPoints.ssf$SW <- SW.full.step.cov
AllPoints.ssf$AG <- AG.full.step.cov
AllPoints.ssf$DEV <- DEV.full.step.cov
# #RSF
# rsf.allpoints.sp$BA <- BA.full.rsf.cov
# rsf.allpoints.sp$DtE <- DtE.full.rsf.cov
# rsf.allpoints.sp$Ht <- Ht.full.rsf.cov
# rsf.allpoints.sp$SW <- SW.full.rsf.cov
# rsf.allpoints.sp$AG <- AG.full.rsf.cov
# rsf.allpoints.sp$DEV <- DEV.full.rsf.cov

##########################
### Weather Covariates ###
##########################
weather.raw <- read.csv("WeatherVariables.csv") %>%
  mutate(DATE = as.Date(DATE, format = "%m/%d/%Y")) %>%
  filter(yday(DATE) < 75  | yday(DATE) > 364) %>%
  select(-STATION, -NAME, -DATE, -WT05, -WT03, -WESF, -WDF5, -WDF2,
         -WT01, -WT02, -WT04, -WT06, -WT08, -WT09, -TMAX, -TAVG, -WSF5, -WSF2,
         -WESD, -TSUN, -PSUN, -PGTM, -DAPR, -MDPR)
# require(missMDA)
# nb <- estim_ncpPCA(weather.raw, method.cv = "Kfold", verbose = F)
# nb$ncp
# weather.comp <- as.matrix(imputePCA(weather.raw, ncp = nb$ncp))
# weather.pca <- prcomp(t(weather.comp))
# weather.pca <- prcomp((weather.comp))



###########################
### Final Data Clean Up ###
###########################
#Filter down the extra points you generated to avoid BA == NA for roost locations
AllPoints.ssf.final <- as.data.frame(AllPoints.ssf@data)  %>%
  group_by(ID, case_, t2_) %>%
  filter(hour(t2_) %notin% 4:7 | BA != 0) %>%
  slice(1:10)

#Make sure you got 10 available for each used points
AllPoints.ssf.test <- AllPoints.ssf.final  %>%
  group_by(ID, t2_) %>%
  summarize(Total = n()) %>%
  filter(Total != 11)




#Create Roost Site Objects
#SSF
RoostPoints.ssf.final <- AllPoints.ssf.final  %>%
  mutate(Roost = ifelse(hour(t2_) %in% 4:7, 1, 0)) %>%
  filter(Roost == 1)
# #RSF
# rsf.roostpoints.sp <- rsf.allpoints.sp %>%
#   filter(Roost == 1)

#Z-Standardize Covariates
#SSF
BA.z.full.step.cov <- as.numeric(scale(BA.full.step.cov, center = T, scale = T))
DtE.z.full.step.cov <- as.numeric(scale(DtE.full.step.cov, center = T, scale = T))
Ht.z.full.step.cov <- as.numeric(scale(Ht.full.step.cov, center = T, scale = T))
SW.z.full.step.cov <- as.numeric(scale(SW.full.step.cov, center = T, scale = T))
AG.z.full.step.cov <- as.numeric(scale(AG.full.step.cov, center = T, scale = T))
DEV.z.full.step.cov <- as.numeric(scale(DEV.full.step.cov, center = T, scale = T))
# #RSF
# BA.z.full.rsf.cov <- as.numeric(scale(BA.full.rsf.cov, center = T, scale = T))
# DtE.z.full.rsf.cov <- as.numeric(scale(DtE.full.rsf.cov, center = T, scale = T))
# Ht.z.full.rsf.cov <- as.numeric(scale(Ht.full.rsf.cov, center = T, scale = T))
# SW.z.full.rsf.cov <- as.numeric(scale(SW.full.rsf.cov, center = T, scale = T))
# AG.z.full.rsf.cov <- as.numeric(scale(AG.full.rsf.cov, center = T, scale = T))
# DEV.z.full.rsf.cov <- as.numeric(scale(DEV.full.rsf.cov, center = T, scale = T))

########################
### SSF ROOST MODELS ###
########################
require(survival)
R.SSF.null <- clogit(case_ ~ 1, method = "approximate", data = RoostPoints.ssf)
R.SSF.1 <- clogit(case_ ~ ID|ID, method = "approximate", data = RoostPoints.ssf)
