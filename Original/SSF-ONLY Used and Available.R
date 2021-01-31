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

### Loop to download points, create dBBMM home range, and create df of all locations
### Only Looking at data between January 1 through March 15
rm(full_all)#Clear objects so you don't duplicate data (for if not starting new session)

'%notin%' <- Negate('%in%') #Need the opposite of %in%

for(i in 1:nrow(birdlist)){
  BirdID <- as.character(birdlist[i,1])
  Year <- birdlist[i,2]
  StartDate <- paste(Year, "0101000000000", sep = "")
  EndDate <- paste(Year, "0316000000000", sep = "")
  
  #Get all points from Jan1 through March15
  winterpoints_ind <- getMovebankData(study = "Eastern Wild Turkey, Gonnerman, Maine", login = login, 
                                      animalName = BirdID, timestamp_start = StartDate, timestamp_end = EndDate)
  
  full_ind <- winterpoints_ind@data %>%  
    mutate(Time = hour(timestamp)) %>%
    mutate(BirdID = BirdID)
  
  if(i==1){ #rbind ind bird data to create one large df
    full_all <- full_ind
  }else{
    full_all <- rbind(full_all, full_ind)
  }
}

write.csv(full_all, "AllUsedPoints.csv", row.names = F) #All points 


#############################################
### Generate Available Points for SSF/RSF ###
#############################################
# https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html#5_ssf_with_multiple_animals
pcks <- list("sp", "sf", "dplyr", "raster", "rgdal", "lubridate", "amt")
sapply(pcks, require, char = TRUE)

#Load used points
fulllocations.obs.raw <- read.csv("AllUsedPoints.csv")

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


step.lengths <- unlist(map(nested.full$track, step_lengths))
quantile(step.lengths, c(.5, .75, .9, .95), na.rm = T)
#363m
hist(step.lengths, breaks = 20, main="Histogram of step lengths", xlab="Step lengths")
###################################
### Load and Prepare Covariates ###
###################################
require(raster)

#Rasters were prepared in ArcGIS
#Load Rasters
BA.rast <- raster("./Covariate Data/BA_final_0.tif")
DtFE.rast <- raster("./Covariate Data/DtFE_final.tif")
Ht.rast <- raster("./Covariate Data/Ht_final.tif")
SW.rast <- raster("./Covariate Data/SW_final_0.tif")
Ag.rast <- raster("./Covariate Data/AG_final.tif")
Dev.rast <- raster("./Covariate Data/DEV_final.tif")
Aspect.rast <- raster("./Covariate Data/Aspect_final.tif")
NLCD.rast <- raster("./Covariate Data/NLCDClass_final.tif")
PropAg.rast <- raster("./Covariate Data/PropAg_final.tif")
PropDev.rast <- raster("./Covariate Data/PropDev_final.tif")
PropSW.rast <- raster("./Covariate Data/PropSW_final.tif")

#Extract Spatial Covariates
BA.full.step.cov <- raster::extract(BA.rast, AllPoints.ssf)
DtFE.full.step.cov <- raster::extract(DtFE.rast, AllPoints.ssf)
Ht.full.step.cov <- raster::extract(Ht.rast, AllPoints.ssf)
SW.full.step.cov <- raster::extract(SW.rast, AllPoints.ssf)
AG.full.step.cov <- raster::extract(Ag.rast, AllPoints.ssf)
DEV.full.step.cov <- raster::extract(Dev.rast, AllPoints.ssf)
Aspect.full.step.cov <- raster::extract(Aspect.rast, AllPoints.ssf)
NLCDclass.full.step.cov <- raster::extract(NLCD.rast, AllPoints.ssf)
PropAg.full.step.cov <- raster::extract(PropAg.rast, AllPoints.ssf)
PropDev.full.step.cov <- raster::extract(PropDev.rast, AllPoints.ssf)
PropSW.full.step.cov <- raster::extract(PropSW.rast, AllPoints.ssf)

#Add covariates to SP Object as attributes
AllPoints.ssf$BA <- BA.full.step.cov
AllPoints.ssf$DtFE <- DtFE.full.step.cov
AllPoints.ssf$Ht <- Ht.full.step.cov
AllPoints.ssf$SW <- SW.full.step.cov
AllPoints.ssf$AG <- AG.full.step.cov
AllPoints.ssf$DEV <- DEV.full.step.cov
AllPoints.ssf$Aspect <- Aspect.full.step.cov
AllPoints.ssf$NLCDClass <- as.factor(NLCDclass.full.step.cov)
AllPoints.ssf$PropAg <- PropAg.full.step.cov #360m
AllPoints.ssf$PropDev <- PropDev.full.step.cov #360m
AllPoints.ssf$PropSW <- PropSW.full.step.cov #360m
AllPoints.ssf$PropFoodSub <- PropAg.full.step.cov + PropDev.full.step.cov


###############################
### Prep Weather Covariates ###
###############################
#Daily weather summaries downloaded for Bangor Airport weather station
weather.raw <- read.csv("WeatherVariables.csv") %>%
  mutate(DATE = as.Date(DATE, format = "%m/%d/%Y")) %>%
  mutate(AWND = ifelse(is.na(AWND), 0, AWND)) %>% #If AWND is NA, that means no wind
  dplyr::select(Date = DATE, TMIN, AWND, WDF2, SD = SNWD) %>%
  mutate(WC = 13.12 + (.6215*TMIN)-(11.37*(AWND^0.16))+(.3965*TMIN*(AWND^0.16))) %>% #calculate windchill
  mutate(WC_prev = lag(WC, 1)) #want column with previous days wind chill


###########################
### Final Data Clean Up ###
###########################
AllPoints.ssf.df <- AllPoints.ssf@data %>%
  mutate(Date = as.Date(t2_)) %>%
  mutate(Use = as.numeric(case_)) %>%
  mutate(ID = as.factor(ID))

AllPoints.ssf.ogpoints <- cbind(AllPoints.ssf.df, AllPoints.ssf@coords) #in case you want the actual location points
AllPoints.ssf.ogpoints$ConditionalID <- paste(AllPoints.ssf.ogpoints$ID, AllPoints.ssf.ogpoints$step_id_, sep="_")


##Merge Weather
#Add windchill and snow depth by date
AllPoints.ssf.weather <- merge(AllPoints.ssf.ogpoints, weather.raw, by = "Date", all.x = T)


##Aspect
#Calculate difference in Aspect and wind direction to come up with wind exposure metric
AllPoints.ssf.final <- AllPoints.ssf.weather %>%
  mutate(Wind_Exp = ifelse(Aspect == -1, 90,
                           pmin(pmax(Aspect, WDF2) - pmin(Aspect, WDF2), 360 + pmin(Aspect, WDF2) - pmax(Aspect, WDF2)))) %>%
  mutate(Wind_Exp = ifelse(is.na(WDF2), 0, Wind_Exp))

##NLCD
#Need to change the categories for NLCD classification of point
AllPoints.ssf.final <- AllPoints.ssf.final %>%
  filter(NLCDClass != 11) %>%
  mutate(NLCD_Point = as.factor(ifelse(NLCDClass %in% c("21", "22", "23", "24"), "Dev", 
                             ifelse(NLCDClass %in% c("41", "42", "43"), "For", 
                                    ifelse(NLCDClass %in% c("71", "81", "82", "52", "31"), "Open", 
                                           ifelse(NLCDClass %in% c("90", "95"), "Wet", NLCDClass))))))


##Multiple Data Objects needed
#1st scale = all daytime points, any location
DayAll.final <- AllPoints.ssf.final  %>%
  group_by(ID, case_, t2_) %>%
  filter(hour(t2_) %notin% 4:7) %>%
  slice(1:10) %>%
  as.data.frame()

#2nd scale = daytime points, only forest
#need to first identify which used points were in forest
#then limit the available for the used forest points to also be in forest
DayForest.used <- AllPoints.ssf.final  %>%
  filter(hour(t2_) %notin% 4:7) %>% #Only daytime
  filter(BA > 0) %>% #Only in forest
  filter(case_ == T) %>% #Only used
  mutate(MergeID = paste(ID, step_id_, sep="_"))
DayForest.avail <- AllPoints.ssf.final%>%
  group_by(ID, case_, t2_) %>%
  filter(hour(t2_) %notin% 4:7) %>% #Only daytime
  filter(case_ == F) %>% #only available
  filter(BA > 0) %>% #only in forest
  mutate(MergeID = paste(ID, step_id_, sep="_")) %>%
  filter(MergeID %in% unique(DayForest.used$MergeID)) %>% #Only associated with used forest sites
  slice(1:10) %>%
  ungroup()
DayForest.final <- rbind(DayForest.used, DayForest.avail) %>%
  arrange(ID, t2_, case_) %>%
  select(-MergeID) %>%
  as.data.frame()

#3rd scale = only roost points
#Available must match criteria of used
RoostAll.final <- AllPoints.ssf.final  %>%
  group_by(ID, case_, t2_) %>%
  filter(hour(t2_) %in% 4:7) %>%
  filter(BA > 0) %>%
  slice(1:10) %>%
  as.data.frame()

##Z-Standardize Covariates
#All Daytime Points
DayAll.final$BA.Z <- as.numeric(scale(DayAll.final$BA, center = T, scale = T))
DayAll.final$DtFE.Z <- as.numeric(scale(DayAll.final$DtFE, center = T, scale = T))
DayAll.final$Ht.Z <- as.numeric(scale(DayAll.final$Ht, center = T, scale = T))
DayAll.final$SW.Z <- as.numeric(scale(DayAll.final$SW, center = T, scale = T))
DayAll.final$PropAg.Z <- as.numeric(scale(DayAll.final$PropAg, center = T, scale = T))
DayAll.final$PropDev.Z <- as.numeric(scale(DayAll.final$PropDev, center = T, scale = T))
DayAll.final$PropSW.Z <- as.numeric(scale(DayAll.final$PropSW, center = T, scale = T))
DayAll.final$SD.Z <- as.numeric(scale(DayAll.final$SD, center = T, scale = T))
DayAll.final$WC.Z <- as.numeric(scale(DayAll.final$WC, center = T, scale = T))
DayAll.final$WC_prev.Z <- as.numeric(scale(DayAll.final$WC_prev, center = T, scale = T))
DayAll.final$Wind.Exp.Z <- as.numeric(scale(DayAll.final$Wind_Exp, center = T, scale = T))
DayAll.final$StepLength.Z <- as.numeric(scale(DayAll.final$sl_, center = T, scale = T))
DayAll.final$PropFoodSub.Z <- as.numeric(scale(DayAll.final$PropFoodSub, center = T, scale = T))

#Forested Daytime Points
DayForest.final$BA.Z <- as.numeric(scale(DayForest.final$BA, center = T, scale = T))
DayForest.final$DtFE.Z <- as.numeric(scale(DayForest.final$DtFE, center = T, scale = T))
DayForest.final$Ht.Z <- as.numeric(scale(DayForest.final$Ht, center = T, scale = T))
DayForest.final$SW.Z <- as.numeric(scale(DayForest.final$SW, center = T, scale = T))
DayForest.final$PropAg.Z <- as.numeric(scale(DayForest.final$PropAg, center = T, scale = T))
DayForest.final$PropDev.Z <- as.numeric(scale(DayForest.final$PropDev, center = T, scale = T))
DayForest.final$PropSW.Z <- as.numeric(scale(DayForest.final$PropSW, center = T, scale = T))
DayForest.final$SD.Z <- as.numeric(scale(DayForest.final$SD, center = T, scale = T))
DayForest.final$WC.Z <- as.numeric(scale(DayForest.final$WC, center = T, scale = T))
DayForest.final$WC_prev.Z <- as.numeric(scale(DayForest.final$WC_prev, center = T, scale = T))
DayForest.final$Wind.Exp.Z <- as.numeric(scale(DayForest.final$Wind_Exp, center = T, scale = T))
DayForest.final$StepLength.Z <- as.numeric(scale(DayForest.final$sl_, center = T, scale = T))
DayForest.final$PropFoodSub.Z <- as.numeric(scale(DayForest.final$PropFoodSub, center = T, scale = T))

#Roost Points
RoostAll.final$BA.Z <- as.numeric(scale(RoostAll.final$BA, center = T, scale = T))
RoostAll.final$DtFE.Z <- as.numeric(scale(RoostAll.final$DtFE, center = T, scale = T))
RoostAll.final$Ht.Z <- as.numeric(scale(RoostAll.final$Ht, center = T, scale = T))
RoostAll.final$SW.Z <- as.numeric(scale(RoostAll.final$SW, center = T, scale = T))
RoostAll.final$PropAg.Z <- as.numeric(scale(RoostAll.final$PropAg, center = T, scale = T))
RoostAll.final$PropDev.Z <- as.numeric(scale(RoostAll.final$PropDev, center = T, scale = T))
RoostAll.final$PropSW.Z <- as.numeric(scale(RoostAll.final$PropSW, center = T, scale = T))
RoostAll.final$SD.Z <- as.numeric(scale(RoostAll.final$SD, center = T, scale = T))
RoostAll.final$WC.Z <- as.numeric(scale(RoostAll.final$WC, center = T, scale = T))
RoostAll.final$WC_prev.Z <- as.numeric(scale(RoostAll.final$WC_prev, center = T, scale = T))
RoostAll.final$Wind.Exp.Z <- as.numeric(scale(RoostAll.final$Wind_Exp, center = T, scale = T))
RoostAll.final$StepLength.Z <- as.numeric(scale(RoostAll.final$sl_, center = T, scale = T))
RoostAll.final$PropFoodSub.Z <- as.numeric(scale(RoostAll.final$PropFoodSub, center = T, scale = T))

#Create Stratum ID (may already have one but this is numeric)
RoostAll.final <- RoostAll.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))

DayForest.final <- DayForest.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))

DayAll.final <- DayAll.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))


write.csv(RoostAll.final, "RoostAll.csv", row.names = F)
write.csv(DayForest.final, "DayForest.csv", row.names = F)
write.csv(DayAll.final, "DayAll.csv", row.names = F)

