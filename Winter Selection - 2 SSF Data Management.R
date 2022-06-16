#### This is to consolidate all the code for the SSF/RSF Prep
#### Need Used and available points for a roost specific and all points analysis
#### For roost specific, get two sets of available points for 2 scales of analysis
#### 1)Available based on roost to roost, 2) Available based on point previous to roost
#### UPDATE: Could not get option 1 to work, and since Erik wanted 2 anyways, didnt push through that issue yet

# #Install Packages
# install.packages("move", "ggmap", "mapproj", "sf", "devtools")
# library(devtools)
# install_github("bacollier/moveud")

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
    mutate(BirdID = BirdID) %>%
    mutate(timestamp = with_tz(timestamp, tzone = "America/New_York"))
  
  if(i==1){ #rbind ind bird data to create one large df
    full_all <- full_ind
  }else{
    full_all <- rbind(full_all, full_ind)
  }
}

#####################################
### Combine SSF Used and HMM Data ###
#####################################
hmm_data.raw <- read.csv("Results/HMMBehavioralStates_output.csv") %>%
  select(ID, Timestamp, State) %>%
  mutate(BirdID =  sub('^[^X]*X(\\d+).*', '\\1', ID)) %>%
  rename(timestamp = Timestamp) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "EST"))

#Merge by ID and nearest date/time
require(data.table)
setDT(full_all)
setDT(hmm_data.raw)
setkey(full_all, BirdID, timestamp)[,dateMatch:=timestamp]
setkey(hmm_data.raw, BirdID, timestamp)[,dateMatch:=timestamp]
full_all_hmm2 <- hmm_data.raw[full_all, roll = 'nearest'] %>%
  filter(!is.na(sensor_type_id)) %>%
  select(ID, BirdID, location_lat, location_long, timestamp, State) %>% 
  arrange(ID, timestamp) %>%
  group_by(ID) %>%
  mutate(sincelast = difftime(timestamp, lag(timestamp), units = "secs")) %>%
  mutate(tillnext = difftime(lead(timestamp),timestamp, units = "secs")) %>%
  mutate(keep = ifelse(abs(sincelast) > 120 & abs(tillnext) > 120, 1,
                       ifelse(minute(timestamp) != 0 & second(timestamp) != 0, 1, 0))) %>%
  ungroup() %>% 
  arrange(ID, timestamp) %>%
  filter(keep == 1) %>%
  dplyr::select(-keep, -sincelast, -tillnext)

write.csv(full_all_hmm2, "AllUsedPoints.csv", row.names = F) #All points 

# hmm_diff_check <- hmm_data.raw[full_all, roll = 'nearest'] %>%
#   filter(!is.na(sensor_type_id)) %>%
#   mutate(Time_Diff = timestamp - dateMatch)


#############################################
### Generate Available Points for SSF/RSF ###
#############################################
# https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html#5_ssf_with_multiple_animals
pcks <- list("sp", "sf", "dplyr", "raster", "rgdal", "lubridate", "amt")
sapply(pcks, require, char = TRUE)

#Load used points
fulllocations.obs.raw <- read.csv("AllUsedPoints.csv") %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")) %>%
  rename(Lat = location_lat, Long = location_long) %>%
  mutate(State = as.factor(State))

#Save them as a shapefile for easier viewing
usedlocations.sp <- fulllocations.obs.raw
coordinates(usedlocations.sp) <- ~Long+Lat
proj4string(usedlocations.sp) <- CRS("+proj=longlat +datum=WGS84")
usedlocations.sp <- spTransform(usedlocations.sp, CRS = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(usedlocations.sp, ".", layer = "usedlocations", driver = "ESRI Shapefile", overwrite_layer = T)


### All Locations SSF Available Points ###
nested.full <- fulllocations.obs.raw %>%
  nest(-ID) %>%
  mutate(track = map(data, ~ mk_track(., Long, Lat, timestamp, State, BirdID, crs= CRS("+proj=longlat +datum=WGS84"))))
full.true.steps <- map(nested.full$track, steps, keep_cols = 'end')

#Separate stationary and mobile points to create separate distributions for random points
s.true.steps <- lapply(full.true.steps, function(x) filter(x, State == 2))
m.true.steps <- lapply(full.true.steps, function(x) filter(x, State == 3))

#Need to rewrap roosts to get new step lenght/turning anlge values
nested.r <- fulllocations.obs.raw %>%
  filter(State == 1) %>%
  nest(-ID) %>%
  mutate(track = map(data, ~ mk_track(., Long, Lat, timestamp, State, BirdID, crs= CRS("+proj=longlat +datum=WGS84"))))
r.true.steps <- map(nested.r$track, steps, keep_cols = 'end')

r.dist <- do.call("rbind", r.true.steps) 
s.dist <- do.call("rbind", s.true.steps) 
m.dist <- do.call("rbind", m.true.steps) 

full.random.steps.r <- map(r.true.steps, random_steps, n=200, include_observed = T,
                           sl_distr = fit_distr(r.dist$sl_, "gamma"), ta_distr = fit_distr(r.dist$ta_, "vonmises"))
full.random.steps.s <- map(s.true.steps, random_steps, n=200, include_observed = T,
                           sl_distr = fit_distr(s.dist$sl_, "gamma"), ta_distr = fit_distr(s.dist$ta_, "vonmises"))
full.random.steps.m <- map(m.true.steps, random_steps, n=200, include_observed = T,
                           sl_distr = fit_distr(m.dist$sl_, "gamma"), ta_distr = fit_distr(m.dist$ta_, "vonmises"))

full.steps.r <- bind_rows(full.random.steps.r, .id="ID")
full.steps.s <- bind_rows(full.random.steps.s, .id="ID")
full.steps.m <- bind_rows(full.random.steps.m, .id="ID")

full.steps <- rbind(full.steps.r, full.steps.m, full.steps.s) %>%
  arrange(BirdID, t1_, desc(case_))

AllPoints.ssf <- full.steps %>%
  dplyr::select(-dt_)
coordinates(AllPoints.ssf) <- ~x2_+y2_
proj4string(AllPoints.ssf) <- CRS("+proj=longlat +datum=WGS84")
AllPoints.ssf <- spTransform(AllPoints.ssf, CRS = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# writeOGR(AllPoints.ssf, ".", layer = "alllocations", driver = "ESRI Shapefile", overwrite_layer = T)
# step.lengths <- unlist(map(nested.full$track, step_lengths))
# quantile(step.lengths, c(.5, .75, .9, .95), na.rm = T)
# #363m
# hist(step.lengths, breaks = 20, main="Histogram of step lengths", xlab="Step lengths")


### Separate Availability for Roost-to-Roost analysis

# roostonly.sp <- fulllocations.obs.raw %>%
#   filter(State == 1) %>%
#   mutate(Hour = hour(timestamp)) %>% #Just to check 
#   nest(-ID) %>%
#   mutate(track = map(data, ~ mk_track(., Long, Lat, timestamp, State, BirdID, crs= CRS("+proj=longlat +datum=WGS84"))))
# roostonly.true.steps <- map(roostonly.sp$track, steps, keep_cols = 'end')
# roostonly.random.steps <- map(roostonly.true.steps, random_steps, n=200, include_observed = T)
# roostonly.steps <- bind_rows(roostonly.random.steps, .id="ID")
# RoostOnlyPoints.ssf <- roostonly.steps %>%
#   dplyr::select(-dt_)
# coordinates(RoostOnlyPoints.ssf) <- ~x2_+y2_
# proj4string(RoostOnlyPoints.ssf) <- CRS("+proj=longlat +datum=WGS84")
# RoostOnlyPoints.ssf <- spTransform(RoostOnlyPoints.ssf, CRS = "+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")




###################################
### Load and Prepare Covariates ###
###################################
require(raster)

#Load Rasters
#Rasters were prepared in ArcGIS
#For continuous covariates where cell size was less than GPS error (BA, DtFE, Ht, SW)
#Performed moving window and averaged the cell values for 1 cell out in any direction
#Clip raster to smaller size to make easier to work with (snap to OG raster, extent/proj to same raster)
#Focal statistic (mean) for prop, 3x3 rect for NCLD (cell size 30), 9x9 for cell size  centered on point.

BA.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/BA_final_0.tif")
DtFE.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/DtFE_final.tif")
Ht.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/Ht_final.tif")
SW.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/SW_final_0.tif")
# Ag.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/AG_final.tif")
# Dev.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/DEV_final.tif")
Aspect.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/Aspect_final.tif")
NLCD.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/NLCDClass_final.tif")
PropAg.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/PropAg_final.tif")
PropDev.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/PropDev_final.tif")
PropSW.rast <- raster("E:/Maine Drive/Analysis/Kaj Thesis/Covariate Data/PropSW_final.tif")

#Extract Spatial Covariates
BA.full.step.cov <- raster::extract(BA.rast, AllPoints.ssf)
DtFE.full.step.cov <- raster::extract(DtFE.rast, AllPoints.ssf)
Ht.full.step.cov <- raster::extract(Ht.rast, AllPoints.ssf)
SW.full.step.cov <- raster::extract(SW.rast, AllPoints.ssf)
# AG.full.step.cov <- raster::extract(Ag.rast, AllPoints.ssf)
# DEV.full.step.cov <- raster::extract(Dev.rast, AllPoints.ssf)
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
# AllPoints.ssf$AG <- AG.full.step.cov
# AllPoints.ssf$DEV <- DEV.full.step.cov
AllPoints.ssf$Aspect <- Aspect.full.step.cov
AllPoints.ssf$NLCDClass <- as.factor(NLCDclass.full.step.cov)
AllPoints.ssf$PropAg <- PropAg.full.step.cov #1 cell out
AllPoints.ssf$PropDev <- PropDev.full.step.cov #1 cell out
AllPoints.ssf$PropSW <- PropSW.full.step.cov #1 cell out
AllPoints.ssf$PropFoodSub <- PropAg.full.step.cov + PropDev.full.step.cov


###############################
### Prep Weather Covariates ###
###############################
#Daily weather summaries downloaded for Bangor Airport weather station
weather.raw <- read.csv("./Data/WeatherVariables.csv") %>%
  mutate(WC_prev = lag(WC, 1)) %>% #want column with previous days wind chill
  mutate(WC_prev = ifelse(is.na(WC_prev), lead(WC_prev), WC_prev),
         BirdID = stringr::str_extract(ID, "(?<=X).*?(?=\\_)"),
         Hour = lubridate::hour(Timestamp),
         Date = as.Date(Timestamp)) %>%
  select(BirdID, Date, Hour, SD, WC, WC_prev, WDF2 = WindDir)

###########################
### Final Data Clean Up ###
###########################
AllPoints.ssf.df <- AllPoints.ssf@data %>%
  mutate(Date = as.Date(t1_),
         Use = as.numeric(case_),
         ID = as.factor(ID),
         Hour = lubridate::hour(t1_))


AllPoints.ssf.ogpoints <- cbind(AllPoints.ssf.df, AllPoints.ssf@coords) #in case you want the actual location points
AllPoints.ssf.ogpoints$ConditionalID <- paste(AllPoints.ssf.ogpoints$ID, AllPoints.ssf.ogpoints$step_id_, sep="_")


##Merge Weather
#Add windchill and snow depth by date
AllPoints.ssf.weather <- merge(AllPoints.ssf.ogpoints, weather.raw, by = c("BirdID", "Date", "Hour"), all.x = T)


##Aspect
#Calculate difference in Aspect and wind direction to come up with wind exposure metric
AllPoints.ssf.final <- AllPoints.ssf.weather %>%
  mutate(Wind_Exp = ifelse(Aspect == -1, 90,
                           pmin(pmax(Aspect, WDF2) - pmin(Aspect, WDF2), 360 + pmin(Aspect, WDF2) - pmax(Aspect, WDF2)))) %>%
  mutate(Wind_Exp = ifelse(is.na(WDF2), 0, Wind_Exp))

##NLCD
#Need to change the categories for NLCD classification of point
AllPoints.ssf.final <- AllPoints.ssf.final %>%
  mutate(NLCD_Point = as.factor(ifelse(NLCDClass %in% c("21", "22", "23", "24"), "Dev", 
                                       ifelse(NLCDClass %in% c("41", "42", "43"), "For", 
                                              ifelse(NLCDClass %in% c("71", "81", "82", "52", "31"), "Open", 
                                                     ifelse(NLCDClass %in% c("90", "95"), "Wet", 
                                                            ifelse(NLCDClass == "11", "Water", NLCDClass)))))))

##Multiple Data Objects needed
Foraging.final <- AllPoints.ssf.final  %>%
  filter(State == 3) %>%
  group_by(ID, step_id_, case_) %>%
  slice(1:10) %>%
  as.data.frame()

foraging.use <- Foraging.final %>% filter(Use == 1)
foraging.use.vec <- foraging.use$ConditionalID
Foraging.final <- Foraging.final %>% filter(ConditionalID %in% foraging.use.vec)

Loafing.final <- AllPoints.ssf.final  %>%
  filter(State == 2) %>%
  group_by(ID, step_id_, case_) %>%
  slice(1:10) %>%
  as.data.frame()

loafing.use <- Loafing.final %>% filter(Use == 1)
loafing.use.vec <- loafing.use$ConditionalID
Loafing.final <- Loafing.final %>% filter(ConditionalID %in% loafing.use.vec)

Roost.final <- AllPoints.ssf.final  %>%
  filter(State == 1) %>%
  filter(BA > 0) %>%
  group_by(ID, step_id_, case_) %>%
  slice(1:10) %>%
  as.data.frame()

roost.use <- Roost.final %>% filter(Use == 1)
roost.use.vec <- roost.use$ConditionalID
Roost.final <- Roost.final %>% filter(ConditionalID %in% roost.use.vec)

AllPoints.final <- AllPoints.ssf.final  %>%
  group_by(ID, step_id_, case_) %>%
  slice(1:10) %>%
  as.data.frame()

all.use <- AllPoints.final %>% filter(Use == 1)
all.use.vec <- all.use$ConditionalID
AllPoints.final <- AllPoints.final %>% filter(ConditionalID %in% all.use.vec)


##Z-Standardize Covariates
#All Daytime Points
AllPoints.final$BA.Z <- as.numeric(scale(AllPoints.final$BA, center = T, scale = T))
AllPoints.final$DtFE.Z <- as.numeric(scale(AllPoints.final$DtFE, center = T, scale = T))
AllPoints.final$Ht.Z <- as.numeric(scale(AllPoints.final$Ht, center = T, scale = T))
AllPoints.final$SW.Z <- as.numeric(scale(AllPoints.final$SW, center = T, scale = T))
AllPoints.final$PropAg.Z <- as.numeric(scale(AllPoints.final$PropAg, center = T, scale = T))
AllPoints.final$PropDev.Z <- as.numeric(scale(AllPoints.final$PropDev, center = T, scale = T))
AllPoints.final$PropSW.Z <- as.numeric(scale(AllPoints.final$PropSW, center = T, scale = T))
AllPoints.final$SD.Z <- as.numeric(scale(AllPoints.final$SD, center = T, scale = T))
AllPoints.final$WC.Z <- as.numeric(scale(AllPoints.final$WC, center = T, scale = T))
AllPoints.final$WC_prev.Z <- as.numeric(scale(AllPoints.final$WC_prev, center = T, scale = T))
AllPoints.final$Wind.Exp.Z <- as.numeric(scale(AllPoints.final$Wind_Exp, center = T, scale = T))
AllPoints.final$StepLength.Z <- as.numeric(scale(AllPoints.final$sl_, center = T, scale = T))
AllPoints.final$PropFoodSub.Z <- as.numeric(scale(AllPoints.final$PropFoodSub, center = T, scale = T))


#All Daytime Points
Foraging.final$BA.Z <- as.numeric(scale(Foraging.final$BA, center = T, scale = T))
Foraging.final$DtFE.Z <- as.numeric(scale(Foraging.final$DtFE, center = T, scale = T))
Foraging.final$Ht.Z <- as.numeric(scale(Foraging.final$Ht, center = T, scale = T))
Foraging.final$SW.Z <- as.numeric(scale(Foraging.final$SW, center = T, scale = T))
Foraging.final$PropAg.Z <- as.numeric(scale(Foraging.final$PropAg, center = T, scale = T))
Foraging.final$PropDev.Z <- as.numeric(scale(Foraging.final$PropDev, center = T, scale = T))
Foraging.final$PropSW.Z <- as.numeric(scale(Foraging.final$PropSW, center = T, scale = T))
Foraging.final$SD.Z <- as.numeric(scale(Foraging.final$SD, center = T, scale = T))
Foraging.final$WC.Z <- as.numeric(scale(Foraging.final$WC, center = T, scale = T))
Foraging.final$WC_prev.Z <- as.numeric(scale(Foraging.final$WC_prev, center = T, scale = T))
Foraging.final$Wind.Exp.Z <- as.numeric(scale(Foraging.final$Wind_Exp, center = T, scale = T))
Foraging.final$StepLength.Z <- as.numeric(scale(Foraging.final$sl_, center = T, scale = T))
Foraging.final$PropFoodSub.Z <- as.numeric(scale(Foraging.final$PropFoodSub, center = T, scale = T))

#Forested Daytime Points
Loafing.final$BA.Z <- as.numeric(scale(Loafing.final$BA, center = T, scale = T))
Loafing.final$DtFE.Z <- as.numeric(scale(Loafing.final$DtFE, center = T, scale = T))
Loafing.final$Ht.Z <- as.numeric(scale(Loafing.final$Ht, center = T, scale = T))
Loafing.final$SW.Z <- as.numeric(scale(Loafing.final$SW, center = T, scale = T))
Loafing.final$PropAg.Z <- as.numeric(scale(Loafing.final$PropAg, center = T, scale = T))
Loafing.final$PropDev.Z <- as.numeric(scale(Loafing.final$PropDev, center = T, scale = T))
Loafing.final$PropSW.Z <- as.numeric(scale(Loafing.final$PropSW, center = T, scale = T))
Loafing.final$SD.Z <- as.numeric(scale(Loafing.final$SD, center = T, scale = T))
Loafing.final$WC.Z <- as.numeric(scale(Loafing.final$WC, center = T, scale = T))
Loafing.final$WC_prev.Z <- as.numeric(scale(Loafing.final$WC_prev, center = T, scale = T))
Loafing.final$Wind.Exp.Z <- as.numeric(scale(Loafing.final$Wind_Exp, center = T, scale = T))
Loafing.final$StepLength.Z <- as.numeric(scale(Loafing.final$sl_, center = T, scale = T))
Loafing.final$PropFoodSub.Z <- as.numeric(scale(Loafing.final$PropFoodSub, center = T, scale = T))

#Roost Points
Roost.final$BA.Z <- as.numeric(scale(Roost.final$BA, center = T, scale = T))
Roost.final$DtFE.Z <- as.numeric(scale(Roost.final$DtFE, center = T, scale = T))
Roost.final$Ht.Z <- as.numeric(scale(Roost.final$Ht, center = T, scale = T))
Roost.final$SW.Z <- as.numeric(scale(Roost.final$SW, center = T, scale = T))
Roost.final$PropAg.Z <- as.numeric(scale(Roost.final$PropAg, center = T, scale = T))
Roost.final$PropDev.Z <- as.numeric(scale(Roost.final$PropDev, center = T, scale = T))
Roost.final$PropSW.Z <- as.numeric(scale(Roost.final$PropSW, center = T, scale = T))
Roost.final$SD.Z <- as.numeric(scale(Roost.final$SD, center = T, scale = T))
Roost.final$WC.Z <- as.numeric(scale(Roost.final$WC, center = T, scale = T))
Roost.final$WC_prev.Z <- as.numeric(scale(Roost.final$WC_prev, center = T, scale = T))
Roost.final$Wind.Exp.Z <- as.numeric(scale(Roost.final$Wind_Exp, center = T, scale = T))
Roost.final$StepLength.Z <- as.numeric(scale(Roost.final$sl_, center = T, scale = T))
Roost.final$PropFoodSub.Z <- as.numeric(scale(Roost.final$PropFoodSub, center = T, scale = T))

#Create Stratum ID (may already have one but this is numeric)
Roost.final <- Roost.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))

Loafing.final <- Loafing.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))

Foraging.final <- Foraging.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))

AllPoints.final <- AllPoints.final %>%
  mutate(StratID = as.numeric(paste(ID, "0000", step_id_, sep ="")))


write.csv(Roost.final, "Roostfinal_HMM.csv", row.names = F)
write.csv(Loafing.final, "Stationaryfinal_HMM.csv", row.names = F)
write.csv(Foraging.final, "Mobilefinal_HMM.csv", row.names = F)
write.csv(AllPoints.final, "AllPoints.csv", row.names = F)