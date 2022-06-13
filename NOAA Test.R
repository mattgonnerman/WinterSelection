require(sf)
require(dplyr)
require(ggplot2)
require(tidyr)

fullmaine.weather.raw <- read.csv("MaineWeather.csv") %>%
  mutate(Date = as.Date(DATE, "%m/%d/%Y")) %>%
  filter(lubridate::month(Date) %in% c(1:3,12))

maine.weather.sf <- fullmaine.weather.raw %>%
  dplyr::rename(Station = STATION, Location = NAME, lat = LATITUDE, long = LONGITUDE) %>%
  dplyr::select(-ELEVATION, -DATE) %>%
  filter(!is.na(SNWD)) %>%
  st_as_sf(., coords = c("long", "lat"), dim = "XY")


fulllocations.obs.raw <- read.csv("AllUsedPoints.csv") %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")) %>%
  rename(lat = location_lat, long = location_long) %>%
  mutate(State = as.factor(State)) %>%
  st_as_sf(., coords = c("long", "lat"), dim = "XY")



ggplot(maine.weather.sf) +
  geom_sf() +
  geom_sf(data = fulllocations.obs.raw, size = .2, color = "red")

nearest <- st_nearest_points(fulllocations.obs.raw,maine.weather.sf)


### Check Correlation among sites

cor.check <- fullmaine.weather.raw %>%
  dplyr::rename(Station = STATION, Location = NAME, lat = LATITUDE, long = LONGITUDE) %>%
  dplyr::select(-ELEVATION, -DATE) %>%
  filter(!is.na(SNWD)) %>%
  dplyr::select(Station, SNWD, Date) %>%
  pivot_wider(names_from = Station, values_from = SNWD) %>%
  select_if(~ !any(is.na(.)))

cor(cor.check[2:ncol(cor.check)])



download_SNODAS(as.Date("2020-02-02"), out_dir = "C:/Users/mgonn/Downloads/Raw/")

dates.2018 <- as.Date(as.Date("2018-01-01"):as.Date("2018-03-15"), origin = "1970-01-01")
dates.2019 <- as.Date(as.Date("2018-12-01"):as.Date("2019-03-15"), origin = "1970-01-01")
dates.2020 <- as.Date(as.Date("2019-12-01"):as.Date("2020-03-15"), origin = "1970-01-01")

sapply(dates.2018, FUN = download_SNODAS, out_dir = "C:/Users/mgonn/Downloads/Raw/")
sapply(dates.2019, FUN = download_SNODAS, out_dir = "C:/Users/mgonn/Downloads/Raw/")
sapply(dates.2020, FUN = download_SNODAS, out_dir = "C:/Users/mgonn/Downloads/Raw/")


unpack_SNODAS()
rasterize_SNODAS(data_dir = "C:/Users/mgonn/Downloads/Unpack/", 
                 out_dir = "C:/Users/mgonn/Downloads/Rasters/", 
                 format = "GTiff")



### Loop through points
source("./Weather Functions.R")
options(prism.path = "C:/Users/mgonn/Downloads/prismtmp")
lapply(c("sf", "raster", "prism", "rWind", "gdistance"), require, character.only = T)

dates <- unique(sort(as.Date(turkeyData.zm1$Timestamp, tz = "America/New_York")))

for(i in 1:length(dates)){
  date <- dates[i]
  
  #Used Locations
  points <- turkeyData.zm1 %>% filter(as.Date(Timestamp, tz = "America/New_York") == date) %>%
    dplyr::select(ID, long = location_long, lat = location_lat, Timestamp) %>%
    mutate(Date = as.Date(Timestamp, tz = "America/New_York")) %>%
    filter(!is.na(long)) %>%
    st_as_sf(., coords = c("long", "lat"), crs =  4326)
  
  #Snow Raster
  download_SNODAS(date, out_dir = "C:/Users/mgonn/Downloads/Raw/")
  unpack_SNODAS()
  rasterize_SNODAS(data_dir = "C:/Users/mgonn/Downloads/Unpack/", 
                   out_dir = "C:/Users/mgonn/Downloads/Rasters/", 
                   format = "GTiff")
  snow.raster <- raster(list.files("C:/Users/mgonn/Downloads/Rasters/", full.names = T))
  
  #Wind Raster
  wind.df <- wind.dl(yyyy = year(date), mm = month(date), dd = day(date), 10,
                     lon1 = floor(min(st_coordinates(points)[,1])-1), lon2 = ceiling(max(st_coordinates(points)[,1])+1),
                     lat1 = floor(min(st_coordinates(points)[,2])-1), lat2 = ceiling(max(st_coordinates(points)[,2])+1))
  wind.raster <- wind2raster(wind.df)
  
  #Temp Raster
  get_prism_dailys(type = "tmin", minDate = date, maxDate = date, keepZip = F)
  prismdata.list <- ls_prism_data()
  prismlist.pos <- which(grepl(gsub("-", x = date, ""), ls_prism_data()[,1]))
  temp.raster <- raster(ls_prism_data(absPath=T)[prismlist.pos,2])
    
  #Extract Values
  points$SnowDepth <- raster::extract(snow.raster, points)
  points$Temp <- raster::extract(temp.raster, points)
    
  #Cleanup Files
  
}

x <- wind.dl(2015,2,12,0,0,10,35,45)
