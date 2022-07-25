# group data into regions of interest  

# Packages ----------------------------------------------------------------
library(plyr)
library(tidyverse)
library(lubridate)
library(reshape)
library(pdp)
library(changepoint)
library(changepoint.np)
library(broom)
library(sf)

# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


# data --------------------------------------------------------------------
load("data/nyc-crime-data.Rdata")
head(nyc_crime_data)
glimpse(nyc_crime_data)
# identify regions of interest (NYC) --------------------------------------

# SORUCES
# https://www.icphusa.org/wp-content/uploads/2016/04/Poverty.pdf
# https://www.cubesmart.com/blog/city-guides/nyc/5-of-the-best-neighborhoods-to-live-in-manhattan/
# https://homevestorsfranchise.com/blog/northeast/2020/03/the-worst-neighborhoods-in-queens-that-you-should-be-buying-investment-houses-in-now/


#East Harlem 40.79472 -73.94250
nycRegion1 <- data.frame(lon = -73.94250, lat = 40.79472) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# Bronx (south and central) 40.84665 -73.87859
nycRegion2 <- data.frame(lon = -73.87859, lat = 40.84665) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# Morrisania (district 17 ) 40.82927 -73.90653 is part of the bronx; poor neighborhood 40% of residents below poverty level 
nycRegion3 <- data.frame(lon = -73.90653, lat = 40.82927) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

#Tottenville (district 51) 40.51122 -74.24931; wealthy neighborhood 6% of residents below poverty level
nycRegion4 <- data.frame(lon = -74.24931, lat = 40.51122) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# Upper East Side 40.77370 -73.96412
nycRegion5 <- data.frame(lon = -73.96412, lat = 40.77370) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# Washington heights 40.84020 -73.94022
nycRegion6 <- data.frame(lon = -73.94022, lat = 40.84020) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

# brownsville (Brooklyn) 40.66452 -73.91183
nycRegion7 <- data.frame(lon = -73.91183, lat = 40.66452) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

#Ozone park (Queens) 40.67677 -73.84375
nycRegion8 <- data.frame(lon = -73.84375, lat = 40.67677) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)


# prepare data for framework ----------------------------------------------

nyc_data1 <- nyc_crime_data %>% 
  filter(str_detect(PD_DESC,"ASSAULT|SEXUAL ABUSE|RAPE")) %>%
  rename_all(tolower) 
  

nyc_data1 %>%
  rename_all(tolower)%>%
  group_by(pd_desc)%>%
  count()

nyc_data1 %>% count()

nyc_data1 <- nyc_crime_data %>% 
  filter(str_detect(PD_DESC,"ASSAULT|SEXUAL ABUSE|RAPE|BURGLAR|THEFT")) %>%
  rename_all(tolower) %>% 
  mutate(
    arrest_date  = mdy(arrest_date), 
    crime_type   = case_when(str_detect(pd_desc,"ASSAULT|SEXUAL ABUSE|RAPE") ~ "assault",
                             str_detect(pd_desc,"BURGLAR|THEFT") ~ "theft"),
    per_day      = arrest_date %>% floor_date("day"),
    per_week     = arrest_date %>% floor_date("week"),
    per_month    = arrest_date %>% floor_date("month"),
    per_biannual = arrest_date %>% floor_date("halfyear"),
    per_annual   = arrest_date %>% floor_date("year"),
    boro         = arrest_boro
  ) 


head(nyc_data1)
nrow(nyc_data1)
max(nyc_data1$arrest_date)
min(nyc_data1$arrest_date)

nyc_spatial <- nyc_data1 %>%
  select(crime_type,starts_with("per_"),longitude,latitude, boro) %>% 
  gather(date_group, date, starts_with("per_"),-boro,-crime_type, -longitude, -latitude) %>%
  drop_na %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    mutate(
      region1 = as.numeric(st_distance(nycRegion1, ., by_element = T)),
      region2 = as.numeric(st_distance(nycRegion2, ., by_element = T)),
      region3 = as.numeric(st_distance(nycRegion3, ., by_element = T)),
      region4 = as.numeric(st_distance(nycRegion4, ., by_element = T)),
      region5 = as.numeric(st_distance(nycRegion5, ., by_element = T)),
      region6 = as.numeric(st_distance(nycRegion6, ., by_element = T)),
      region7 = as.numeric(st_distance(nycRegion7, ., by_element = T)),
      region8 = as.numeric(st_distance(nycRegion8, ., by_element = T))
    )

# save data ---------------------------------------------------------------
#save(nyc_spatial, file = "data/nyc-spatialNew.Rdata")

