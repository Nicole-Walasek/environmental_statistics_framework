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
library(raster)
library(RColorBrewer)

# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


# explore shape data ------------------------------------------------------

# shapefiles from 2000
nycCencusTracks <- st_read("data/2000 Census Tracts/geo_export_c605b1b8-6f2a-44ae-8099-4357c3b8f555.shp")


pumaNAMES <- c(3701:3710, 3801:3810,3901:3903,4001:4018,4101:4114)
#districtNames <- c(208,212,210,211,263,207,205,204,209,221,112,109,110,111,108,107,154,106,103,121,503,502,501,301,304,303,302,306,308,316,305,318,317,309,307,310,312,314,315,311,313,401,403,407,411,413,408,404,406,402,405,409,412,410,414)


# exclude irrelavnt PUMAs 
nycCensusTracksReduced <- nycCencusTracks %>%
  filter((puma %in% pumaNAMES))


nycCensusTracksReducedBORO <- nycCensusTracksReduced %>%
  group_by(boro_name)%>%
  dplyr::summarise(across(geometry, ~ sf::st_union(.)), .groups = "keep")

plot(nycCensusTracksReducedBORO)


nycCensusTracksReducedPUMA <- nycCensusTracksReduced %>%
  group_by(puma)%>%
  dplyr::summarise(across(geometry, ~ sf::st_union(.)), .groups = "keep")

plot(nycCensusTracksReducedPUMA)



# data --------------------------------------------------------------------
load("data/nyc-crime-data.Rdata")
head(nyc_crime_data)
glimpse(nyc_crime_data)
# identify regions of interest (NYC) --------------------------------------

# SORUCES
# https://www.icphusa.org/wp-content/uploads/2016/04/Poverty.pdf
# https://www.cubesmart.com/blog/city-guides/nyc/5-of-the-best-neighborhoods-to-live-in-manhattan/
# https://homevestorsfranchise.com/blog/northeast/2020/03/the-worst-neighborhoods-in-queens-that-you-should-be-buying-investment-houses-in-now/


# prepare data for framework ----------------------------------------------

nyc_data1 <- nyc_crime_data %>% 
  filter(str_detect(PD_DESC,"ASSAULT|SEXUAL ABUSE|RAPE")) %>%
  rename_all(tolower) 
  

nyc_data1 %>%
  rename_all(tolower)%>%
  group_by(pd_desc)%>%
  count()


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


# need to make sure the NYC shape data and the crime data have the same
# coordinate system 
nycCensusTracksReducedTransPUMA <- st_transform(nycCensusTracksReducedPUMA, 4326)


nyc_spatial <- nyc_data1 %>%
  dplyr::select(crime_type,starts_with("per_"),longitude,latitude, boro) %>% 
  gather(date_group, date, starts_with("per_"),-boro,-crime_type, -longitude, -latitude) %>%
  drop_na %>%
  mutate(
    longitude_c = as.character(longitude),
    latitude_c = as.character(latitude)
  )%>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  mutate(
    intersection = as.integer(st_intersects(geometry, nycCensusTracksReducedTransPUMA))
    , district = if_else(is.na(intersection), '', str_wrap(nycCensusTracksReducedTransPUMA$puma[intersection]))
  )


# change the district labels to match those of the other batches
# nyc_spatial$district <- as.factor(nyc_spatial$district)
# names(nyc_spatial$district)
# 
# 
# pumaTOcd <- setNames(as.numeric(districtNames), as.numeric(pumaNAMES))
# nyc_spatial$district <- unname(pumaTOcd[as.character(nyc_spatial$district)])


# save data ---------------------------------------------------------------
save(nyc_spatial, file = "data/nyc-spatialCommunityDistricts2000.Rdata")

