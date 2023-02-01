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

# I am using shapefiles from 2010 because the PUMS microdata of the census
# bureau are also based on community districts from 2010 


nycCommunityDistricts <- st_read("data/nycd_10cav/nycd_10cav/nycd.shp")


nycDistricts <- st_read("data/nypuma2010_22c/nypuma2010_22c/nypuma2010.shp")


nycBoroughs <- st_read("data/nybb_10cav/nybb_10cav/nybb.shp")

#plot the shape file

regionNames <- sort(unique(nycBoroughs$BoroName))
colorNames <- brewer.pal(n = length(regionNames), name = 'BrBG') #RdYlBu
colorNames[3] <- '#a9a9a9'
names(colorNames) = regionNames

# 
# Brooklyn     Manhattan        Queens Staten Island     The Bronx 
# "#A6611A"     "#DFC27D"     "#a9a9a9"     "#80CDC1"     "#018571" 

# add labels for plotting
nycDistricts <- nycDistricts %>% mutate(PUMALabels = substr(PUMA, nchar(PUMA)-2+1, nchar(PUMA)))
nycDistricts$PUMALabels <- round(as.numeric(nycDistricts$PUMALabels))

p_boroughs <- ggplot() + 
  geom_sf(data = nycBoroughs, size = 0,aes(fill = BoroName)) + 
  geom_sf(data = nycDistricts, linewidth = 0.8, color = "white", alpha = 0) + 
  #geom_sf(data = nycCommunityDistricts, linewidth = 0.6, color = "black", alpha = 0) +
  scale_fill_manual(values = colorNames) +
  geom_sf_text(data = nycDistricts,aes(label = str_wrap(PUMALabels, 1)), size = 6)+ 
  ggtitle("NYC boroughs and PUMAs") + coord_sf()+ theme_bw()+
  labs(y = 'Longitude\n', x = '\nLatitude', fill = 'Borough')+
  theme(
    plot.title = element_text(face = "bold", size = 35) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.border = element_rect(colour = 'grey70', size = 1),
    #legend.position = 'None',
    text = element_text(size = 35)
  ) + theme(strip.background = element_blank(), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) + #scale_y_continuous(expand = c(0.01, 0.01))
  scale_x_continuous(expand=c(0,0))

p_boroughs


ggsave("boroughAnalysis/boroughsAndDistrictsMap2.png", units="in", width=15, height=15, dpi=600)


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
nycDistrictsTrans <- st_transform(nycDistrictsReduced, 4326)


nyc_spatial <- nyc_data1 %>%
  dplyr::select(crime_type,starts_with("per_"),longitude,latitude, boro) %>% 
  gather(date_group, date, starts_with("per_"),-boro,-crime_type, -longitude, -latitude) %>%
  drop_na %>%
  mutate(
    longitude_c = as.character(longitude),
    latitude_c = as.character(latitude)
  )%>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  mutate(
    intersection = as.integer(st_intersects(geometry, nycDistrictsTrans))
    , district = if_else(is.na(intersection), '', str_wrap(nycDistrictsTrans$PUMA[intersection]))
  )
  

# save data ---------------------------------------------------------------

save(nyc_spatial, file = "data/nyc-spatialCommunityDistricts.Rdata")

