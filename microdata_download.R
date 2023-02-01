
# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


# association microdata at the district level 
# https://walker-data.com/tidycensus/articles/pums-data.html
library(tidyverse)
library(tidycensus)

# provide API key for entire session
census_api_key("d4ede4df89da2f188e0a25513aecbc97aaf41474", install = T)
#d4ede4df89da2f188e0a25513aecbc97aaf41474


#getting the data 2016-2020
# analysis with these should be quite straightforqard

vt_pums <- get_pums(
  variables = c("PUMA", "SEX", "AGEP", "SCHL", "SCHG","POVPIP", "WAGP", "HINCP", "ESR"),
  state = "NY",
  survey = "acs5",
  year = 2020,
  rep_weights = "person"
)

save(vt_pums, file = "data/nyc-pumsCommunityDistricts2016_2020.Rdata")


#getting the data 2012-2016
# analysis with these should be quite straightforqard

vt_pums <- get_pums(
  variables = c("PUMA", "SEX", "AGEP", "SCHL", "SCHG","POVPIP", "WAGP", "HINCP", "ESR"),
  state = "NY",
  survey = "acs5",
  year = 2016,
  rep_weights = "person"
)

save(vt_pums, file = "data/nyc-pumsCommunityDistricts2012_2016.Rdata")



#getting the data 2006-2010
# 
# in combination with the 2000 census data?
# shapefiles for 2000 : https://data.cityofnewyork.us/City-Government/2000-Census-Tracts/ysjj-vb9j
# show this only in the appendix or commend on at least on the limitations 


vt_pums <- get_pums(
  variables = c("PUMA", "SEX", "AGEP", "SCHL", "SCHG","POVPIP", "WAGP", "HINCP", "ESR"),
  state = "NY",
  survey = "acs5",
  year = 2010,
  rep_weights = "person"
)

save(vt_pums, file = "data/nyc-pumsCommunityDistricts2006_2010.Rdata")



# get_pums() also always returns SERIALNO, SPORDER, WGTP, PWGTP, and ST. 
# SERIALNO and SPORDER are the variables that uniquely identify observations, 
# WGTP and PWGTP are the housing-unit and person weights, and ST is the state code.


# "SEX" : 1 male , 2 female 
# "AGEP": age i years
# "SCHL": educational attainment 
# "SCHG" : Grade level attending 15,16 covers university and graduate level 
# "POVPIP": income-to-poverty-ratio-recode 
# "WAGP": wages or salary income past  12 months 
# "HINCP" household income past 12 months 
# "ESR": employment status recode 

