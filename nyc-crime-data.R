# Packages ----------------------------------------------------------------
library(tidyverse)
library(lubridate)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


# Links -------------------------------------------------------------------
nyc_links <- list(
  crime     = "https://data.cityofnewyork.us/api/views/8h9b-rp9u/rows.csv?accessType=DOWNLOAD",
  shootings = "https://data.cityofnewyork.us/api/views/833y-fsy8/rows.csv?accessType=DOWNLOAD"
)

# Read Data ---------------------------------------------------------------
nyc_crime_data <- read_csv(nyc_links$crime)
nyc_shooting_data <- read_csv(nyc_links$shootings)

# Save Data ---------------------------------------------------------------
save(nyc_shooting_data, file = "data/nyc-shooting-data.Rdata")
save(nyc_crime_data,file="data/nyc-crime-data.Rdata")


