# simulate realistic move trajectories based on NYC dataset 

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
library(markovchain)
library(rlang)
library(ggtext)
library(gghighlight)
library(ggmap)
library(cowplot)
library(gtable)
library(grid)
library(RColorBrewer)
library(xlsx)

# helper function 
element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}


# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()


# aggregate assaults across time scales--------------------------

load("data/nyc-spatialCommunityDistricts2000.Rdata")
head(nyc_spatial)
glimpse(nyc_spatial)


dataUnfiltered <- nyc_spatial %>% 
  as_tibble() %>% 
  group_by(date_group, date, crime_type, boro, district) %>% 
  drop_na() %>% 
  dplyr::summarize(crimes = n()) %>% 
  ungroup() %>% 
  mutate(
    date_group = factor(date_group, 
                        levels = c("per_day","per_week","per_month","per_biannual","per_annual"),
                        labels = c("Daily","Weekly","Monthly","Bianual","Annual"))
  )%>%
  filter(date<ymd("2021-01-01"))%>%
  drop_na() %>%   
  spread(crime_type,crimes) %>% 
  mutate(assault = ifelse(is.na(assault),0,assault),
         theft   = ifelse(is.na(theft),0,theft)) 


# add correction for population density per borough and per district 

# first, prepare the population information per borough and district
# information is based on: nyc_decennialcensusdata_2000_2010_change
boro <- c("M","B","K","Q","S")
population <- (c(1585873,1385108,2504700,2230722,468730)+c(1537195,1332650,2465326,2229379,443728))/2

nyc_pop_boro <- tibble(
  boro        = boro,
  populationB  = population
)

nyc_pop_district <- read.csv2("data/district_popPUMA.csv", header=TRUE)

names(nyc_pop_district) <- c("district", "p_2010")
nyc_pop_district$populationD <- (as.numeric(nyc_pop_district$p_2010)) 


#exclusionList <- c(164,226,227,228,355,356,480,481,482,483,484,595)

#nyc_pop_district <- nyc_pop_district %>%
#  filter(!(district %in% exclusionList))

nyc_pop_district <- as_tibble(nyc_pop_district)
nyc_pop_district$district <- as.factor(nyc_pop_district$district)



#####

# 263 : 203 & 206, 221: 201 & 202, 154: 104 & 105, 121: 101 & 102, 

# nyc_pop_district$pums <- nyc_pop_district$district
# nyc_pop_district$pums <-
#   recode(
#     nyc_pop_district$pums,
#     "203" = "263",
#     "206" = "263",
#     "201" = "221",
#     "202" = "221",
#     "104" = "154",
#     "105" = "154",
#     "101" = "121",
#     "102" = "121",
#     .default = levels(nyc_pop_district$pums)
#   ) 
# 
# nyc_pop_pums <- nyc_pop_district %>%
#   group_by(pums) %>%
#   dplyr::summarize(populationD = sum(populationD)) %>%
#   dplyr::rename(district = pums)

#####

# correct for borough population
dataUnfiltered$district <- as.factor(dataUnfiltered$district)

dataUnfilteredPopCorrected <- dataUnfiltered %>%
  inner_join(nyc_pop_boro, "boro")%>%
  mutate(
    assaultAdjustedB = (assault/populationB)*100000,
    theftAdjustedB = (theft/populationB)*100000
  )

# correct for district population
dataUnfilteredPopCorrected <- dataUnfilteredPopCorrected %>%
  inner_join(nyc_pop_district, "district")%>%
  mutate(
    assaultAdjustedD = (assault/populationD)*100000,
    theftAdjustedD = (theft/populationD)*100000
  )

# data across all temporal resolutions per borough
dataBorough_unfiltered <- dataUnfilteredPopCorrected %>%
  group_by(date_group,date,boro)%>%
  dplyr::summarize(assaultAdjustedB = sum(assaultAdjustedB),
                   theftAdjustedB = sum(theftAdjustedB),
                   assaultB = sum(assault),
                   theftB = sum(theft))

# data across all temporal resolutions per district
dataDistrict_unfiltered <- dataUnfilteredPopCorrected %>%
  group_by(date_group,date,district)%>%
  dplyr::summarize(assaultAdjustedD = sum(assaultAdjustedD),
                   theftAdjustedD = sum(theftAdjustedD),
                   assaultD = sum(assault),
                   theftD = sum(theft))

# monthly data per borough
dataBorough_monthly <- dataUnfilteredPopCorrected %>% 
  filter(date_group == "Monthly") %>% 
  group_by(date, boro) %>% 
  dplyr::summarize(assaultAdjustedB = sum(assaultAdjustedB),
            theftAdjustedB = sum(theftAdjustedB),
            assaultB = sum(assault),
            theftB = sum(theft))

# monthly data per distric
dataDistrict_monthly <- dataUnfilteredPopCorrected %>% 
  filter(date_group == "Monthly") %>% 
  group_by(date, district) %>% 
  dplyr::summarize(assaultAdjustedD = sum(assaultAdjustedD),
                   theftAdjustedD = sum(theftAdjustedD),
                   assaultD = sum(assault),
                   theftD = sum(theft))

# monthly data per borough and district (need this for plotting)
dataBoroughDistrict_monthly <- dataUnfilteredPopCorrected %>% 
  filter(date_group == "Monthly") %>% 
  group_by(date, district,boro) %>% 
  dplyr::summarize(assaultAdjustedD = sum(assaultAdjustedD),
                   theftAdjustedD = sum(theftAdjustedD),
                   assaultD = sum(assault),
                   theftD = sum(theft))


#plot data per region with correction
dataBorough_unfiltered %>% ggplot(aes(x = date, y = assaultAdjustedB)) +
  geom_line() +
  scale_x_date(minor_breaks = "3 months")+
  facet_grid(date_group~boro, scales = "free_y")+
  theme(legend.position = "none") + theme_bw() + labs(y = 'Number of assaults (per 100,000)') +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = c(0.06, 0.75), 
    text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  theme(
    strip.background = element_blank(),
    strip.text = element_textbox_highlight(
      size = 16, face = "bold",
      fill = "white", box.color = "grey70", color = "black",
      halign = .5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
      padding = margin(5, 0, 3, 0), margin = margin(0, 1, 3, 1)
    ))

ggsave("boroughAnalysis//boroughs_acrossTime2000.png", units="in", width=20, height=15, dpi=600)

# add PUMS labels ---------------------------------------------------------
# 
# 
# # 263 : 203 & 206, 221: 201 & 202, 154: 104 & 105, 121: 101 & 102, 
# dataDistrict_monthly$pums <-
#   as.factor(dataDistrict_monthly$district)
# 
# 
# 
# dataBoroughDistrict_monthly$district <-
#   as.factor(dataBoroughDistrict_monthly$district)
# 
# dataBoroughDistrict_monthly$pums <- dataBoroughDistrict_monthly$district

# save data ---------------------------------------------------------------
save(dataBorough_monthly, file = "data/monthly_borough2000.Rdata")
save(dataDistrict_monthly, file = "data/monthly_districts2000.Rdata")
save(dataBoroughDistrict_monthly, file = "data/monthly_borough_districts2000.Rdata")




