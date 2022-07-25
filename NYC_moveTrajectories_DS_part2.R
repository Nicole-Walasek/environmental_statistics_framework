# simulate relatic move trajectories based on NYC dataset 

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


# divide data according to regions of interest --------------------------

load("data/nyc-spatialNew.Rdata")
head(nyc_spatial)
glimpse(nyc_spatial)


dataUnfiltered <- nyc_spatial %>% 
  as_tibble() %>% 
  select(date_group, date, crime_type,boro,starts_with("region")) %>% 
  gather(region, distance, -date_group, -date, -crime_type,-boro) %>% 
  mutate(
    region = case_when(region=="region3" & distance<5000 ~ "Morrisania (Bronx, Harlem)",
                       region=="region4" & distance<5000 ~ "Tottenville",
                       region=="region5" & distance<5000 ~ "Upper East Side",
                       region=="region7" & distance<5000 ~ "Brownsville",
                       region=="region8" & distance<5000 ~ "Ozone Park",
                       T ~ NA_character_)
  ) %>% 
  group_by(date_group, date, region,crime_type, boro) %>% 
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



# add correction for population density per boro
nyc_pop_2010 <- tibble(
  boro        = c("B","K","M","Q","S"),
  population  = c(1385108, 2504700,1585873,2230722,468730)
)



dataUnfilteredCorrected <- dataUnfiltered %>%
  inner_join(nyc_pop_2010, "boro") %>%
  group_by(date_group,date,region)%>%
  dplyr::summarize(assaultAdjusted = (sum(assault)/sum(population))*100000,
            theftAdjusted = (sum(theft)/sum(population))*100000,
            assault = sum(assault),
            theft = sum(theft))


dataCorrected <- dataUnfiltered %>% 
  filter(date_group == "Monthly") %>% 
  inner_join(nyc_pop_2010, "boro") %>% 
  group_by(date, region) %>% 
  dplyr::summarize(assaultAdjusted = (sum(assault)/sum(population))*100000,
            theftAdjusted = (sum(theft)/sum(population))*100000,
            assault = sum(assault),
            theft = sum(theft))

#plot data per region 
dataUnfilteredCorrected %>% ggplot(aes(x = date, y = assault)) +
  geom_line() +
  scale_x_date(minor_breaks = "3 months")+
  facet_grid(date_group~region, scales = "free_y")+ theme(legend.position = "none") + theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = c(0.06, 0.75), text = element_text(size = 20))


#plot data per region with correction
dataUnfilteredCorrected %>% ggplot(aes(x = date, y = assaultAdjusted)) +
  geom_line() +
  scale_x_date(minor_breaks = "3 months")+
  facet_grid(date_group~region, scales = "free_y")+
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

#ggsave("ROI2.png", units="in", width=20, height=15, dpi=600)


# for publication: make new plot ------------------------------------------

# plot only the monthly data 

regionNames <- unique(dataCorrected$region)
colorNames <- brewer.pal(n = length(regionNames), name = 'RdYlBu')
colorNames[3] <- '#ac7339'
names(colorNames) = regionNames


pTemporal2 <-
  ggplot(data = dataCorrected, aes(x = date, y = assaultAdjusted, group = region)) +
  geom_line(aes(color = region), size = 0.8) +
  annotate(
    "rect",
    xmin = as.Date("2013-01-01"),
    xmax = as.Date("2014-01-01"),
    ymin = 0,
    ymax = 43,
    alpha = .2
  ) +scale_color_manual(values=colorNames)+
  scale_x_date(expand = c(0.01,0.01) ,limits = as.Date(c('2006-01-01','2021-01-01')),breaks = seq.Date(as.Date("2006-01-01"), as.Date("2021-01-01"), by = "1 year"),
    labels = c(2006:2021)
  ) +
  theme(legend.position = "none") + theme_bw() + labs(y = 'Number of assaults (per 100,000)\n', x = 'Year') +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = 'None',
    text = element_text(size = 35),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )) + theme(strip.background = element_blank()) + scale_y_continuous(expand = c(0.01,0.01)) 

print(pTemporal2)
plot_filtered <- pTemporal2

make_circles <- function(centers, radius, nPoints = 100){
  # centers: the data frame of centers with ID
  # radius: radius measured in kilometer
  #
  meanLat <- mean(centers$lattitudes)
  # length per longitude changes with lattitude, so need correction
  radiusLon <- radius /111 / cos(meanLat/57.3) 
  radiusLat <- radius / 111
  circleDF <- data.frame(ID = rep(centers$region, each = nPoints))
  angle <- seq(0,2*pi,length.out = nPoints)
  
  circleDF$lon <- unlist(lapply(centers$longitudes, function(x) x + radiusLon * cos(angle)))
  circleDF$lat <- unlist(lapply(centers$lattitudes, function(x) x + radiusLat * sin(angle)))
  return(circleDF)
}

# regions of interest 
longitudes <- c(-73.90653,-74.24931,-73.96412,-73.91183, -73.84375)
lattitudes <- c(40.82927,40.51122, 40.77370, 40.66452, 40.67677)
region <- c("M", "T", "UES", "B", "OP")

ROI_df <- data.frame(longitudes, lattitudes, region)
ROI_df$region <- as.factor(ROI_df$region)

## get map and plot station locations 
register_google(key = "INSERT YOUR KEY")

newyork.map <-  get_map(location = c(lon = mean(ROI_df$longitudes), 
                                     lat = mean(ROI_df$lattitudes)), 
                        maptype = "roadmap", color = 'bw',zoom = 10, scale =4)#


# change opacity of basemap
mapatt <- attributes(newyork.map)
map_transparent <- matrix(adjustcolor(newyork.map, alpha.f = 0.95), nrow = nrow(newyork.map))
attributes(map_transparent) <- mapatt

# here is the data frame for all circles
myCircles <- make_circles(ROI_df, 5)
regionNames <- levels(myCircles$ID)
names(colorNames) = regionNames


p <-
  ggmap(map_transparent, inherit.aes = FALSE) + geom_polygon(
    inherit.aes = FALSE,
    data = myCircles,
    aes(lon, lat, group = ID, fill = ID),
    alpha = 0.5
  ) +
  geom_point(inherit.aes = FALSE,
    data = ROI_df,
    aes(x = longitudes, y = lattitudes),
    shape = region,
    size = 12,
    color = 'black'
  ) +
  labs(y = 'Longitude\n', x = '\nLatitude') + xlim(-74.35, -73.7) + ylim(40.45, 40.9) + scale_fill_manual(values = colorNames) +
  theme(legend.position = "None", text = element_text(size = 35))
print(p) 



# add crime data for one month 

dfPlotting <- nyc_spatial %>% 
  filter(date >= as.Date("2013-01-01") & date < as.Date("2014-01-01")) %>%
  filter(crime_type == 'assault')

crimes <- data.frame(st_coordinates(dfPlotting$geometry))

pNew <- p + geom_point(data = crimes, aes(x= X, y = Y),shape =20, size = 1,alpha = 0.02, color = '#196F3D')+
  theme(text = element_text(size = 35)) 

print(pNew)

# now put them both together with cowplot 

plot_grid(pNew, plot_filtered, scale = 0.9,labels = "AUTO",label_size = 35, 
          ncol =1, rel_widths = c(1, 3), rel_heights = c(2.5,1),
          label_x = 0, label_y = 0.95,
          hjust = -0.8, vjust = -1)


ggsave("ROINew3Test.png", units="in", width=21, height=22, dpi=600)


# additional data cleaning ------------------------------------------------
data<- dataCorrected
data %>%
  group_by(region) %>%
  dplyr::count(region)


data_by_district <-
  data %>% group_by(region)

data_by_district_split <- group_split(data_by_district)
missingDates = !(unique(data_by_district_split[[5]]$date) %in% unique(data_by_district_split[[4]]$date))
missingDates <- unique(data_by_district_split[[5]]$date)[missingDates]
#
# # to only work with complete cases run the following
# dataNew <- data[!(data$date %in% missingDates),]
# data <- dataNew
#

# to add values to the data frame for the missing measures run the following
missingData <- data.frame(missingDates,rep("Tottenville",length(missingDates)), rep(0, length(missingDates)), rep(0,length(missingDates)),rep(0, length(missingDates)), rep(0,length(missingDates)) )
names(missingData) <- names(data)
missingData
data <- bind_rows(data, missingData)
data <- data[order(data$date),]


# save regionwise data 
head(data)
save(data, file = "data/crimeDataRegion_monthly.Rdata")


# simulate participants  --------------------------------------------------

glimpse(data)
head(data)
tail(data)
# what is the goal? 
# simulate a dataset of participants who have moved between these different locations

# treat n unique locations as states in a Markov Model and define transition probabilities between states
# base case: probability of x for staying in the same state and probability of (1-x)/(n-1) to
# any other state
# randomly select starting point
# simulate uniqueTimes many datapoints; this is tracjteory for one participany

uniqueTimes <- unique(data$date)
nDates <- length(uniqueTimes) #180 months, 15 years
uniqueRegions <- unique(data$region)
nRegions <- length(uniqueRegions) #6 regions 

noMoveProb <- 0.98
moveProb <- (1-noMoveProb)/(nRegions-1)
numParticipants <- 500

pMat <- matrix(0,nRegions, nRegions, byrow = TRUE)
diag(pMat) <- rep(noMoveProb,nRegions)
pMat[lower.tri(pMat)] <- moveProb
pMat[upper.tri(pMat)] <- moveProb


mc <- new("markovchain", states = uniqueRegions, transitionMatrix = pMat)
curves <- matrix(nrow = 5, ncol = 5, 0.1)
plot(mc, package= "diagram", box.size = 0.05, curve = curves,
     name = c('B', 'M', 'OP', 'T', 'UES'), lwd = 1,arr.lcol = 'grey',arr.col = 'black', 
     shadow.size = 0, dtext = 0.7, self.lwd = 1, self.cex = 0.9,arr.pos = 0.6, cex.txt = 0.7,
self.shiftx = c(0.1,0.1,0,-0.1,-0.1,0), 
self.shifty = c(0,0,-0.1,0,0,0.1))

show(mc)

regionsData <- c()
# this is data for one participant; apply this n times as many pp you want and store in long format? 
for(idx in 1:numParticipants){
  regionsData <- c(regionsData,rmarkovchain(nDates,mc, t0 = sample(uniqueRegions,1),  include.t0 = FALSE))
}

# make the pp and time vector in long format 
ppVec <- as.factor(rep(1:numParticipants, each = nDates))
timeVec <- rep(uniqueTimes, numParticipants)

#select data from dataframe based on specific combinations of columns values 

resIdx <- numeric(length(regionsData))
for(idx in 1:length(regionsData)){
  resIdx[idx] = which(data$date == timeVec[idx] & data$region == regionsData[idx])  
}


crimeDataSim <- data[resIdx,]
crimeDataSim$ppID <- ppVec
crimeDataSim

#next save for analysis 
save(crimeDataSim, file = "data/crimeDataSim_098_monthly_500ppCorrectedNew.Rdata")
