# Packages ----------------------------------------------------------------
library(plyr)
library(tidyverse)
library(lubridate)
library(reshape)
library(pdp)
library(changepoint)
library(changepoint.np)
library(broom)
library(ggplot2)
library(ggthemes)
library(gsubfn)
library(ggpubr)
library(pracma)
library(corrplot)
library(tseries)
library(formattable)
library(kableExtra)
library(knitr)
library(rlang)
library(ggtext)
library(lemon)
library(gtable)
library(grid)



# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

# data --------------------------------------------------------------------
load("data/monthly_borough_districts.Rdata")
data <- dataBoroughDistrict_monthly
originalData <- data


# only inlcude full cases -------------------------------------------------

# only use districts with 180 data points
dataByBorough <- data %>%
  group_by(boro) %>%
  group_split()

# only keep districts with full number of data points
for(idx in 1:length(dataByBorough)){
  data <- dataByBorough[[idx]]
  countDf <- data %>%
    count(district)
  districtIDX <- countDf[countDf$n == 180,]$district
  data <- data %>% filter(district %in% districtIDX)
  dataByBorough[[idx]] <- data
}


# load framework components -----------------------------------------------
source("explore.R")
source("preprocess.R")
source("framework.R")

# prepare data for framework ----------------------------------------------

# create a time variable based on the dates
dataByBoroughNew <- list()

for(idx in 1:length(dataByBorough)){
  data <- dataByBorough[[idx]]
  time <- as.factor(data$date)
  data$time <- plyr::mapvalues(time, from = as.factor(unique(data$date)), to = c(1:length(unique(data$date))), warn_missing = F)
  data$time <- as.numeric(data$time)
  
  # user needs to specify their variable names
  data_x <- "time"
  data_y <- "assaultAdjustedD"
  
  # assign the variable names that are neceesary for the subsequent functions
  names(data)[names(data) == data_x] <- "time"
  names(data)[names(data) == data_y] <- "y"
  
  #relabel partipant ID factor levels for subsequent analysis
  data$region <- droplevels(as.factor(data$district))
  data$ppID <- droplevels(as.factor(data$region))
  
  levels(data$ppID) <- c(1:length(unique(data$ppID)))
  dataByBoroughNew[[idx]] <- data
}



# ANALYZE -----------------------------------------------------------------
# if interested in polynomial effects 


modelResultsByBorough <- list()

for(idx in 1:length(dataByBoroughNew)){
  data <- dataByBoroughNew[[idx]]
  data <- data
  data$time_c <- scale(data$time, center = TRUE, scale = FALSE)
  data$time2 <- data$time_c^2
  data$time3 <- data$time_c^3
  
  # could also add an arima model in here if that was desired 
  #look add jebb paper and add an example with arima models and auto arima models  
  # maybe provide one example here 
  meanModel <- function(data){
    #specify the model of the mean here 
    lm(y ~ time, data = data)
  }
  
  # if a model of the variance is specified use .var and residuals = FALSE instead
  varModel <- function(data){
    #specify the model of the variance here 
    lm(yVar ~ time, data = data)
  }
  
  
  acLAG = 12
  freq = 12
  freqUnit = "Year"
  yLabel = "Number of assaults (per 100,000)"
  minseglen = 6 # minimum length of segments for change point detection 
  percentile = 0.95 # upper bound for the chronicity score 
  
  
  # parameters for stationarity:
  diffDegree = 1
  removeSeason = FALSE
  logTransform = FALSE
  type = 'additive'
  
  
  modelResults <-
    framework.applyModel(data, acLAG,minseglen,percentile, meanModel, varModel, 
                         residuals = FALSE, sqaured = TRUE,diffDegree = diffDegree,
                         removeSeason = removeSeason,
                         logTransform = logTransform, freq = freq,type = type)
  
  modelResultsByBorough[[idx]] <- modelResults
  
}


modelResultsByBorough

# different measures of unpredictability ----------------------------------------------

# autocorrelation lag 1 low autocorrelation means high unpred
# entropy high means high
# sd high means high
# change points mean high means high
# change points variance high mean high


##### Table 1 with unpredictability statistics #####

# load the boroughwise stats
load(file = "data/boroughStats.Rdata")
boroughStats <- corrMatVar
boroughStats <- boroughStats %>%
  dplyr::rename(borough = region)

boroughStats = subset(boroughStats, select = -c(npNumberMBIC,stat_ac1) )


which.min(abs(boroughStats$varNumberAIC))
which.max(abs(boroughStats$varNumberAIC))


res_df <- data.frame()


for(idx in 1:length(modelResultsByBorough)){
  
  modelResults <- modelResultsByBorough[[idx]]
  
  dfVis <- bind_cols(modelResults$Statistics_Mean,modelResults$Statistics_Residuals, modelResults$Statistics_Stationary,
                     modelResults$Statistics_ChangePointsMBIC,modelResults$Statistics_ChangePointsAIC, modelResults$Compressed_Statistics, .name_repair = 'minimal')
  
  corrMatVar <- dfVis %>%
    select(m_sd,
           m_spectralCoef,
           m_apprEntropy,
           m_ac1,
           #stat_ac1,
           meanNumberMBIC,
           varNumberAIC) 
  
  # save the average and sd for each stat
  data_long <- gather(corrMatVar, factor_key=TRUE, key = 'statistic')
  
  df <- data_long%>% group_by(statistic)%>%
    summarise(mean= mean(value), sd= sd(value), max = max(value),min = min(value)) 
  
  
  boro <- dataByBoroughNew[[idx]]$boro[1]
  df$borough <- boro
  res_df <- bind_rows(df,res_df)

  
}

res_df$borough <- as.factor(res_df$borough)

test <- res_df %>%
  group_by(statistic) %>%
  group_split()

which.min(abs(test[[4]]$mean))
which.max(abs(test[[3]]$mean))

# bar plot  ---------------------------------------------------------------

newNames <-
  c(
    'Standard deviation',
    'Color of noise', # random component of decomposed TS
    'Entropy',
    'Autocorrelation (lag 1)' ,
    #'Autocorrelation (lag 1, stat.)',
    'Changepoints (mean)',
    'Changepoints (var.)'
  )
oldNames <- unique(res_df$statistic)
names(newNames) <- oldNames

regions <- c("Brooklyn", "Manhattan", "Queens", "Staten Island", "The Bronx")
names(regions) <- c('K','M','Q','S','B') 
regionLabeler <- as_labeller(regions)

res_df$borough <- unlist(regionLabeler(res_df$borough))

regionNames <- c('Brooklyn', 'Manhattan', 'Queens', 'Staten Island', 'The Bronx')
colorNames <- c("#A6611A","#DFC27D","#a9a9a9","#80CDC1","#018571" )
names(colorNames) = regionNames


boroughStats$borough <- as.factor(boroughStats$borough)
boroughStats_long <- gather(boroughStats, factor_key=TRUE, key = 'statistic', value = 'value',-borough)



plotDF = merge(x = res_df, y = boroughStats_long, by = c('borough','statistic'))

ggplot(plotDF, aes(x=borough, y=mean, fill=borough))+
  geom_bar(stat='identity')+  scale_fill_manual(values = colorNames) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin = value, ymax = value), color = 'black',size =1)+
  facet_wrap(~statistic, scales="free_y", labeller = as_labeller(newNames))+ 
  labs(x = "Borough", y = "Value")+
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white', size = 0),
    legend.position = 'None',
    text = element_text(size = 35),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) 

ggsave("boroughAnalysis/boroughAcrossDistrictsStats2.png", units="in", width=20, height=13, dpi=600)

# V2

sdLow <- plotDF$sd

for (idx in 1:nrow(plotDF)){
  mean <- plotDF$mean[idx]
  sd <- plotDF$sd[idx]
  if ((mean >0) & ((mean-sd) < 0)){
    sdLow[idx] <- 0.01
  }
  
}

plotDF$sdLow <- sdLow
ggplot(plotDF, aes(x = borough, y = mean, fill = borough)) +
  geom_bar(stat = 'identity',
           width = 0.4,
           position = position_nudge(x = -0.21)) +  scale_fill_manual(values = colorNames) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = .2,
                position = position_nudge(x = -0.21)) +
  geom_bar(
    aes(x = borough, y = value, fill = borough),
    inherit.aes = F,
    color = 'black',
    stat = 'identity',
    alpha = .3,
    width = 0.4,
    position = position_nudge(x = 0.21)
  ) +
  #geom_errorbar(aes(ymin = value, ymax = value), color = 'black',size =1)+
  facet_wrap(~ statistic, scales = "free_y", labeller = as_labeller(newNames)) +
  labs(x = "Borough", y = "Value") +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white', size = 0),
    legend.position = 'None',
    text = element_text(size = 40),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )
  )


ggsave("boroughAnalysis/boroughAcrossDistrictsStats_V22.png", units="in", width=20, height=11, dpi=600)




