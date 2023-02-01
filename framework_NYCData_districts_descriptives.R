

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
library(pracma) # appr. and sample entropy
library(corrplot)
library(tseries)
library(rlang)
library(ggtext)
library(lemon)
library(gtable)
library(grid)
library(formattable)
library(kableExtra)
library(sjPlot)
library(stargazer)
library(gridGraphics)
library(gridExtra)
library(Hmisc)
library(brms)
library(betareg)
library(tidybayes)
library(emmeans)
library(parallel)
library(RColorBrewer)
library(brms)
library(emmeans)
library(GGally)
# label for saving the results --------------------------------------------

resLabel <- "2016_2020"

# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

# data --------------------------------------------------------------------

if (grepl("2006", resLabel, fixed = TRUE)) {
  load("data/monthly_districts2000.Rdata")
} else{
  load("data/monthly_districts.Rdata")
}

head(dataDistrict_monthly)
data <- dataDistrict_monthly

# load framework components -----------------------------------------------
source("explore.R")
source("preprocess.R")
source("framework.R")




# in case you want to filter the data by year ----------------------
year1 <- str_split(resLabel, "_")[[1]][1]
year2 <- str_split(resLabel, "_")[[1]][2]
year2 <- as.character((as.numeric(year2) + 1))

data <- data %>%
  filter(date >= as.Date(sprintf("%s-01-01", year1)) &
           date < as.Date(sprintf("%s-01-01", year2)))


# prepare data for framework ----------------------------------------------

# create a time variable based on the dates
time <- as.factor(data$date)
data$time <-
  plyr::mapvalues(
    time,
    from = as.factor(unique(data$date)),
    to = c(1:length(unique(data$date))),
    warn_missing = F
  )
data$time <- as.numeric(data$time)

# user needs to specify their variable names
data_x <- "time"
data_y <- "assaultAdjustedD"

# assign the variable names that are neceesary for the subsequent functions
names(data)[names(data) == data_x] <- "time"
names(data)[names(data) == data_y] <- "y"

#relabel partipant ID factor levels for subsequent analysis
data$region <- as.factor(data$district)
data$ppID <- as.factor(data$region)
levels(data$ppID) <- c(1:length(unique(data$ppID)))
region_names <- setNames(levels(data$region), levels(data$ppID))

origData <- data # in case decide to use this dataset after all
data

yLabel = "Number of assaults (per 100,000)"
freq = 12
freqUnit = "Year"

# get sample size and number of repeated measures
nSample = length(unique(data$ppID))
nObs <- nrow(data) / nSample


# explore stationarity
lagMax = 20
diffDegree = 1
removeSeason = F
logTransform = F
type = "additive"
freq = 12

# distribution of dickey_fuller values in raw data and after attempting stationarity
plotStat <-
  explore.visualize_stationarity(
    data = data,
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    type = type,
    freq = freq
  )

plotStat
# ANALYZE -----------------------------------------------------------------
# if interested in polynomial effects

data$time_c <- scale(data$time, center = TRUE, scale = FALSE)
data$time2 <- data$time_c ^ 2
data$time3 <- data$time_c ^ 3
data$time_s <- scale(data$time)
data$time_s2 <- data$time_s**2
data$y_s <- scale(data$y)

# could also add an arima model in here if that was desired
#look add jebb paper and add an example with arima models and auto arima models
# maybe provide one example here
meanModel <- function(data) {
  #specify the model of the mean here
  lm(y_s ~ time_s, data = data)
}

# if a model of the variance is specified use .var and residuals = FALSE instead
varModel <- function(data) {
  #specify the model of the variance here
  lm(yVar ~ time, data = data)
}


acLAG = 4
freq = 12
freqUnit = "Year"
yLabel = "Number of assaults (per 100,000)"
minseglen = 6 # minimum length of segments for change point detection
percentile = 0.95 # upper bound for the chronicity score
r = 0.35 # distance parameter for estimating the sample entropy

# parameters for stationarity:
diffDegree = 1
removeSeason = FALSE
logTransform = FALSE
type = 'additive'

modelResults <-
  framework.applyModel(
    data,
    acLAG,
    minseglen,
    percentile,
    meanModel,
    varModel,
    residuals = FALSE,
    sqaured = TRUE,
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    freq = freq,
    type = type,
    r
  )



# min and max plots for different statistics
dataModel <- modelResults$data
time <- modelResults$data[[1]]$time

#slope
framework.plotMinMax(time,
                     abs(modelResults$Statistics_Mean$m_slope1),
                     dataModel,
                     'mean slope')


# variance
framework.plotMinMax(time,
                     modelResults$Statistics_Variance$v_slope1,
                     dataModel,
                     'variance slope')


#change points
framework.plotMinMax(
  time,
  modelResults$Statistics_ChangePointsMBIC$meanNumberMBIC,
  dataModel,
  'changepoints mean'
)
framework.plotMinMax(
  time,
  modelResults$Statistics_ChangePointsMBIC$varNumberMBIC,
  dataModel,
  'changepoints variance'
)



# autocorrelation stationary distribution
framework.plotMinMax(
  time,
  abs(modelResults$Statistics_Stationary$stat_ac1),
  dataModel,
  'autocorrelation (stationary time series)'
)


# autocorrelation raw time series
framework.plotMinMax(
  time,
  modelResults$Statistics_Mean$m_ac1,
  dataModel,
  'autocorrelation (raw time series)'
)

framework.plotMinMax(time,
                     modelResults$Statistics_Mean$m_apprEntropy,
                     dataModel,
                     'entropy data')

framework.plotMinMax(
  time,
  modelResults$Statistics_Stationary$stat_apprEntropy,
  dataModel,
  'entropy stationary'
)


# for publication: correlation among measures ----------------------------------------------



# to what extent do measures of unpredictability correlate?
# autocorrelation lag 1 low autocorrelation means high unpred
# entropy high means high
# sd high means high
# change points mean high means high
# change points variance high mean high

dfVis <-
  bind_cols(
    modelResults$Statistics_Mean,
    modelResults$Statistics_Residuals,
    modelResults$Statistics_Stationary,
    modelResults$Statistics_ChangePointsMBIC,
    modelResults$Statistics_ChangePointsAIC,
    modelResults$Compressed_Statistics,
    .name_repair = 'minimal'
  )

glimpse(dfVis)


corrMatVar <- dfVis %>%
  select(
    m_sd,
    m_spectralCoef,
    m_apprEntropy,
    m_ac1,
    stat_ac1,
    meanNumberMBIC,
    varNumberAIC,
    #npNumberMBIC,
    m_mean,
    #m_slope1
  )



names(corrMatVar) <-
  c(
    'SD',
    'C.o.N',
    # random component of decomposed TS
    'Entropy',
    'AC' ,
    'AC (stat.)',
    'CP (mean)',
    'CP (var.)',
    #'CP (both)',
    'Mean'
    #'Slope'
  )


ggpairs(
  corrMatVar,
  upper = 'blank',
  lower = list(continuous = "points"),
  diag =  list(continuous = wrap("barDiag", bins =30))
) + theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'black', size = 1),
    axis.title.x= element_text(size = 18),
    axis.title.y= element_text(size = 18),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13,
                               angle = 45,
                               vjust = 1,
                               hjust = 1)
  )


# ggsave(
#   sprintf("districtAnalysis/descriptives%s.png", resLabel),
#   units = "in",
#   width = 15,
#   height = 15,
#   dpi = 600
# )


# for publication: association with district wise data --------------------

load(file = sprintf("data/nyc-microDatByDistrict%s.Rdata", resLabel))
head(microDat)
#Continue here
# transform this into a nice table for the appendix with all the info


microDat <- microDat %>% 
  mutate(District = (str_remove(District, "0")))

names(microDat) <-
  c(
    'District',
    'Bachelor or above',
    'Bachelor or above (SE)',
    'Bachelor or above (CV)',
    'Bachelor',
    'Prop. bachelor or above (SE)',
    'Prop. bachelor or above (CV)',
    'Unemployed',
    'Unemployed (SE)',
    'Unemployed (CV)',
    'Unempl.',
    'Prop. unemployed (SE)',
    'Prop. unemployed (CV)',
    'Income-poverty below 100',
    'Income-poverty below 100 (SE)',
    'Income-poverty below 100 (CV)',
    'I-pov.<100',
    'Prop. income-poverty below 100 (SE)',
    'Prop. income-poverty below 100 (CV)',
    'Income-poverty above 500',
    'Income-poverty above 500 (SE)',
    'Income-poverty above 500 (CV)',
    'I-pov.>500',
    'Prop. income-poverty above 500 (SE)',
    'Prop. income-poverty above 500 (CV)',
    'Gini'
  )



# compute correlations between unpredictability statistics and community data
communityDat <- microDat %>%
  select(
    'District',
    'Bachelor',
    'Unempl.',
    'I-pov.<100',
    'I-pov.>500'
  )
communityDat$District <- as.factor(communityDat$District)

corrMatVar$District <- region_names

# merge the two data frames
community_unpredictability_data <-
  inner_join(communityDat, corrMatVar, by = 'District')

#maybe turn this into a function?
community_unpredictability_data <-
  column_to_rownames(community_unpredictability_data, var = 'District')

#visually inspect the data

ggplot(gather(community_unpredictability_data), aes(value)) + 
  geom_histogram(bins = 20) + labs(x = 'Value', y = 'Frequency') +
  facet_wrap(~key, scales = 'free_x') +theme_bw()+ theme(text = element_text(size = 20))

ggsave(
  sprintf("districtAnalysis/histograms%s.png", resLabel),
  units = "in",
  width = 12,
  height = 10,
  dpi = 600
)

ggpairs(
  community_unpredictability_data,
  upper = 'blank',
  lower = list(continuous = "points"),
  diag =  'blank')+
theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'black', size = 1),
    axis.title.x= element_text(size = 18),
    axis.title.y= element_text(size = 18),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14,
                               angle = 60,
                               vjust = 1,
                               hjust = 1.1)
  )


ggsave(
  sprintf("districtAnalysis/descriptives%s.png", resLabel),
  units = "in",
  width = 18,
  height = 18,
  dpi = 600
)


