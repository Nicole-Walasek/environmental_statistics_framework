# apply framework to semi-real data

# this dataset contains simulated move trajectories of 100 participants based
# on real data

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
library(rlang)
library(ggtext)
library(lemon)
library(gtable)
library(grid)


# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

# data --------------------------------------------------------------------
load("data/crimeDataSim_098_monthly_500ppCorrectedNew.Rdata")

# load framework components -----------------------------------------------
source("explore.R")
source("preprocess.R")
source("framework.R")


data <- crimeDataSim

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
data_ppID <- "ppID"
data_x <- "time"
data_y <- "assaultAdjusted"

# assign the variable names that are neceesary for the subsequent functions
names(data)[names(data) == data_ppID] <- "ppID"
names(data)[names(data) == data_x] <- "time"
names(data)[names(data) == data_y] <- "y"
origData <- data # in case decide to use this dataset after all
data



# EXPLORE -----------------------------------------------------------------
yLabel = "Number of assaults (per 100,000)"
freq = 12
freqUnit = "Year"

# get sample size and number of repeated measures
nSample = length(unique(data$ppID))
nObs <- nrow(data) / nSample
numParticipants <- 4 # number of participants to show in the plot
nSample


pOrg <-
  explore.plotTimeSeries(
    origData,
    numParticipants,
    nSample,
    nObs,
    "None",
    "None",
    yLabel,
    'Raw Data',
    yBreak = 5,
    freq,
    freqUnit,
    MEAN = FALSE,
    SD = FALSE,
    SE = FALSE
  )
pOrg

pOrgSeparate <- pOrg + facet_wrap( ~ ppID, ncol = 2)

plot1 <- pOrgSeparate



# time series decomposition

# If the seasonality and residual components are independent of the trend,
# then you have an additive series. If the seasonality and residual components
# are in fact dependent, meaning they fluctuate on trend, then you have a multiplicative series

type <- 'additive'

# The additive decomposition model is most appropriate when the magnitude of the
# trend-cycle and seasonal components remain constant over the course of the
# series. However, when the magnitude of these components varies but still appears
# proportional over time (i.e., it changes by a multiplicative factor),
# the series may be better represented by the multiplicative decomposition model

numParticipantsDecompose <- 1 #how many lines to plot

# this plots the different components of a time series for the
# specified number of participants
# TS_decompose contains one time series object per participant

# across districts
outcome <-
  explore.decompose_TS(data, numParticipantsDecompose, type, freq, freqUnit)

pTS <- outcome$p
dataTS <- outcome$data_TS_decomposed


#data_TS_seasonal <-
#  data.frame((sapply(dataTS["seasonal", ], c)))
# Do we want to do something with that? Do we want to report the range between
# seasonal peaks and valleys?


pTS
plot2 <- pTS


# PART 3
# autocorrelation and partial autocorrelation at different lags

lagMax = 20
diffDegree = 1
removeSeason = FALSE
logTransform = FALSE
type = "additive"
freq = 12

# distribution of dickey_fuller values in raw data and after attempting stationarity
plot3 <-
  explore.visualize_stationarity(
    data = data,
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    type = type,
    freq = freq
  )



outcome <-
  explore.visualize_autocorr(
    data = data,
    lagMax,
    nSample,
    typeAutoCorr = 'raw',
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    type = type,
    freq = freq,
    regionsFlag = FALSE
  )


plot4 <- outcome


plot5 <-
  outcome <-
  explore.visualize_autocorr(
    data = data,
    lagMax,
    nSample,
    typeAutoCorr = 'stationary',
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    type = type,
    freq = freq
  )


# explore change points
p_cp <-
  explore.changePoints(data = data,
                       minseglen = 12,
                       yLab = 'participant',
                       freqUnit,
                       freq) #x-axis needs some adjustment here
plot6 <- p_cp

# turn all exploratory plots into a single plot
library(ggpubr)

plotList = list(plot1,
                plot4$p_acf,
                plot2 ,
                plot4$p_pacf,
                plot3,
                plot5$p_acf,
                plot6,
                plot5$p_pacf)

plotList = list(
  as.grob(plot1),
  as.grob(plot4$p_acf),
  as.grob(plot2) ,
  as.grob(plot4$p_pacf),
  as.grob(plot3),
  as.grob(plot5$p_acf),
  as.grob(plot6),
  as.grob(plot5$p_pacf)
)


ggarrange(
  plot1,
  plot4$p_acf,
  plot2,
  plot4$p_pacf,
  plot3,
  plot5$p_acf,
  plot6,
  plot5$p_pacf,
  nrow = 4,
  ncol = 2,
  widths = c(1.5, 1),
  heights = c(1, 1.1, 1, 1.1)
)#save with 1800x1600




# PREPROCESS --------------------------------------------------------------
# various options are offered to the user

#1a. seasonally adjust the time series

type <- 'additive'

data_SeasonalaAdj <- preprocess.adjust_season(dataTS, type)
#inspect the adjusted data
outcome <-
  explore.decompose_TS(data_SeasonalaAdj ,
                       numParticipantsDecompose,
                       type,
                       freq,
                       freqUnit)
pTS <- outcome$p
dataTSNew <- outcome$data_TS_decomposed
pTS

# 1b. remove the trend in the data
data_TrendaAdj <- preprocess.adjust_trend(dataTS, type)
#inspect the adjusted data
outcome <-
  explore.decompose_TS(data_TrendaAdj , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTSNew <- outcome$data_TS_decomposed
pTS





#2. difference the time series to achieve stationarity in mean or apply
#   a logarithmic transformation to achieve stationarity in variance

# difference once
data_diffAdj <- preprocess.difference(data, degree = 1, FALSE)
outcome <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq, freqUnit)


pTS <- outcome$p
dataTSNew <- outcome$data_TS_decomposed
pTS


p_cp <-
  explore.changePoints(data = data_diffAdj,
                       minseglen = 4,
                       yLab = 'district',
                       freqUnit,
                       freq) #x-axis needs some adjustment here
p_cp


# logarithmic transformation and no differencing; reduces variance
data_diffAdj <- preprocess.difference(data, degree = 0, TRUE)
head(data_diffAdj)
list[pTS, dataTSNew] <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq, freqUnit)
pTS

# logarithmic transformation and one round of differencing
data_diffAdj <- preprocess.difference(data, degree = 1, TRUE)
head(data_diffAdj)
outcome <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTSNew <- outcome$data_TS_decomposed
pTS

# add seasonal adjustment
data_SeasonalaAdj <- preprocess.adjust_season(dataTSNew, type)
#inspect the adjusted data
outcome <-
  explore.decompose_TS(data_SeasonalaAdj ,
                       numParticipantsDecompose,
                       type,
                       freq,
                       freqUnit)
pTS <- outcome$p
dataTSNew <- outcome$data_TS_decomposed
pTS


# partition data based on different time indices
splitValues <- c(50, 100, 125)
dataSplitted <- preprocess.splitData(data = data, splitValues)
# these individual data sets may be used as input to the framework function
dataEarly <- dataSplitted[[1]]
dataMiddle <- dataSplitted[[2]]
dataLate <- dataSplitted[[3]]



# ANALYZE -----------------------------------------------------------------
# if interested in polynomial effects

data <- data
data$time_c <- scale(data$time, center = TRUE, scale = FALSE)
data$time2 <- data$time_c ^ 2
data$time3 <- data$time_c ^ 3

# could also add an arima model in here if that was desired
#look add jebb paper and add an example with arima models and auto arima models
# maybe provide one example here
meanModel <- function(data) {
  #specify the model of the mean here
  lm(y ~ time, data = data)
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


# parameters for stationarity:
diffDegree = 1
removeSeason = FALSE
logTransform = TRUE
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
    type = type
  )

# to request the "statistics of the environment"
modelResults$Statistics_Mean
names(modelResults$Statistics_Mean)

modelResults$Statistics_Variance
names(modelResults$Statistics_Variance)

modelResults$Statistics_Residuals
names(modelResults$Statistics_Residuals)

modelResults$Statistics_Stationary
names(modelResults$Statistics_Stationary)

modelResults$Statistics_ChangePointsMBIC
names(modelResults$Statistics_ChangePointsMBIC)

modelResults$Statistics_ChangePointsAIC
names(modelResults$Statistics_ChangePointsAIC)

modelResults$Compressed_Statistics # not sure if thats what we want ; add the percentile as a variable
names(modelResults$Compressed_Statistics)

# season
modelResults$dataMaxPeriod


# plotting the models of the mean and variance
outcome <-
  framework.plotModel(modelResults, 5, "mean", "mean", freq, freqUnit, yLabel, segment = FALSE)
pModel <- outcome$pModel
pResiduals <- outcome$pRes

pModel
pModel_Sep <- pModel + facet_wrap( ~ ppID, ncol = 2)
pModel_Sep

pResiduals
pResiduals_Sep <- pResiduals + facet_wrap( ~ ppID, ncol = 2)
pResiduals_Sep

outcome <-
  framework.plotModel(modelResults,
                      5,
                      "variance",
                      "variance",
                      freq,
                      freqUnit,
                      yLabel,
                      segment = FALSE)
pModel <- outcome$pModel
pResiduals <- outcome$pRes

pModel
pModel_Sep <- pModel + facet_wrap( ~ ppID, ncol = 2)
pModel_Sep

pResiduals
pResiduals_Sep <- pResiduals + facet_wrap( ~ ppID, ncol = 2)
pResiduals_Sep


#plot color of noise distribution
# function outputs the same output for a function with and without trend
pColorOfNoise <-
  framework.plotColorOfNoise(modelResults$Statistics_Mean$m_spectralCoef, binwidth = 0.01)
pColorOfNoise


#visualize change points in mean and variance where both are estimated individually
pNew <-
  framework.plotCP(
    modelResults,
    5,
    "Changepoints in mean and variance",
    "MBIC",
    "individual" ,
    freq,
    freqUnit,
    yLabel,
    sampleList = c()
  )
pNew_Sep <- pNew + facet_wrap( ~ ppID, ncol = 2)
pNew_Sep

# where both estimated together
pNew <-
  framework.plotCP(
    modelResults,
    5,
    "Changepoints in mean and variance",
    "MBIC",
    "simulateneous" ,
    freq,
    freqUnit,
    yLabel
  )
pNew_Sep <- pNew + facet_wrap( ~ ppID, ncol = 2)
pNew_Sep



# min and max plots for different statistics
dataModel <- modelResults$data
time <- modelResults$data[[1]]$time

framework.plotMinMax(
  time,
  modelResults$Compressed_Statistics$chronicity_score,
  dataModel,
  'chronicity score'
)
framework.plotMinMax(
  time,
  modelResults$Compressed_Statistics$intensity_score,
  dataModel,
  'intensity score'
)
framework.plotMinMax(time,
                     modelResults$Compressed_Statistics$scaled_mean,
                     dataModel,
                     'scaled mean')

framework.plotMinMax(time,
                     modelResults$Statistics_Mean$m_slope1,
                     dataModel,
                     'mean slope')
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

# spectral coefficient
framework.plotMinMax(time,
                     abs(modelResults$Statistics_Mean$m_spectralCoef),
                     dataModel,
                     'spectral coef')
framework.plotMinMax(time,
                     modelResults$Statistics_Mean$m_spectralCoef,
                     dataModel,
                     'spectral coef')


# mean and sd of raw data
framework.plotMinMax(time,
                     abs(modelResults$Statistics_Mean$m_mean),
                     dataModel,
                     'mean raw data')
framework.plotMinMax(time,
                     modelResults$Statistics_Mean$m_sd,
                     dataModel,
                     'sd raw data')


# mean and sd of var data
framework.plotMinMax(time,
                     abs(modelResults$Statistics_Variance$v_mean),
                     dataModel,
                     'mean variance data')
framework.plotMinMax(time,
                     modelResults$Statistics_Variance$v_sd,
                     dataModel,
                     'sd variance data')



# autocorrelation stationary distribution
framework.plotMinMax(
  time,
  abs(modelResults$Statistics_Stationary$stat_ac1),
  dataModel,
  'autocorrelation (stationary time series)'
)

modelResults$Statistics_Stationary$stat_ac1[317]
modelResults$Statistics_Stationary$stat_ac1[15]


framework.plotMinMax(
  time,
  abs(modelResults$Statistics_Stationary$statVar_ac1),
  dataModel,
  'autocorrelation variance '
)

# autocorrelation raw time series
framework.plotMinMax(time,
                     abs(modelResults$Statistics_Mean$m_ac1),
                     dataModel,
                     'autocorrelation (raw time series)')

modelResults$Statistics_Mean$m_ac1[305]
modelResults$Statistics_Mean$m_ac1[87]


framework.plotMinMax(
  time,
  abs(modelResults$Statistics_Variance$v_ac1),
  dataModel,
  'autocorrelation variance '
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
framework.plotMinMax(
  time,
  modelResults$Statistics_Variance$v_apprEntropy,
  dataModel,
  'entropy variance'
)
framework.plotMinMax(
  time,
  modelResults$Statistics_Residuals$r_apprEntropy,
  dataModel,
  'entropy residuals'
)


# histograms of statistics
histStat <-
  ggplot(data = modelResults$Statistics_Mean, aes(m_apprEntropy)) +
  geom_histogram(binwidth = 0.02) +
  theme(legend.position = "none") + theme_bw() + xlab('Entropy') + labs(title = 'Distribution of entropy values') +
  theme(
    plot.title = element_text(face = "bold", size = 8) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = c(0.06, 0.75)
  )

histStat


# visualize ac of variance
dataVar <- modelResults %>%
  unnest(data) %>%
  select(ppID, yVar)
dataVar$y <- dataVar$yVar


outcome <-
  explore.visualize_autocorr(
    data = dataVar,
    lagMax = 20,
    nSample,
    typeAutoCorr = 'raw',
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    type = type,
    freq = freq
  )

outcome <-
  explore.visualize_autocorr(
    data = dataVar,
    lagMax = 20,
    nSample,
    typeAutoCorr = 'stationary',
    diffDegree = diffDegree,
    removeSeason = removeSeason,
    logTransform = logTransform,
    type = type,
    freq = freq
  )


# correlation among measures ----------------------------------------------

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
    modelResults$Compressed_Statistics,
    .name_repair = 'minimal'
  )


glimpse(dfVis)

corrMatVar <- dfVis %>%
  select(m_sd,
         stat_sd,
         r_sd,
         m_spectralCoef,
         m_apprEntropy,
         stat_apprEntropy,
         stat_ac1,
         m_ac1)

names(corrMatVar) <-
  c(
    "sd (raw)",
    "sd (stationary)",
    "sd (residuals)",
    "spectral exponent",
    "entropy",
    "entropy (stationary)",
    "AC (lag 1, stat)",
    "AC (lag 1, raw)"
  )
M <- cor(corrMatVar)
head(round(M, 2))

col <-
  colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(
  M,
  method = "color",
  col = col(200),
  type = 'upper',
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  diag = F,
  order = "hclust"
)


# compare different measures of intensity / chronicity
corrMatVar <- dfVis %>%
  select(m_mean, chronicity_score, intensity_score, scaled_mean)
M <- cor(corrMatVar)
head(round(M, 2))

col <-
  colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(
  M,
  method = "color",
  col = col(200),
  type = 'upper',
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  diag = F
)


# for publication ---------------------------------------------------------


# same everwhere 265,69,230
# V1       V2 V3
# 69   69 3.052263  7
# 230 230 3.086424  7
# 265 265 3.048399  7



#80,382
# V1       V2 V3
# 80   80 4.749210 17
# 382 382 4.764534 10

sampleList <- c(80,382)


source("framework.R")

unitOffset = 2006
statsPositions = c(27,80,145)

yLabel = "Number of assaults (per 100,000)"
pNew <-
  framework.plotPublicationTest(
    modelResults,
    length(sampleList),
    "title",
    "AIC",
    "individual" ,
    freq,
    freqUnit,
    yLabel,
    sampleList = sampleList, regionFlag = FALSE, 
    unitOffset = unitOffset, 
    stats_positions = statsPositions
  )

pNew


ggsave("SameSD.png", units="in", width=10, height=7, dpi=600)



# plot highlighting statistics of stability, systematic change, unsystematic change 
source("framework.R")
#pp: 105,55
# pp 179,55
yLabel = "Number of assaults (per 100,000)"
sampleList <- c(105,55) #124
pNew <-
  framework.plotPublication2(
    modelResults,
    length(sampleList),
    "title",
    "AIC",
    "individual" ,
    freq,
    freqUnit,
    yLabel,
    sampleList = sampleList,
    type = 'mean', 
    unitOffset = unitOffset
  )

pNew


ggsave("stabilityAndChange_rawNew.png", units="in", width=10, height=7, dpi=600)

yLabel = 'Variance in number of assaults (per 100,000)'
pNew <-
  framework.plotPublication3(
    modelResults,
    length(sampleList),
    "title",
    "AIC",
    "individual" ,
    freq,
    freqUnit,
    yLabel,
    sampleList = sampleList,
    type = 'variance',
    unitOffset = unitOffset
  )

pNew


ggsave("stabilityAndChange_varianceNew.png", units="in", width=10, height=7, dpi=600)





