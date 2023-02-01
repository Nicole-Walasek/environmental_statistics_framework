
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
# label for saving the results --------------------------------------------

resLabel <- "2006_2010"

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


# EXPLORE -----------------------------------------------------------------
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

pOrgSeparate <- pOrg + facet_wrap(~ ppID, ncol = 2)

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
pModel_Sep <- pModel + facet_wrap(~ ppID, ncol = 2)
pModel_Sep

pResiduals
pResiduals_Sep <- pResiduals + facet_wrap(~ ppID, ncol = 2)
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
pModel_Sep <- pModel + facet_wrap(~ ppID, ncol = 2)
pModel_Sep

pResiduals
pResiduals_Sep <- pResiduals + facet_wrap(~ ppID, ncol = 2)
pResiduals_Sep


# plot color of noise distribution
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
pNew_Sep <- pNew + facet_wrap(~ ppID, ncol = 2)
pNew_Sep

# where both estimated together
pNew <-
  framework.plotCP(
    modelResults,
    5,
    "Changepoints in mean and variance",
    "AIC",
    "simulateneous" ,
    freq,
    freqUnit,
    yLabel
  )
pNew_Sep <- pNew + facet_wrap(~ ppID, ncol = 2)
pNew_Sep



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


# for the appendix: statistics of all districts ----------------------------------

#### Table 1 with unpredictability statistics ####
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
  select(m_sd,
         m_spectralCoef,
         m_apprEntropy,
         m_ac1,
         meanNumberMBIC,
         varNumberAIC)

corrMatVar$region <- region_names


names(corrMatVar) <-
  c(
    'Standard deviation',
    'Color of noise',
    'Entropy',
    'Autocorrelation (lag 1)' ,
    'Changepoints (mean)',
    'Changepoints (var.)',
    'District'
  )
corrMatVar <- remove_rownames(corrMatVar)
corrMatVar <- column_to_rownames(corrMatVar, var = 'District')

# continue here and transform this into a table
# then make a new table with different ranks for each variable
corrMatVar <- round(corrMatVar, 2)
corrMatVar <- rownames_to_column(corrMatVar, var = "District")

# use package formattable
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"
customGrey = "#CACFD2"

# helper functions for creating the bars
myNormalize <- function(x, absFlag, highest) {
  if (absFlag) {
    x <- abs(x)
  }
  if (highest) {
    normVal <- (x - min(x)) / (max(x) - min(x))
    
  } else{
    normVal <- (x - min(x)) / (max(x) - min(x))
    normVal <- 1 - normVal
  }
  return(normVal)
}



f.color_bar <-
  function (color = "lightgray", fun = "proportion", ...){
    fun <- match.fun(fun)
    formattable::formatter(
      "span",
      style = function(x)
        style(
          display = "inline-block",
          direction = "ltr",
          `border-radius` = "4px",
          `padding-right` = "2px",
          `background-color` = csscolor(color),
          width = percent(fun(as.numeric(x), ...))
        )
    )
  }


#integrate formattable and kableExtra

corrMatVar %>%
  dplyr::mutate(
    District = formatter("span", style = ~ style(color = "black"))(District),
    `Standard_deviation` = f.color_bar(customGrey, function(x)
      myNormalize(x, TRUE, TRUE))(`Standard deviation`),
    `Color of noise` = f.color_bar(customGrey, function(x)
      myNormalize(x, TRUE, FALSE))(`Color of noise`),
    `Entropy` = f.color_bar(customGrey, function(x)
      myNormalize(x, TRUE, TRUE))(`Entropy`),
    `Autocorrelation (lag 1)` = f.color_bar(customGrey, function(x)
      myNormalize(x, TRUE, FALSE))(`Autocorrelation (lag 1)`),
    `Changepoints (mean)` = f.color_bar(customGrey, function(x)
      myNormalize(x, TRUE, TRUE))(`Changepoints (mean)`),
    `Changepoints (var.)` = f.color_bar(customGrey, function(x)
      myNormalize(x, TRUE, TRUE))(`Changepoints (var.)`)
  ) %>%
  kable("html", escape = F) %>%
  kable_styling("hover", full_width = F, font_size = 20) %>%
  row_spec(0, color = "black", bold = T) %>%
  column_spec(1,
              width = "4cm",
              color = 'black',
              italic = T) %>%
  column_spec(2, width = "3cm", color = 'black') %>%
  column_spec(3, width = "3cm", color = 'black') %>%
  column_spec(4, width = "3cm", color = 'black') %>%
  column_spec(5, width = "3cm", color = 'black') %>%
  column_spec(6, width = "3cm", color = 'black') %>%
  column_spec(7, width = "3cm", color = 'black') %>%
  column_spec(8, width = "3cm", color = 'black') %>%
  #column_spec (1:7,border_left = T, border_right = T)%>%
  kable_classic() %>%
  save_kable(
    sprintf("districtAnalysis/districtPUMSStats%s.png", resLabel),
    density = 600,
    zoom = 2
  )


# for publication: visualizing unpredictability stats per dist---------------------------------------------------------

sortedMean <- order(modelResults$Statistics_Mean$m_mean)

minVal <- sortedMean[1:3]
maxVal <- sortedMean[(length(sortedMean)-2):length(sortedMean)]
sampleList <- c(minVal, maxVal)#c(2, 13, 21, 39, 53)

source("framework.R")

unitOffset = as.numeric(year1)
statsPositions = c(10, 30, 50)#c(27,80,145)

yLabel = "Number of assaults (per 100,000)"
pNew <-
  framework.plotPublicationTest(
    modelResults,
    length(sampleList),
    "title",
    "MBIC", # mean penalty
    "AIC", #var penalty
    "individual" ,
    freq,
    freqUnit,
    yLabel,
    sampleList = sampleList,
    regionFlag = TRUE,
    unitOffset = unitOffset,
    stats_positions = statsPositions
  )

pNew

ggsave(
  sprintf("districtAnalysis/oneDistrictPerBorrough%s.png", resLabel),
  units = "in",
  width = 10,
  height = 12,
  dpi = 600
)


# for publication: correlation among measures ----------------------------------------------


# function for correlations

getLowerTri <- function(cor.mat) {
  cor.mat[upper.tri(cor.mat, diag = T)] <- NA
  return(cor.mat)
}


cors <- function(df) {
  # turn all three matrices (r, n, and P into a data frame)
  M <- Hmisc::rcorr(as.matrix(df), type = "spearman") #spearman
  # return the three data frames in a list return(Mdf)
  Mdf <- map(M, ~ data.frame(.x, check.names = F))
  Mdf$r <- getLowerTri(Mdf$r)
  return(Mdf)
}


formatted_cors <- function(df) {
  cors(df) %>%
    map( ~ rownames_to_column(.x, var = "measure1")) %>%
    map( ~ pivot_longer(.x,-measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(
      sig_p = ifelse(P < .05, T, F),
      p_if_sig = ifelse(P < .05, P, NA),
      r_if_sig = ifelse(P < .05, r, NA)
    )
}



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
    m_mean
  )


corrMatVar <- corrMatVar %>%
  mutate(
    m_spectralCoef = abs(m_spectralCoef),
    m_ac1 = abs(m_ac1),
    stat_ac1 = abs(stat_ac1),
    #m_slope1 = abs(m_slope1)
  )

names(corrMatVar) <-
  c(
    'Standard deviation',
    'Color of noise',
    # random component of decomposed TS
    'Entropy',
    'Autocorrelation (lag 1)' ,
    'Autocorrelation (lag 1, stationary)',
    'Changepoints (mean)',
    'Changepoints (var.)',
    #'Changepoints (both)',
    'Mean'
    
  )


plotDf <- formatted_cors(corrMatVar)
plotDf$measure1 <-
  factor(plotDf$measure1, levels = unique(plotDf$measure1))
plotDf$measure2 <-
  factor(plotDf$measure2, levels = unique(plotDf$measure2))
plotDf <- droplevels(plotDf[!(is.na(plotDf$r)), ])

limits <- c(-1, 1)
plotDf %>%
  ggplot(aes(measure1, measure2)) +
  geom_tile(aes(fill=r_if_sig),color = 'black')  + geom_text(aes(label=round(r,2)),size =5)+
  labs(x = NULL, y = NULL, fill = "Correlation") +
  scale_fill_distiller(
    palette = "BrBG",
    na.value = 'white',
    limits = limits,
    #breaks = breaks,
    guide = guide_colourbar(fill = guide_colourbar(
      barheight = 10, label.position = "left"
    ))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'white', size = 0),
    text = element_text(size = 20),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )
  )


ggsave(
  sprintf("districtAnalysis/correlationMatPUMS%s.png", resLabel),
  units = "in",
  width = 12,
  height = 10,
  dpi = 600
)


# for publication: association with district wise data --------------------

load(file = sprintf("data/nyc-microDatByDistrict%s.Rdata", resLabel))
head(microDat)
range(microDat$giniCoef)
#Continue here
# transform this into a nice table for the appendix with all the info
names(microDat) <-
  c(
    'District',
    'Bachelor or above',
    'Bachelor or above (SE)',
    'Bachelor or above (CV)',
    'Prop. bachelor or above',
    'Prop. bachelor or above (SE)',
    'Prop. bachelor or above (CV)',
    'Unemployed',
    'Unemployed (SE)',
    'Unemployed (CV)',
    'Prop. unemployed',
    'Prop. unemployed (SE)',
    'Prop. unemployed (CV)',
    'Income-poverty below 1',
    'Income-poverty below 1 (SE)',
    'Income-poverty below 1 (CV)',
    'Prop. income-poverty below 1',
    'Prop. income-poverty below 1 (SE)',
    'Prop. income-poverty below 1 (CV)',
    'Income-poverty above 5',
    'Income-poverty above 5 (SE)',
    'Income-poverty above 5 (CV)',
    'Prop. income-poverty above 5',
    'Prop. income-poverty above 5 (SE)',
    'Prop. income-poverty above 5 (CV)',
    'Gini'
  )

if (year1 == '2016'){
  microDat <- microDat %>% 
    mutate(District = (str_remove(District, "0")))
}

# microDat %>%
#   select(District, contains('achelor or above')) %>%
#   tab_df(file = sprintf("districtAnalysis/microDat_BA_des%s.doc", resLabel))
# 
# microDat %>%
#   select(District, contains('employed')) %>%
#   tab_df(file = sprintf("districtAnalysis/microDat_unemployed_des%s.doc", resLabel))
# 
# microDat %>%
#   select(District, contains('1')) %>%
#   tab_df(file = sprintf("districtAnalysis/microDat_poverty100_des%s.doc", resLabel))
# 
# microDat %>%
#   select(District, contains('5')) %>%
#   tab_df(file = sprintf("districtAnalysis/microDat_poverty500_des%s.doc", resLabel))


# compute correlations between unpredictability statistics and community data
communityDat <- microDat %>%
  select(
    'District',
    'Prop. bachelor or above',
    'Prop. unemployed',
    'Prop. income-poverty below 1',
    'Prop. income-poverty above 5',
    'Gini'
  )
communityDat$District <- as.factor(communityDat$District)

corrMatVar$District <- region_names

# merge the two data frames
community_unpredictability_data <-
  inner_join(communityDat, corrMatVar, by = 'District')

#maybe turn this into a function?
community_unpredictability_data <-
  column_to_rownames(community_unpredictability_data, var = 'District')

plotDf <- formatted_cors(community_unpredictability_data)
plotDf$measure1 <-
  factor(plotDf$measure1, levels = unique(plotDf$measure1))
plotDf$measure2 <-
  factor(plotDf$measure2, levels = unique(plotDf$measure2))


plotDf <- droplevels(plotDf[!(is.na(plotDf$r)), ])


limits = c(-1, 1)
plotDf %>%
  ggplot(aes(measure1, measure2)) +
  geom_tile(aes(fill=r_if_sig),color = '#878787',size = 0.5)  + geom_text(aes(label=round(r,2)),size =5.5, angle = 45)+
  geom_segment(aes(y = 4.5, yend = 4.5, x = 4.5, xend = 11.5), color = 'black', linewidth = 2.2)+
  labs(x = NULL, y = NULL, fill = "Correlation") +
  scale_fill_distiller(
    palette = "BrBG",
    na.value = 'white',
    limits = limits,
    #breaks = breaks,
    guide = guide_colourbar(fill = guide_colourbar(
      barheight = 10, label.position = "left"
    ))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'white', size = 0),
    text = element_text(size = 25),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

ggsave(
  sprintf(
    "districtAnalysis/correlationMatPUMS_Unpredictability%s.png",
    resLabel
  ),
  units = "in",
  width = 12,
  height = 10,
  dpi = 600
)

# making a nice correlation mat
tab_corr(
  community_unpredictability_data,
  file = sprintf(
    "districtAnalysis/corrMatVar2PUMS_Unpredictability%s.doc",
    resLabel
  ),
  triangle = "lower",
  p.numeric = T,
  corr.method = c("spearman")
)



# beta regression ---------------------------------------------------------


# next to beta regression models between each outcome and each unpredictability predictor, controlling for the mean

# first standardize the variables
varList <- names(community_unpredictability_data)
varList <- varList[5:length(varList)]

community_unpredictability_data_S <-
  community_unpredictability_data %>%
  mutate(across(eval(varList), scale))

# rename for bayesian regression
names(community_unpredictability_data_S) <-
  c(
    'Prop._BA_above',
    'Prop._unemployed',
    'Prop._poverty100',
    'Prop._poverty500',
    'Standard_deviation',
    'Color_of_noise',# random component of decomposed TS
    'Entropy',
    'Autocorrelation_lag1' ,
    'Autocorrelation_lag1_stationary',
    'Changepoints_mean',
    'Changepoints_var',
    'Mean'
  )


# here is where the loop starts

yVars <- c('Prop._BA_above',
           'Prop._unemployed',
           'Prop._poverty100',
           'Prop._poverty500')
xVars <- c(
  'Standard_deviation',
  'Color_of_noise',
  # random component of decomposed TS
  'Entropy',
  'Autocorrelation_lag1' ,
  'Autocorrelation_lag1_stationary',
  'Changepoints_mean',
  'Changepoints_var'
)


resultsBetaRef <-
  data.frame(matrix(NA, nrow = length(xVars), ncol = length(yVars)))
rownames(resultsBetaRef) <- xVars
colnames(resultsBetaRef) <- yVars

lHPDBetaRef <- resultsBetaRef
hHPDBetaRef <- resultsBetaRef

for (idx in 1:length(yVars)) {
  y <- yVars[idx]
  resultModel <- list()
  lHPD <- list()
  hHPD <- list()
  for (idxx in 1:length(xVars)) {
    x <- xVars[idxx]
    
    # subset the dataframe
    dfSubset <-
      community_unpredictability_data_S[names(community_unpredictability_data_S) %in% c(x, y, 'Mean')]
    names(dfSubset)[names(dfSubset) == x] <- "x"
    names(dfSubset)[names(dfSubset) == y] <- "y"
    
    if (x == "Standard_deviation"){
      model_beta_bayes <- brms::brm(
        bf(y ~ x),
        data = dfSubset,
        family = Beta(),
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        seed = 1234,
        backend = "cmdstanr"
      )
      
    }else{
      model_beta_bayes <- brms::brm(
        bf(y ~ x+ Mean),
        data = dfSubset,
        family = Beta(),
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        seed = 1234,
        backend = "cmdstanr"
      )
      
      
    }
    
    margEffect <-
      emtrends(model_beta_bayes, ~ 1, var = 'x', regrid = "response")
    margEffectDf <- summary(margEffect)
    resultModel[idxx] <- margEffectDf[1, 2]
    lHPD[idxx] <- margEffectDf$lower.HPD
    hHPD[idxx] <- margEffectDf$upper.HPD
    
  }
  resultsBetaRef[, idx] <- unlist(resultModel)
  lHPDBetaRef[, idx] <- unlist(lHPD)
  hHPDBetaRef[, idx] <- unlist(hHPD)
  
}



# plot the results
resultsBetaRef$predictor <- rownames(resultsBetaRef)
resultsBetaRef_long <-
  gather(
    resultsBetaRef,
    factor_key = TRUE,
    key = 'outcome',
    value = 'estimate',
    -predictor
  )

lHPDBetaRef$predictor <- rownames(lHPDBetaRef)
lHPDBetaRef_long <-
  gather(
    lHPDBetaRef,
    factor_key = TRUE,
    key = 'outcome',
    value = 'lHDP',
    -predictor
  )

hHPDBetaRef$predictor <- rownames(hHPDBetaRef)
hHPDBetaRef_long <-
  gather(
    hHPDBetaRef,
    factor_key = TRUE,
    key = 'outcome',
    value = 'hHDP',
    -predictor
  )

# merge the dataframes

plottingDF <-
  inner_join(resultsBetaRef_long,
             lHPDBetaRef_long,
             by = c('predictor', 'outcome'))
plottingDF <-
  inner_join(plottingDF, hHPDBetaRef_long, by = c('predictor', 'outcome'))


newNames <-
  c(
    "Prop. bachelor or above",
    "Prop. unemployed",
    "Prop. income-poverty below 100",
    "Prop. income-poverty above 500"
  )
oldNames <- unique(plottingDF$outcome)
names(newNames) <- oldNames


predictorNames <- unique(plottingDF$predictor)
newPredictorNames <- c(
  'Standard deviation',
  'Color of noise',
  # random component of decomposed TS
  'Entropy',
  'Autocorrelation (lag 1)' ,
  'Autocorrelation (lag 1, stat.)',
  'Changepoints (mean)',
  'Changepoints (var.)',
  'Changepoints (both)'
)
names(newPredictorNames) <- predictorNames
plottingDF$predictor <- newPredictorNames[plottingDF$predictor]


ggplot(plottingDF, aes(x = predictor, y = 100*estimate, fill = predictor)) +
  geom_bar(stat = 'identity') +  scale_fill_brewer(palette = "BrBG") +
  geom_errorbar(aes(ymin = 100*lHDP, ymax = 100*hHDP),
                width = .2,
                position = position_dodge(.9)) +
  facet_wrap( ~ outcome, scales = "free_y", labeller = as_labeller(newNames)) +
  labs(x = "Predictor", y = "Estimate (in %)") +
  theme(
    #plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white', size = 0),
    legend.position = 'None',
    text = element_text(size = 30),
    axis.text.x = element_text(
      angle = 60,
      vjust = 1,
      hjust = 1.1
    )
  )
ggsave(
  sprintf("districtAnalysis/betaRegressionResult_%s.png", resLabel),
  units = "in",
  width = 30,
  height = 16,
  dpi = 600
)


# plot biavariate association between raw autocorrelation
# and pov 100 and umemployment 

# #visualizing the posterior
# posterior_beta <- model_beta_bayes %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(component = ifelse(str_detect(.variable, "phi_"), "Precision", "Mean"),
#          intercept = str_detect(.variable, "Intercept"))
#
# ggplot(posterior_beta, aes(x = .value, y = fct_rev(.variable), fill = component)) +
#   geom_vline(xintercept = 0) +
#   stat_halfeye(aes(slab_alpha = intercept),
#                .width = c(0.8, 0.95), point_interval = "median_hdi") +
#   scale_fill_viridis_d(option = "viridis", end = 0.6) +
#   scale_slab_alpha_discrete(range = c(1, 0.4)) +
#   guides(fill = "none", slab_alpha = "none") +
#   labs(x = "Coefficient", y = "Variable",
#        caption = "80% and 95% credible intervals shown in black") +
#   facet_wrap(vars(component), ncol = 1, scales = "free_y") +
#   theme_clean()
#
