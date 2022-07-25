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
load("data/crimeDataRegion_monthly.Rdata")
# load framework components -----------------------------------------------
source("explore.R")
source("preprocess.R")
source("framework.R")

# prepare data for framework ----------------------------------------------

# create a time variable based on the dates 
time <- as.factor(data$date)
data$time <- plyr::mapvalues(time, from = as.factor(unique(data$date)), to = c(1:length(unique(data$date))), warn_missing = F)
data$time <- as.numeric(data$time)

# user needs to specify their variable names
data_x <- "time"
data_y <- "assaultAdjusted"

# assign the variable names that are neceesary for the subsequent functions
names(data)[names(data) == data_x] <- "time"
names(data)[names(data) == data_y] <- "y"

#relabel partipant ID factor levels for subsequent analysis
data$region <- as.factor(data$region)
data$ppID <- as.factor(data$region)

levels(data$ppID) <- c(1:length(unique(data$ppID)))

#look at the data
data
origData <- data # in case decide to use this dataset after all 

# fix spelling error
head(data)
data$region <- revalue(data$region,c("Morrisania (Bronx, Harelm)"="Morrisania"))
data$region <- revalue(data$region,c("Upper East Side"="Upper\nEast Side"))



# little helper for plotting region wise statistics
region_names <- setNames(levels(data$region), levels(data$ppID))

# EXPLORE -----------------------------------------------------------------
yLabel = "Number of assaults (per 100,000)"
freq = 12
freqUnit = "Year"

# get sample size and number of repeated measures
nSample = length(unique(data$ppID))
nObs <- nrow(data)/nSample
numParticipants <- 5 # number of participants to show in the plot
nSample


nObs

pOrg <-
  explore.plotTimeSeries(origData,
                         numParticipants,
                         nSample,
                         nObs,
                         "None",
                         "None", yLabel,'raw data',yBreak = 100, freq, freqUnit,MEAN = FALSE, SD= FALSE, SE = FALSE)
pOrg

pOrgSeparate <- pOrg + facet_wrap(~ppID, ncol = 2, labeller=as_labeller(region_names)) 
pOrgSeparate

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

data_by_district <- 
  data %>% group_by(ppID)

data_by_district_split <- group_split(data_by_district)

outcome <- explore.decompose_TS(data_by_district_split[[1]] , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTS <- outcome$data_TS_decomposed
pTS + labs(title = unique(data$region)[1])


outcome <- explore.decompose_TS(data_by_district_split[[2]] , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTS <- outcome$data_TS_decomposed
pTS + labs(title = unique(data$region)[2])


outcome <- explore.decompose_TS(data_by_district_split[[3]] , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTS <- outcome$data_TS_decomposed
pTS + labs(title = unique(data$region)[3])


outcome <- explore.decompose_TS(data_by_district_split[[4]] , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTS <- outcome$data_TS_decomposed
pTS + labs(title = unique(data$region)[4])

outcome <- explore.decompose_TS(data_by_district_split[[5]] , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p
dataTS <- outcome$data_TS_decomposed
pTS + labs(title = unique(data$region)[5])


# PART 3
# autocorrelation and partial autocorrelation at different lags

lagMax = 24
diffDegree = 1
removeSeason = TRUE
logTransform = FALSE 
type = "additive" 
freq = 12

# distribution of dickey_fuller values in raw data and after attempting stationarity 
explore.visualize_stationarity(data = data,diffDegree = diffDegree,
                               removeSeason = removeSeason,
                               logTransform = logTransform, type = type, freq = freq)



outcome <- explore.visualize_autocorr(data = data, lagMax, nSample,typeAutoCorr = 'raw',
                                      diffDegree = diffDegree,
                                      removeSeason = removeSeason,
                                      logTransform = logTransform, type = type, freq = freq, regionsFlag = TRUE)

ggarrange(outcome$p_acf, outcome$p_pacf, nrow = 1)

#ggsave("autocorrelation(raw).png", units="in", width=20, height=5, dpi=600)


outcome <- explore.visualize_autocorr(data = data, lagMax, nSample,typeAutoCorr = 'stationary',
                                      diffDegree = diffDegree,
                                      removeSeason = removeSeason,
                                      logTransform = logTransform, type = type, freq = freq, regionsFlag = TRUE)


ggarrange(outcome$p_acf, outcome$p_pacf, nrow = 1)

#ggsave("autocorrelation(stationary).png", units="in", width=20, height=5, dpi=600)


# explore change points
p_cp <- explore.changePoints(data = data, minseglen = 12, yLab = 'participant', freqUnit, freq) #x-axis needs some adjustment here 
p_cp

# need to make those same adjustments for the actual analysis 

# PREPROCESS --------------------------------------------------------------
# various options are offered to the user

#1. seasonally adjust the time series

type <- 'additive'

data_SeasonalaAdj <- preprocess.adjust_season(dataTS, type)
#inspect the adjusted data
outcome <-
  explore.decompose_TS(data_SeasonalaAdj , numParticipantsDecompose, type, freq,freqUnit)
pTS <- outcome$p 
dataTSNew <- outcome$data_TS_decomposed
pTS



#2. difference the time series to achieve stationarity in mean or apply 
#   a logarithmic transformation to achieve stationarity in variance

# difference once
data_diffAdj <- preprocess.difference(data, degree =1, FALSE)
head(data_diffAdj)

outcome <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq,freqUnit)

pTS <- outcome$p 
dataTSNew <- outcome$data_TS_decomposed
pTS


p_cp <- explore.changePoints(data = data_diffAdj, minseglen = 4, yLab = 'district', freqUnit, freq) #x-axis needs some adjustment here 
p_cp


# logarithmic transformation and no differencing; reduces variance
data_diffAdj <- preprocess.difference(data, degree =0, TRUE)
head(data_diffAdj)
list[pTS,dataTSNew] <-
  explore.decompose_TS(data_diffAdj , numParticipantsDecompose, type, freq,freqUnit)
pTS

# logarithmic transformation and one round of differencing
data_diffAdj <- preprocess.difference(data, degree =1, TRUE)
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
  explore.decompose_TS(data_SeasonalaAdj , numParticipantsDecompose, type, freq, freqUnit)
pTS <- outcome$p 
dataTSNew <- outcome$data_TS_decomposed
pTS


# partition data based on different time indices  
splitValues <- c(50,100,125)
dataSplitted <- preprocess.splitData(data = data,splitValues)
# these individual data sets may be used as input to the framework function
dataEarly <- dataSplitted[[1]]
dataMiddle <- dataSplitted[[2]]
dataLate <- dataSplitted[[3]]



# ANALYZE -----------------------------------------------------------------
# if interested in polynomial effects 

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
  lm(yVar ~ time_c+time2, data = data)
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
logTransform = TRUE
type = 'additive'


modelResults <-
  framework.applyModel(data, acLAG,minseglen,percentile, meanModel, varModel, 
                       residuals = FALSE, sqaured = TRUE,diffDegree = diffDegree,
                       removeSeason = removeSeason,
                       logTransform = logTransform, freq = freq,type = type)



# to request the "statistics of the environment"
modelResults$Statistics_Mean
modelResults$Statistics_Variance
modelResults$Statistics_Residuals
modelResults$Statistics_Stationary

modelResults$Statistics_ChangePointsMBIC
glimpse(modelResults$Statistics_ChangePointsMBIC)
modelResults$Statistics_ChangePointsAIC
modelResults$Compressed_Statistics 

# returns the length of a period
# returns nonsensical values if there is no period
modelResults$dataMaxPeriod 


# plotting the models of the mean and variance
outcome <- framework.plotModel(modelResults, 5, "mean", "mean", freq, freqUnit, yLabel, segment = FALSE)
pModel <- outcome$pModel
pResiduals <- outcome$pRes

pModel
pModel_Sep <- pModel + facet_wrap(~ppID, ncol = 2, labeller=as_labeller(region_names))
pModel_Sep

pResiduals
pResiduals_Sep <- pResiduals + facet_wrap(~ppID, ncol = 2, labeller=as_labeller(region_names))
pResiduals_Sep

outcome <- framework.plotModel(modelResults, 5, "variance", "variance", freq, freqUnit, yLabel, segment = FALSE)
pModel <- outcome$pModel
pResiduals <- outcome$pRes

pModel
pModel_Sep <- pModel + facet_wrap(~ppID, ncol = 2, labeller=as_labeller(region_names))
pModel_Sep

pResiduals
pResiduals_Sep <- pResiduals + facet_wrap(~ppID, ncol = 2, labeller=as_labeller(region_names))
pResiduals_Sep


#plot color of noise distribution 
# function outputs the same output for a function with and without trend 
pColorOfNoise <- framework.plotColorOfNoise(modelResults$Statistics_Mean$m_spectralCoef, binwidth = 0.01)
pColorOfNoise


#visualize change points in mean and variance where both are estimated individually
pNew <-framework.plotCP(modelResults, 5, "Changepoints in mean and variance","AIC","individual" ,freq, freqUnit, yLabel)
pNew_Sep <- pNew + facet_wrap(~ppID, ncol = 2,labeller=as_labeller(region_names))
pNew_Sep

# where both estimated together 
pNew <-framework.plotCP(modelResults, 5, "Changepoints in mean and variance","MBIC","together" ,freq, freqUnit, yLabel)
pNew_Sep <- pNew + facet_wrap(~ppID, ncol = 2, labeller=as_labeller(region_names))
pNew_Sep

# min and max plots for different statistics 
dataModel <- modelResults$data
time <- modelResults$data[[1]]$time

framework.plotMinMax(time, modelResults$Compressed_Statistics$chronicity_score, dataModel, 'chronicity score')
framework.plotMinMax(time, modelResults$Compressed_Statistics$intensity_score, dataModel, 'intensity score')
framework.plotMinMax(time, modelResults$Compressed_Statistics$scaled_mean, dataModel, 'scaled mean')

framework.plotMinMax(time, modelResults$Statistics_Mean$m_slope1, dataModel, 'mean slope')
framework.plotMinMax(time, modelResults$Statistics_Variance$v_slope1, dataModel, 'variance slope')

framework.plotMinMax(time, modelResults$Statistics_ChangePointsMBIC$meanNumberMBIC, dataModel, 'changepoints mean')
framework.plotMinMax(time, modelResults$Statistics_ChangePointsMBIC$varNumberMBIC, dataModel, 'changepoints variance')


framework.plotMinMax(time, modelResults$Statistics_ChangePointsAIC$meanNumberAIC, dataModel, 'changepoints mean')
framework.plotMinMax(time, modelResults$Statistics_ChangePointsAIC$varNumberAIC, dataModel, 'changepoints variance')


framework.plotMinMax(time, abs(modelResults$Statistics_Mean$m_spectralCoef), dataModel, 'spectral coef')

framework.plotMinMax(time, abs(modelResults$Statistics_Stationary$stat_ac1), dataModel, 'autocorrelation data')
framework.plotMinMax(time, abs(modelResults$Statistics_Stationary$statVar_ac1), dataModel, 'autocorrelation variance ')


framework.plotMinMax(time, modelResults$Statistics_Mean$m_apprEntropy, dataModel, 'entropy data')
framework.plotMinMax(time, modelResults$Statistics_Stationary$stat_apprEntropy, dataModel, 'entropy stationary')
framework.plotMinMax(time, modelResults$Statistics_Variance$v_apprEntropy, dataModel, 'entropy variance')
framework.plotMinMax(time, modelResults$Statistics_Residuals$r_apprEntropy, dataModel, 'entropy residuals')


# histograms of statistics 
histStat <- ggplot(data = modelResults$Statistics_Mean, aes(m_apprEntropy))+
  geom_histogram(binwidth = 0.02) +
  theme(legend.position = "none") + theme_bw() + xlab('Entropy')+ labs(title = 'Distribution of entropy values')+
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = c(0.06, 0.75))

histStat


# visualize ac of variance
dataVar <- modelResults%>% 
     unnest(data) %>%
  select(ppID, yVar)
dataVar$y <- dataVar$yVar


dataVar$region <- as.factor(region_names[dataVar$ppID])



outcome <- explore.visualize_autocorr(data = dataVar, lagMax, nSample,typeAutoCorr = 'raw',
                                      diffDegree = diffDegree,
                                      removeSeason = removeSeason,
                                      logTransform = logTransform, type = type, freq = freq, regionsFlag = TRUE)

outcome <- explore.visualize_autocorr(data = dataVar, lagMax, nSample,typeAutoCorr = 'stationary',
                                      diffDegree = diffDegree,
                                      removeSeason = removeSeason,
                                      logTransform = logTransform, type = type, freq = freq, regionsFlag = TRUE)


# correlation among measures ----------------------------------------------

# to what extent do measures of unpredictability correlate?
# autocorrelation lag 1 low autocorrelation means high unpred
# entropy high means high
# sd high means high
# change points mean high means high
# change points variance high mean high


dfVis <- bind_cols(modelResults$Statistics_Mean,modelResults$Statistics_Residuals, modelResults$Statistics_Stationary,
                   modelResults$Statistics_ChangePointsAIC, modelResults$Compressed_Statistics, .name_repair = 'minimal')

glimpse(dfVis)


corrMatVar <- dfVis %>%
  select(m_sd,m_spectralCoef,m_apprEntropy,m_ac1, stat_ac1,meanNumberAIC,varNumberAIC,npNumberAIC )

names(corrMatVar) <- c('sd (raw)', 'color of noise', 'entropy', 'ac (raw)', 'ac (stationary)','CP (mean)', 'CP (variance)', 'CP (both)')
corrMatVar$region <- region_names


corrMatVar <- remove_rownames(corrMatVar) 
corrMatVar <- column_to_rownames(corrMatVar, var = 'region')

#continue here and transform this into a nice table 
# then make a new table with different ranks for each variable
corrMatVar


M <- cor(corrMatVar)
head(round(M,2))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method = "color", col=col(200),type = 'upper',
         tl.col="black", tl.srt=45,addCoef.col = "black", diag = F, order="hclust")


# compare different measures of intensity / chronicity 
corrMatVar <- dfVis %>%
  select(m_mean, chronicity_score, intensity_score, scaled_mean)


M <- cor(corrMatVar)
head(round(M,2))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method = "color", col=col(200),type = 'upper',
         tl.col="black", tl.srt=45,addCoef.col = "black", diag = F)



# for publication ---------------------------------------------------------



##### Table 1 with unpredictability statistics #####
dfVis <- bind_cols(modelResults$Statistics_Mean,modelResults$Statistics_Residuals, modelResults$Statistics_Stationary,
                   modelResults$Statistics_ChangePointsAIC, modelResults$Compressed_Statistics, .name_repair = 'minimal')

glimpse(dfVis)

corrMatVar <- dfVis %>%
  select(m_sd,
         m_spectralCoef,
         m_apprEntropy,
         m_ac1,
         meanNumberAIC,
         varNumberAIC) #stat_pac6,

names(corrMatVar) <-
  c(
    'Standard deviation',
    'Color of noise',
    'Entropy',
    'Autocorrelation (lag 1)' ,
    'Changepoints (mean)',
    'Changepoints (var.)'
  ) #'AC (lag 6)',
corrMatVar$region <- region_names

corrMatVar <- remove_rownames(corrMatVar)
corrMatVar <- column_to_rownames(corrMatVar, var = 'region')

# continue here and transform this into a table
# then make a new table with different ranks for each variable
corrMatVar <- round(corrMatVar, 2)
corrMatVar <- rownames_to_column(corrMatVar, var = "Region")

# use package formattable
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"
customGrey = "#CACFD2"

# helper functions for creating the bars
myNormalize <- function(x,absFlag, highest){
  if(absFlag){
    x <- abs(x)
  }
  if(highest){
    normVal <- (x-min(x))/(max(x)-min(x))
    
  }else{
    normVal <- (x-min(x))/(max(x)-min(x))
    normVal <- 1-normVal
  }
  return(normVal)
}



f.color_bar <- function (color = "lightgray", fun = "proportion", ...) 
{
  fun <- match.fun(fun)
  formattable::formatter("span", style = function(x) style(display = "inline-block", 
                                                           direction = "ltr", `border-radius` = "4px", 
                                                           `padding-right` = "2px", `background-color` = csscolor(color), 
                                                           width = percent(fun(as.numeric(x), ...))))
}


#integrate formattable and kableExtra 
corrMatVar %>%
  mutate(
    Region = formatter("span", style = ~ style(color = "black"))(Region),
    `Standard deviation` = f.color_bar(customGrey, function(x) myNormalize(x,TRUE,TRUE))(`Standard deviation`),
    `Color of noise` = f.color_bar(customGrey,function(x) myNormalize(x,TRUE,FALSE))(`Color of noise`),
    `Entropy` = f.color_bar(customGrey,function(x) myNormalize(x,TRUE,TRUE))(`Entropy`),
    `Autocorrelation (lag 1)` = f.color_bar(customGrey, function(x) myNormalize(x,TRUE,FALSE))(`Autocorrelation (lag 1)`),
    `Changepoints (mean)` = f.color_bar(customGrey,function(x) myNormalize(x,TRUE,TRUE))(`Changepoints (mean)`),
    `Changepoints (var.)` = f.color_bar(customGrey,function(x) myNormalize(x,TRUE,TRUE))(`Changepoints (var.)`)
  ) %>%
  kable("html", escape = F) %>%
  kable_styling("hover", full_width = F, font_size = 20) %>%
  row_spec(0, color = "black", bold =T) %>%
  column_spec(1, width = "4cm",color = 'black', italic = T)%>%
  column_spec(2, width = "3cm",color = 'black')%>%
  column_spec(3, width = "3cm",color = 'black')%>%
  column_spec(4, width = "3cm",color = 'black')%>%
  column_spec(5, width = "3cm",color = 'black')%>%
  column_spec(6, width = "3cm",color = 'black')%>%
  column_spec(7, width = "3cm",color = 'black')%>%
  #column_spec(8, width = "3cm",color = 'black')%>%
  #column_spec (1:7,border_left = T, border_right = T)%>%
  kable_classic()

#%>%
#  save_kable("tableNew.png",density = 600, zoom = 2)
  

# plot regionwise change points

sampleList <- c(1,2,3,4,5)
source("framework.R")
yLabel = "Number of assaults (per 100,000)"
unitOffset <- 2006
statsPositions <- c(27,80,145)

pNew <-
  framework.plotPublicationTest(
    modelResults,
    5,
    "Regions of interest",
    "AIC",
    "individual" ,
    freq,
    freqUnit,
    yLabel,
    sampleList = sampleList, regionFlag = TRUE,
    unitOffset = unitOffset,
    stats_positions = statsPositions)

pNew

ggsave("regionCPNew2Test2.png", units="in", width=10, height=10, dpi=600)


