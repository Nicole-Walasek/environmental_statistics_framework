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

resLabel <- "entire"

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

data$borough <-round(as.numeric(as.character(data$district))/100)

boroughs <- c("M", "B", "K","Q", "S")
names(boroughs) <- c(38,37,40,41,39)


data$borough <- as.factor(unname(boroughs[as.character(data$borough)]))

# in case you want to filter the data by year ----------------------
year1 <- str_split(resLabel, "_")[[1]][1]
year2 <- str_split(resLabel, "_")[[1]][2]
year2 <- as.character((as.numeric(year2) + 1))

data <- data %>%
  filter(date >= as.Date(sprintf("%s-01-01", year1)) &
           date < as.Date(sprintf("%s-01-01", year2)))



# compute descriptives

head(data)


result <- data %>% 
  group_by(borough) %>%
  dplyr::summarize(n = length(unique(district)),
                   averageCrime = mean(assaultD),
                   sdCrime = sd(assaultD))

result %>%
  tab_df(file = sprintf("districtAnalysis/descriptives%s.doc", resLabel))

# next descriptives for the microdata
# average number of value, standard deviation and coeficient of variation 



load(file = sprintf("data/nyc-microDatByDistrict%s.Rdata", resLabel))
head(microDat)
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
    'Income-poverty below 100',
    'Income-poverty below 100 (SE)',
    'Income-poverty below 100 (CV)',
    'Prop. income-poverty below 100',
    'Prop. income-poverty below 100 (SE)',
    'Prop. income-poverty below 100 (CV)',
    'Income-poverty above 500',
    'Income-poverty above 500 (SE)',
    'Income-poverty above 500 (CV)',
    'Prop. income-poverty above 500',
    'Prop. income-poverty above 500 (SE)',
    'Prop. income-poverty above 500 (CV)',
    'Gini'
  )

if (year1 == '2016'){
  microDat <- microDat %>% 
    mutate(District = (str_remove(District, "0")))
}


resultMicroDat <- microDat %>%
  select(c(`Prop. bachelor or above`, `Prop. bachelor or above (CV)`,
           `Prop. unemployed`,`Prop. unemployed (CV)`,
           `Prop. income-poverty below 100`, `Prop. income-poverty below 100 (CV)`,
           `Prop. income-poverty above 500`, `Prop. income-poverty above 500 (CV)`))%>%
  gather(varType, value)
  

resultMicroDat$varType <- as.factor(resultMicroDat$varType)

resultMicroDatFinal <- resultMicroDat %>%
  group_by(varType) %>%
  dplyr::summarise(m = mean(value),
                   sd = sd(value))



resultMicroDatFinal %>%
  tab_df(file = sprintf("districtAnalysis/descriptivesMicroDat%s.doc", resLabel))














