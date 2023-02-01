# set working directory ---------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

# load packages -----------------------------------------------------------
library(tidyverse)
library(tidycensus)
library(srvyr, warn.conflicts = FALSE)
library(dineq)




# loading & prepping data --------------------------------------------------------------------

savLabel <- "2016-2020"

if (!file.exists(sprintf("data/nyc-pumsCommunityDistrictsReduced%s.Rdata", savLabel))) {
  
  # file does not exist yet, prepare it with this code
  load(sprintf("data/nyc-pumsCommunityDistricts%s.Rdata", savLabel))
  
  # prep: recode puma codes to match community district codes
  # based on:
  # https://www.census.gov/geographies/reference-maps/2010/geo/2010-pumas/new-york.html
  
  # 263 : 203 & 206, 221: 201 & 202, 154: 104 & 105, 121: 101 & 102,
  
  
  pumaNAMES <- c(3701:3710, 3801:3810,3901:3903,4001:4018,4101:4114)
  #districtNames <- c(208,212,210,211,263,207,205,204,209,221,112,109,110,111,108,107,154,106,103,121,503,502,501,301,304,303,302,306,308,316,305,318,317,309,307,310,312,314,315,311,313,401,403,407,411,413,408,404,406,402,405,409,412,410,414)
  
  #pumaTOcd <- setNames(as.numeric(districtNames), as.numeric(pumaNAMES))
  
  
  vt_pums$PUMA_num <- as.numeric(vt_pums$PUMA)
  
  # next exclude geographical areas that we are not interested in
  # we are only looking at community districts / PUMAs in NYC
  vt_pumsReduced <- vt_pums %>%
    filter((PUMA_num %in% pumaNAMES))
  
  vt_pumsReduced$districts <- vt_pums$PUMA_num#unname(pumaTOcd[as.character(vt_pumsReduced$PUMA_num)])
  
  # save the data
  save(vt_pumsReduced, file = sprintf("data/nyc-pumsCommunityDistrictsReduced%s.Rdata", savLabel))
  
  
  
} else { # file exists just load it
  
  load(file = sprintf("data/nyc-pumsCommunityDistrictsReduced%s.Rdata", savLabel))
  
}


# convert into survey format ----------------------------------------------

# see for variable selection: https://data.census.gov/mdat/#/

# convert data into a survey object to also get standard errors 
vt_pumsReduced$districts <- as.factor(vt_pumsReduced$PUMA) 
vt_survey_design <- to_survey(vt_pumsReduced)

# select variables we are interested in
# educational attainment BA or higher if 25 years or older
# unemployment (if 25 years or older); need to check this
# poverty to income ratio 


# first: educational attainment 

unique(vt_pumsReduced$SCHL)

if (grepl("2006", savLabel, fixed =TRUE)){
  eduList <- c("13", "14", "15", "16")
} else {
  eduList <-c("21", "22", "23", "24")
}



eduDat <- vt_survey_design %>% 
  mutate(ba_above = SCHL %in% eduList) %>% 
  filter(as.numeric(AGEP) >= 25) %>% 
  group_by(districts) %>% 
  summarize(
    ba_above_n = survey_total(ba_above, vartype = c("se","cv")),
    ba_above_pct = survey_mean(ba_above, vartype = c("se","cv"))
  )

# second: unemployment 
unemplDat <- vt_survey_design %>% 
  mutate(unemployed = ESR %in% c("3")) %>% 
  filter(as.numeric(AGEP) >= 16 & SCHG %in% c("bb","0")) %>%
  group_by(districts) %>% 
  summarize(
    unemployed_n = survey_total(unemployed, vartype = c("se", "cv")),
    unemployed_n_pct = survey_mean(unemployed, vartype = c("se", "cv"))
  )


# third: income to poverty ratio  
# make variable for the percentage of individuals under 100 and 200 percent 
# of the poverty line


povDat100 <-vt_survey_design %>% 
  mutate(pov_100 = (as.numeric(POVPIP) >= 0 & as.numeric(POVPIP) < 100)) %>% 
  filter(as.numeric(POVPIP) >= 0)%>%
  group_by(districts) %>% 
  summarize(
    pov_100n = survey_total(pov_100,vartype = c("se", "cv")),
    pov_100_n_pct = survey_mean(pov_100, vartype = c("se", "cv"))
  )


povDat500 <-vt_survey_design %>% 
  mutate(pov_500 = (as.numeric(POVPIP) > 500)) %>% 
  filter(as.numeric(POVPIP) >= 0)%>%
  group_by(districts) %>% 
  summarize(
    pov_500n = survey_total(pov_500,vartype = c("se", "cv")),
    pov_500_n_pct = survey_mean(pov_500, vartype = c("se", "cv"))
  )

povDat <- inner_join(povDat100,povDat500, by = 'districts')

povDatGini <-vt_survey_design %>% 
  filter(as.numeric(POVPIP) >= 0)%>%
  group_by(districts) %>%
  summarize(
    giniCoef = gini.wtd(as.numeric(POVPIP))
    )

povDat <- inner_join(povDat,povDatGini, by = 'districts')

# next put it all together in one dataframe and use the district labels instead. 
microDat <- inner_join(eduDat, unemplDat, by = 'districts')
microDat <- inner_join(microDat, povDat, by = 'districts')
head(microDat)

save(microDat, file = sprintf("data/nyc-microDatByDistrict%s.Rdata", savLabel) )








