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

load("data/nyc-spatialCommunityDistricts.Rdata")
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
# information is based on: nyc_decennialcensusdata_2010_2020_change
boro <- c("M","B","K","Q","S")
population <- (c(1585873,1385108,2504700,2230722,468730)+c(1694251,1472654,2736074,2405464,495747))/2

nyc_pop_boro <- tibble(
  boro        = boro,
  populationB  = population
)

nyc_pop_district <- read.csv2("data/district_popPUMA.csv", header=TRUE)

names(nyc_pop_district) <- c("district", "p_2010")
nyc_pop_district$populationD <- (as.numeric(nyc_pop_district$p_2010)) 


nyc_pop_district <- as_tibble(nyc_pop_district)
nyc_pop_district$district <- as.factor(nyc_pop_district$district)

# correct for borough population
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

ggsave("boroughAnalysis/boroughs_acrossTime.png", units="in", width=20, height=15, dpi=600)

# for each borough plot the monthly ts and 2 random districts 

# Brooklyn     Manhattan        Queens Staten Island     The Bronx 
# "#A6611A"     "#DFC27D"     "#a9a9a9"     "#80CDC1"     "#018571" 
regionNames <- c('K', 'M', 'Q', 'S', 'B')
colorNames <- c("#A6611A","#DFC27D","#a9a9a9","#80CDC1","#018571")
names(colorNames) = regionNames


p_boroughs_monthly <-
  ggplot(data = dataBorough_monthly, aes(x = date, y = assaultAdjustedB, group = boro)) +
  geom_line(aes(color = boro), linewidth = 1.3) +
  theme_bw() + labs(y = 'Number of assaults (per 100,000)\n', x = 'Year', color = 'Borough') +
  scale_color_manual(values = colorNames) +
  scale_x_date(
    expand = c(0.01, 0.01) ,
    limits = as.Date(c('2006-01-01', '2021-01-01')),
    breaks = seq.Date(as.Date("2006-01-01"), as.Date("2021-01-01"), by = "1 year"),
    labels = c(2006:2021)
  ) +  #guides(colour = guide_legend(override.aes = list(size=10, linewidth =3)))+
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = 'None',
    text = element_text(size = 40),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )
  ) + theme(strip.background = element_blank()) + scale_y_continuous(expand = c(0.01, 0.01))

print(p_boroughs_monthly)

# next two random districts per borough

selectionDf <- dataBoroughDistrict_monthly %>%
  group_by(boro)%>%
  dplyr::count(district) %>%
  filter(n ==180) %>% group_by(boro)%>%
  slice_sample(n= 2)

plt1 <- droplevels(selectionDf[seq_along(selectionDf$district) %% 2 > 0,])
plt2 <- droplevels(selectionDf[seq_along(selectionDf$district) %% 2 == 0,])


plt2$district[plt2$boro == 'Q'] <- '4102' 

plt1_df <- dataBoroughDistrict_monthly %>%
  filter(boro == plt1$boro[1] & district == plt1$district[1])


for (idx in 2:nrow(plt1)){
  currBoro <- plt1$boro[idx]
  currDistrict <- plt1$district[idx]
  currDF <- dataBoroughDistrict_monthly %>%
    filter(boro == currBoro & district == currDistrict)
  plt1_df <- full_join(plt1_df, currDF)
}



plt2_df <- dataBoroughDistrict_monthly %>%
  filter(boro == plt2$boro[1] & district == plt2$district[1])


for (idx in 2:nrow(plt2)){
  currBoro <- plt2$boro[idx]
  currDistrict <- plt2$district[idx]
  currDF <- dataBoroughDistrict_monthly %>%
    filter(boro == currBoro & district == currDistrict)
  plt2_df <- full_join(plt2_df, currDF)
}



test <- unique(plt1_df$district)
names(test) <- unique(plt1_df$boro)
colorNamesPlt1 <- colorNames
names(colorNamesPlt1) <- test[names(colorNamesPlt1)]



# number 2

test <- unique(plt2_df$district)
names(test) <- unique(plt2_df$boro)
colorNamesPlt2 <- colorNames
names(colorNamesPlt2) <- test[names(colorNamesPlt2)]


# try a different plot 
plt_df <- rbind(plt1_df,plt2_df)


colorNames_pltDistricts <- c(colorNamesPlt1, colorNamesPlt2)



regions <- c("Brooklyn", "Manhattan", "Queens", "Staten Island", "The Bronx")
names(regions) <- regionNames 

p_districts <-
  ggplot(data = plt_df, aes(x = date, y = assaultAdjustedD, group = district)) +
  geom_line(aes(color = district), size = 0.7) +
  theme(legend.position = "none") + theme_bw() + labs(y = '', x = '', color = 'District') +
  scale_color_manual(values = colorNames_pltDistricts) +
  facet_wrap(vars(boro),nrow = 1, ncol = length(unique(plt_df$boro)), labeller = as_labeller(regions))+
  scale_x_date(
    expand = c(0.01, 0.01) ,
    limits = as.Date(c('2006-01-01', '2021-01-01')),
    breaks = c(),#seq.Date(as.Date("2006-01-01"), as.Date("2021-01-01"), by = "1 year"),
    labels = c()
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15) ,
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'grey70', size = 1),
    legend.position = 'None',
    text = element_text(size = 40),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )
  ) + theme(strip.background = element_blank()) + scale_y_continuous(expand = c(0.01, 0.01))

p_districts

# putting them all together 


# now put them both together with cowplot 


plot_grid(p_boroughs_monthly, p_districts, scale = 0.9,labels = "AUTO",label_size = 35, 
          ncol =1, rel_widths = c(1, 3), rel_heights = c(2,1),
          label_x = 0, label_y = 0.94,
          hjust = -0.8, vjust = -1)



ggsave("boroughAnalysis/boroughsAndDistricts.png", units="in", width=25, height=20, dpi=600)



# save data ---------------------------------------------------------------
save(dataBorough_monthly, file = "data/monthly_borough.Rdata")
save(dataDistrict_monthly, file = "data/monthly_districts.Rdata")
save(dataBoroughDistrict_monthly, file = "data/monthly_borough_districts.Rdata")




