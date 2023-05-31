###### MSWEP weather data for NPKD Net #######

library(tidyverse)
library(lubridate)
library(zoo)

clim <- read_csv('/Users/wilf0020/Library/CloudStorage/Dropbox/IDE/data_raw/climate/mswep/mswep_ppt_daily_all-sites_2022-11-20.csv')
comb <- read_csv('/Users/wilf0020/Library/CloudStorage/Dropbox/NPKDNet/comb-by-plot_2023-05-07.csv')
bio <- read_csv('/Users/wilf0020/Library/CloudStorage/Dropbox/NPKDNet/full-biomass_2023-05-07.csv')
trt <- read_csv('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NPKD Net/npkdnet/data/npkd-trt-start-dates.csv')
drtnet_clim <- read_csv('/Users/wilf0020/Library/CloudStorage/Dropbox/Drought-Net/IDE MS_Single year extreme/Data/precip/worldclim_interannual_cv.csv')
elev <- read_csv('/Users/wilf0020/Library/CloudStorage/Dropbox/IDE/data_processed/Site_Elev-Disturb.csv')

unique(clim$site_code)

clim_npkd <- filter(clim,site_code %in% unique(comb$site_code)) %>% mutate(year = substr(date,1,4))
unique(clim_npkd$site_code)

annual_ppt <- clim_npkd %>% mutate(year = substr(date,1,4)) %>% group_by(site_code,year) %>% summarize(ppt = sum(precip))
mappt <- annual_ppt %>% summarize(map = mean(ppt), map_sd = sd(ppt))
mappt

ppt90 <- 

# does worldclim MAP line up with MSWEP?

ppt_comp <- left_join(mappt,
                      #drtnet_clim %>% select(site_code,map) %>% distinct(site_code,.keep_all = TRUE),
                      comb %>% select(site_code,MAP_v2) %>% distinct(site_code,.keep_all = TRUE),
                      by= 'site_code',suffix=c('_mswep','_wc'))

ppt_comp
plot(ppt_comp$map,ppt_comp$MAP_v2, xlab='MAP MSWEP',ylab = 'MAP WorldClim')
#plot(ppt_comp$map,ppt_comp$MAP_v2)
abline(a=0,b=1)
# looks super

#ppt_comp %>% filter(map_mswep > 700 & map_wc < 600)

#### calculate percentage of precip received in each year previous to biomass sample dates

unique(bio$max_date) - years(1)

clim_npkd <- left_join(clim_npkd,trt, by = 'site_code') %>% 
  group_by(site_code) %>% 
  mutate(drt_ppt = if_else(date < first_treatment_date, precip, precip * ((100-drought_trt)/100)))


bio %>% distinct(site_code,year,max_date) %>% print(n=150)

#need to fix SQL view, but for now...
date_max <- bio %>%
  group_by(site_code, year) %>%
  mutate(date = max(max_date),n_nutrient_days= max(n_nutrient_days),n_drought_days= max(n_drought_days)) %>%
  distinct(site_code,date,n_nutrient_days,n_drought_days)

harvest_date <- 
  date_max %>% filter(n_drought_days > 0) %>%
  group_by(site_code) %>% 
  filter(n_nutrient_days == min(n_nutrient_days[n_nutrient_days > 0])) %>%
  mutate(md = substr(date,6,10)) %>% 
  distinct(site_code,md)

ppt90 <-  left_join(harvest_date,
                    clim_npkd %>%
    group_by(site_code) %>%
    mutate(precip_sum_90 = rollsum(precip, k = 90, fill = NA, align = "right"),
           md = substr(date,6,10)),
    by = c('site_code','md')) %>% 
  group_by(site_code) %>% 
  summarize(mean90 = mean(precip_sum_90,na.rm=T),
            sd90 = sd(precip_sum_90,na.rm = T))


avg_ppt <- left_join(mappt,ppt90, by = 'site_code')
  
date_ppt <- left_join(date_max,
  clim_npkd %>%
  group_by(site_code) %>%
  mutate(precip_sum_365 = rollsum(precip, k = 365, fill = NA, align = "right"),
         precip_sum_365_drt = rollsum(drt_ppt, k = 365, fill = NA, align = "right")) %>%
    mutate(precip_sum_90 = rollsum(precip, k = 90, fill = NA, align = "right"),
           precip_sum_90_drt = rollsum(drt_ppt, k = 90, fill = NA, align = "right")) %>%
  ungroup() %>% 
    select(site_code, date, first_treatment_date,first_nutrient_date, drought_trt,precip_sum_365,precip_sum_365_drt,
           precip_sum_90,precip_sum_90_drt),
  by = c('site_code','date'='date')) %>% 
  left_join(.,avg_ppt,by='site_code') %>% 
  mutate(ctrl_percentile_365 = (pnorm(precip_sum_365,mean=mean(map,na.rm=T),sd = map_sd))*100,
         drt_percentile_365 = (pnorm(precip_sum_365_drt,mean=mean(map,na.rm=T),sd = map_sd))*100,
         ctrl_percentile_90 = (pnorm(precip_sum_90,mean=mean(mean90,na.rm=T),sd = sd90))*100,
         drt_percentile_90 = (pnorm(precip_sum_90_drt,mean=mean(mean90,na.rm=T),sd = sd90))*100) 

setwd('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NPKD Net/npkdnet/data')
# 
# write_csv(date_ppt,
#           'MSWEP-precipitation-by-treatment-date_2023-05-11.csv')

# explore long-term fidelity using Cedar Creek data

cdcr <- read_csv('/Users/wilf0020/Library/Mobile Documents/com~apple~CloudDocs/Documents/NPKD Net/Data/cdcr-e080_Daily climate summary.csv')

cdcr_sub <- cdcr %>% 
  mutate(Date = mdy(Date), cdcr_precip = `Precip(inches)` * 25.4) %>% 
  filter(Date %in% unique(clim$date)) %>% 
  select(Date,cdcr_precip) %>% 
  left_join(., clim %>% filter(site_code == 'cedarsav.us'),
            by = c('Date' = 'date'))

summary(lm(cdcr_precip ~ precip, data = cdcr_sub))
plot(cdcr_sub$cdcr_precip ~ cdcr_sub$precip)
abline(a=0, b= 1)

cdcr_mon <- cdcr_sub %>% 
  mutate(yd = substr(Date,1,7)) %>% 
  group_by(yd) %>% 
  summarize(month_cdcr = sum(cdcr_precip), month_mswep = sum(precip))

summary(lm(month_cdcr ~ month_mswep, data = cdcr_mon))
plot(cdcr_mon$month_cdcr ~ cdcr_mon$month_mswep)
abline(a=0, b= 1)

cdcr_yr <- cdcr_sub %>% 
  mutate(y = substr(Date,1,4)) %>% 
  group_by(y) %>% 
  summarize(year_cdcr = sum(cdcr_precip), year_mswep = sum(precip))

summary(lm(year_cdcr ~ year_mswep, data = cdcr_yr))
plot(cdcr_yr$year_cdcr ~ cdcr_yr$year_mswep)
abline(a=0, b= 1)

cdcr_sub %>% filter(cdcr_precip > 80 | precip > 80)


### explore site wide fidelity using weather station data 

site_clim <- read_csv('/Users/wilf0020/Library/CloudStorage/Dropbox/Drought-Net/IDE Meeting_Oct2019/data/precip/submitted_daily_weather_bad_vals_removed_2022-04-19.csv')

unique(substr(site_clim$date,1,4))
unique(site_clim$note_weather)


site_swep <- left_join(left_join(clim_npkd,
                       site_clim %>% select(site_code,date,precip,note_weather),
                       by = c('site_code','date'),
                       suffix = c('_mswep','_station')) %>% 
  filter(!is.na(precip_station)) %>% 
  mutate(y = substr(date,1,4), ym = substr(date,1,7)) %>% 
  mutate(note = if_else(grepl("(?i)monthly|weekly", note_weather),"aggregated", 'daily')),
  elev %>% select(site_code,elev) %>% distinct(site_code,elev)) %>%
  filter(y >= 2010)

unique(site_swep$note_weather)
unique(site_swep$note)

### site wide summary and graphs

summary(lm(precip_mswep ~ precip_station, data = site_swep))
# ggplot(site_swep, aes(x = precip_station, y = precip_mswep, col = note)) + geom_point() +
#   geom_abline(slope=1, intercept = 0)


site_mon <- site_swep %>% 
  group_by(site_code,ym) %>% 
  summarize(month_station = sum(precip_station), month_mswep = sum(precip_mswep),elevation = max(elev))

summary(lm(month_station ~ month_mswep, data = site_mon))
ggplot(site_mon, aes(x = month_station, y = month_mswep, col = as.numeric(elevation))) + geom_point() +
  geom_abline(slope=1, intercept = 0)
abline(a=0, b= 1)

site_yr <- site_swep %>% 
  group_by(site_code,y) %>% 
  summarize(year_station = sum(precip_station), year_mswep = sum(precip_mswep),elevation = max(elev)) %>% 
  group_by(site_code) %>% 
  arrange(desc(year_station)) %>% 
  mutate(rank_station = row_number()) %>% 
  ungroup() %>% 
  group_by(site_code) %>% 
  arrange(desc(year_mswep)) %>% 
  mutate(rank_mswep = row_number())


summary(lm(year_station ~ year_mswep, data = site_yr))
ggplot(site_yr, aes(x = year_station, y = year_mswep, col = as.numeric(elevation))) + geom_point() +
  geom_abline(slope=1, intercept = 0) +  
  #geom_abline(slope=0.5, intercept = 0, linetype = 'dashed') +  
  #geom_abline(slope=2, intercept = 0, linetype = 'dashed') +  geom_abline(slope=0.8, intercept = 0, linetype = 'dotted') +  
  #geom_abline(slope=1.2, intercept = 0, linetype = 'dashed') + 
  xlab('Precipitation (mm/yr) - site provided') + ylab('Precipitation (mm/yr) - MSWEP')
abline(a=0, b= 1)

summary(lm(rank_station ~ rank_mswep, data = site_yr))
ggplot(site_yr, aes(x = rank_station, y = rank_mswep, col = as.numeric(elevation))) + geom_jitter() +
  geom_abline(slope=1, intercept = 0) +  geom_abline(slope=0.5, intercept = 0, linetype = 'dashed') +  
  geom_abline(slope=2, intercept = 0, linetype = 'dashed') +  geom_abline(slope=0.8, intercept = 0, linetype = 'dotted') +  
  geom_abline(slope=1.2, intercept = 0, linetype = 'dashed')
abline(a=0, b= 1)



site_yr %>% filter(year_station / year_mswep < 0.5) %>% print(n=50)
