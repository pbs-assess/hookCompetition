# Script for paper to estimate the relative abundance of Pacific Halibut and Yelloweye

# Load packages
library(plyr)
library(rgeos)
library(rgdal)
library(maptools)
library(sp)
library(spdep)
library(INLA)
library(inlabru)
library(ggplot2)
library(tidyverse)
library(mgcv)
library(GGally)
library(boot)
library(rerddap)
library(raster)
library(brms)
library(tidymv)

setwd("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final")

fit_mods <- F
#Again - we advise getting a pardiso license
#inla.setOption(pardiso.license="pardiso.lic")
#inla.pardiso.check()

# Load the bootstrap indices of abundance from Anderson et al 2019
pacific_halibut_boot <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Andy Indices Scripts/all_species_results/all-species-results/pacific-halibut-results.rds")[[2]][[1]]
yelloweye_rockfish_boot <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Andy Indices Scripts/all_species_results/all-species-results/yelloweye-rockfish-results.rds")[[2]][[1]]

# Load the censored Poisson indices
combined_indices_halibut_INLA_Alldata <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final/combined_indices_halibut_INLA_Alldata.rds")
combined_indices_yelloweye_INLA_Alldata <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final/combined_indices_yelloweye_INLA_Alldata.rds")

# Load the catch data
load('WS_Modelling_Alldata.RData')

# Set the seed
seed <- 06082021 # today's date

# Load the function that fits the censored Poisson models
# Load all the functions required for analysis
source('../Deliverable_2_Final/All_Functions.R')

All_data_sp@data$twentyhooks=ifelse(
  is.na(as.numeric(as.matrix(All_data_sp@data[,paste0('N_itAll_halibut')])[,1])),
  1,0)

Reduced_data_sp <- 
  All_data_sp[which(All_data_sp$station %in% (names(table(All_data_sp$station))[table(All_data_sp$station)>1])),]

Reduced_data_sp <- Reduced_data_sp[which(!(Reduced_data_sp$year %in% c(2012))),]

# NOTE THAT 2013 USED ONLY FIRST 20 HOOKS BUT THE OBSERVED NUMBER OF HOOKS IS NA. 
# MAP 1997'S EFFECTIVE SKATE VALUES TO 2013 VALUES TO OBTAIN THE NHOOKS.
Reduced_data_sp[which(Reduced_data_sp$year==2013),]$obsHooksPerSet
mean_effskate_perhook_1997 <-
  by(Reduced_data_sp[which(Reduced_data_sp$year==1997),]$effSkateIPHC,
     INDICES = as.factor(
       Reduced_data_sp[which(Reduced_data_sp$year==1997),]$obsHooksPerSet),
     FUN=mean
  )
# map 2013 values of eff skate to the closest matching n hooks
Reduced_data_sp@data[which(Reduced_data_sp$year==2013),]$obsHooksPerSet <-
  c(as.numeric(names(mean_effskate_perhook_1997))[
    apply(
      outer(Reduced_data_sp[which(Reduced_data_sp$year==2013),]$effSkateIPHC,mean_effskate_perhook_1997,'-')^2,
      1, which.min
    )
  ])
# update prop_removed variable
Reduced_data_sp@data[which(Reduced_data_sp$year==2013),]$prop_removed <-
  (Reduced_data_sp@data[which(Reduced_data_sp$year==2013),]$obsHooksPerSet -
     Reduced_data_sp@data[which(Reduced_data_sp$year==2013),]$N_it_hook) / 
  Reduced_data_sp@data[which(Reduced_data_sp$year==2013),]$obsHooksPerSet

## Exploratory analysis - investigate a reasonable value for cprop
binomial_smooth <- function(...) {
  geom_smooth(method='gam',formula = y ~ -1 + x + s(x, bs='cs'),method.args=list(family='quasibinomial'), ...)
}

All_data_sp@data %>%
  filter(region_INLA==2) %>%
  ggplot(aes(x=prop_removed, y=N_it_halibut/(obsHooksPerSet-N_it_hook), group=factor(region_INLA),colour=factor(region_INLA))) +
  binomial_smooth(aes(weight=obsHooksPerSet-N_it_hook), fullrange=T) +
  xlim(c(0.7,1)) + 
  ylab('Average catch proportion of halibut') +
  xlab('Proportion of baited hooks removed') +
  guides(fill='none',colour='none')
# Sharp changepoint after 95% of baits removed

All_data_sp@data %>%
  filter(region_INLA==2) %>%
  ggplot(aes(x=prop_removed, y=N_it_yelloweye/(obsHooksPerSet-N_it_hook), group=factor(region_INLA),colour=factor(region_INLA))) +
  binomial_smooth(aes(weight=obsHooksPerSet-N_it_hook), fullrange=T) +
  xlim(c(0.7,1)) + 
  ylab('Average catch proportion of yelloweye') +
  xlab('Proportion of baited hooks removed') +
  guides(fill='none',colour='none')
# gradual changepoint after 90% of baits removed

All_data_sp@data %>%
  filter(region_INLA==2) %>%
  ggplot(aes(x=prop_removed, y=N_it_halibut/effSkateIPHC, group=factor(region_INLA),colour=factor(region_INLA), weight=effSkateIPHC)) +
  geom_smooth() +
  xlim(c(0.7,1)) + 
  ylab('Average CPUE of halibut') +
  xlab('Proportion of baited hooks removed') +
  guides(fill='none',colour='none')
# Sharp changepoint after 95% of baits removed

All_data_sp@data %>%
  filter(region_INLA==2) %>%
  ggplot(aes(x=prop_removed, y=N_it_yelloweye/effSkateIPHC, group=factor(region_INLA),colour=factor(region_INLA), weight=effSkateIPHC)) +
  geom_smooth() +
  xlim(c(0.7,1)) + 
  ylab('Average CPUE of yelloweye') +
  xlab('Proportion of baited hooks removed') +
  guides(fill='none',colour='none')
# no change

# Based on plots, fit the censored Poisson model with cprop = 0.95 and upper quantile = 0.8
#halibut_cenPoisson <- modified_censored_index_fun_Alldata(Reduced_data_sp, species='halibut',M=1000,cprop=0.95, upper_bound_quantile2 = 0.8, upper_bound_quantile = 0.8)
#halibut_cenPoisson$pred_poisson
#halibut_cenPoisson <- final_censored_index_fun(Reduced_data_sp, species='halibut',M=1000,cprop=0.95, upper_bound_quantile2 = 1, use_upper_bound = F)


# Repeat for yelloweye rockfish with cprop = 0.975 and upper quantile = 0.6
# Note that these values were required for convergence!
#yelloweye_cenPoisson <- modified_censored_index_fun_Alldata(Reduced_data_sp, species='yelloweye',M=1000,cprop=0.975, upper_bound_quantile2 = 0.6, upper_bound_quantile = 0.6)

mean(Reduced_data_sp$prop_removed>=0.95, na.rm=T) # 49% of fishing events return >= 95% hook saturation
mean(Reduced_data_sp$prop_removed>=0.975, na.rm=T) # 36% of fishing events return >= 97.5% hook saturation

# Try instead to use brms
brms_dat <- Reduced_data_sp@data[which(!is.na(Reduced_data_sp$prop_removed)),]
brms_dat$censored_txt_95 <- rep('none',dim(brms_dat)[1])
brms_dat$censored_txt_95[which(brms_dat$prop_removed>=0.95)] <- 'right'
brms_dat$censored_txt_80 <- rep('none',dim(brms_dat)[1])
brms_dat$censored_txt_80[which(brms_dat$prop_removed>=0.80)] <- 'right'
brms_dat$year <- brms_dat$year - min(brms_dat$year) + 1
brms_dat$station_ID <- as.factor(brms_dat$station)
brms_dat$event_ID <- factor(1:length(brms_dat$year))
brms_dat$region <- factor(brms_dat$region_INLA)
nyear <- length(unique(brms_dat$year))
brms_dat <- brms_dat[!is.na(brms_dat$region),]

Exploratory_halibut <-
  gam(N_it_halibut ~ region + twentyhooks + s(year, by = region) + s(station_ID, bs='re') + s(prop_removed),
      data=brms_dat, family = 'nb', offset = log(brms_dat$effSkateIPHC*brms_dat$prop_removed))
plot(Exploratory_halibut, select=6, scale=0, 
     ylab='Estimated Effect on Log Scale', 
     xlab='Proportion of Baits Removed', 
     main='Estimated Change in Log Halibut Catch Count \nvs. Proportion of Baits Removed ')

Exploratory_yelloweye <-
  gam(N_it_yelloweye ~ region + twentyhooks + s(year, by = region) + s(station_ID, bs='re') + s(prop_removed),
      data=brms_dat, family = 'nb', offset = log(brms_dat$effSkateIPHC*brms_dat$prop_removed))
plot(Exploratory_yelloweye, select=6, scale=0, 
     ylab='Estimated Effect on Log Scale', 
     xlab='Proportion of Baits Removed', 
     main='Estimated Change in Log Yelloweye Catch Count \nvs. Proportion of Baits Removed ')

# Based on the plots, both species appear to suffer from gear saturation when >95%
# of hooks are removed.

# What is the trend with effective skate
Exploratory_halibut2 <-
  gam(N_it_halibut ~ region + twentyhooks + s(year, by = region) + s(station_ID, bs='re') + s(I(log(effSkateIPHC))),
      data=brms_dat, family = 'nb', offset = log(brms_dat$effSkateIPHC))
plot(Exploratory_halibut2, select=6, scale=0, 
     ylab='Estimated Nonlinear Effect on Log Scale', 
     xlab='Log Effective Skate',
     main='Estimated Change in Log Halibut Catch Count \nvs. Log Effective Skate ')

Exploratory_yelloweye2 <-
  gam(N_it_yelloweye ~ region + twentyhooks + s(year, by = region) + s(station_ID, bs='re') + s(I(log(effSkateIPHC))),
      data=brms_dat, family = 'nb', offset = log(brms_dat$effSkateIPHC))
plot(Exploratory_yelloweye2, select=6, scale=0, 
     ylab='Estimated Nonlinear Effect on Log Scale', 
     xlab='Log Effective Skate', 
     main='Estimated Nonlinear Change in Log Yelloweye Catch Count \nvs. Log Effective Skate')

# Plot the average station locations and mark those that are only online for 2 years
all_stations <- Reduced_data_sp[which(!is.na(Reduced_data_sp$region_INLA)),]
count <- 1
for(i in unique(all_stations$station))
{
  station <- gCentroid(all_stations[which(all_stations$station==i),])
  station$n_obs <- length(which(all_stations$station==i))
  station$region <- all_stations@data$region_INLA[which(all_stations$station==i)[1]]
  if(count==1)
  {
    stations_sp <- station
  }
  if(count>1)
  {
    stations_sp <- rbind(stations_sp,station)
  }
  count <- count + 1
}

Coast <- readOGR('~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final/Shapefiles',
                 layer='Coastline_Medres')
Coast <- spTransform(Coast, Reduced_data_sp@proj4string)
Coast <- gSimplify(Coast, tol = 1)
table(stations_sp$n_obs) # 114 / 279 visited twice. The rest visited for 18-22 years

table(all_stations$year[all_stations$station %in% unique(all_stations$station)[ stations_sp$n_obs==2]])
# all these temporary stations were in 1996 and 1997
sum(all_stations$year %in% c(1996, 1997))
# These were the ONLY stations online!

ggplot() + 
  gg(stations_sp, colour=stations_sp$region+1, shape=stations_sp$n_obs==2) +
  gg(Coast) + xlab('Eastings (km)') + ylab('Northings (km)') +
  coord_fixed() +
  ggtitle('The locations of all 279 IPHC fishing stations', 
          subtitle = 'The colour denotes the region, a circle plotting character denotes the \nstation was temporary, with a square denoting a permanent station')

# Note that the yelloweye empirical count distribution is very right skew
hist(Reduced_data_sp$N_it_yelloweye)
hist(Reduced_data_sp$N_it_halibut)

sd(Reduced_data_sp$N_it_yelloweye, na.rm=T)/mean(Reduced_data_sp$N_it_yelloweye, na.rm=T)
# CV is 2.942808 and sd(IID) in censored model is 1.07
sd(Reduced_data_sp$N_it_halibut, na.rm=T)/mean(Reduced_data_sp$N_it_halibut, na.rm=T)
# CV is 1.118502 and sd(IID) in censored model is 0.66

setwd("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Paper")

INLA_dat <- Reduced_data_sp[which(!(is.na(Reduced_data_sp$region_INLA) |
                                      is.na(Reduced_data_sp$prop_removed))),]

if(fit_mods)
{
  brms_mod_halibut <-
    final_censored_index_fun(data=INLA_dat,
                             species = 'halibut',
                             cprop=1.1, 
                             use_upper_bound = F, 
                             upper_bound_quantile = 1)
}
if(!fit_mods)
{
  brms_mod_halibut <- readRDS('halibut_mod_final_INLA.rds')
}

brms_mod_halibut_pred <- brms_mod_halibut$pred_overdisp
brms_mod_halibut_pred$model = 'CPUE-based'

# Now fit the hook competition scale factor-adjusted model
INLA_dat$N_it_halibut_compfactor <-
  round(INLA_dat$N_it_halibut * 
          comp_factor_fun(INLA_dat$prop_removed, INLA_dat$obsHooksPerSet))
INLA_dat2 <- INLA_dat
INLA_dat2$N_it_halibut <- INLA_dat2$N_it_halibut_compfactor

if(fit_mods)
{
  brms_mod_halibut_compfactor <-
    final_censored_index_fun(data=INLA_dat2,
                             species = 'halibut',
                             cprop=1.1, 
                             use_upper_bound = F, 
                             upper_bound_quantile = 1)
}
if(!fit_mods)
{
  brms_mod_halibut_compfactor <- readRDS('halibut_compfactor_mod_final.rds')
}
summary(brms_mod_halibut_compfactor)

brms_mod_halibut_compfactor_pred <- brms_mod_halibut_compfactor$pred_overdisp
brms_mod_halibut_compfactor_pred$model = 'Instantaneous catch-rate adjusted'

# Now fit the censored model
if(fit_mods)
{
  brms_mod_halibut_cens <-
    final_censored_index_fun(data=INLA_dat,
                             species = 'halibut',
                             cprop=0.95, 
                             use_upper_bound = T, 
                             upper_bound_quantile = 0.8)
}
if(!fit_mods)
{
  brms_mod_halibut_cens <- readRDS('halibut_cens_mod_final.rds')
}

brms_mod_halibut_cens_pred <- brms_mod_halibut_cens$pred_overdisp
brms_mod_halibut_cens_pred$model = 'Censored Poisson'

boot_halibut_pred <- brms_mod_halibut_pred
# recode factor to match brm mod output
# combined_indices_halibut_INLA_Alldata$combined_ts_all$region <-
#   forcats::fct_recode(combined_indices_halibut_INLA_Alldata$combined_ts_all$region,
#                       `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
#                       `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI',
#                       `All` = 'All')
# boot_halibut_pred$region <-
#   forcats::fct_recode(boot_halibut_pred$region,
#                       `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
#                       `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')
# # calibrate the year variable
# combined_indices_halibut_INLA_Alldata$combined_ts_all$year <- 
#   combined_indices_halibut_INLA_Alldata$combined_ts_all$year 
# match by year and region
boot_halibut_pred <- inner_join(
  boot_halibut_pred, combined_indices_halibut_INLA_Alldata$combined_ts_all[
    combined_indices_halibut_INLA_Alldata$combined_ts_all$model=='boot',
  ],
  by = c('region','year')
)
boot_halibut_pred$mean <- boot_halibut_pred$mean.y
boot_halibut_pred$q0.025 <- NA#boot_halibut_pred$LCL
boot_halibut_pred$q0.975 <- NA#boot_halibut_pred$UCL
boot_halibut_pred$model <- boot_halibut_pred$model.y
boot_halibut_pred$Est.Error <- NA

# keep only the relevant variables
boot_halibut_pred <- boot_halibut_pred[,names(brms_mod_halibut_pred)]

# Add the hook saturation to the dataframe
hook_sat_df <-
  Reduced_data_sp@data %>%
  group_by(year, region_INLA) %>%
  summarise(Estimate2 = mean(prop_removed[!is.nan(prop_removed)]>=0.95)) %>%
  mutate(model='hook saturation level') %>%
  filter(!is.na(Estimate2))

#hook_sat_df2$Species <- 'Proportion of Fishing Events with >95% Baits Removed'
hook_sat_df$region <- forcats::fct_recode(factor(hook_sat_df$region_INLA),
                                          `HS` = '1', `QCS` = '2',
                                          `WCHG` = '3', `WCVI` = '4')

# predict the linear trend of hook saturation
hook_sat_df <- hook_sat_df[!is.na(hook_sat_df$region),]
hook_sat_df$mean <- 
  predict(
    lm(data=hook_sat_df, 
       Estimate2 ~ region * year)
  )
hook_sat_df$q0.025 <- hook_sat_df$mean - 
  1.96* predict(
    lm(data=hook_sat_df, 
       Estimate2 ~ region * year),
    se=T
  )$se.fit
hook_sat_df$q0.975 <- hook_sat_df$mean + 
  1.96* predict(
    lm(data=hook_sat_df, 
       Estimate2 ~ region * year),
    se=T
  )$se.fit

halibut_pred <- rbind(brms_mod_halibut_pred,brms_mod_halibut_compfactor_pred,brms_mod_halibut_cens_pred,boot_halibut_pred)

# halibut_pred <- rbind(brms_mod_halibut_pred,brms_mod_halibut_compfactor_pred,brms_mod_halibut_cens_pred)
#   
halibut_pred <- halibut_pred[halibut_pred$region!='All',]

#halibut_pred$year <- halibut_pred$year + 1995
# remove the predictions from WCVI pre 1999
halibut_pred[halibut_pred$region=='WCVI' & 
               halibut_pred$year < 1999, c('mean','Q2.5','Q97.5')] <- NA 

halibut_pred <- 
  full_join(halibut_pred,
            hook_sat_df[,c("year", "mean", "Estimate2", "model", "region","q0.025",'q0.975')], 
            suffix=c('',''))

halibut_pred$region <-
  forcats::fct_recode(halibut_pred$region,
                      `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
                      `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')

####

halibut_pred %>%
  filter(!is.na(region)) %>%
  filter(year > 1997) %>%
  ggplot(aes(x=year, y=mean, ymax=q0.975, ymin=q0.025,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Pacific Halibut Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Baits Removed') +
  geom_hline(yintercept = 0)

halibut_pred %>%
  filter(!is.na(region)) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q97.5, ymin=Q2.5,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Pacific Halibut Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Baits Removed') +
  geom_hline(yintercept = 0)

# What could be causing the lack of difference in censored Poisson indices?
by(INLA_dat$N_it_halibut, 
   INLA_dat$prop_removed>=0.95, 
   quantile, probs=c(0.95,0.96,0.97,0.98,0.99,1))
# The largest catch counts correspond to censored observations! We do not observe upper tail!

# Observe an increasing upper tail with time
by(INLA_dat$N_it_halibut[INLA_dat$prop_removed>=0.95],
   INLA_dat$year[INLA_dat$prop_removed>=0.95], quantile, probs=c(0.95,0.96,0.97,0.98,0.99,1))

library(MASS)
INLA_dat$censored_txt_95 <-
  ifelse(INLA_dat$prop_removed>=0.95, 'right','none')
INLA_dat %>%
  group_by(censored_txt_95, year, region) %>%
  mutate(region=
           fct_recode(factor(region),
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
  )%>%
  rename(censored=censored_txt_95) %>%
  summarise(Q80 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.8),
            Q98 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.98),
            Q99 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.99),
            QMax = quantile(N_it_yelloweye/effSkateIPHC, probs=1)) %>%
  pivot_longer(cols = c('Q80','Q98','Q99','QMax'), names_to = 'Quantile') %>%
  ggplot(aes(x=year, y=value, group=censored, colour=censored, fill=censored)) +
  geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=F) +
  facet_wrap(~Quantile+region, scales='free') +
  ylab('Catch Counts') +
  ggtitle('Pacific Halibut Catch Count Percentiles vs. Year', 
          subtitle = 'The 80th, 98th, 99th, and 100th percentiles are shown across the rows \nthe regions are shown as columns')

# There is evidence that we HAVE observed the upper tail for all regions except WCHG
# Perhaps this explains the increase in the censored Poisson index at the end

brms_dat %>%
  group_by(censored_txt_95, year, region) %>%
  mutate(region=
           fct_recode(factor(region),
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
  )%>%
  rename(censored=censored_txt_95) %>%
  mutate(Q80 = quantile(N_it_halibut/effSkateIPHC, probs=0.8),
         Q98 = quantile(N_it_halibut/effSkateIPHC, probs=0.98),
         Q99 = quantile(N_it_halibut/effSkateIPHC, probs=0.99),
         QMax = quantile(N_it_halibut/effSkateIPHC, probs=1)) %>%
  filter(N_it_halibut/effSkateIPHC >= QMax) %>%
  ggplot(aes(x=year+1995, y=N_it_halibut/effSkateIPHC, group=censored, colour=censored, fill=censored)) +
  # geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
  #                                                              data,
  #                                                              weights=weight,
  #                                                              method="MM"),
  #             fullrange=F) +
  geom_smooth(method='lm', fullrange=T) +
  facet_wrap(~region, scales='free') +
  geom_point() +
  ylab('CPUE') +
  ggtitle('The Maximum Pacific Halibut CPUE Value vs. Year', 
          subtitle = 'The maximum values are shown for each region and each censorship status')


cor_results <- 
  halibut_pred %>% 
  dplyr::select(region, year, Estimate, model) %>%
  pivot_wider(names_from=model, values_from = c(Estimate)) %>% 
  group_by(region) %>% 
  summarise(Cor_CPUE_ICR = cor(`CPUE-based`,
                               `Instantaneous catch-rate adjusted`,
                               method = 'spearman', use='complete.obs'),
            Cor_CPUE_Cen = cor(`CPUE-based`,
                               `Censored Poisson`,
                               method = 'spearman', use='complete.obs'),
            Cor_Cen_ICR = cor(`Censored Poisson`,
                              `Instantaneous catch-rate adjusted`,
                              method = 'spearman', use='complete.obs'),
            Cor_HookSat_CPUE = cor(`CPUE-based`,
                                   `hook saturation level`,
                                   method = 'spearman', use='complete.obs'),
            Cor_HookSat_ICR = cor(`hook saturation level`,
                                  `Instantaneous catch-rate adjusted`,
                                  method = 'spearman', use='complete.obs'),
            Cor_HookSat_Cen = cor(`Censored Poisson`,
                                  `hook saturation level`,
                                  method = 'spearman', use='complete.obs'))
cor_results

#### Repeat for yelloweye rockfish

brms_dat$censored_txt_975 <- rep('none',dim(brms_dat)[1])
brms_dat$censored_txt_975[which(brms_dat$prop_removed>=0.975)] <- 'right'

if(fit_mods)
{
  brms_mod_yelloweye <-
    brm(N_it_yelloweye | rate(effSkateIPHC) ~
          region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8) , 
        data = brms_dat,
        family='Poisson',
        iter=4000, warmup = 3000,
        cores = 4, chains = 4,
        control = list(adapt_delta = 0.99, max_treedepth = 15),
        seed=8102021,
        prior = c(prior(student_t(7, 0, 2), class = sd, group = event_ID, coef = Intercept)))
}
if(!fit_mods)
{
  brms_mod_yelloweye <- readRDS('yelloweye_mod_final.rds')
}
summary(brms_mod_yelloweye)

# Manually predict on the linear predictor scale and subtract the sampled intercept terms
newdata <- inlabru:::cprod(
  data.frame(effSkateIPHC=1, region=factor(c('1','2','3','4')), twentyhooks=0, station_ID='2004', event_ID='1'),
  data.frame(year=min(brms_dat$year):max(brms_dat$year))
)

brms_mod_yelloweye_pred <- brms_trend_plotter(brms_mod_yelloweye, newdata = newdata)
brms_mod_yelloweye_pred$model = 'CPUE-based'

# Now fit the hook competition scale factor-adjusted model
brms_dat$N_it_yelloweye_compfactor <-
  round(brms_dat$N_it_yelloweye * 
          comp_factor_fun(brms_dat$prop_removed, brms_dat$obsHooksPerSet))

if(fit_mods)
{
  brms_mod_yelloweye_compfactor <-
    brm(N_it_yelloweye_compfactor | rate(effSkateIPHC) ~
          region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8) , 
        data = brms_dat,
        family='Poisson',
        iter=4000, warmup = 3000,
        cores = 4, chains = 4,
        control = list(adapt_delta = 0.99, max_treedepth = 15),
        seed=8102021,
        prior = c(prior(student_t(7, 0, 2), class = sd, group = event_ID, coef = Intercept)))
}
if(!fit_mods)
{
  brms_mod_yelloweye_compfactor <- readRDS('yelloweye_compfactor_mod_final.rds')
}
summary(brms_mod_yelloweye_compfactor)

brms_mod_yelloweye_compfactor_pred <- brms_trend_plotter(brms_mod_yelloweye_compfactor, newdata = newdata)
brms_mod_yelloweye_compfactor_pred$model = 'Instantaneous catch-rate adjusted'

# Now fit the censored model
if(fit_mods)
{
  brms_mod_yelloweye_cens <-
    brm(N_it_yelloweye | cens(censored_txt_95) + rate(effSkateIPHC) ~
          region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8) , 
        data = brms_dat,
        family='Poisson',
        iter=4000, warmup = 3000,
        cores = 4, chains = 4,
        control = list(adapt_delta = 0.99, max_treedepth = 15),
        seed=8102021,
        prior = c(prior(student_t(7, 0, 2), class = sd, group = event_ID, coef = Intercept)))
}
if(!fit_mods)
{
  brms_mod_yelloweye_cens <- readRDS('yelloweye_cens_mod_final.rds')
}
summary(brms_mod_yelloweye_cens)

brms_mod_yelloweye_cens_pred <- brms_trend_plotter(brms_mod_yelloweye_cens, newdata = newdata)
brms_mod_yelloweye_cens_pred$model = 'Censored Poisson'

## New
# Now fit the censored model with the upper bound in place

# Define the upper bound
# upper_bound <- rep(0, length(brms_dat$prop_removed))
# 
# scale_fac <- rep(0, length(brms_dat$prop_removed))
# scale_fac[brms_dat$prop_removed>0.95] <- 
#     comp_factor_fun(signif((brms_dat[brms_dat$prop_removed>0.95,]$prop_removed-0.95)/(1-0.95),5),
#                     round((1-0.95)*brms_dat$obsHooksPerSet[brms_dat$prop_removed>0.95]))
#   
# upper_bound[brms_dat$prop_removed>0.95] <- round(
#     (brms_dat$prop_removed[brms_dat$prop_removed>0.95]-0.95)*brms_dat$obsHooksPerSet[brms_dat$prop_removed>0.95]*
#       scale_fac[brms_dat$prop_removed>0.95])

brms_dat$N_it_yelloweye_cens <- brms_dat$N_it_yelloweye

# brms_dat$upper_yelloweye <- 
#   brms_dat$N_it_yelloweye + upper_bound

brms_dat$int_censored_txt_95 <- brms_dat$censored_txt_95
# define interval censorship variable
# brms censorship interval is open on left - cannot handle 0's so need left
# brms_dat$int_censored_txt_95[brms_dat$N_it_yelloweye == 0 &
#                              brms_dat$int_censored_txt_95 == 'right'] <- 'left'
# 
# brms_dat$int_censored_txt_95[brms_dat$N_it_yelloweye > 0 &
#                                brms_dat$int_censored_txt_95 == 'right'] <- 'interval'
# 
# brms_dat$N_it_yelloweye_cens[brms_dat$int_censored_txt_95 == 'interval'] <- 
#   brms_dat$N_it_yelloweye[brms_dat$int_censored_txt_95 == 'interval'] -1 #remember open left interval
# define as upper bound
# brms_dat$N_it_yelloweye_cens[brms_dat$int_censored_txt_95 == 'left'] <-
#   brms_dat$upper_yelloweye[brms_dat$int_censored_txt_95 == 'left']

# Remove censorship from upper 2% of values per region
year_map_fun <- function(ind)
{
  return(pmatch(ind, sort(unique(brms_dat$year)), duplicates.ok = T))
}

quant_regions <- 
  matrix(by(brms_dat$N_it_yelloweye, 
            list(brms_dat$region_INLA,brms_dat$year),
            FUN = function(x){quantile(x,0.8, na.rm=T)}, simplify = T), nrow=4, ncol=length(unique(brms_dat$year)))

# Remove censorship
brms_dat$int_censored_txt_95[
  brms_dat$N_it_yelloweye >= quant_regions[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))]
] <- 'none'

# revert value back to original
# brms_dat$N_it_yelloweye_cens[
#   brms_dat$N_it_yelloweye >= quant_regions[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))]
# ] <- 
#   brms_dat$N_it_yelloweye[
#   brms_dat$N_it_yelloweye >= quant_regions[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))]
# ]

if(fit_mods)
{
  brms_mod_yelloweye_cens_upperbound <-
    brm(N_it_yelloweye_cens | cens(censored_txt_95) + rate(effSkateIPHC) ~
          region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8) , 
        data = brms_dat,
        family='Poisson',
        iter=4000, warmup = 3000,
        cores = 4, chains = 4,
        control = list(adapt_delta = 0.99, max_treedepth = 15),
        #24092021 seed=8102021,
        prior = c(prior(student_t(7, 0, 2), class = sd, group = event_ID, coef = Intercept)))
}
if(!fit_mods)
{
  brms_mod_yelloweye_cens_80 <- readRDS('yelloweye_cens_80_mod_final.rds')
}
summary(brms_mod_yelloweye_cens_80)

brms_mod_yelloweye_cens_80_pred <- brms_trend_plotter(brms_mod_yelloweye_cens_80, newdata = newdata)
brms_mod_yelloweye_cens_80_pred$model = 'Censored Poisson 80'


boot_yelloweye_pred <- brms_mod_yelloweye_pred
# recode factor to match brm mod output
combined_indices_yelloweye_INLA_Alldata$combined_ts_all$region <-
  forcats::fct_recode(combined_indices_yelloweye_INLA_Alldata$combined_ts_all$region,
                      `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
                      `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')
boot_yelloweye_pred$region <-
  forcats::fct_recode(boot_yelloweye_pred$region,
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
# calibrate the year variable
combined_indices_yelloweye_INLA_Alldata$combined_ts_all$year <- 
  combined_indices_yelloweye_INLA_Alldata$combined_ts_all$year - 1995
# match by year and region
boot_yelloweye_pred <- inner_join(
  boot_yelloweye_pred, combined_indices_yelloweye_INLA_Alldata$combined_ts_all[
    combined_indices_yelloweye_INLA_Alldata$combined_ts_all$model=='boot',
  ],
  by = c('region','year')
)
boot_yelloweye_pred$Estimate <- boot_yelloweye_pred$mean
boot_yelloweye_pred$Q2.5 <- NA#boot_yelloweye_pred$LCL
boot_yelloweye_pred$Q97.5 <- NA#boot_yelloweye_pred$UCL
boot_yelloweye_pred$model <- boot_yelloweye_pred$model.y
boot_yelloweye_pred$Est.Error <- NA

# keep only the relevant variables
boot_yelloweye_pred <- boot_yelloweye_pred[,names(brms_mod_yelloweye_pred)]

yelloweye_pred <- rbind(brms_mod_yelloweye_pred,brms_mod_yelloweye_compfactor_pred,brms_mod_yelloweye_cens_pred,boot_yelloweye_pred)

yelloweye_pred$region <-
  forcats::fct_recode(yelloweye_pred$region,
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
yelloweye_pred$year <- yelloweye_pred$year + 1995
# remove the predictions from WCVI pre 1999
yelloweye_pred[yelloweye_pred$region=='West Coast Van Island' & 
                 yelloweye_pred$year < 1999, c('Estimate','Q2.5','Q97.5')] <- NA 

yelloweye_pred <- 
  full_join(yelloweye_pred,
            hook_sat_df[,c("year", "Estimate", "Estimate2", "model", "region","Q2.5",'Q97.5')], 
            suffix=c('',''))

# Are these results consistent with INLA?
INLA_mod_yelloweye <- 
  final_censored_index_fun(data=Reduced_data_sp[which(!(is.na(Reduced_data_sp$region_INLA) |
                                                          is.na(Reduced_data_sp$prop_removed))),],
                           species = 'yelloweye',
                           cprop=0.95, 
                           use_upper_bound = T, 
                           upper_bound_quantile2 = 0.98)

# rename variables to make it ready for merging with brms indices
INLA_yelloweyepred <-
  INLA_mod_yelloweye$pred_overdisp %>%
  rename(Q2.5=q0.025,
         Q97.5=q0.975,
         Estimate = mean) %>%
  mutate(model='Adjusted Censored Poisson') %>%
  filter(region != 'All')

INLA_yelloweyepred$region <-
  forcats::fct_recode(INLA_yelloweyepred$region,
                      `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
                      `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')

yelloweye_pred <-
  full_join(yelloweye_pred,INLA_yelloweyepred)

yelloweye_pred$model <- 
  factor(yelloweye_pred$model,
         levels=c('hook saturation level',"boot","CPUE-based","Instantaneous catch-rate adjusted", "Adjusted Censored Poisson", "Censored Poisson"),
         ordered = T)

yelloweye_pred %>%
  filter(!is.na(region)) %>%
  filter(year > 1997) %>%
  filter(!is.na(model) & !(model %in% c('hook saturation level','boot', 'CPUE-based'))) %>%#,'Instantaneous catch-rate adjusted') )) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q97.5, ymin=Q2.5,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Yelloweye Rockfish Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Baits Removed') +
  geom_hline(yintercept = 0)

yelloweye_pred %>%
  filter(!is.na(region) & model!='Adjusted Censored Poisson') %>%
  #filter(year > 1997) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q97.5, ymin=Q2.5,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Yelloweye Rockfish Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Baits Removed') +
  geom_hline(yintercept = 0)

# What could be causing the difference in censored Poisson indices?
by(brms_dat$N_it_yelloweye, 
   brms_dat$censored_txt_95=='right', 
   quantile, probs=c(0.95,0.96,0.97,0.98,0.99,1))
# The largest catch counts correspond to censored observations! We do not observe upper tail!

# Observe an increasing upper tail with time
by(brms_dat$N_it_yelloweye[brms_dat$censored_txt_95=='right'],
   brms_dat$year[brms_dat$censored_txt_95=='right'], quantile, probs=c(0.95,0.96,0.97,0.98,0.99,1))

library(MASS)
brms_dat %>%
  group_by(censored_txt_95, year, region) %>%
  mutate(region=
           fct_recode(factor(region),
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
  )%>%
  rename(censored=censored_txt_95) %>%
  summarise(Q80 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.8),
            Q98 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.98),
            Q99 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.99),
            QMax = quantile(N_it_yelloweye/effSkateIPHC, probs=1)) %>%
  pivot_longer(cols = c('Q80','Q98','Q99','QMax'), names_to = 'Quantile') %>%
  ggplot(aes(x=year+1995, y=value, group=censored, colour=censored, fill=censored)) +
  geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=F) +
  facet_wrap(~Quantile+region, scales='free') +
  ylab('Catch Counts') +
  ggtitle('Yelloweye Rockfish Catch Count Percentiles vs. Year', 
          subtitle = 'The 80th, 98th, 99th, and 100th percentiles are shown across the rows \nthe regions are shown as columns')


# The upper 2% of values observed under censored fishing events are increasing with time
# The maximum goes doubles in WCHG and over quadruples in region 2. No change in region 4 or region 1
# Could this explain the dramatic change seen in regions 3 and 2?

brms_dat %>%
  group_by(censored_txt_95, year, region) %>%
  mutate(region=
           fct_recode(factor(region),
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
  )%>%
  rename(censored=censored_txt_95) %>%
  mutate(Q80 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.8),
         Q98 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.98),
         Q99 = quantile(N_it_yelloweye/effSkateIPHC, probs=0.99),
         QMax = quantile(N_it_yelloweye/effSkateIPHC, probs=1)) %>%
  filter(N_it_yelloweye/effSkateIPHC >= QMax) %>%
  ggplot(aes(x=year+1995, y=N_it_yelloweye/effSkateIPHC, group=censored, colour=censored, fill=censored)) +
  # geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
  #                                                              data,
  #                                                              weights=weight,
  #                                                              method="MM"),
  #             fullrange=F) +
  geom_smooth(method='lm', fullrange=T) +
  facet_wrap(~region, scales='free') +
  geom_point() +
  ylab('CPUE') +
  ggtitle('The Maximum Yelloweye Rockfish CPUE Value vs. Year', 
          subtitle = 'The maximum values are shown for each region and each censorship status')

brms_dat %>%
  group_by(censored_txt_95, year, region) %>%
  mutate(region=
           fct_recode(factor(region),
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
  )%>%
  rename(censored=censored_txt_95) %>%
  ggplot(aes(x=year+1995, y=N_it_yelloweye/effSkateIPHC, group=censored, colour=censored, fill=censored)) +
  geom_smooth(fullrange=T) +
  facet_wrap(~region, scales='free') +
  #geom_point() +
  ylab('CPUE') +
  ggtitle('Mean Yelloweye Rockfish CPUE Values vs. Year', 
          subtitle = 'Values are shown for each region and each censorship status')

# What happens if we perform only a complete case analysis?
# Remove all censored zeros and consider nonzero censored observations as observed
CC_dat <- Reduced_data_sp[which(!(is.na(Reduced_data_sp$region_INLA) |
                                    is.na(Reduced_data_sp$prop_removed))),] #|
#(Reduced_data_sp$prop_removed >= 0.95 &
#   Reduced_data_sp$N_it_yelloweye==0) )),]


INLA_mod_yelloweye_CC <- 
  final_censored_index_fun(data=CC_dat,
                           species = 'yelloweye',
                           cprop=0.95, 
                           use_upper_bound = F, 
                           upper_bound_quantile2 = 1)

# rename variables to make it ready for merging with brms indices
INLA_yelloweyepred_CC <-
  INLA_mod_yelloweye_CC$pred_overdisp %>%
  rename(Q2.5=q0.025,
         Q97.5=q0.975,
         Estimate = mean) %>%
  mutate(model='CC Censored Poisson') %>%
  filter(region != 'All')

INLA_yelloweyepred_CC$region <-
  forcats::fct_recode(INLA_yelloweyepred_CC$region,
                      `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
                      `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')


cor_results_yelloweye <- 
  yelloweye_pred %>% 
  dplyr::select(region, year, Estimate, model) %>%
  pivot_wider(names_from=model, values_from = c(Estimate)) %>% 
  group_by(region) %>% 
  summarise(Cor_CPUE_ICR = cor(`CPUE-based`,
                               `Instantaneous catch-rate adjusted`,
                               method = 'spearman', use='complete.obs'),
            Cor_CPUE_Cen = cor(`CPUE-based`,
                               `Censored Poisson`,
                               method = 'spearman', use='complete.obs'),
            Cor_Cen_ICR = cor(`Censored Poisson`,
                              `Instantaneous catch-rate adjusted`,
                              method = 'spearman', use='complete.obs'),
            Cor_HookSat_CPUE = cor(`CPUE-based`,
                                   `hook saturation level`,
                                   method = 'spearman', use='complete.obs'),
            Cor_HookSat_ICR = cor(`hook saturation level`,
                                  `Instantaneous catch-rate adjusted`,
                                  method = 'spearman', use='complete.obs'),
            Cor_HookSat_Cen = cor(`Censored Poisson`,
                                  `hook saturation level`,
                                  method = 'spearman', use='complete.obs'))
cor_results_yelloweye

Reduced_data_sp@data %>% 
  filter(!is.na(region_INLA)) %>% 
  group_by(year, region_INLA) %>% 
  summarise(prop_saturation = mean(prop_removed >= 0.95, na.rm=T)) %>% 
  ggplot(aes(x=year, y=prop_saturation, group=region_INLA)) +  
  geom_line() +
  geom_smooth(method = 'lm') + 
  facet_wrap(~region_INLA)
# Notice that the 2 regions (QCS and WCHG) with the negative correlations experience the biggest
# increases in >=95% hook saturation events across time. 
# Hook saturation does not increase in the other two regions  

#save.image('Paper_Case_Study.RData')

###
hook_sat_df <-
  Reduced_data_sp@data %>%
  group_by(year, region_INLA) %>%
  summarise(prop_saturation = mean(prop_removed[!is.nan(prop_removed)], na.rm=T)) %>%
  rename(year=year) %>%
  slice(rep(1:n(), each=length(unique(yelloweye_pred$model)))) %>%
  mutate(model=rep(unique(yelloweye_pred$model),times=n()/length(unique(yelloweye_pred$model))))

#hook_sat_df$Species <- 'Proportion of Baits Removed'
hook_sat_df$region <- forcats::fct_recode(factor(hook_sat_df$region_INLA),
                                          `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                          `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
#hook_sat_df$Q2.5 <- NA
#hook_sat_df$Q97.5 <- NA

hook_sat_df2 <-
  Reduced_data_sp@data %>%
  group_by(year, region_INLA) %>%
  summarise(prop_saturation = mean(prop_removed[!is.nan(prop_removed)]>=0.95)) %>%
  rename(year=year) %>%
  slice(rep(1:n(), each=length(unique(yelloweye_pred$model)))) %>%
  mutate(Model=rep(unique(yelloweye_pred$model),times=n()/length(unique(yelloweye_pred$model)))) %>%
  filter(!is.na(prop_saturation))

#hook_sat_df2$Species <- 'Proportion of Fishing Events with >95% Baits Removed'
hook_sat_df2$region <- forcats::fct_recode(factor(hook_sat_df2$region_INLA),
                                           `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                           `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
#hook_sat_df2$Q2.5 <- NA
#hook_sat_df2$Q97.5 <- NA

left_join(yelloweye_pred,
          hook_sat_df, suffix=c('','')) %>%
  filter(year < 2020) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q2.5, ymin=Q97.5, group=model, colour=model)) +
  geom_line() +
  facet_grid(model~region)

yelloweye_pred$species = 'yelloweye rockfish'
halibut_pred$species = 'pacific halibut'

Estimator_Properties <-
  rbind(yelloweye_pred, halibut_pred) %>%
  filter(year < 2020 & year!=2002) %>%
  group_by(model, species) %>%
  summarise(Range = diff(range(Estimate, na.rm=T)), 
            Uncertainty = median(Q97.5-Q2.5, na.rm=T))
Estimator_Properties$Uncertainty[Estimator_Properties$model=='boot' & 
                                   Estimator_Properties$species=='yelloweye rockfish'] <- 
  median(apply(combined_indices_yelloweye_INLA_Alldata$combined_ts_all[
    combined_indices_yelloweye_INLA_Alldata$combined_ts_all$model=='boot',
    c('LCL','UCL')],1,diff), na.rm=T)
Estimator_Properties$Uncertainty[Estimator_Properties$model=='boot' & 
                                   Estimator_Properties$species=='pacific halibut'] <- 
  median(apply(combined_indices_halibut_INLA_Alldata$combined_ts_all[
    combined_indices_halibut_INLA_Alldata$combined_ts_all$model=='boot',
    c('LCL','UCL')],1,diff), na.rm=T)

Estimator_Properties$Uncertainty[Estimator_Properties$species=='yelloweye rockfish']#/Estimator_Properties$Uncertainty[Estimator_Properties$species=='yelloweye rockfish'][4]
Estimator_Properties$Uncertainty[Estimator_Properties$species=='pacific halibut']#/Estimator_Properties$Uncertainty[Estimator_Properties$species=='pacific halibut'][4]
# Credible interval width is only 63% and 80% as wide for yelloweye and 3% and 24% wider for pacific halibut
# Above results comparing CPUE-based and instantaneous catch-rate adjusted indices of relative abundance
Estimator_Properties$Range[Estimator_Properties$species=='yelloweye rockfish']#/Estimator_Properties$Range[Estimator_Properties$Species=='yelloweye rockfish'][4]
Estimator_Properties$Range[Estimator_Properties$species=='pacific halibut']#/Estimator_Properties$Range[Estimator_Properties$Species=='pacific halibut'][4]
# Range of index values are only 77% and 89% as large for yelloweye and only 92% as high and 59% higher for pacific halibut
# Above results comparing CPUE-based and instantaneous catch-rate adjusted indices of relative abundance

left_join(yelloweye_pred,
          hook_sat_df, suffix=c('','')) %>%
  filter(year < 2020) %>%
  group_by(region, model) %>%
  ggplot(aes(x=prop_saturation, y=Estimate, colour=region)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(region~model) 

library(MASS)
left_join(yelloweye_pred,
          hook_sat_df2, suffix=c('','')) %>%
  filter(year < 2020) %>%
  group_by(region, model) %>%
  ggplot(aes(x=prop_saturation, y=Estimate, colour=region)) +
  geom_point() +
  geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=F) +
  facet_grid(region~model, scales = 'free_y') +
  xlab('Annual Proportion of Fishing Events Returning <5% baits') +
  ylab('Relative Abundance') +
  ggtitle('Relative Yelloweye Rockfish Abundance vs. Degree of Hook Saturation',
          subtitle = 'Plots shown separately for each region and method')

left_join(halibut_pred,
          hook_sat_df2, suffix=c('','')) %>%
  filter(year < 2020) %>%
  group_by(region, model) %>%
  ggplot(aes(x=prop_saturation, y=Estimate, colour=region)) +
  geom_point() +
  geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=F) +
  facet_grid(region~model, scales = 'free_y')  +
  xlab('Annual Proportion of Fishing Events Returning <5% baits') +
  ylab('Relative Abundance') +
  ggtitle('Relative Pacific Halibut Abundance vs. Degree of Hook Saturation',
          subtitle = 'Plots shown separately for each region and method')

# Plot the indices vs the proportion of saturation events > 95%
left_join(yelloweye_pred,
          hook_sat_df2, suffix=c('','')) %>%
  filter(year < 2020) %>%
  group_by(region, model) %>%
  mutate(Diff_Index = c(NA,diff(Estimate)-mean(diff(Estimate),na.rm=T)),
         Diff_Sat = c(NA,diff(prop_saturation)-mean(diff(prop_saturation),na.rm=T))) %>%
  ungroup() %>%
  ggplot(aes(y=Diff_Index, x=Diff_Sat, group=model, colour=model)) +
  geom_point() +
  geom_smooth(method='gam',
              fullrange=F, formula=y~s(x, bs="cs", k=4)) +
  # stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
  #                                                              data,
  #                                                              weights=weight,
  #                                                              method="MM"),
  #             fullrange=TRUE) +
  facet_grid(model~region, scales = 'free') +
  #geom_hline(yintercept = 0) +
  ylab('Detrended Change in Annual Normalized Relative Abundance Index') + 
  xlab('Change in Annual Proportion of Fishing Events Returning <5% baits') +
  guides(colour='none')

left_join(yelloweye_pred,
          hook_sat_df, suffix=c('','')) %>%
  filter(year < 2020) %>%
  mutate(delta_saturation = c(NA,diff(prop_saturation)),
         delta_abundance = c(NA, diff(Estimate))) %>%
  group_by(region, model) %>%
  ggplot(aes(x=delta_saturation, y=delta_abundance, colour=region)) +
  geom_point() +
  geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=F) +
  facet_grid(region~model) 

left_join(halibut_pred,
          hook_sat_df, suffix=c('','')) %>%
  filter(year < 2020) %>%
  mutate(delta_saturation = c(NA,diff(prop_saturation)),
         delta_abundance = c(NA, diff(Estimate))) %>%
  group_by(region, model) %>%
  ggplot(aes(x=delta_saturation, y=delta_abundance, colour=region)) +
  geom_point() +
  geom_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=F) +
  facet_grid(region~model) 

### Posterior Predictive Checks
library(bayesplot)
bayesplot::ppc_intervals_grouped(y=brms_mod_halibut$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'], 
                                 yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 group = brms_mod_halibut$data$region[brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 x = brms_mod_halibut$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none'])
bayesplot::ppc_intervals_grouped(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'], 
                                 yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 group = brms_mod_halibut_cens$data$region[brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 x = brms_mod_halibut_cens$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none'])

bayesplot::ppc_ecdf_overlay(y=brms_mod_halibut$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'], 
                            yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_mod_halibut_cens$data$censored_txt_95=='none'])
bayesplot::ppc_ecdf_overlay(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'], 
                            yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none'])

bayesplot::ppc_scatter_avg_grouped(y=brms_mod_halibut$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'], 
                                   yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                   group=forcats::fct_recode(brms_mod_halibut$data$region,
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')[brms_mod_halibut_cens$data$censored_txt_95=='none']) + 
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue')
bayesplot::ppc_scatter_avg_grouped(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'], 
                                   yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                   group=forcats::fct_recode(brms_mod_halibut_cens$data$region,
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')[brms_mod_halibut_cens$data$censored_txt_95=='none']) + 
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue')

# Again for yelloweye
bayesplot::ppc_intervals_grouped(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'], 
                                 yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 group = brms_mod_yelloweye$data$region[brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 x = brms_mod_yelloweye$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none'])
bayesplot::ppc_intervals_grouped(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'], 
                                 yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 group = brms_mod_yelloweye_cens$data$region[brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 x = brms_mod_yelloweye_cens$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none'])

library(loo)
log_lik <- log_lik(brms_mod_yelloweye_cens, merge_chains = FALSE)
rel_n_eff <- relative_eff(exp(log_lik), chain_id = rep(1:4, each=1000), cores = 10)
loo(log_lik[,which(brms_mod_yelloweye_cens$data$censored_txt_95=='none') ], r_eff = rel_n_eff[which(brms_mod_yelloweye_cens$data$censored_txt_95=='none')], cores = 10)
loo(log_lik, r_eff = rel_n_eff, cores = 10, save_psis = T)

bayesplot::ppc_loo_intervals(y=brms_mod_yelloweye_cens$data$N_it_yelloweye, 
                             yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens)),
                             psis_object = psis(log_lik, r_eff = rel_n_eff, cores = 10, save_psis = T),
                             subset = which(brms_mod_yelloweye_cens$data$censored_txt_95=='none'))

bayesplot::ppc_ecdf_overlay(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'], 
                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'])
bayesplot::ppc_ecdf_overlay(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'], 
                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'])

bayesplot::ppc_scatter_avg_grouped(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'], 
                                   yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                   group=forcats::fct_recode(brms_mod_yelloweye$data$region,
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')[brms_mod_yelloweye_cens$data$censored_txt_95=='none']) + 
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue')
bayesplot::ppc_scatter_avg_grouped(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'], 
                                   yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                   group=forcats::fct_recode(brms_mod_yelloweye_cens$data$region,
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')[brms_mod_yelloweye_cens$data$censored_txt_95=='none']) + 
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue')

plot_list_halibut <- list()
plot_list_halibut2 <- list()
plot_list_halibut4 <- list()
plot_list_yelloweye <- list()
plot_list_yelloweye2 <- list()
plot_list_yelloweye4 <- list()
set.seed(21092021)
samp_ind <- sample.int(n=4000, size=100)
yrange_halibut <- c(0,0)
yrange_halibut2 <- c(0,0)
xrange_halibut <- range(c(0,brms_dat$N_it_halibut[brms_dat$censored_txt_95=='right']))
yrange_yelloweye <- c(0,0)
yrange_yelloweye2 <- c(0,0)
xrange_yelloweye <- range(c(0,brms_dat$N_it_yelloweye[brms_dat$censored_txt_95=='right']))
for(j in 1:4)
{
  plot_list_halibut[[1+(j-1)*3]] <- 
    bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_halibut$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                 brms_mod_halibut_cens$data$region==j], 
                                          yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                                brms_mod_halibut_cens$data$region==j],
                                          x=1995+brms_mod_halibut_cens$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                   brms_mod_halibut_cens$data$region==j]) + 
    geom_smooth(colour='red', method='gam') +
    geom_hline(yintercept = 0) + ggtitle('Prediction Error CPUE-Based', subtitle = paste0('Restricted to fishing events with <95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + xlab('Year')
  
  plot_list_halibut[[2+(j-1)*3]] <- 
    bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_halibut_compfactor$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                            brms_mod_halibut_cens$data$region==j], 
                                          yrep=as.matrix(posterior_predict(brms_mod_halibut_compfactor))[,brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                                           brms_mod_halibut_cens$data$region==j],
                                          x=1995+brms_mod_halibut_cens$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                   brms_mod_halibut_cens$data$region==j]) + 
    geom_smooth(colour='red', method='gam') +
    geom_hline(yintercept = 0) + ggtitle('Prediction Error Comp-Factor Adjust', subtitle = paste0('Restricted to fishing events with <95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + xlab('Year')
  
  
  plot_list_halibut[[3+(j-1)*3]] <- 
    bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                      brms_mod_halibut_cens$data$region==j], 
                                          yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                                     brms_mod_halibut_cens$data$region==j],
                                          x=1995+brms_mod_halibut_cens$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                   brms_mod_halibut_cens$data$region==j]) + 
    geom_smooth(colour='red', method='gam') +
    geom_hline(yintercept = 0) + ggtitle('Prediction Error Censored', subtitle = paste0('Restricted to fishing events with <95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + xlab('Year')
  
  post_pred_df <-
    data.frame(Censored_pred = as.numeric(posterior_predict(brms_mod_halibut_cens,draw_ids = samp_ind)[,brms_mod_halibut_cens$data$censored_txt_95=='right' &
                                                                                                         brms_mod_halibut_cens$data$region==j]),
               CPUE_pred = as.numeric(posterior_predict(brms_mod_halibut,draw_ids = samp_ind)[,brms_mod_halibut_cens$data$censored_txt_95=='right' &
                                                                                                brms_mod_halibut_cens$data$region==j]),
               Adj_pred = as.numeric(posterior_predict(brms_mod_halibut_compfactor,draw_ids = samp_ind)[,brms_mod_halibut_cens$data$censored_txt_95=='right' &
                                                                                                          brms_mod_halibut_cens$data$region==j]),
               Observed = rep(as.numeric(brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='right' &
                                                                                   brms_mod_halibut_cens$data$region==j]),
                              each=100))
  
  plot_list_halibut2[[1+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=Censored_pred)) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range(post_pred_df)) + ggtitle('Censored Posterior Predictions vs. Observed Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed Halibut Catch Counts') + ylab('Posterior Predicted Halibut Counts')
  
  plot_list_halibut2[[2+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=Adj_pred)) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range(post_pred_df)) + ggtitle('Scale Factor Adjusted Posterior Predictions vs. Observed Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed Halibut Catch Counts') + ylab('Posterior Predicted Halibut Counts')
  
  # Now plot the variance
  plot_list_halibut4[[1+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=(Censored_pred-Observed)^2 )) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range((post_pred_df[,c(1,3)]-post_pred_df[,c(4)])^2)) + ggtitle('Censored Posterior Squared Prediction Difference vs. Observed Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed Halibut Catch Counts') + ylab('Posterior Squared Prediction Difference')
  
  plot_list_halibut4[[2+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=(Adj_pred-Observed)^2)) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range((post_pred_df[,c(1,3)]-post_pred_df[,c(4)])^2)) + ggtitle('Scale Factor Adjusted Posterior Squared Prediction Difference vs. Observed Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed Halibut Catch Counts') + ylab('Posterior Squared Prediction Difference')
  
  # Update y-axis scale for later plots
  yrange_halibut <- range(c(yrange_halibut,range(post_pred_df)))
  yrange_halibut2 <- range(c(yrange_halibut2,range((post_pred_df[,c(1,3)]-post_pred_df[,c(4)])^2)))
  
  # Repeat for yelloweye
  plot_list_yelloweye[[1+(j-1)*3]] <- 
    bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                     brms_mod_yelloweye_cens$data$region==j], 
                                          yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                                  brms_mod_yelloweye_cens$data$region==j],
                                          x=1995+brms_mod_yelloweye_cens$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                     brms_mod_yelloweye_cens$data$region==j]) + 
    geom_smooth(colour='red', method='gam') +
    geom_hline(yintercept = 0) + ggtitle('Prediction Error CPUE-Based', subtitle = paste0('Restricted to fishing events with <95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + xlab('Year')
  
  plot_list_yelloweye[[2+(j-1)*3]] <- 
    bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_yelloweye_compfactor$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                                brms_mod_yelloweye_cens$data$region==j], 
                                          yrep=as.matrix(posterior_predict(brms_mod_yelloweye_compfactor))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                                             brms_mod_yelloweye_cens$data$region==j],
                                          x=1995+brms_mod_yelloweye_cens$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                     brms_mod_yelloweye_cens$data$region==j]) + 
    geom_smooth(colour='red', method='gam') +
    geom_hline(yintercept = 0) + ggtitle('Prediction Error Comp-Factor Adjust', subtitle = paste0('Restricted to fishing events with <95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + xlab('Year')
  
  
  plot_list_yelloweye[[3+(j-1)*3]] <- 
    bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                          brms_mod_yelloweye_cens$data$region==j], 
                                          yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                                       brms_mod_yelloweye_cens$data$region==j],
                                          x=1995+brms_mod_yelloweye_cens$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                     brms_mod_yelloweye_cens$data$region==j]) + 
    geom_smooth(colour='red', method='gam') +
    geom_hline(yintercept = 0) + ggtitle('Prediction Error Censored', subtitle = paste0('Restricted to fishing events with <95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + xlab('Year')
  
  post_pred_df <-
    data.frame(Censored_pred = as.numeric(posterior_predict(brms_mod_yelloweye_cens,draw_ids = samp_ind)[,brms_mod_yelloweye_cens$data$censored_txt_95=='right' &
                                                                                                           brms_mod_yelloweye_cens$data$region==j]),
               CPUE_pred = as.numeric(posterior_predict(brms_mod_yelloweye,draw_ids = samp_ind)[,brms_mod_yelloweye_cens$data$censored_txt_95=='right' &
                                                                                                  brms_mod_yelloweye_cens$data$region==j]),
               Adj_pred = as.numeric(posterior_predict(brms_mod_yelloweye_compfactor,draw_ids = samp_ind)[,brms_mod_yelloweye_cens$data$censored_txt_95=='right' &
                                                                                                            brms_mod_yelloweye_cens$data$region==j]),
               Observed = rep(as.numeric(brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='right' &
                                                                                       brms_mod_yelloweye_cens$data$region==j]),
                              each=100))
  
  plot_list_yelloweye2[[1+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=Censored_pred)) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range(post_pred_df)) + ggtitle('Censored Posterior Predictions vs. Observed yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed yelloweye Catch Counts') + ylab('Posterior Predicted yelloweye Counts')
  
  plot_list_yelloweye2[[2+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=Adj_pred)) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range(post_pred_df)) + ggtitle('Scale Factor Adjusted Posterior Predictions vs. Observed yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed yelloweye Catch Counts') + ylab('Posterior Predicted yelloweye Counts')
  
  plot_list_yelloweye4[[1+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=(Censored_pred-Observed)^2 )) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range((post_pred_df[,c(1,3)]-post_pred_df[,c(4)])^2)) + ggtitle('Censored Posterior Squared Prediction Difference vs. Observed Yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed Yelloweye Catch Counts') + ylab('Posterior Squared Prediction Difference')
  
  plot_list_yelloweye4[[2+(j-1)*2]] <- 
    ggplot(data=post_pred_df,
           aes(x=Observed, y=(Adj_pred-Observed)^2)) +
    geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
    ylim(range((post_pred_df[,c(1,3)]-post_pred_df[,c(4)])^2)) + ggtitle('Scale Factor Adjusted Posterior Squared Prediction Difference vs. Observed Yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% hook saturation in ' ,c('Hectate Strait', 'Queen Charlotte Sound','West Coast Haida Gwaii', 'West Coast Van Island')[j])) + 
    xlab('Observed Yelloweye Catch Counts') + ylab('Posterior Squared Prediction Difference')
  
  # Update y-axis scale for later plots
  yrange_yelloweye <- range(c(yrange_yelloweye,range(post_pred_df)))
  yrange_yelloweye2 <- range(c(yrange_yelloweye2,range((post_pred_df[,c(1,3)]-post_pred_df[,c(4)])^2)))
  
}

plot_list_halibut3 <- plot_list_halibut2
plot_list_yelloweye3 <- plot_list_yelloweye2
plot_list_halibut5 <- plot_list_halibut4
plot_list_yelloweye5 <- plot_list_yelloweye4
for(i in 1:length(plot_list_halibut2))
{
  plot_list_halibut2[[i]] <- plot_list_halibut2[[i]] + ylim(yrange_halibut) + xlim(xrange_halibut)
  plot_list_yelloweye2[[i]] <- plot_list_yelloweye2[[i]] + ylim(yrange_yelloweye) + xlim(xrange_yelloweye)
  
  plot_list_halibut3[[i]] <- plot_list_halibut2[[i]] + coord_cartesian(ylim=c(0,1500))
  # remove the geom_point layer
  plot_list_halibut3[[i]]$layers[[1]] <- NULL
  plot_list_yelloweye3[[i]] <- plot_list_yelloweye2[[i]] + coord_cartesian(ylim=c(0,1400))
  # remove the geom_point layer
  plot_list_yelloweye3[[i]]$layers[[1]] <- NULL
  
  plot_list_halibut4[[i]] <- plot_list_halibut4[[i]] + ylim(yrange_halibut2) + xlim(xrange_halibut)
  plot_list_yelloweye4[[i]] <- plot_list_yelloweye4[[i]] + ylim(yrange_yelloweye2) + xlim(xrange_yelloweye)
  
  plot_list_halibut5[[i]] <- plot_list_halibut4[[i]] + coord_cartesian(ylim=c(0,rep(c(120000,1200000),times=4)[i]) )
  # remove the geom_point layer
  plot_list_halibut5[[i]]$layers[[1]] <- NULL
  plot_list_yelloweye5[[i]] <- plot_list_yelloweye4[[i]] + coord_cartesian(ylim=c(0,rep(c(1000000,1000000),times=4)[i]) )
  # remove the geom_point layer
  plot_list_yelloweye5[[i]]$layers[[1]] <- NULL
}

#library(inlabru)
inlabru::multiplot(plotlist = plot_list_halibut,
                   layout=matrix(c(1:12), nrow=4, ncol=3, byrow=T))
inlabru::multiplot(plotlist = plot_list_halibut2,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))
inlabru::multiplot(plotlist = plot_list_halibut3,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))
inlabru::multiplot(plotlist = plot_list_halibut4,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))
inlabru::multiplot(plotlist = plot_list_halibut5,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))

inlabru::multiplot(plotlist = plot_list_yelloweye,
                   layout=matrix(c(1:12), nrow=4, ncol=3, byrow=T))
inlabru::multiplot(plotlist = plot_list_yelloweye2,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))
inlabru::multiplot(plotlist = plot_list_yelloweye3,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))
inlabru::multiplot(plotlist = plot_list_yelloweye4,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))
inlabru::multiplot(plotlist = plot_list_yelloweye5,
                   layout=matrix(c(1:8), nrow=4, ncol=2, byrow=T))

# What could be causing the higher variance in the censored observations?
# Compute the posterior SD of the IID event random effects for both the censored and uncensored events
tmp <- apply(as_draws_matrix(brms_mod_halibut, variable='r_event_ID'),2, FUN = function(x){sd(x)})
tmp <- data.frame(SD_IID = tmp,
                  Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')),
                  Species='Halibut')

tmp2 <- apply(as_draws_matrix(brms_mod_halibut_cens, variable='r_event_ID'),2, FUN = function(x){sd(x)})
tmp2 <- data.frame(SD_IID = tmp2,
                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')),
                   Species='Halibut')

tmp3 <- apply(as_draws_matrix(brms_mod_halibut_compfactor, variable='r_event_ID'),2, FUN = function(x){sd(x)})
tmp3 <- data.frame(SD_IID = tmp3,
                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')),
                   Species='Halibut')

tmp4 <- apply(as_draws_matrix(brms_mod_yelloweye, variable='r_event_ID'),2, FUN = function(x){sd(x)})
tmp4 <- data.frame(SD_IID = tmp4,
                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')),
                   Species='Yelloweye')

tmp5 <- apply(as_draws_matrix(brms_mod_yelloweye_cens, variable='r_event_ID'),2, FUN = function(x){sd(x)})
tmp5 <- data.frame(SD_IID = tmp5,
                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')),
                   Species='Yelloweye')

tmp6 <- apply(as_draws_matrix(brms_mod_yelloweye_compfactor, variable='r_event_ID'),2, FUN = function(x){sd(x)})
tmp6 <- data.frame(SD_IID = tmp6,
                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')),
                   Species='Yelloweye')

multiplot(
  ggplot(data=tmp, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
    geom_density(alpha=0.2) + coord_cartesian(xlim=c(0.18,0.7)) +
    ggtitle('Density Curve of the SDs of the IID Per-Event Random Effects',
            subtitle = 'From the CPUE-based Halibut model and separated into the \nfishing events with hook saturation <95% and >=95%') ,
  ggplot(data=tmp2, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
    geom_density(alpha=0.2) + coord_cartesian(xlim=c(0.18,0.7)) + 
    ggtitle('Density Curve of the SDs of the IID Per-Event Random Effects',
            subtitle = 'From the Censored Halibut model and separated into the \nfishing events with hook saturation <95% and >=95%') ,
  ggplot(data=tmp3, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
    geom_density(alpha=0.2) + coord_cartesian(xlim=c(0.18,0.7)) + 
    ggtitle('Density Curve of the SDs of the IID Per-Event Random Effects',
            subtitle = 'From the Scale Factor-Adjusted Halibut model and separated into the \nfishing events with hook saturation <95% and >=95%'),
  ggplot(data=tmp4, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
    geom_density(alpha=0.2) + coord_cartesian(xlim=c(0.21,1.35)) + 
    ggtitle('Density Curve of the SDs of the IID Per-Event Random Effects',
            subtitle = 'From the CPUE-based Yelloweye model and separated into the \nfishing events with hook saturation <95% and >=95%') ,
  ggplot(data=tmp5, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
    geom_density(alpha=0.2) + coord_cartesian(xlim=c(0.21,1.35)) + 
    ggtitle('Density Curve of the SDs of the IID Per-Event Random Effects',
            subtitle = 'From the Censored Yelloweye model and separated into the \nfishing events with hook saturation <95% and >=95%') ,
  ggplot(data=tmp6, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
    geom_density(alpha=0.2) + coord_cartesian(xlim=c(0.21,1.35)) +
    ggtitle('Density Curve of the SDs of the IID Per-Event Random Effects',
            subtitle = 'From the Scale Factor-Adjusted Yelloweye model and separated into the \nfishing events with hook saturation <95% and >=95%'),
  layout=matrix(c(1:6), nrow=3, ncol=2, byrow=F))


# 
# 
# tmp <- apply(as_draws_matrix(brms_mod_yelloweye, variable='r_event_ID'),2, FUN = function(x){sd(x)})
# tmp <- data.frame(SD_IID = tmp,
#                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')))
# ggplot(data=tmp, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
#   geom_density(alpha=0.2) 
# 
# tmp <- apply(as_draws_matrix(brms_mod_yelloweye_cens, variable='r_event_ID'),2, FUN = function(x){sd(x)})
# tmp <- data.frame(SD_IID = tmp,
#                   Censored = factor(ifelse(brms_dat$censored_txt_95=='right','yes','no')))
# ggplot(data=tmp, aes(x=SD_IID, group=Censored, colour=Censored, fill=Censored)) + 
#   geom_density(alpha=0.2) 
# 
# 
# 
# tmp <- as_draws(brms_mod_yelloweye, variable='r_event_ID')


######
# Incorporate the upper bound and quantile adjustment
# brms_dat$censored_txt_95 <- rep('none',dim(brms_dat)[1])
# brms_dat$censored_txt_975 <- rep('none',dim(brms_dat)[1])
# 
# upper_bound_95 <- rep(0, length(brms_dat$prop_removed))
# upper_bound_975 <- rep(0, length(brms_dat$prop_removed))
# 
# scale_fac <- rep(0, length(brms_dat$prop_removed))
# scale_fac[brms_dat$prop_removed>0.95] <- 
#   comp_factor_fun(signif((brms_dat[brms_dat$prop_removed>0.95,]$prop_removed-0.95)/(1-0.95),5),
#                   round((1-0.95)*brms_dat$obsHooksPerSet[brms_dat$prop_removed>0.95]))
# 
# upper_bound_95[brms_dat$prop_removed>0.95] <- round(
#   (brms_dat$prop_removed[brms_dat$prop_removed>0.95]-0.95)*brms_dat$obsHooksPerSet[brms_dat$prop_removed>0.95]*
#     scale_fac[brms_dat$prop_removed>0.95])
# 
# upper_bound_975[brms_dat$prop_removed>0.975] <- round(
#   (brms_dat$prop_removed[brms_dat$prop_removed>0.975]-0.975)*brms_dat$obsHooksPerSet[brms_dat$prop_removed>0.975]*
#     scale_fac[brms_dat$prop_removed>0.975])
# 
# # Incorporate the quantile adjustment
# year_map_fun <- function(ind)
# {
#   return(pmatch(ind, sort(unique(brms_dat$year)), duplicates.ok = T))
# }
# quant_regions_halibut <- 
#   matrix(by(brms_dat$N_it_halibut, 
#             list(brms_dat$region_INLA,brms_dat$year),
#             FUN = function(x){quantile(x,0.8, na.rm=T)}, simplify = T), nrow=4, ncol=nyear)
# 
# quant_regions_yelloweye <- 
#   matrix(by(brms_dat$N_it_yelloweye, 
#             list(brms_dat$region_INLA,brms_dat$year),
#             FUN = function(x){quantile(x,0.8, na.rm=T)}, simplify = T), nrow=4, ncol=nyear)
# 
# brms_dat$high_halibut <- brms_dat$N_it_halibut
# brms_dat$high_halibut[which(brms_dat$prop_removed >= 0.95 & 
#                               0 < upper_bound_95 &
#                               brms_dat$N_it_halibut < quant_regions_halibut[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])] <- 
#   brms_dat$N_it_halibut[which(brms_dat$prop_removed >= 0.95 & 
#                                 0 < upper_bound_95 &
#                                 brms_dat$N_it_halibut < quant_regions_halibut[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])] +
#   upper_bound_95[which(brms_dat$prop_removed >= 0.95 & 
#                          0 < upper_bound_95 &
#                          brms_dat$N_it_halibut < quant_regions_halibut[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])]
# 
# brms_dat$censored_txt_95[which(brms_dat$prop_removed >= 0.95 & 
#                                  0 < upper_bound_95 &
#                                  brms_dat$N_it_halibut < quant_regions_halibut[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])] <- 'interval'
# 
# # Need to subtract 1 for censored interval observations when using brms 
# brms_dat$N_it_halibut[which(brms_dat$censored_txt_95 == 'interval')] <- pmax(0,
#                                                                              brms_dat$N_it_halibut[which(brms_dat$censored_txt_95 == 'interval')] - 1)
# 
# brms_dat$high_yelloweye <- brms_dat$N_it_yelloweye
# brms_dat$high_yelloweye[which(brms_dat$prop_removed >= 0.975 & 
#                                 0 < upper_bound_975 &
#                                 brms_dat$N_it_yelloweye < quant_regions_yelloweye[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])] <- 
#   brms_dat$N_it_yelloweye[which(brms_dat$prop_removed >= 0.975 & 
#                                   0 < upper_bound_975 &
#                                   brms_dat$N_it_yelloweye < quant_regions_yelloweye[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])] +
#   upper_bound_975[which(brms_dat$prop_removed >= 0.975 & 
#                           0 < upper_bound_975 &
#                           brms_dat$N_it_yelloweye < quant_regions_yelloweye[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])]
# 
# brms_dat$censored_txt_975[which(brms_dat$prop_removed >= 0.975 & 
#                                   0 < upper_bound_975 &
#                                   brms_dat$N_it_yelloweye < quant_regions_yelloweye[cbind(brms_dat$region_INLA,year_map_fun(brms_dat$year))])] <- 'interval'
# 
# brms_dat$N_it_yelloweye[which(brms_dat$censored_txt_975 == 'interval')] <- pmax(0,
#                                                                                 brms_dat$N_it_yelloweye[which(brms_dat$censored_txt_975 == 'interval')] - 1)
# 
# brms_mod_halibut2 <-
#   brm(N_it_halibut | cens(censored_txt_95, high_halibut) + rate(effSkateIPHC) ~ 
#         region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8), 
#       data = brms_dat,
#       family='Poisson',
#       cores=1)
# summary(brms_mod_halibut2)
# 
# brms_mod_yelloweye2 <-
#   brm(N_it_yelloweye | cens(censored_txt_95, high_yelloweye) + rate(effSkateIPHC) ~ 
#         region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8), 
#       data = brms_dat,
#       family='Poisson')
# summary(brms_mod_yelloweye2)
