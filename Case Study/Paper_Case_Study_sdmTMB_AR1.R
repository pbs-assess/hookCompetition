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
library(sdmTMB)

setwd("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final")

fit_mods <- T
#Again - we advise getting a pardiso license
#inla.setOption(pardiso.license="pardiso.lic")
#inla.pardiso.check()

# Load the bootstrap indices of abundance from Anderson et al 2019
#pacific_halibut_boot <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Andy Indices Scripts/all_species_results/all-species-results/pacific-halibut-results.rds")[[2]][[1]]
#yelloweye_rockfish_boot <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Andy Indices Scripts/all_species_results/all-species-results/yelloweye-rockfish-results.rds")[[2]][[1]]

# Load the censored Poisson indices
#combined_indices_halibut_INLA_Alldata <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final/combined_indices_halibut_INLA_Alldata.rds")
#combined_indices_yelloweye_INLA_Alldata <- readRDS("~/OneDrive - The University Of British Columbia/DFO IPHC Project 2021/Inlabru Approach/Deliverable_1_Final/combined_indices_yelloweye_INLA_Alldata.rds")

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

## Use sdmTMB to estimate trends
year_map_fun <- function(ind)
{
  return(pmatch(ind, sort(unique(sdmTMB_dat$year)), duplicates.ok = T))
}

# Create a list for storing the results
results_list <- vector('list', length=length(simple_species_vec))
names(results_list) <- data_species_vec
count <- 1
sdmTMB_dat <- Reduced_data_sp@data[which(!is.na(Reduced_data_sp$prop_removed)),]
#sdmTMB_dat <- sdmTMB_dat %>% filter(year > 1997)
sdmTMB_dat$year <- sdmTMB_dat$year - min(sdmTMB_dat$year) + 1
sdmTMB_dat$region <- factor(sdmTMB_dat$region_INLA)
sdmTMB_dat <- sdmTMB_dat[which(!is.na(sdmTMB_dat$region)),]
sdmTMB_dat$offset <- log(sdmTMB_dat$effSkateIPHC)
sdmTMB_dat$event_ID <- factor(1:length(sdmTMB_dat$year))
sdmTMB_dat$station_ID <- as.factor(sdmTMB_dat$station)
# Store a record of the type of censored Poisson index required
cens_poisson_type <- rep('right-censored no modification',
                         times=length(simple_species_vec))
# Loop through species and estimate relative abundance index each time
for(i in simple_species_vec)
{
  lower <- as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))
  upper <- as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))
  upper[which(sdmTMB_dat$prop_removed>=0.95)] <- NA

  # check no NA's
  #summary(sdmTMB_dat[,c('region','year','station_ID','event_ID')])

  mod_sdmTMB_poisson <-
    sdmTMB(
      data=sdmTMB_dat,
      formula=as.formula(paste0(paste0('N_it_',i),' ~ -1 + offset + (1|station_ID) + (1|event_ID)')),
      mesh=make_mesh(sdmTMB_dat,
                     xy_cols=c('lon','lat'),
                     n_knots=3),
      family=poisson(link = 'log'),
      #experimental = list(lwr=lower_halibut,upr=upper_halibut),
      time_varying = ~ -1 + region,
      time = 'year',
      spatial = 'off',
      spatiotemporal = 'off'
    )
  #summary(mod_sdmTMB_poisson)

  ## Now fit the scale-factor adjusted model
  # Multiply by scale factor
  sdmTMB_dat$N_it_ICR <-
    round(as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) *
            comp_factor_fun(sdmTMB_dat$prop_removed, sdmTMB_dat$obsHooksPerSet))

  mod_sdmTMB_ICR <-
    sdmTMB(
      data=sdmTMB_dat,
      formula=N_it_ICR ~ -1 + offset + (1|station_ID) + (1|event_ID),
      mesh=make_mesh(sdmTMB_dat,
                     xy_cols=c('lon','lat'),
                     n_knots=3),
      family=poisson(link = 'log'),
      #experimental = list(lwr=lower_halibut,upr=upper_halibut),
      time_varying = ~ -1 + region,
      time = 'year',
      spatial = 'off',
      spatiotemporal = 'off'
    )
  #summary(mod_sdmTMB_ICR)

  ## Now fit the censored-Poisson model
  mod_sdmTMB_censored <-
    sdmTMB(
          data=sdmTMB_dat,
          formula=as.formula(paste0(paste0('N_it_',i),' ~ -1 + offset + (1|station_ID) + (1|event_ID)')),
          mesh=make_mesh(sdmTMB_dat,
                         xy_cols=c('lon','lat'),
                         n_knots=3),
          family=censored_poisson(link = 'log'),
          experimental = list(lwr=lower,upr=upper),
          time_varying = ~ -1 + region,
          time = 'year',
          spatial = 'off',
          spatiotemporal = 'off',
          silent = T,
          previous_fit = mod_sdmTMB_poisson
        )

  if(mod_sdmTMB_censored$sd_report$pdHess == FALSE)
  {
    # If fails to converge, change to interval censoring and remove
    # censoring from upper 5% of data
    upper <- as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))

    quant_regions <-
      matrix(by(as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])),
                list(sdmTMB_dat$region,sdmTMB_dat$year),
                FUN = function(x){quantile(x,0.95, na.rm=T)}, simplify = T),
             nrow=4,
             ncol=length(unique(sdmTMB_dat$year)))

    scale_fac <- rep(0, length(sdmTMB_dat$prop_removed))
    scale_fac[sdmTMB_dat$prop_removed>0.95] <-
      comp_factor_fun(signif((sdmTMB_dat[sdmTMB_dat$prop_removed>0.95,]$prop_removed-0.95)/(1-0.95),5),
                      round((1-0.95)*sdmTMB_dat$obsHooksPerSet[sdmTMB_dat$prop_removed>0.95]))

    upper[sdmTMB_dat$prop_removed>0.95 &
            as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]] <-
      as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))[which(sdmTMB_dat$prop_removed >0.95 &
                                                                    as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])] +
      round(
      (sdmTMB_dat$prop_removed[sdmTMB_dat$prop_removed>0.95 &
                                 as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]]-0.95)*
        sdmTMB_dat$obsHooksPerSet[sdmTMB_dat$prop_removed>0.95 &
                                    as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]]*
        scale_fac[sdmTMB_dat$prop_removed>0.95 &
                    as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]])

    # upper[which(sdmTMB_dat$prop_removed >0.95 &
    #               as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])] <-
    #   as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))[which(sdmTMB_dat$prop_removed >0.95 &
    #                                                                 as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])] +
    #   upper[which(sdmTMB_dat$prop_removed >0.95 &
    #                 as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])]

    mod_sdmTMB_censored <-
    sdmTMB(
        data=sdmTMB_dat,
        formula=as.formula(paste0(paste0('N_it_',i),' ~ -1 + offset + (1|station_ID) + (1|event_ID)')),
        mesh=make_mesh(sdmTMB_dat,
                       xy_cols=c('lon','lat'),
                       n_knots=3),
        family=censored_poisson(link = 'log'),
        experimental = list(lwr=lower,upr=upper),
        time_varying = ~ -1 + region,
        time = 'year',
        spatial = 'off',
        spatiotemporal = 'off',
        silent = T,
        previous_fit = mod_sdmTMB_poisson
      )
    cens_poisson_type[count] <- 'interval-censored and upper 5% of data uncensored'
  }
  if(mod_sdmTMB_censored$sd_report$pdHess == FALSE)
    {
      # If fails to converge again, change to interval censoring and remove
      # censoring from upper 15% of data
      upper <- as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))

      quant_regions <-
        matrix(by(as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])),
                  list(sdmTMB_dat$region,sdmTMB_dat$year),
                  FUN = function(x){quantile(x,0.85, na.rm=T)}, simplify = T),
               nrow=4,
               ncol=length(unique(sdmTMB_dat$year)))

      scale_fac <- rep(0, length(sdmTMB_dat$prop_removed))
      scale_fac[sdmTMB_dat$prop_removed>0.95] <-
        comp_factor_fun(signif((sdmTMB_dat[sdmTMB_dat$prop_removed>0.95,]$prop_removed-0.95)/(1-0.95),5),
                        round((1-0.95)*sdmTMB_dat$obsHooksPerSet[sdmTMB_dat$prop_removed>0.95]))

      upper[sdmTMB_dat$prop_removed>0.95 &
              as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]] <-
        as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))[which(sdmTMB_dat$prop_removed >0.95 &
                                                                      as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])] +
        round(
          (sdmTMB_dat$prop_removed[sdmTMB_dat$prop_removed>0.95 &
                                     as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]]-0.95)*
            sdmTMB_dat$obsHooksPerSet[sdmTMB_dat$prop_removed>0.95 &
                                        as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]]*
            scale_fac[sdmTMB_dat$prop_removed>0.95 &
                        as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))]])

      # upper[which(sdmTMB_dat$prop_removed >0.95 &
      #               as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])] <-
      #   as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)]))[which(sdmTMB_dat$prop_removed >0.95 &
      #                                                                 as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])] +
      #   upper[which(sdmTMB_dat$prop_removed >0.95 &
      #                 as.numeric(as.matrix(sdmTMB_dat[,paste0('N_it_',i)])) < quant_regions[cbind(sdmTMB_dat$region_INLA,year_map_fun(sdmTMB_dat$year))])]

      mod_sdmTMB_censored <-
        sdmTMB(
          data=sdmTMB_dat,
          formula=as.formula(paste0(paste0('N_it_',i),' ~ -1 + offset + (1|station_ID) + (1|event_ID)')),
          mesh=make_mesh(sdmTMB_dat,
                         xy_cols=c('lon','lat'),
                         n_knots=3),
          family=censored_poisson(link = 'log'),
          experimental = list(lwr=lower,upr=upper),
          time_varying = ~ -1 + region,
          time = 'year',
          spatial = 'off',
          spatiotemporal = 'off',
          silent = T,
          previous_fit = mod_sdmTMB_poisson
        )
      cens_poisson_type[count] <- 'interval-censored and upper 15% of data uncensored'
  }
  if(mod_sdmTMB_censored$sd_report$pdHess == FALSE)
  {
    cens_poisson_type[count] <- 'Failed to Converge'
  }
  #summary(mod_sdmTMB_censored)

  # Create relative indices by simulating from model
  index_poisson <- NULL
  if(mod_sdmTMB_poisson$sd_report$pdHess)
  {
    index_poisson <-
      expand.grid(Region=levels(sdmTMB_dat$region),
                  Year=as.numeric(sort(unique(sdmTMB_dat$year))),
                  event_ID=sdmTMB_dat$event_ID[1],
                  station_ID=sdmTMB_dat$station_ID[1],
                  lat=1, lon=1,
                  Method='CPUE-based',
                  Species=data_species_vec[count]) %>%
      mutate(preds =
               predict(mod_sdmTMB_poisson,
                       newdata=expand.grid(region=levels(sdmTMB_dat$region),
                                           year=as.numeric(sort(unique(sdmTMB_dat$year))),
                                           event_ID=sdmTMB_dat$event_ID[1],
                                           station_ID=sdmTMB_dat$station_ID[1],
                                           lat=1, lon=1,
                                           #twentyhooks=0,
                                           offset=1),
                       sims=1000,
                       re_form=NULL,
                       re_form_iid=NA)) %>%
      group_by(Region, Year, Method, Species) %>%
      summarise(across(.cols=starts_with("preds"),
                       list(Mean = function(x){mean(exp(x),na.rm=T)},
                            UCL = function(x){quantile(exp(x),probs=c(0.975), na.rm=T)},
                            LCL = function(x){quantile(exp(x),probs=c(0.025), na.rm=T)},
                            SD = function(x){sd(exp(x), na.rm=T)}))) %>%
      group_by(Region, Method, Species) %>%
      mutate(preds_UCL = preds_UCL/mean(preds_Mean),
             preds_LCL = preds_LCL/mean(preds_Mean),
             preds_SD = preds_SD/mean(preds_Mean),
             preds_Mean = preds_Mean/mean(preds_Mean))
  }
  index_ICR <- NULL
  if(mod_sdmTMB_ICR$sd_report$pdHess)
  {
  index_ICR <-
    expand.grid(Region=levels(sdmTMB_dat$region),
                Year=as.numeric(sort(unique(sdmTMB_dat$year))),
                event_ID=sdmTMB_dat$event_ID[1],
                station_ID=sdmTMB_dat$station_ID[1],
                lat=1, lon=1,
                Method='ICR-based',
                Species=data_species_vec[count]) %>%
    mutate(preds =
             predict(mod_sdmTMB_ICR,
                     newdata=expand.grid(region=levels(sdmTMB_dat$region),
                                         year=as.numeric(sort(unique(sdmTMB_dat$year))),
                                         event_ID=sdmTMB_dat$event_ID[1],
                                         station_ID=sdmTMB_dat$station_ID[1],
                                         lat=1, lon=1,
                                         #twentyhooks=0,
                                         offset=1),
                     sims=1000,
                     re_form=NULL,
                     re_form_iid=NA)) %>%
    group_by(Region, Year, Method, Species) %>%
    summarise(across(.cols=starts_with("preds"),
                     list(Mean = function(x){mean(exp(x),na.rm=T)},
                          UCL = function(x){quantile(exp(x),probs=c(0.975), na.rm=T)},
                          LCL = function(x){quantile(exp(x),probs=c(0.025), na.rm=T)},
                          SD = function(x){sd(exp(x), na.rm=T)}))) %>%
    group_by(Region, Method, Species) %>%
    mutate(preds_UCL = preds_UCL/mean(preds_Mean),
           preds_LCL = preds_LCL/mean(preds_Mean),
           preds_SD = preds_SD/mean(preds_Mean),
           preds_Mean = preds_Mean/mean(preds_Mean))
  }
  index_censored <- NULL
  if(mod_sdmTMB_censored$sd_report$pdHess)
  {
  index_censored <-
    expand.grid(Region=levels(sdmTMB_dat$region),
                Year=as.numeric(sort(unique(sdmTMB_dat$year))),
                event_ID=sdmTMB_dat$event_ID[1],
                station_ID=sdmTMB_dat$station_ID[1],
                lat=1, lon=1,
                Method='Censored Poisson',
                Species=data_species_vec[count]) %>%
    mutate(preds =
             predict(mod_sdmTMB_censored,
                     newdata=expand.grid(region=levels(sdmTMB_dat$region),
                                         year=as.numeric(sort(unique(sdmTMB_dat$year))),
                                         event_ID=sdmTMB_dat$event_ID[1],
                                         station_ID=sdmTMB_dat$station_ID[1],
                                         lat=1, lon=1,
                                         #twentyhooks=0,
                                         offset=1),
                     sims=1000,
                     re_form=NULL,
                     re_form_iid=NA)) %>%
    group_by(Region, Year, Method, Species) %>%
    summarise(across(.cols=starts_with("preds"),
                     list(Mean = function(x){mean(exp(x),na.rm=T)},
                          UCL = function(x){quantile(exp(x),probs=c(0.975), na.rm=T)},
                          LCL = function(x){quantile(exp(x),probs=c(0.025), na.rm=T)},
                          SD = function(x){sd(exp(x), na.rm=T)}))) %>%
    group_by(Region, Method, Species) %>%
    mutate(preds_UCL = preds_UCL/mean(preds_Mean),
           preds_LCL = preds_LCL/mean(preds_Mean),
           preds_SD = preds_SD/mean(preds_Mean),
           preds_Mean = preds_Mean/mean(preds_Mean))
  }

  # Store the indices
  results_list[[count]] <- rbind(index_poisson, index_ICR, index_censored)
  print(paste0('finished processing species number ',count,' out of ',length(results_list)))
  count <- count + 1
}

# Plot the indices
results_list[[1]] %>%
  ggplot(aes(x=Year, y=preds_Mean,
             ymax=preds_UCL, ymin=preds_LCL,
             colour=Region, fill=Region)) +
  geom_line() + geom_ribbon(alpha=0.2) +
  facet_grid(Region ~ Method, scales='free_y')

results_list[[1]] %>%
  ggplot(aes(x=Year, y=preds_Mean,
             ymax=preds_UCL, ymin=preds_LCL,
             colour=Method, fill=Method)) +
  geom_line() + geom_ribbon(alpha=0.2) +
  facet_grid(~Region, scales='free_y')

# Compute metrics comparing estimators
results <- do.call('rbind',results_list)
results_summary <-
results %>%
  filter(Year>=3) %>%
  group_by(Species, Region) %>%
  summarize(Correlation_PI = ifelse(sum(Method=='CPUE-based') > 0 &
                                      sum(Method=='ICR-based') > 0,
                                    cor(preds_Mean[Method=='CPUE-based'],
                                      preds_Mean[Method=='ICR-based'],
                                      use='pairwise.complete.obs',
                                      method='spearman'),
                                    NA),
            Correlation_PC = ifelse(sum(Method=='CPUE-based') > 0 &
                                      sum(Method=='Censored Poisson') > 0,
                                    cor(preds_Mean[Method=='CPUE-based'],
                                        preds_Mean[Method=='Censored Poisson'],
                                        use='pairwise.complete.obs',
                                        method='spearman'),
                                    NA),
            Correlation_IC = ifelse(sum(Method=='Censored Poisson') > 0 &
                                      sum(Method=='ICR-based') > 0,
                                    cor(preds_Mean[Method=='Censored Poisson'],
                                        preds_Mean[Method=='ICR-based'],
                                        use='pairwise.complete.obs',
                                        method='spearman'),
                                    NA),
            Percentage_Overlap_PI = ifelse(sum(Method=='CPUE-based') > 0 &
                                      sum(Method=='ICR-based') > 0,
                                    mean(preds_LCL[Method=='CPUE-based'] >
                                         preds_UCL[Method=='ICR-based'] |
                                           preds_UCL[Method=='CPUE-based'] <
                                            preds_UCL[Method=='ICR-based']),
                                    NA),
            Percentage_Overlap_PC = ifelse(sum(Method=='CPUE-based') > 0 &
                                             sum(Method=='Censored Poisson') > 0,
                                           mean(preds_LCL[Method=='CPUE-based'] >
                                                  preds_UCL[Method=='Censored Poisson'] |
                                                  preds_UCL[Method=='CPUE-based'] <
                                                  preds_UCL[Method=='Censored Poisson']),
                                           NA),
            Percentage_Overlap_IC = ifelse(sum(Method=='Censored Poisson') > 0 &
                                             sum(Method=='ICR-based') > 0,
                                           mean(preds_LCL[Method=='Censored Poisson'] >
                                                  preds_UCL[Method=='ICR-based'] |
                                                  preds_UCL[Method=='Censored Poisson'] <
                                                  preds_UCL[Method=='ICR-based']),
                                           NA),
            Interval_Width_Poisson = ifelse(sum(Method=='CPUE-based') > 0,
                                    median(preds_UCL[Method=='CPUE-based']/
                                      preds_LCL[Method=='CPUE-based']),
                                    NA),
            Interval_Width_ICR = ifelse(sum(Method=='ICR-based') > 0,
                                        median(preds_UCL[Method=='ICR-based']/
                                           preds_LCL[Method=='ICR-based']),
                                        NA),
            Interval_Width_Censored = ifelse(sum(Method=='Censored Poisson') > 0,
                                        median(preds_UCL[Method=='Censored Poisson']/
                                             preds_LCL[Method=='Censored Poisson']),
                                        NA),
            MAD_Poisson = ifelse(sum(Method=='CPUE-based') > 0,
                                mad(preds_Mean[Method=='CPUE-based']),
                                    NA),
            MAD_ICR = ifelse(sum(Method=='ICR-based') > 0,
                            mad(preds_Mean[Method=='ICR-based']),
                            NA),
            MAD_Censored = ifelse(sum(Method=='Censored Poisson') > 0,
                                 mad(preds_Mean[Method=='Censored Poisson']),
                                 NA),
            Range_Poisson = ifelse(sum(Method=='CPUE-based') > 0,
                                 range(preds_Mean[Method=='CPUE-based'])[2]/
                                   range(preds_Mean[Method=='CPUE-based'])[1],
                                 NA),
            Range_ICR = ifelse(sum(Method=='ICR-based') > 0,
                             range(preds_Mean[Method=='ICR-based'])[2]/
                               range(preds_Mean[Method=='ICR-based'])[1],
                             NA),
            Range_Censored = ifelse(sum(Method=='Censored Poisson') > 0,
                                  range(preds_Mean[Method=='Censored Poisson'])[2]/
                                    range(preds_Mean[Method=='Censored Poisson'])[1],
                                  NA)
  ) %>%
  ungroup() %>%
  summarise(Comparison=rep(c('CPUE-based vs ICR-based','CPUE-based vs Censored','ICR-based vs Censored'),5),
            Measure=rep(c('Correlation','Percentage Overlap','Interval Width','MAD','Range'),each=3),
            Median=c(median(Correlation_PI, na.rm=T),
                    median(Correlation_PC, na.rm=T),
                    median(Correlation_IC, na.rm=T),
                    median(Percentage_Overlap_PI, na.rm=T),
                    median(Percentage_Overlap_PC, na.rm=T),
                    median(Percentage_Overlap_IC, na.rm=T),
                    median(Interval_Width_Poisson-Interval_Width_ICR, na.rm=T),
                    median(Interval_Width_Poisson-Interval_Width_Censored, na.rm=T),
                    median(Interval_Width_ICR-Interval_Width_Censored, na.rm=T),
                    median(MAD_Poisson-MAD_ICR, na.rm=T),
                    median(MAD_Poisson-MAD_Censored, na.rm=T),
                    median(MAD_ICR-MAD_Censored, na.rm=T),
                    median(Range_Poisson-Range_ICR, na.rm=T),
                    median(Range_Poisson-Range_Censored, na.rm=T),
                    median(Range_ICR-Range_Censored, na.rm=T)),
            MAD=c(mad(Correlation_PI, na.rm=T),
                  mad(Correlation_PC, na.rm=T),
                  mad(Correlation_IC, na.rm=T),
                  mad(Percentage_Overlap_PI, na.rm=T),
                  mad(Percentage_Overlap_PC, na.rm=T),
                  mad(Percentage_Overlap_IC, na.rm=T),
                  mad(Interval_Width_Poisson-Interval_Width_ICR, na.rm=T),
                  mad(Interval_Width_Poisson-Interval_Width_Censored, na.rm=T),
                  mad(Interval_Width_ICR-Interval_Width_Censored, na.rm=T),
                  mad(MAD_Poisson-MAD_ICR, na.rm=T),
                  mad(MAD_Poisson-MAD_Censored, na.rm=T),
                  mad(MAD_ICR-MAD_Censored, na.rm=T),
                  mad(Range_Poisson-Range_ICR, na.rm=T),
                  mad(Range_Poisson-Range_Censored, na.rm=T),
                  mad(Range_ICR-Range_Censored, na.rm=T)),
            N_comparisons=c(sum(!is.na(Correlation_PI)),
                            sum(!is.na(Correlation_PC)),
                            sum(!is.na(Correlation_IC)),
                            sum(!is.na(Percentage_Overlap_PI)),
                            sum(!is.na(Percentage_Overlap_PC)),
                            sum(!is.na(Percentage_Overlap_IC)),
                            sum(!is.na(Interval_Width_Poisson-Interval_Width_ICR)),
                            sum(!is.na(Interval_Width_Poisson-Interval_Width_Censored)),
                            sum(!is.na(Interval_Width_ICR-Interval_Width_Censored)),
                            sum(!is.na(MAD_Poisson-MAD_ICR)),
                            sum(!is.na(MAD_Poisson-MAD_Censored)),
                            sum(!is.na(MAD_ICR-MAD_Censored)),
                            sum(!is.na(Range_Poisson-Range_ICR)),
                            sum(!is.na(Range_Poisson-Range_Censored)),
                            sum(!is.na(Range_ICR-Range_Censored))))

results_summary$Comparison <- factor(results_summary$Comparison,
                                     levels=c('CPUE-based vs ICR-based',
                                              'CPUE-based vs Censored',
                                              'ICR-based vs Censored'),
                                     ordered = T)
results_summary$Measure <- factor(results_summary$Measure,
                                     levels=c('Correlation','Percentage Overlap','Interval Width','MAD','Range'),
                                     ordered = T)

results_summary %>%
  ggplot(aes(x=Comparison, y=Median,
         ymax=Median+1.96*(MAD/sqrt(N_comparisons)),
         ymin=Median-1.96*(MAD/sqrt(N_comparisons)),
         colour=Comparison, group=Comparison)) +
  geom_errorbar() + geom_point() +
  geom_hline(mapping=aes(yintercept=ifelse(Measure%in%c('Correlation','Percentage Overlap'),1,0))) +
  facet_wrap(~Measure, scales='free_y', nrow=3, ncol=2) +
  ylab('Median Correlation or Median Difference') +
  xlab('Estimators Being Compared')









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

if(fit_mods)
{
  brms_mod_halibut <-
    brm(N_it_halibut | rate(effSkateIPHC) ~
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
  brms_mod_halibut <- readRDS('halibut_mod_final.rds')
}
summary(brms_mod_halibut)
# This function plots the posterior prediction on the log scale
#conditional_smooths(brms_mod_halibut, spaghetti = T)

# Unfortunately, brms does not allow the smooths to be plotted on the original scale
# without incorporating the (large) uncertainty associated with the intercept
# e.g. see conditional_effects(brms_mod_halibut)
# Manually predict on the linear predictor scale and subtract the sampled intercept terms
newdata <- inlabru:::cprod(
  data.frame(effSkateIPHC=1, region=factor(c('1','2','3','4')), twentyhooks=0, station_ID='2004', event_ID='1'),
  data.frame(year=min(brms_dat$year):max(brms_dat$year))
)
# Write a function to perform this
brms_trend_plotter <- function(mod, newdata)
{
  # sample from the posterior
  tmp2 <- posterior_linpred(mod,
                            newdata = newdata,
                            summary=F,
                            re_formula = NA,
                            transform = F
  )
  tmp3 <- posterior_samples(mod,"^b")
  # remove the intercept
  tmp2 <- tmp2 - tmp3$b_Intercept
  # remove the region 2 effect
  tmp2[,which(newdata$region==2)] <- tmp2[,which(newdata$region==2)] - tmp3$b_region2
  tmp2[,which(newdata$region==3)] <- tmp2[,which(newdata$region==3)] - tmp3$b_region3
  tmp2[,which(newdata$region==4)] <- tmp2[,which(newdata$region==4)] - tmp3$b_region4

  return(cbind(newdata,posterior_summary(exp(tmp2))))
}
brms_mod_halibut_pred <- brms_trend_plotter(brms_mod_halibut, newdata = newdata)
brms_mod_halibut_pred$model = 'CPUE-based'

# Now fit the hook competition scale factor-adjusted model
brms_dat$N_it_halibut_compfactor <-
  round(brms_dat$N_it_halibut *
  comp_factor_fun(brms_dat$prop_removed, brms_dat$obsHooksPerSet))

if(fit_mods)
{
brms_mod_halibut_compfactor <-
  brm(N_it_halibut_compfactor | rate(effSkateIPHC) ~
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
  brms_mod_halibut_compfactor <- readRDS('halibut_compfactor_mod_final.rds')
}
summary(brms_mod_halibut_compfactor)

brms_mod_halibut_compfactor_pred <- brms_trend_plotter(brms_mod_halibut_compfactor, newdata = newdata)
brms_mod_halibut_compfactor_pred$model = 'ICR-based'

# Now fit the censored model
# Due to qwerk of brms, need to redefine count variable (censorship interval open on left)
brms_dat_halibut <- brms_dat[which(!(brms_dat$censored_txt_95 == 'right' &
                                brms_dat$N_it_halibut == 0)),]
brms_dat_halibut$N_it_halibut[which(brms_dat_halibut$censored_txt_95 == 'right')] <-
  brms_dat_halibut$N_it_halibut[which(brms_dat_halibut$censored_txt_95 == 'right')] - 1

brms_dat_yelloweye <- brms_dat[which(!(brms_dat$censored_txt_95 == 'right' &
                                       brms_dat$N_it_yelloweye == 0)),]
brms_dat_yelloweye$N_it_yelloweye[which(brms_dat_yelloweye$censored_txt_95 == 'right')] <-
  brms_dat_yelloweye$N_it_yelloweye[which(brms_dat_yelloweye$censored_txt_95 == 'right')] - 1

if(fit_mods)
{
brms_mod_halibut_cens <-
  brm(N_it_halibut | cens(censored_txt_95) + rate(effSkateIPHC) ~
        region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8) ,
      data = brms_dat_halibut,
      family='Poisson',
      iter=4000, warmup = 3000,
      cores = 4, chains = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      seed=8102021,
      prior = c(prior(student_t(7, 0, 2), class = sd, group = event_ID, coef = Intercept)))
}
if(!fit_mods)
{
  brms_mod_halibut_cens <- readRDS('halibut_cens_mod_final.rds')
}
summary(brms_mod_halibut_cens)

brms_mod_halibut_cens_pred <- brms_trend_plotter(brms_mod_halibut_cens, newdata = newdata)
brms_mod_halibut_cens_pred$model = 'Censored Poisson'

boot_halibut_pred <- brms_mod_halibut_pred
# recode factor to match brm mod output
combined_indices_halibut_INLA_Alldata$combined_ts_all$region <-
  forcats::fct_recode(combined_indices_halibut_INLA_Alldata$combined_ts_all$region,
                      `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
                      `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')
boot_halibut_pred$region <-
  forcats::fct_recode(boot_halibut_pred$region,
                      `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                      `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
# calibrate the year variable
combined_indices_halibut_INLA_Alldata$combined_ts_all$year <-
  combined_indices_halibut_INLA_Alldata$combined_ts_all$year - 1995
# match by year and region
boot_halibut_pred <- inner_join(
  boot_halibut_pred, combined_indices_halibut_INLA_Alldata$combined_ts_all[
    combined_indices_halibut_INLA_Alldata$combined_ts_all$model=='boot',
  ],
  by = c('region','year')
)
boot_halibut_pred$Estimate <- boot_halibut_pred$mean
boot_halibut_pred$Q2.5 <- NA#boot_halibut_pred$LCL
boot_halibut_pred$Q97.5 <- NA#boot_halibut_pred$UCL
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
                                          `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                          `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')

# predict the linear trend of hook saturation
hook_sat_df <- hook_sat_df[!is.na(hook_sat_df$region),]
hook_sat_df$Estimate <-
  predict(
    lm(data=hook_sat_df,
       Estimate2 ~ region * year)
  )
hook_sat_df$Q2.5 <- hook_sat_df$Estimate -
  1.96* predict(
    lm(data=hook_sat_df,
       Estimate2 ~ region * year),
    se=T
  )$se.fit
hook_sat_df$Q97.5 <- hook_sat_df$Estimate +
  1.96* predict(
    lm(data=hook_sat_df,
       Estimate2 ~ region * year),
    se=T
  )$se.fit

halibut_pred <- rbind(brms_mod_halibut_pred,brms_mod_halibut_compfactor_pred,brms_mod_halibut_cens_pred,boot_halibut_pred)

# halibut_pred <- rbind(brms_mod_halibut_pred,brms_mod_halibut_compfactor_pred,brms_mod_halibut_cens_pred)
#
halibut_pred$region <-
   forcats::fct_recode(halibut_pred$region,
                       `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                       `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
halibut_pred$year <- halibut_pred$year + 1995
# remove the predictions from WCVI pre 1999
halibut_pred[halibut_pred$region=='West Coast Van Island' &
               halibut_pred$year < 1999, c('Estimate','Q2.5','Q97.5')] <- NA

halibut_pred <-
  full_join(halibut_pred,
            hook_sat_df[,c("year", "Estimate", "Estimate2", "model", "region","Q2.5",'Q97.5')],
            suffix=c('',''))

####
halibut_pred <-
full_join(halibut_pred,
          hook_sat_df[,c("year", "Estimate", "Estimate2", "model", "region","Q2.5",'Q97.5')],
          suffix=c('',''))

# Are these results consistent with INLA?
INLA_mod_halibut <-
  final_censored_index_fun(data=Reduced_data_sp[which(!(is.na(Reduced_data_sp$region_INLA) |
                                                          is.na(Reduced_data_sp$prop_removed))),],
                           species = 'halibut',
                           cprop=0.95,
                           use_upper_bound = T,
                           upper_bound_quantile = 0.8,
                           init_vals = c(1.287, 1.042, 1.409, -1.288),
                           n_knots = 8)

# rename variables to make it ready for merging with brms indices
INLA_halibutpred <-
  INLA_mod_halibut$pred_overdisp %>%
  rename(Q2.5=q0.025,
         Q97.5=q0.975,
         Estimate = mean) %>%
  mutate(model='Adjusted Censored Poisson') %>%
  filter(region != 'All')

INLA_halibutpred$region <-
  forcats::fct_recode(INLA_halibutpred$region,
                      `Hectate Strait` = 'HS', `Queen Charlotte Sound` = 'QCS',
                      `West Coast Haida Gwaii` = 'WCHG', `West Coast Van Island` = 'WCVI')

halibut_pred <-
  full_join(halibut_pred,INLA_halibutpred)

halibut_pred$model <-
  factor(halibut_pred$model,
         levels=c('hook saturation level',"boot","CPUE-based","ICR-based", "Adjusted Censored Poisson", "Censored Poisson"),
         ordered = T)

halibut_pred %>%
  filter(!is.na(region) & model!='Adjusted Censored Poisson') %>%
  #filter(year > 1997) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q97.5, ymin=Q2.5,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Pacific Halibut Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Saturated Events') +
  geom_hline(yintercept = 0)

halibut_pred %>%
  filter(!is.na(region) & (model %in% c('Adjusted Censored Poisson','ICR-based','Censored Poisson') )) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q97.5, ymin=Q2.5,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Pacific Halibut Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Censored Observations') +
  geom_hline(yintercept = 0)

# What could be causing the lack of difference in censored Poisson indices?
by(brms_dat$N_it_halibut,
   brms_dat$censored_txt_95=='right',
   quantile, probs=c(0.95,0.96,0.97,0.98,0.99,1))
# The largest catch counts correspond to censored observations! We do not observe upper tail!

# Observe an increasing upper tail with time
by(brms_dat$N_it_halibut[brms_dat$censored_txt_95=='right'],
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
                               `ICR-based`,
                               method = 'spearman', use='complete.obs'),
            Cor_CPUE_Cen = cor(`CPUE-based`,
                               `Censored Poisson`,
                               method = 'spearman', use='complete.obs'),
            Cor_Cen_ICR = cor(`Censored Poisson`,
                               `ICR-based`,
                              method = 'spearman', use='complete.obs'),
            Cor_HookSat_CPUE = cor(`CPUE-based`,
                               `hook saturation level`,
                               method = 'spearman', use='complete.obs'),
            Cor_HookSat_ICR = cor(`hook saturation level`,
                                  `ICR-based`,
                                  method = 'spearman', use='complete.obs'),
            Cor_HookSat_Cen = cor(`Censored Poisson`,
                              `hook saturation level`,
                              method = 'spearman', use='complete.obs'))
cor_results
xtable::xtable(cor_results)

# Correlations with the adjusted method?
cor_results_sens <-
  halibut_pred %>%
  dplyr::select(region, year, Estimate, model) %>%
  pivot_wider(names_from=model, values_from = c(Estimate)) %>%
  group_by(region) %>%
  summarise(Cor_CPUE_Adjust = cor(`CPUE-based`,
                               `Adjusted Censored Poisson`,
                               method = 'spearman', use='complete.obs'),
            Cor_ICR_Adjust = cor(`ICR-based`,
                               `Adjusted Censored Poisson`,
                               method = 'spearman', use='complete.obs'),
            Cor_Cen_Adjust = cor(`Censored Poisson`,
                              `Adjusted Censored Poisson`,
                              method = 'spearman', use='complete.obs'))
cor_results_sens
xtable::xtable(t(cor_results_sens))
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
brms_mod_yelloweye_compfactor_pred$model = 'ICR-based'

# Now fit the censored model

if(fit_mods)
{
brms_mod_yelloweye_cens <-
  brm(N_it_yelloweye | cens(censored_txt_95) + rate(effSkateIPHC) ~
        region + twentyhooks + (1 | station_ID) + (1 | event_ID) + s(year, by=region, k=8) ,
      data = brms_dat_yelloweye,
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
by(brms_dat$N_it_yelloweye, brms_dat$region_INLA, quantile, probs=c(0.8, 0.85, 0.9))
# Note that 80th percentile is zero! Sensitivity analysis at 0.85 instead else it will
# be the CPUE-based estimator

INLA_mod_yelloweye <-
  final_censored_index_fun(data=Reduced_data_sp[which(!(is.na(Reduced_data_sp$region_INLA) |
                                                          is.na(Reduced_data_sp$prop_removed))),],
                           species = 'yelloweye',
                           cprop=0.95,
                           use_upper_bound = T,
                           upper_bound_quantile = 0.85,
                           n_knots = 8,
                           init_vals = c(0.615, -1.985, 2.654, -0.677))

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
         levels=c('hook saturation level',"boot","CPUE-based","ICR-based", "Adjusted Censored Poisson", "Censored Poisson"),
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
  ylab('Estimated Relative Abundance or Proportion of Censored Observations') +
  geom_hline(yintercept = 0)

yelloweye_pred %>%
  filter(!is.na(region) & model!='Adjusted Censored Poisson') %>%
  filter(year > 1997) %>%
  ggplot(aes(x=year, y=Estimate, ymax=Q97.5, ymin=Q2.5,colour=region, group=region, fill=region)) +
  geom_ribbon(alpha=0.2) + geom_line() + facet_grid(model~region, scales='free_y') +
  ggtitle('Yelloweye Rockfish Relative Abundance Indices',
          subtitle = 'Three different indices are shown across 4 different regions') +
  geom_hline(yintercept = 1, linetype='dotted') +
  geom_point(aes(x=year, y=Estimate2), colour='grey') +
  scale_x_continuous('Year', c(1995, 2000, 2005, 2010, 2015), labels=c(1995, 2000, 2005, 2010, 2015)) +
  ylab('Estimated Relative Abundance or Proportion of Censored Observations') +
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

cor_results_yelloweye <-
  yelloweye_pred %>%
  dplyr::select(region, year, Estimate, model) %>%
  pivot_wider(names_from=model, values_from = c(Estimate)) %>%
  group_by(region) %>%
  summarise(Cor_CPUE_ICR = cor(`CPUE-based`,
                               `ICR-based`,
                               method = 'spearman', use='complete.obs'),
            Cor_CPUE_Cen = cor(`CPUE-based`,
                               `Censored Poisson`,
                               method = 'spearman', use='complete.obs'),
            Cor_Cen_ICR = cor(`Censored Poisson`,
                              `ICR-based`,
                              method = 'spearman', use='complete.obs'),
            Cor_HookSat_CPUE = cor(`CPUE-based`,
                                   `hook saturation level`,
                                   method = 'spearman', use='complete.obs'),
            Cor_HookSat_ICR = cor(`hook saturation level`,
                                  `ICR-based`,
                                  method = 'spearman', use='complete.obs'),
            Cor_HookSat_Cen = cor(`Censored Poisson`,
                                  `hook saturation level`,
                                  method = 'spearman', use='complete.obs'))
cor_results_yelloweye
xtable::xtable(t(cor_results_yelloweye))

cor_results_yelloweye_sens <-
  yelloweye_pred %>%
  dplyr::select(region, year, Estimate, model) %>%
  pivot_wider(names_from=model, values_from = c(Estimate)) %>%
  group_by(region) %>%
  summarise(Cor_CPUE_Adjust = cor(`CPUE-based`,
                                  `Adjusted Censored Poisson`,
                                  method = 'spearman', use='complete.obs'),
            Cor_ICR_Adjust = cor(`ICR-based`,
                                 `Adjusted Censored Poisson`,
                                 method = 'spearman', use='complete.obs'),
            Cor_Cen_Adjust = cor(`Censored Poisson`,
                                 `Adjusted Censored Poisson`,
                                 method = 'spearman', use='complete.obs'))
cor_results_yelloweye_sens
xtable::xtable(t(cor_results_yelloweye_sens))

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
  group_by(model, species, region) %>%
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

xtable::xtable(Estimator_Properties[!(Estimator_Properties$model %in% c('hook saturation level','boot')) &
                                      Estimator_Properties$species=='pacific halibut', ])

xtable::xtable(Estimator_Properties[!(Estimator_Properties$model %in% c('hook saturation level','boot')) &
                                      Estimator_Properties$species=='yelloweye rockfish', ])

left_join(yelloweye_pred,
          hook_sat_df, suffix=c('','')) %>%
  filter(year < 2020) %>%
  group_by(region, model) %>%
  ggplot(aes(x=prop_saturation, y=Estimate, colour=region)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(region~model)

### Posterior Predictive Checks
library(bayesplot)
bayesplot::ppc_intervals_grouped(y=brms_mod_halibut$data$N_it_halibut[brms_dat$prop_removed >= 0.95],
                            yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_dat$prop_removed >= 0.95],
                            group = forcats::fct_recode(factor(brms_mod_halibut$data$region[brms_dat$prop_removed >= 0.95]),
                                                        `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                        `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4'),
                            x = 1995+brms_mod_halibut$data$year[brms_dat$prop_removed >= 0.95]) + xlab('Year') + ylab('Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Halibut CPUE-based Model')
bayesplot::ppc_intervals_grouped(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 group = forcats::fct_recode(factor(brms_mod_halibut_cens$data$region[brms_mod_halibut_cens$data$censored_txt_95=='none']),
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4'),
                                 x = 1995+brms_mod_halibut_cens$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none']) + xlab('Year') + ylab('Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Halibut Censored Poisson Model')

bayesplot::ppc_ecdf_overlay(y=brms_mod_halibut$data$N_it_halibut[brms_dat$prop_removed >= 0.95],
                            yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_dat$prop_removed >= 0.95]) + xlab('Year') + ylab('Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Halibut CPUE-based Model')

bayesplot::ppc_ecdf_overlay(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'],
                            yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none']) + xlab('Count') + ylab('Empirical Cumulative Density Function') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Halibut Censored Poisson Model')

bayesplot::ppc_scatter_avg_grouped(y=brms_mod_halibut$data$N_it_halibut[brms_dat$prop_removed >= 0.95],
                                   yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_dat$prop_removed >= 0.95],
                                   group=forcats::fct_recode(brms_mod_halibut$data$region,
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')[brms_dat$prop_removed >= 0.95]) +
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue')
bayesplot::ppc_scatter_avg_grouped(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none'],
                                 group=forcats::fct_recode(brms_mod_halibut_cens$data$region,
                                                           `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                           `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')[brms_mod_halibut_cens$data$censored_txt_95=='none']) +
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue') + xlab('Average Sampled Catch Count') + ylab('Observed Catch Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Halibut Censored Poisson Model')

# Again for yelloweye
bayesplot::ppc_intervals_grouped(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 group = forcats::fct_recode(factor(brms_mod_yelloweye$data$region[brms_mod_yelloweye_cens$data$censored_txt_95=='none']),
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4'),
                                 x = 1995+brms_mod_yelloweye$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none']) + xlab('Year') + ylab('Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Yelloweye CPUE-based Model')
bayesplot::ppc_intervals_grouped(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                                 group = forcats::fct_recode(factor(brms_mod_yelloweye_cens$data$region[brms_mod_yelloweye_cens$data$censored_txt_95=='none']),
                                                             `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
                                                             `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4'),
                                 x = 1995+brms_mod_yelloweye_cens$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none']) + xlab('Year') + ylab('Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Yelloweye Censored Poisson Model')

bayesplot::ppc_ecdf_overlay(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none']) + xlab('Year') + ylab('Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Yelloweye CPUE-based Model')

bayesplot::ppc_ecdf_overlay(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none'],
                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none']) + xlab('Count') + ylab('Empirical Cumulative Density Function') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Yelloweye Censored Poisson Model')

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
  geom_smooth(colour='red') + geom_abline(intercept = 0, slope=1, colour='blue') + xlab('Average Sampled Catch Count') + ylab('Observed Catch Count') +
  ggtitle('Posterior Predictive Checks',subtitle = 'Yelloweye Censored Poisson Model')


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
      bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_halibut$data$N_it_halibut[brms_dat$prop_removed < 0.95 &
                                                                              brms_dat$region==j],
                                         yrep=as.matrix(posterior_predict(brms_mod_halibut))[,brms_dat$prop_removed < 0.95 &
                                                                                               brms_dat$region==j],
                                         x=1995+brms_mod_halibut$data$year[brms_dat$prop_removed < 0.95 &
                                                                             brms_dat$region==j]) +
      geom_smooth(colour='red', method='gam') +
      geom_hline(yintercept = 0) + ggtitle('Prediction Error CPUE-Based', subtitle = paste0('Restricted to fishing events with <95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) + xlab('Year')

    plot_list_halibut[[2+(j-1)*3]] <-
      bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_halibut_compfactor$data$N_it_halibut[brms_dat$prop_removed < 0.95 &
                                                                                              brms_dat$region==j],
                                         yrep=as.matrix(posterior_predict(brms_mod_halibut_compfactor))[,brms_dat$prop_removed < 0.95 &
                                                                                                          brms_dat$region==j],
                                         x=1995+brms_mod_halibut_compfactor$data$year[brms_dat$prop_removed < 0.95 &
                                                                                        brms_dat$region==j]) +
      geom_smooth(colour='red', method='gam') +
      geom_hline(yintercept = 0) + ggtitle('Prediction Error Comp-Factor Adjust', subtitle = paste0('Restricted to fishing events with <95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) + xlab('Year')


    plot_list_halibut[[3+(j-1)*3]] <-
      bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                        brms_mod_halibut_cens$data$region==j],
                                         yrep=as.matrix(posterior_predict(brms_mod_halibut_cens))[,brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                                                    brms_mod_halibut_cens$data$region==j],
                                         x=1995+brms_mod_halibut_cens$data$year[brms_mod_halibut_cens$data$censored_txt_95=='none' &
                                                                             brms_mod_halibut_cens$data$region==j]) +
      geom_smooth(colour='red', method='gam') +
      geom_hline(yintercept = 0) + ggtitle('Prediction Error Censored', subtitle = paste0('Restricted to fishing events with <95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) + xlab('Year')

    post_pred_df <-
      data.frame(Censored_pred = as.numeric(posterior_predict(brms_mod_halibut_cens,draw_ids = samp_ind)[,brms_mod_halibut_cens$data$censored_txt_95=='right' &
                                                                                                           brms_mod_halibut_cens$data$region==j]),
                 CPUE_pred = as.numeric(posterior_predict(brms_mod_halibut,draw_ids = samp_ind)[,brms_dat$prop_removed >= 0.95 &
                                                                                                  brms_dat$region==j &
                                                                                                  brms_dat$N_it_halibut != 0]),
                 Adj_pred = as.numeric(posterior_predict(brms_mod_halibut_compfactor,draw_ids = samp_ind)[,brms_dat$prop_removed >= 0.95 &
                                                                                                            brms_dat$region==j &
                                                                                                            brms_dat$N_it_halibut != 0]),
                 Observed = rep(as.numeric(brms_mod_halibut_cens$data$N_it_halibut[brms_mod_halibut_cens$data$censored_txt_95=='right' &
                                                                                     brms_mod_halibut_cens$data$region==j]),
                                each=100))

    plot_list_halibut2[[1+(j-1)*2]] <-
    ggplot(data=post_pred_df,
           aes(x=Observed, y=Censored_pred)) +
      geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(post_pred_df)) + ggtitle('Censored Posterior Predictions vs. \nObserved Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed Halibut Catch Counts') + ylab('Posterior Predicted Halibut Counts')

    plot_list_halibut2[[2+(j-1)*2]] <-
    ggplot(data=post_pred_df,
           aes(x=Observed, y=Adj_pred)) +
      geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(post_pred_df)) + ggtitle('Scale Factor Adjusted Posterior Predictions vs. \nObserved Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed Halibut Catch Counts') + ylab('Posterior Predicted Halibut Counts')

    # Now plot the variance
    plot_list_halibut4[[1+(j-1)*2]] <-
      ggplot(data=post_pred_df,
             aes(x=Observed, y=abs(Censored_pred-Observed) )) +
      geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(abs(post_pred_df[,c(1,3)]-post_pred_df[,c(4)]))) + ggtitle('Censored Posterior Absolute Prediction Difference vs. \nObserved Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed Halibut Catch Counts') + ylab('Posterior Absolute Prediction Difference')

    plot_list_halibut4[[2+(j-1)*2]] <-
      ggplot(data=post_pred_df,
             aes(x=Observed, y=abs(Adj_pred-Observed))) +
      geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(abs(post_pred_df[,c(1,3)]-post_pred_df[,c(4)]))) + ggtitle('Scale Factor Adjusted Posterior Absolute Prediction Difference vs. \nObserved Halibut Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed Halibut Catch Counts') + ylab('Posterior Absolute Prediction Difference')

    # Update y-axis scale for later plots
    yrange_halibut <- range(c(yrange_halibut,range(post_pred_df)))
    yrange_halibut2 <- range(c(yrange_halibut2,range(abs(post_pred_df[,c(1,3)]-post_pred_df[,c(4)]))))

    # Repeat for yelloweye
    plot_list_yelloweye[[1+(j-1)*3]] <-
      bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_yelloweye$data$N_it_yelloweye[brms_dat$prop_removed < 0.95 &
                                                                                       brms_dat$region==j],
                                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye))[,brms_dat$prop_removed < 0.95 &
                                                                                                    brms_dat$region==j],
                                            x=1995+brms_mod_yelloweye$data$year[brms_dat$prop_removed < 0.95 &
                                                                                       brms_dat$region==j]) +
      geom_smooth(colour='red', method='gam') +
      geom_hline(yintercept = 0) + ggtitle('Prediction Error CPUE-Based', subtitle = paste0('Restricted to fishing events with <95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) + xlab('Year')

    plot_list_yelloweye[[2+(j-1)*3]] <-
      bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_yelloweye_compfactor$data$N_it_yelloweye[brms_dat$prop_removed < 0.95 &
                                                                                                  brms_dat$region==j],
                                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye_compfactor))[,brms_dat$prop_removed < 0.95 &
                                                                                                               brms_dat$region==j],
                                            x=1995+brms_mod_yelloweye_compfactor$data$year[brms_dat$prop_removed < 0.95 &
                                                                                       brms_dat$region==j]) +
      geom_smooth(colour='red', method='gam') +
      geom_hline(yintercept = 0) + ggtitle('Prediction Error Comp-Factor Adjust', subtitle = paste0('Restricted to fishing events with <95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) + xlab('Year')


    plot_list_yelloweye[[3+(j-1)*3]] <-
      bayesplot::ppc_error_scatter_avg_vs_x(y=brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                        brms_mod_yelloweye_cens$data$region==j],
                                            yrep=as.matrix(posterior_predict(brms_mod_yelloweye_cens))[,brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                                       brms_mod_yelloweye_cens$data$region==j],
                                            x=1995+brms_mod_yelloweye_cens$data$year[brms_mod_yelloweye_cens$data$censored_txt_95=='none' &
                                                                                     brms_mod_yelloweye_cens$data$region==j]) +
      geom_smooth(colour='red', method='gam') +
      geom_hline(yintercept = 0) + ggtitle('Prediction Error Censored', subtitle = paste0('Restricted to fishing events with <95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) + xlab('Year')

    post_pred_df <-
      data.frame(Censored_pred = as.numeric(posterior_predict(brms_mod_yelloweye_cens,draw_ids = samp_ind)[,brms_mod_yelloweye_cens$data$censored_txt_95=='right' &
                                                                                                           brms_mod_yelloweye_cens$data$region==j]),
                 CPUE_pred = as.numeric(posterior_predict(brms_mod_yelloweye,draw_ids = samp_ind)[,brms_dat$prop_removed >= 0.95 &
                                                                                                    brms_dat$region==j &
                                                                                                    brms_dat$N_it_yelloweye != 0]),
                 Adj_pred = as.numeric(posterior_predict(brms_mod_yelloweye_compfactor,draw_ids = samp_ind)[,brms_dat$prop_removed >= 0.95 &
                                                                                                              brms_dat$region==j &
                                                                                                              brms_dat$N_it_yelloweye != 0]),
                 Observed = rep(as.numeric(brms_mod_yelloweye_cens$data$N_it_yelloweye[brms_dat$prop_removed >= 0.95 &
                                                                                         brms_dat$region==j &
                                                                                         brms_dat$N_it_yelloweye != 0]),
                                each=100))

    plot_list_yelloweye2[[1+(j-1)*2]] <-
      ggplot(data=post_pred_df,
             aes(x=Observed, y=Censored_pred)) +
      geom_point() + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(post_pred_df)) + ggtitle('Censored Posterior Predictions vs. \nObserved yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed yelloweye Catch Counts') + ylab('Posterior Predicted yelloweye Counts')

    plot_list_yelloweye2[[2+(j-1)*2]] <-
      ggplot(data=post_pred_df,
             aes(x=Observed, y=Adj_pred)) +
      geom_point() + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(post_pred_df)) + ggtitle('Scale Factor Adjusted Posterior Predictions vs. \nObserved yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed yelloweye Catch Counts') + ylab('Posterior Predicted yelloweye Counts')

    plot_list_yelloweye4[[1+(j-1)*2]] <-
      ggplot(data=post_pred_df,
             aes(x=Observed, y=abs(Censored_pred-Observed) )) +
      geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(abs(post_pred_df[,c(1,3)]-post_pred_df[,c(4)]))) + ggtitle('Censored Posterior Absolute Prediction Difference vs. \nObserved Yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed Yelloweye Catch Counts') + ylab('Posterior Absolute Prediction Difference')

    plot_list_yelloweye4[[2+(j-1)*2]] <-
      ggplot(data=post_pred_df,
             aes(x=Observed, y=abs(Adj_pred-Observed))) +
      geom_point() + geom_smooth(colour='red', method='gam') + geom_abline(intercept = 0, slope = 1, colour='blue') +
      ylim(range(abs(post_pred_df[,c(1,3)]-post_pred_df[,c(4)]))) + ggtitle('Scale Factor Adjusted Posterior Absolute Prediction Difference vs. \nObserved Yelloweye Catch Counts', subtitle = paste0('Restricted to fishing events with >95% \nhook saturation in ' ,c('HS', 'QCS','WCHG', 'WCVI')[j])) +
      xlab('Observed Yelloweye Catch Counts') + ylab('Posterior Absolute Prediction Difference')

    # Update y-axis scale for later plots
    yrange_yelloweye <- range(c(yrange_yelloweye,range(post_pred_df)))
    yrange_yelloweye2 <- range(c(yrange_yelloweye2,range(abs(post_pred_df[,c(1,3)]-post_pred_df[,c(4)]))))

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

  plot_list_halibut5[[i]] <- plot_list_halibut4[[i]] + coord_cartesian(ylim=c(0,rep(sqrt(c(120000,1200000)),times=4)[i]) )
  # remove the geom_point layer
  plot_list_halibut5[[i]]$layers[[1]] <- NULL
  plot_list_yelloweye5[[i]] <- plot_list_yelloweye4[[i]] + coord_cartesian(ylim=c(0,rep(sqrt(c(1000000,1000000)),times=4)[i]) )
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
                   Censored = factor(ifelse(brms_dat_halibut$censored_txt_95=='right','yes','no')),
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
                   Censored = factor(ifelse(brms_dat_yelloweye$censored_txt_95=='right','yes','no')),
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

multiplot(
brms_dat %>%
  group_by(region) %>%
  ggplot(aes(x=N_it_halibut, group=region)) +
  geom_density() +
  xlim(c(0,310)) +
  xlab('Catch Count Pacific Halibut') +
  facet_grid(~region),

brms_dat %>%
  group_by(region) %>%
  ggplot(aes(x=N_it_yelloweye, group=region)) +
  geom_density(trim=T) +
  coord_cartesian(xlim=c(0,310)) +
  xlab('Catch Count Yelloweye Rockfish') +
  facet_grid(~region),
layout = matrix(c(1:2),nrow=2,ncol=1)
)

multiplot(
  brms_dat %>%
    mutate(region=
             fct_recode(factor(region),
                        `HS` = '1', `QCS` = '2',
                        `WCHG` = '3', `WCVI` = '4')
    )%>%
    group_by(region, censored_txt_95) %>%
    ggplot(aes(x=N_it_halibut, group=censored_txt_95, colour=censored_txt_95)) +
    geom_density(trim=T) +
    #xlim(c(0,310)) +
    xlab('Catch Count Pacific Halibut') +
    facet_grid(region~censored_txt_95, scales='free_y') +
    theme(legend.position = 0),

  brms_dat %>%
    mutate(region=
             fct_recode(factor(region),
                        `HS` = '1', `QCS` = '2',
                        `WCHG` = '3', `WCVI` = '4')
    )%>%
    group_by(region, censored_txt_95) %>%
    ggplot(aes(x=N_it_yelloweye, group=censored_txt_95, colour=censored_txt_95)) +
    geom_density(trim=T) +
    #xlim(c(0,310)) +
    xlab('Catch Count Yelloweye Rockfish') +
    facet_grid(region~censored_txt_95, scales='free_y') +
    theme(legend.position = 0),

  # brms_dat %>%
  #   mutate(region=
  #            fct_recode(factor(region),
  #                       `Hectate Strait` = '1', `Queen Charlotte Sound` = '2',
  #                       `West Coast Haida Gwaii` = '3', `West Coast Van Island` = '4')
  #   )%>%
  #   group_by(region, censored_txt_95) %>%
  #   ggplot(aes(x=N_it_spiny, group=censored_txt_95, colour=censored_txt_95)) +
  #   geom_density(trim=T) +
  #   #xlim(c(0,310)) +
  #   xlab('Catch Count Spiny Dogfish') +
  #   facet_grid(region~censored_txt_95, scales='free_y') +
  #   theme(legend.position = 0),

  layout = matrix(c(1:2),nrow=2,ncol=1)
)

# How much wider are the credible intervals on average vs CPUE-based?
Estimator_Properties %>%
  group_by(region, species) %>%
  mutate(tmp=100*((Uncertainty-(Uncertainty[model=='CPUE-based']))/(Uncertainty[model=='CPUE-based']))) %>%
  ungroup() %>%
  filter(model=='Censored Poisson') %>%
  summarise(mean(tmp))

Estimator_Properties %>%
  group_by(region, species) %>%
  mutate(tmp=100*((Uncertainty-(Uncertainty[model=='ICR-based']))/(Uncertainty[model=='ICR-based']))) %>%
  ungroup() %>%
  filter(model=='Censored Poisson') %>%
  summarise(mean(tmp))
