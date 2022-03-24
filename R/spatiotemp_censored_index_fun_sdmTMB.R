#'  Compute Spatio-temporal Poisson-lognormal model-based indices of relative abundance using sdmTMB
#'  See the vignette for details
#'
#' @param data a sf points object containing the IPHC data
#' @param species a character name of the species (linked to the variables `N_it` in the dataframe)
#' @param survey_boundaries a sf polygons object containing the survey boundary definitions
#' @param M the number of independent Monte Carlo samples from the posterior used to compute indices
#' @param return logical stating whether or not to return indices as a data.frame
#' @param ICR_adjust A logical determining if the ICR-based scale-factor adjustment should be used
#' @param cprop the minimum proportion of baits needed to be removed to induce censorship. If >1 no censorship.
#' @param keep a logical stating if you want to return the sdmTMB model object. Useful for inspecting DIC, properties etc.,
#' @param use_upper_bound a logical stating if right-censored or interval-censored response (with concervative upper bound derived using Baranov Catch equation) is desired.
#' @param upper_bound_quantile a proportion stating which quantiles of observation to remove censorship from. If greater than 1, do not remove censorship from any observation
#' @param plot logical. Do you want to return a ggplot of the indices along with the data.frame?
#' @param allyears logical determining if upper bound quantile is computed uniquely for each year (if False) or over all years (if TRUE)
#' @param station_effects logical stating if IID station-station random effects wanted
#' @param seed the seed used for Monte Carlo sampling. Note that 0 means the result is non-reproducible, but the code can be significantly faster.
#' @param verbose logical. If TRUE, print INLA diagnostics on the console as the model fits. Can help to diagnose convergence issues.
#' @param n_trajectories integer specifying how many Monte Carlo sampled relative abundance indices to plot on a spaghetti plot. This can help with interpreting uncertainty. Suggested value 10.
#' @param preserve_inter_regional_differences Logical if TRUE estimated inter-regional differences in mean abundance are shown at a cost of higher variance. Does not affect coastwide index.
#' @param mesh an INLA::inla.mesh.2d object
#' @param spde An INLA inla.spde2 or inla.rgeneric object (e.g. `spde_mod` created in `make_spatial_objects()`)
#' @param pixels A sf points (pixels) object used for plotting the maps of relative abundance
#' @param covs Currently not implemented. Please let me know if (spatio-temporal) environmental covariate functionality is wanted in a future version.
#' @param spatiotemporal character string "IID", "AR1", "RW", or "off". Do you want to identify temporal changes in relative abundance uniquely across space? Almost always yes.
#' @param grf_priors a sdmTMBpriors() object specifying pcmatern priors on al grfs.
#' @param prev_fit a previous sdmTMB model fit (probably returned using keep=T). Speeds up model fitting dramatically
#' @param smoothing_spline Logical if TRUE (default) then relative coastwide abundance trend assumed to be smooth - perfect for long-living species with small home range. If FALSE, then a random walk temporal structure is assumed.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import stats
spatiotemp_censored_index_fun_sdmTMB <- function(data, survey_boundaries, species, M=1000, return=T, ICR_adjust=F, cprop=1.1, keep=F, use_upper_bound=FALSE, upper_bound_quantile=1, plot=T, allyears=F, station_effects=T, seed=0, verbose=F, n_trajectories=10, preserve_inter_regional_differences=F, mesh, sdpe, pixels, spatiotemporal='AR1', covs=NULL, grf_priors=sdmTMB::sdmTMBpriors(), prev_fit=NULL, smoothing_spline=T)
{
  # remove geometry features from data
  data$N_dat <- as.numeric(as.matrix(data[,c(paste0('N_it_',species))])[,1])
  # filter out NAs
  data <- data[which(!is.na(data$N_dat) & !is.na(data$effSkateIPHC) & !is.na(data$region_INLA) & !is.na(data$prop_removed)),]

  if(ICR_adjust & cprop<(1+1e-6))
  {
    stop('Combining the ICR_adjustment within a censored likelihood approach has not been tested')
  }

  if(ICR_adjust)
  {
    data$N_dat <- round(data$N_dat*comp_factor_fun(data$prop_removed,data$obsHooksPerSet))
  }
  data$year_INLA <- as.numeric(data$year - min(data$year) + 1)
  data$station_ID <- factor(data$station)
  data$event_ID <- factor(1:length(data$year))#as.numeric(as.factor(data$station))#1:length(data$year)#
  nyear <- length(unique(data$year_INLA))
  nyear_INLA <- max(data$year_INLA)
  min_year <- min(data$year, na.rm=T)
  max_year <- max(data$year, na.rm=T)
  nregion <- max(data$region_INLA, na.rm=T)

  ndata <- length(unique(data$event_ID))
  nstation <- length(unique(data$station_ID))
  data <- data[which(!is.na(data$effSkateIPHC)),]
  data$region_INLA <- factor(data$region_INLA)
  data$offset <- log(data$effSkateIPHC)

  # sdmTMB appears to need a mesh object, even when not in use
  mesh_tmp <- sdmTMB::make_mesh(data=data.frame(sf::st_coordinates(data)),
                                xy_cols = c('X','Y'),
                                mesh = mesh)

  data <- cbind(data,sf::st_coordinates(data))
  data <- sf::st_drop_geometry(data)

  # Define a mapping between the year_INLA value and the index of the data using pmatch
  # Remember 2012 data is missing, so we need to map subsequent years to account for skipped years
  year_map_fun <- function(ind)
  {
    return(pmatch(ind, sort(unique(data$year_INLA)), duplicates.ok = T))
  }

  quant_regions <- matrix(1e9, nrow=nregion, ncol=nyear)
  if(upper_bound_quantile<=1)
  {
    if(allyears)
    {
      quant_regions <-
        matrix( rep(
          as.numeric(by(data$N_dat,
                        data$region_INLA,
                        FUN = function(x){quantile(x,upper_bound_quantile, na.rm=T)}, simplify = T)),
          times=nyear),
          nrow=nregion, ncol=nyear)
    }
    if(!allyears)
    {
      quant_regions <-
        matrix(by(data$N_dat,
                  list(data$region_INLA,data$year_INLA),
                  FUN = function(x){quantile(x,upper_bound_quantile, na.rm=T)}, simplify = T), nrow=nregion, ncol=nyear)
    }
  }

  upper_bound <- rep(0, length(data$prop_removed))
  if(use_upper_bound)
  {
    scale_fac <- rep(0, length(data$prop_removed))
    scale_fac[data$prop_removed>cprop] <-
      comp_factor_fun(signif((data[data$prop_removed>cprop,]$prop_removed-cprop)/(1-cprop),5),
                      round((1-cprop)*data$obsHooksPerSet[data$prop_removed>cprop]))

    upper_bound[data$prop_removed>cprop] <- round(
      (data$prop_removed[data$prop_removed>cprop]-cprop)*data$obsHooksPerSet[data$prop_removed>cprop]*
        scale_fac[data$prop_removed>cprop])
  }

  data$low <- data$N_dat
  data$high <- data$N_dat

  if(use_upper_bound)
  {
    data$high[which(data$prop_removed >= cprop &
                      data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])] <-
      as.numeric(
        data$N_dat[which(data$prop_removed >= cprop &
                           data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])] +
          upper_bound[which(data$prop_removed >= cprop &
                              data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])]
      )
  }
  if(!use_upper_bound)
  {
    data$high[which(data$prop_removed >= cprop &
                      data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])] <-
      NA
  }

  if(smoothing_spline)
  {
    time_varying <- NULL
    if(station_effects)
    {
      # # define formulae
      formulae <- formula('N_dat ~ -1 +
    twentyhooks +
    s(year) + offset +
    (1 | event_ID) + (1 | station_ID)')
    }
    if(!station_effects)
    {
      # # define formulae
      formulae <- formula('N_dat ~ -1 +
    twentyhooks +
    s(year) + offset +
    (1 | event_ID)')
    }
  }
  if(!smoothing_spline)
  {
    time_varying <- formula('~ 1')
    if(station_effects)
    {
      # # define formulae
      formulae <- formula('N_dat ~ -1 +
    twentyhooks +
    offset +
    (1 | event_ID) + (1 | station_ID)')
    }
    if(!station_effects)
    {
      # # define formulae
      formulae <- formula('N_dat ~ -1 +
    twentyhooks +
    offset +
    (1 | event_ID)')
    }
  }

  # sdmTMB only fits the model at years present in the data. Fill in missing years
  missing_years <- NULL
  # Doesn't yet work
   if(sum(!(min_year:max_year %in% data$year)) > 0)
   {
     missing_years <- (min_year:max_year)[!(min_year:max_year %in% data$year)]
   }
  # Manually expand the data.frame to incorporate new weights
  data$weights <- 1
  if(!is.null(missing_years))
  {
    data <- rbind(data[,c('twentyhooks','offset','event_ID','station_ID','region_INLA',
                          'N_dat','weights','year','X','Y','low','high')],
                  expand.grid(twentyhooks=0,
                              offset=1,
                              event_ID=data$event_ID[1],
                              station_ID=data$station_ID[1],
                              region_INLA=unique(data$region_INLA),
                              N_dat=1,
                              weights=0,
                              year=missing_years,
                              X=data$X[1],
                              Y=data$Y[1],
                              low=1,
                              high=1))
    mesh_tmp <- sdmTMB::make_mesh(data=data,
                                  xy_cols = c('X','Y'),
                                  n_knots = 3)
    nyear <- length(unique(data$year))
  }

  # mod <- tryCatch(
  #   {
  #     tmp <- suppressWarnings(

      if(!(cprop <= 1 & !is.null(prev_fit)))
      {
        mod <-  sdmTMB::sdmTMB(formula=formulae,
                               data=data, time='year',
                               family = poisson(link = 'log'),
                               spatial = 'on', spatiotemporal = spatiotemporal,
                               mesh = mesh_tmp,
                               time_varying = time_varying,
                               silent = !verbose,
                               return_tmb_object = TRUE,
                               extra_time = missing_years,
                               priors=grf_priors,
                               previous_fit = prev_fit,
                               share_range = F)
      }
      if(cprop <= 1 & !is.null(prev_fit))
      {
        mod <- prev_fit
      }


      # )

      if(cprop <= 1)
      {
        # Use converged uncensored model as starting values.
        mod <- suppressWarnings(
          sdmTMB::sdmTMB(formula=formulae,
                      data=data, time='year',
                      family = sdmTMB::censored_poisson(link = 'log'),
                      spatial = 'on', spatiotemporal = spatiotemporal,
                      previous_fit = mod,
                      mesh = mesh_tmp,
                      time_varying = time_varying,
                      experimental = list(lwr=as.integer(data$low), upr=as.integer(data$high)),
                      silent = !verbose,
                      return_tmb_object = TRUE,
                      extra_time = missing_years,
                      priors=grf_priors,
                      share_range = F)
        )
       }
  #
  #
  #   },
  #   error=function(cond)
  #   {
  #     return(NULL)
  #   },
  #   finally={
  #
  #   }
  # )

  #browser()

  pred_df<-NULL
  trajectory_samples <- NULL
  trajectory_plot <- NULL
  index_plot <- NULL

  if(!is.null(mod) & !mod$bad_eig  & mod$sd_report$pdHess)
  {
    pred_dat <-
      inlabru::cprod(pixels, data.frame(year = sort(unique(data$year))))

    pred_dat <-
    pred_dat %>%
      dplyr::mutate(station_ID=unique(data$station_ID)[1],
                    event_ID=unique(data$event_ID)[1],
                    twentyhooks=0,
                    offset=1,
                    region_INLA=factor(region_INLA)
      )
    pred_dat <- cbind(pred_dat, sf::st_coordinates(pred_dat))

    pred_mod <- as.matrix(
      predict(mod, newdata=sf::st_drop_geometry(pred_dat), se_fit = T, re_form_iid = NA, re_form = NULL, sims=M)
    )

    samples_overall <-
        data.frame(value=as.numeric(apply(exp(pred_mod), 2, FUN = function(x){by(x, INDICES = pred_dat$year, FUN = mean)})),
                   year=rep(sort(unique(pred_dat$year)),times=M),
                   MC_ind=rep(1:M, each=nyear),
                   region='All') %>%
      dplyr::group_by(MC_ind) %>%
      dplyr::mutate(value = value / gm_mean(value))

    overall_df <-
      samples_overall %>%
      dplyr::group_by(year) %>%
      dplyr::summarize(mean = mean(value),
                       sd = sd(value),
                       q0.025 = quantile(value ,probs=0.025),
                       q0.975 = quantile(value, probs=0.975)
      )
    overall_df$region <- 'All'

    if(!preserve_inter_regional_differences)
    {
      samples_regional <-
        data.frame(value=as.numeric(apply(exp(pred_mod), 2, FUN = function(x){by(x, INDICES = list(pred_dat$year,pred_dat$region_INLA), FUN = mean)})),
                   year=rep(sort(unique(pred_dat$year)),times=M*nregion),
                   MC_ind=rep(1:M, each=nyear*nregion),
                   region=rep(rep(as.character(survey_boundaries$Region), each=nyear), times=M)) %>%
        dplyr::group_by(MC_ind, region) %>%
        dplyr::mutate(value = value / gm_mean(value))

      regional_df <-
        samples_regional %>%
        dplyr::ungroup(MC_ind) %>%
        dplyr::group_by(year,region) %>%
        dplyr::summarize(mean = mean(value),
                         sd = sd(value),
                         q0.025 = quantile(value, probs=0.025),
                         q0.975 = quantile(value, probs=0.975)
        )
    }
    if(preserve_inter_regional_differences)
    {
      samples_regional <-
        data.frame(value=as.numeric(apply(exp(pred_mod), 2, FUN = function(x){by(x, INDICES = list(pred_dat$year,pred_dat$region_INLA), FUN = mean)})),
                   year=rep(sort(unique(pred_dat$year)),times=M*nregion),
                   MC_ind=rep(1:M, each=nyear*nregion),
                   region=rep(rep(as.character(survey_boundaries$Region), each=nyear), times=M)) %>%
        dplyr::group_by(MC_ind) %>%
        dplyr::mutate(value = value / gm_mean(value))

      regional_df <-
        samples_regional %>%
        dplyr::ungroup(MC_ind) %>%
        dplyr::group_by(year,region) %>%
        dplyr::summarize(mean = mean(value),
                         sd = sd(value),
                         q0.025 = quantile(value, probs=0.025),
                         q0.975 = quantile(value, probs=0.975)
        )

    }

    pred_df <- rbind(overall_df[,c('year','region','mean','q0.025','q0.975','sd')],
                     regional_df[,c('year','region','mean','q0.025','q0.975','sd')])

    trajectory_samples <-
      rbind(samples_regional[,c('year','region','value','MC_ind')] %>%
              dplyr::filter(MC_ind <= n_trajectories),
            samples_overall[,c('year','region','value','MC_ind')] %>%
              dplyr::filter(MC_ind <= n_trajectories))

    # Spatial summaries for plotting
    pred_dat$mean <- apply(pred_mod, 1, mean)
    pred_dat$sd <- apply(pred_mod, 1, sd)
    pred_dat$q0.025 <- apply(pred_mod, 1, quantile, probs=0.025)
    pred_dat$q0.975 <- apply(pred_mod, 1, quantile, probs=0.975)
    pred_dat <- pred_dat[,c('X','Y','year','region','mean','sd','q0.025','q0.975')]

    if(plot)
    {
      index_plot <-
        ggplot2::ggplot(pred_df, ggplot2::aes(x=.data$year ,y=.data$mean, ymin=.data$q0.025, ymax=.data$q0.975)) +
        ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::ylab('Catch rate index') + ggplot2::facet_grid(~.data$region) +
        ggplot2::ggtitle(paste0('Spatio-temporal Overdispersed ',ifelse(ICR_adjust,'ICR-Adjusted', ifelse(cprop<(1+1e-6),'Censored','')) ,' Poisson Index ',species),
                         subtitle = ifelse(cprop<(1+1e-6),paste0('censorship proportion ',cprop, ', censorship from data in upper ',max(0,(1-upper_bound_quantile)*100), '% of values removed'),''))
      print(
        index_plot
      )
      if(n_trajectories>0)
      {
        trajectory_plot <-
          ggplot2::ggplot(trajectory_samples, ggplot2::aes(x=.data$year ,y=.data$value, colour=.data$MC_ind, group=.data$MC_ind)) +
          ggplot2::geom_line() + ggplot2::ylab('Catch rate index') + ggplot2::facet_grid(~.data$region) +
          ggplot2::ggtitle(paste0('Spatio-temporal Overdispersed ',ifelse(ICR_adjust,'ICR-Adjusted', ifelse(cprop<(1+1e-6),'Censored','')) ,' Poisson Index ',species),
                           subtitle = ifelse(cprop<(1+1e-6),paste0('censorship proportion ',cprop, ', censorship from data in upper ',max(0,(1-upper_bound_quantile)*100), '% of values removed'),'')) +
          viridis::scale_color_viridis()
        print(
          trajectory_plot
        )
      }
    }

  }

  if(is.null(pred_df) | !mod$sd_report$pdHess)
  {
    print('The model failed to converge. If attempting to fit a censored Poisson model, try increasing cprop, setting use_upper_bound equal TRUE, and/or decreasing upper_bound_quantile')
  }

  if(return & !is.null(pred_df) & mod$sd_report$pdHess)
  {
    if(keep == F)
    {
      return(list(pred_df_plot = pred_dat,
                  pred_spatiotemp=pred_df,
                  trajectory_samples=trajectory_samples,
                  index_plot = index_plot,
                  trajectory_plot = trajectory_plot,
                  preserve_inter_regional_differences=preserve_inter_regional_differences))#, pred_poisson=pred_df_poisson))
    }
    if(keep == T)
    {
      return(list(pred_df_plot = pred_dat,
                  mod=mod,
                  pred_spatiotemp=pred_df,
                  trajectory_samples=trajectory_samples,
                  index_plot = index_plot,
                  trajectory_plot = trajectory_plot,
                  preserve_inter_regional_differences=preserve_inter_regional_differences))#, pred_poisson=pred_df_poisson))
    }
  }

}
