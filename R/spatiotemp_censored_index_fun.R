#' Fit spatio-temporal indices using the GMRF-SPDE approximation of lindgren et al 2011
#'
#' @param data a SpatialPointsDataFrame object containing the IPHC data
#' @param species a character name of the species (linked to the variables `N_it` in the dataframe)
#' @param survey_boundaries a SpatialPolygonsDataFrame object containing the survey boundary definitions
#' @param M the number of independent Monte Carlo samples from the posterior used to compute indices
#' @param return logical stating whether or not to return indices as a data.frame
#' @param ICR_adjust A logical determining if the ICR-based scale-factor adjustment should be used
#' @param cprop the minimum proportion of baits needed to be removed to induce censorship. If >1 no censorship.
#' @param nthreads how many cores do you want to run in parallel when fitting the model?
#' @param keep a logical stating if you want to return the INLA model object. Useful for inspecting DIC, properties etc.,
#' @param use_upper_bound a logical stating if right-censored or interval-censored response (with concervative upper bound derived using Baranov Catch equation) is desired.
#' @param upper_bound_quantile a proportion stating which quantiles of observation to remove censorship from. If greater than 1, do not remove censorship from any observation
#' @param plot logical. Do you want to return a ggplot of the indices along with the data.frame?
#' @param allyears logical determining if upper bound quantile is computed uniquely for each year (if False) or over all years (if TRUE)
#' @param station_effects logical stating if IID station-station random effects wanted
#' @param prior_event a prior distribution for the event-event IID random effect standard deviation (see HT.prior for default half t_7(0,1) distribution)
#' @param prior_station a prior distribution for the station-station IID random effect standard deviation (see HT.prior for default half t_7(0,1) distribution)
#' @param init_vals initial values for fitting model. If NULL let INLA choose them.
#' @param n_knots the number of knots used for the temporal 'spline'. Default is 8. More knots means a 'wigglier' temporal effect is allowed, at a risk of higher variance.
#' @param seed the seed used for Monte Carlo sampling. Note that 0 means the result is non-reproducible, but the code can be significantly faster.
#' @param verbose logical. If TRUE, print INLA diagnostics on the console as the model fits. Can help to diagnose convergence issues.
#' @param n_trajectories integer specifying how many Monte Carlo sampled relative abundance indices to plot on a spaghetti plot. This can help with interpreting uncertainty. Suggested value 10.
#' @param preserve_inter_regional_differences Logical if TRUE estimated inter-regional differences in mean abundance are shown at a cost of higher variance. Does not affect coastwide index.
#' @param mesh An INLA inla.mesh.2d object (e.g. `mesh` created in `make_spatial_objects()`)
#' @param spde An INLA inla.spde2 or inla.rgeneric object (e.g. `spde_mod` created in `make_spatial_objects()`)
#' @param pixels A SpatialPixelsDataFrame object used for plotting the maps of relative abundance
#' @param covs Currently not implemented. Please let me know if (spatio-temporal) environmental covariate functionality is wanted in a future version.
#' @param spatiotemporal Logical. Do you want to identify temporal changes in relative abundance uniquely across space? Almost always yes.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import stats
spatiotemp_censored_index_fun <- function(data, survey_boundaries, species, M=1000, return=T, ICR_adjust=F, cprop=1.1, nthreads=1, keep=F, use_upper_bound=FALSE, upper_bound_quantile=1, plot=T, allyears=F, station_effects=T, prior_event=HT.prior(), prior_station=HT.prior(), init_vals = NULL, n_knots=8, seed=0, verbose=F, covs=NULL, spatiotemporal=F, mesh, spde, pixels, n_trajectories=10, preserve_inter_regional_differences=F)
{
  data$N_dat <- as.matrix(data[,c(paste0('N_it_',species))])[,1]
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
  data$station_ID <- as.numeric(as.factor(data$station))
  data$event_ID <- 1:length(data$year)#as.numeric(as.factor(data$station))#1:length(data$year)#
  nyear <- length(unique(data$year_INLA))
  nyear_INLA <- max(data$year_INLA)
  min_year <- min(data$year, na.rm=T)
  max_year <- max(data$year, na.rm=T)
  nregion <- max(data$region_INLA, na.rm=T)
  tmesh2 <- INLA::inla.mesh.1d(1:max(data$year_INLA))
  A_proj_t2 <- INLA::inla.spde.make.A(tmesh2, loc=data$year_INLA, repl=data$region_INLA, n.repl=nregion)
  t_index2 <- INLA::inla.spde.make.index('yearindex',max(data$year_INLA), n.repl=nregion)
  tmesh <- INLA::inla.mesh.1d(seq(from=1, to=nyear_INLA, length.out=n_knots))
  spde_t <- INLA::inla.spde2.pcmatern(tmesh, constr=T, prior.range = c(4,0.1), prior.sigma = c(1,0.1))
  A_t <- INLA::inla.spde.make.A(tmesh, loc=data$year_INLA, repl=data$region_INLA, n.repl = nregion)
  t_index <- INLA::inla.spde.make.index('yearind',n_knots, n.repl=nregion)

  # Spatial and Spatio-temporal stuff
  s_index <- INLA::inla.spde.make.index('spaceind',mesh$n, n.repl=1)
  A_s <- INLA::inla.spde.make.A(mesh=mesh, loc=sf::st_coordinates(data)[,1:2])
  st_index <- INLA::inla.spde.make.index('spacetimeind',mesh$n, n.group=n_knots)
  A_st <- INLA::inla.spde.make.A(mesh=mesh,
                           loc=sf::st_coordinates(data)[,1:2],
                           group = data$year_INLA,
                           group.mesh = tmesh,
                           n.group = n_knots)


  ndata <- length(unique(data$event_ID))
  nstation <- length(unique(data$station_ID))
  data <- data[which(!is.na(data$effSkateIPHC)),]

  cov_formula <- NULL

  if(!is.null(covs))
  {
    # build the covariate terms for INLA
    for(i in covs)
    {
      cov_formula <- paste0(cov_formula,' + ',i)
    }
  }

  GRF_formula <- '+ f(spaceind, model=spde_mod)'

  if(spatiotemporal==T)
  {
    GRF_formula <- '+ f(spaceind, model=spde_mod) +
                    f(spacetimeind, model=spde_mod,
                      group=spacetimeind.group,
                      control.group=list(model="ar1"))'
  }
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

  data$low <- rep(Inf,dim(data)[1])
  data$high <- rep(Inf,dim(data)[1])

  data$low[which(data$prop_removed >= cprop &
                   data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])] <- as.matrix(data[which(data$prop_removed >= cprop &
                                                                                                                                                    data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))]),
                                                                                                                                            c(paste0('N_it_',species))])[,1]

  if(use_upper_bound)
  {
    data$high[which(data$prop_removed >= cprop &
                      data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])] <-
      data$N_dat[which(data$prop_removed >= cprop &
                         data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])] +
      upper_bound[which(data$prop_removed >= cprop &
                          data$N_dat < quant_regions[cbind(data$region_INLA,year_map_fun(data$year_INLA))])]
  }

  MDAT <- INLA::inla.mdata(cbind(data$N_dat,
                           data$low,
                           data$high))

  year_names <- paste0("year",levels(data$year))

  if(station_effects)
  {
    # # define formulae
    formulae <- formula(paste0('inla.mdata(cbind(N_dat, low, high)) ~ -1 +
           Intercept +
           twentyhooks +
           f(event_ID, n=ndata, model="iid", constr=T, hyper=list(theta=list(prior= prior_event))) +
           f(station_ID, n=nstation, model="iid", constr=T, hyper=list(theta=list(prior=prior_station))) +
           f(yearind, model=spde_t)',
           cov_formula, GRF_formula))
    #f(year, model="ar1", replicate=region_INLA)')

  }
  if(!station_effects)
  {
    # # define formulae
    formulae <- formula(paste0('inla.mdata(cbind(N_dat, low, high)) ~ -1 +
           Intercept +
           twentyhooks +
           f(event_ID, n=ndata, model="iid", constr=T, hyper=list(theta=list(prior=prior_event))) +
           f(yearind, model=spde_t)',
           cov_formula, GRF_formula))
    #f(year, model="ar1", replicate=region_INLA)')

  }
  # remove geometry features from data
  data <- sf::st_drop_geometry(data)

  ind_var <- which(names(data) %in% c('low','high','N_dat'))

  stk <- INLA::inla.stack(
    data = data[,ind_var],
    A = list(A_t, A_s, A_st, 1, 1),
    effects = list(t_index, s_index, st_index, data[,-ind_var], Intercept=rep(1, times=dim(data)[1])),
    remove.unused = F)

  mod <- tryCatch(
    {
      INLA::inla(formula = formulae,
           family='cenpoisson2',
           data= INLA::inla.stack.data(stk),
           E=data$effSkateIPHC,
           control.predictor = list(A = INLA::inla.stack.A(stk), compute=T),
           control.compute=list(config = TRUE, dic=T),
           control.mode = list(theta=init_vals, restart=T),
           #control.inla = list(adapt.hessian.mode=F),
           num.threads = paste0(nthreads,":1"),
           verbose = verbose)
    },
    error=function(cond)
    {
      return(NULL)
    },

    finally={

    }
  )

  if(!is.null(mod) & mod$misc$mode.status==1000)
  {
    print('Numerical approximation not sufficiently accurate - refitting model at previous mode')
    mod <- inla.rerun(mod)
  }

  pred_df<-NULL
  trajectory_samples <- NULL
  index_plot <- NULL
  trajectory_plot <- NULL

  if(!is.null(mod) & !is.nan(mod$dic$dic) & mod$misc$mode.status!=1000)
  {
    df_tmp <- inlabru::cprod(data.frame(year=(c(1:nyear_INLA))),cbind(pixels, data.frame(sf::st_coordinates(pixels))))

    if(seed==0)
    {
      pred_df_samp <- INLA::inla.posterior.sample(result=mod, seed=seed,
                                            n=M, num.threads=paste0(nthreads,":1"))
    }
    if(seed!=0)
    {
      pred_df_samp <- INLA::inla.posterior.sample(result=mod, seed=seed,
                                            n=M, num.threads='1:1')
    }

    A_pred_t <- INLA::inla.spde.make.A(tmesh, loc = df_tmp$year)
    pred_df2 <- A_pred_t %*% INLA::inla.posterior.sample.eval(
      fun='yearind',
      samples=pred_df_samp,
      parallel=ifelse(nthreads==1, FALSE, TRUE))

    # No need to add spatial GRF as it is static across time so will not impact relative abundance indices
    # A_pred_s <- inla.spde.make.A(mesh, loc = cbind(df_tmp$x, df_tmp$y))
    # pred_df2 <- pred_df2 +
    #   A_pred_s %*% inla.posterior.sample.eval(
    #   fun='c(spaceind)',
    #   samples=pred_df_samp)

    if(spatiotemporal)
    {
      # Add spatiotemporal GRF
      A_pred_st <- INLA::inla.spde.make.A(mesh, loc = cbind(df_tmp$X,df_tmp$Y),
                                    group = df_tmp$year, group.mesh = tmesh)
      pred_df2 <- pred_df2 +
        A_pred_st %*% INLA::inla.posterior.sample.eval(
          fun='spacetimeind',
          samples=pred_df_samp,
          parallel=ifelse(nthreads==1, FALSE, TRUE))
    }

    pred_df2 <- exp(pred_df2)
    ## We no longer need to scale by the geometric mean due to the constr=T argument in INLA
    ## To see this, run the code up to here (using browser()) and see that geom_mean_regions ~ 1

    # estimate the geometric mean
     # geom_mean_regions <- inlabru:::bru_summarise(apply(pred_df2, 2, FUN=function(x){
     #   as.numeric(by(x, INDICES = factor(df_tmp$region_INLA), FUN=gm_mean))}))$mean

    # geom_mean_regions <- apply(pred_df2, 2, FUN=function(x){
    #    as.numeric(by(x, INDICES = factor(df_tmp$region_INLA), FUN=gm_mean))})

    # create dataframe for plotting
    pred_df_plot <- cbind(df_tmp,
                          bru_summarise(pred_df2) )
    # update the standard deviation column to avoid mesh aliasing
    pred_df_plot$sd <- as.numeric(sqrt((A_pred_st %*% mod$summary.random$spacetimeind$sd^2) +
                       (A_pred_t %*% mod$summary.random$yearind$sd^2)))
    pred_df_plot <- sf::st_as_sf(
      pred_df_plot,
      coords=c('x','y'),
      crs = sf::st_crs(data)
    )
    pred_df_plot$year <- pred_df_plot$year + min_year - 1

    # Computed coastwide average by computing approximate integral of intensity surface
    samples_overall <- apply(pred_df2, 2, FUN=function(x){
      as.numeric(by(x, INDICES = factor(df_tmp$year), FUN=mean))})

     geom_mean_allregions <- apply(samples_overall, 2, FUN=function(x){
       as.numeric(gm_mean(x))})

    if(!preserve_inter_regional_differences)
    {
      # reweight the samples of overall abundance by geometric mean
      samples_overall <- samples_overall / geom_mean_allregions
    }

    pred_df <- cbind(year=c(min_year:max_year), weight=1,
                     region='All',bru_summarise(samples_overall))

    # reweight the samples of regional abundance by regional geometric mean
    #pred_df2 <- pred_df2 / geom_mean_regions[df_tmp$region_INLA]

    samples_regional <- apply(pred_df2, 2, FUN=function(x){
      as.numeric(by(x, INDICES = list(factor(df_tmp$region_INLA), factor(df_tmp$year)), FUN=mean))})

    geom_mean_regions <- apply(samples_regional, 2, FUN=function(x){
          as.numeric(by(x, INDICES = factor(rep(1:nregion, times=nyear_INLA)), FUN=gm_mean))})

    if(!preserve_inter_regional_differences)
    {
      samples_regional <-
        samples_regional / geom_mean_regions[rep(1:nregion, times=nyear_INLA),]
    }

    pred_df2 <- cbind(inlabru::cprod(data.frame(region=levels(survey_boundaries$Region)),data.frame(year=c(min_year:max_year))),
                      bru_summarise(samples_regional) )

    pred_df <- rbind(pred_df, pred_df2)

    # need to remove the years 1996 and 1997 for WCVI as no data available
    # pred_df <-
    #   pred_df[!(pred_df$year %in% c(1996,1997) &
    #               pred_df$region %in% c('WCVI','All')),]

    if(n_trajectories>0)
    {
      trajectory_samples_regions <- samples_regional[,1:n_trajectories]
      trajectory_samples_All <- samples_overall[,1:n_trajectories]

      trajectory_samples <- rbind(
        data.frame(year=rep(c(min_year:max_year),times=n_trajectories),
                   MC_ind = rep(1:n_trajectories, each=nyear_INLA),
                   weight=1, region='All',
                   mean = as.numeric(trajectory_samples_All)),

        data.frame(
          inlabru::cprod(
                data.frame(region=levels(survey_boundaries$Region)),
                data.frame(year=rep(c(min_year:max_year),times=n_trajectories))
                ),
          MC_ind = rep(1:n_trajectories, each=nyear_INLA*nregion),
          mean = as.numeric(trajectory_samples_regions))
      )

    }

    if(plot)
    {
      index_plot <-
        ggplot2::ggplot(pred_df, ggplot2::aes(x=.data$year ,y=.data$mean, ymin=.data$q0.025, ymax=.data$q0.975)) +
        ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::ylab('Catch rate index') + ggplot2::facet_grid(~.data$region) +
        ggplot2::ggtitle(paste0('Spatiotemporal Overdispersed ',ifelse(ICR_adjust,'ICR-Adjusted', ifelse(cprop<(1+1e-6),'Censored','')) ,' Poisson Index ',species),
                subtitle = ifelse(cprop<(1+1e-6),paste0('censorship proportion ',cprop, ', censorship from data in upper ',max(0,(1-upper_bound_quantile)*100), '% of values removed'),''))
      print(
        index_plot
      )
      if(n_trajectories>0)
      {
        trajectory_plot <-
          ggplot2::ggplot(trajectory_samples, ggplot2::aes(x=.data$year ,y=.data$mean, colour=.data$MC_ind, group=.data$MC_ind)) +
          ggplot2::geom_line() + ggplot2::ylab('Catch rate index') + ggplot2::facet_grid(~.data$region) +
          ggplot2::ggtitle(paste0('Overdispersed ',ifelse(ICR_adjust,'ICR-Adjusted', ifelse(cprop<(1+1e-6),'Censored','')) ,' Poisson Index ',species),
                  subtitle = ifelse(cprop<(1+1e-6),paste0('censorship proportion ',cprop, ', censorship from data in upper ',max(0,(1-upper_bound_quantile)*100), '% of values removed'),'')) +
          viridis::scale_color_viridis()
        print(
          trajectory_plot
        )
      }
    }

  }

  if(is.null(pred_df))
  {
    print('The model failed to converge. If attempting to fit a censored Poisson model, try increasing cprop, setting use_upper_bound equal TRUE, and/or decreasing upper_bound_quantile')
  }

  if(return & !is.null(pred_df))
  {
    if(keep == F)
    {
      return(list(pred_spatiotemp=pred_df,
                  trajectory_samples=trajectory_samples,
                  pred_df_plot=pred_df_plot,
                  posterior_mean = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){INLA::inla.rmarginal(marginal=x, n=1000)}),
                  posterior_sample = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){rpois(lambda=INLA::inla.rmarginal(marginal=x, n=1000), n=1000)}),
                  index_plot = index_plot,
                  trajectory_plot = trajectory_plot,
                  preserve_inter_regional_differences=preserve_inter_regional_differences))#, pred_poisson=pred_df_poisson))
    }
    if(keep == T)
    {
      return(list(mod=mod,
                  pred_spatiotemp=pred_df,
                  trajectory_samples=trajectory_samples,
                  pred_df_plot=pred_df_plot,
                  posterior_mean = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){INLA::inla.rmarginal(marginal=x, n=1000)}),
                  posterior_sample = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){rpois(lambda=INLA::inla.rmarginal(marginal=x, n=1000), n=1000)}),
                  index_plot = index_plot,
                  trajectory_plot = trajectory_plot,
                  preserve_inter_regional_differences=preserve_inter_regional_differences))#, pred_poisson=pred_df_poisson))
    }
  }

}
