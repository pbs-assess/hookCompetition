#' @export
spatiotemp_censored_index_fun <- function(data, survey_boundaries, species, M=1000, return=T, ICR_adjust=F, cprop=1.1, nthreads=1, keep=F, use_upper_bound=FALSE, upper_bound_quantile=1, plot=T, allyears=F, station_effects=T, prior_event=HT.prior, prior_station=HT.prior, init_vals = NULL, n_knots=8, seed=0, verbose=F, covs=NULL, spatiotemporal=F, mesh, spde, pixels)
{
  data$N_dat <- as.matrix(data@data[,c(paste0('N_it_',species))])[,1]
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
  tmesh2 <- inla.mesh.1d(1:max(data$year_INLA))
  A_proj_t2 <- inla.spde.make.A(tmesh2, loc=data$year_INLA, repl=data$region_INLA, n.repl=nregion)
  t_index2 <- inla.spde.make.index('yearindex',max(data$year_INLA), n.repl=nregion)
  tmesh <- inla.mesh.1d(seq(from=1, to=nyear_INLA, length.out=n_knots))
  spde_t <- inla.spde2.pcmatern(tmesh, constr=T, prior.range = c(4,0.1), prior.sigma = c(1,0.1))
  A_t <- inla.spde.make.A(tmesh, loc=data$year_INLA, repl=data$region_INLA, n.repl = nregion)
  t_index <- inla.spde.make.index('yearind',n_knots, n.repl=nregion)

  # Spatial and Spatio-temporal stuff
  s_index <- inla.spde.make.index('spaceind',mesh$n, n.repl=1)
  A_s <- inla.spde.make.A(mesh=mesh, loc=data@coords[,1:2])
  st_index <- inla.spde.make.index('spacetimeind',mesh$n, n.group=n_knots)
  A_st <- inla.spde.make.A(mesh=mesh,
                           loc=data@coords[,1:2],
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
          as.numeric(by(data@data$N_dat,
                        data@data$region_INLA,
                        FUN = function(x){quantile(x,upper_bound_quantile, na.rm=T)}, simplify = T)),
          times=nyear),
          nrow=nregion, ncol=nyear)
    }
    if(!allyears)
    {
      quant_regions <-
        matrix(by(data@data$N_dat,
                  list(data@data$region_INLA,data@data$year_INLA),
                  FUN = function(x){quantile(x,upper_bound_quantile, na.rm=T)}, simplify = T), nrow=nregion, ncol=nyear)
    }
  }

  upper_bound <- rep(0, length(data@data$prop_removed))
  if(use_upper_bound)
  {
    scale_fac <- rep(0, length(data@data$prop_removed))
    scale_fac[data@data$prop_removed>cprop] <-
      comp_factor_fun(signif((data@data[data@data$prop_removed>cprop,]$prop_removed-cprop)/(1-cprop),5),
                      round((1-cprop)*data$obsHooksPerSet[data@data$prop_removed>cprop]))

    upper_bound[data@data$prop_removed>cprop] <- round(
      (data@data$prop_removed[data@data$prop_removed>cprop]-cprop)*data$obsHooksPerSet[data@data$prop_removed>cprop]*
        scale_fac[data@data$prop_removed>cprop])
  }

  data$low <- rep(Inf,dim(data@data)[1])
  data$high <- rep(Inf,dim(data@data)[1])

  data$low[which(data@data$prop_removed >= cprop &
                   data@data$N_dat < quant_regions[cbind(data@data$region_INLA,year_map_fun(data@data$year_INLA))])] <- as.matrix(data@data[which(data@data$prop_removed >= cprop &
                                                                                                                                                    data@data$N_dat < quant_regions[cbind(data@data$region_INLA,year_map_fun(data@data$year_INLA))]),
                                                                                                                                            c(paste0('N_it_',species))])[,1]

  if(use_upper_bound)
  {
    data$high[which(data@data$prop_removed >= cprop &
                      data@data$N_dat < quant_regions[cbind(data@data$region_INLA,year_map_fun(data@data$year_INLA))])] <-
      data$N_dat[which(data@data$prop_removed >= cprop &
                         data@data$N_dat < quant_regions[cbind(data@data$region_INLA,year_map_fun(data@data$year_INLA))])] +
      upper_bound[which(data@data$prop_removed >= cprop &
                          data@data$N_dat < quant_regions[cbind(data@data$region_INLA,year_map_fun(data@data$year_INLA))])]
  }

  MDAT <- inla.mdata(cbind(data$N_dat,
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

  #browser()
  ind_var <- which(names(data@data) %in% c('low','high','N_dat'))

  stk <- inla.stack(
    data = data@data[,ind_var],
    A = list(A_t, A_s, A_st, 1, 1),
    effects = list(t_index, s_index, st_index, data@data[,-ind_var], Intercept=rep(1, times=dim(data)[1])),
    remove.unused = F)

  mod <- tryCatch(
    {
      inla(formula = formulae,
           family='cenpoisson2',
           data= inla.stack.data(stk),
           E=data$effSkateIPHC,
           control.predictor = list(A = inla.stack.A(stk), compute=T),
           control.compute=list(config = TRUE, dic=T),
           control.mode = list(theta=init_vals, restart=T),
           #control.inla = list(adapt.hessian.mode=F),
           num.threads = nthreads,
           verbose = verbose)
    },
    error=function(cond)
    {
      return(NULL)
    },

    finally={

    }
  )

  pred_df<-NULL

  #browser()

  if(!is.null(mod) & !is.nan(mod$dic$dic))
  {
    df_tmp <- cprod(data.frame(year=(c(1:nyear_INLA))),cbind(pixels@data, x=pixels$x, y=pixels$y))

    if(seed==0)
    {
      pred_df_samp <- inla.posterior.sample(result=mod, seed=seed,
                                            n=M, num.threads=nthreads)
    }
    if(seed!=0)
    {
      pred_df_samp <- inla.posterior.sample(result=mod, seed=seed,
                                            n=M, num.threads=1)
    }

    A_pred_t <- inla.spde.make.A(tmesh, loc = df_tmp$year)
    pred_df2 <- A_pred_t %*% inla.posterior.sample.eval(
      fun='c(yearind)',
      samples=pred_df_samp)

    # No need to add spatial GRF as it is static across time so will not impact relative abundance indices
    # A_pred_s <- inla.spde.make.A(mesh, loc = cbind(df_tmp$x, df_tmp$y))
    # pred_df2 <- pred_df2 +
    #   A_pred_s %*% inla.posterior.sample.eval(
    #   fun='c(spaceind)',
    #   samples=pred_df_samp)

    if(spatiotemporal)
    {
      # Add spatiotemporal GRF
      A_pred_st <- inla.spde.make.A(mesh, loc = cbind(df_tmp$x,df_tmp$y),
                                    group = df_tmp$year, group.mesh = tmesh)
      pred_df2 <- pred_df2 +
        A_pred_st %*% inla.posterior.sample.eval(
          fun='c(spacetimeind)',
          samples=pred_df_samp)
    }

    pred_df2 <- exp(pred_df2)
    ## We no longer need to scale by the geometric mean due to the constr=T argument in INLA
    ## To see this, run the code up to here (using browser()) and see that geom_mean_regions ~ 1

    # estimate the geometric mean
    # geom_mean_regions <- inlabru:::bru_summarise(apply(pred_df2, 2, FUN=function(x){
    #   as.numeric(by(x, INDICES = factor(df_tmp$region_INLA), FUN=gm_mean))}))$mean

    # create dataframe for plotting
    pred_df_plot <- cbind(df_tmp,
                          inlabru:::bru_summarise(pred_df2) )
    # update the standard deviation column to avoid mesh aliasing
    pred_df_plot$sd <- as.numeric(sqrt((A_pred_st %*% mod$summary.random$spacetimeind$sd^2) +
                       (A_pred_t %*% mod$summary.random$yearind$sd^2)))
    pred_df_plot <- SpatialPixelsDataFrame(
      points = cbind(pred_df_plot$x,pred_df_plot$y),
      data=pred_df_plot,
      proj4string = data@proj4string
    )
    pred_df_plot$year <- pred_df_plot$year + min_year - 1

    # Computed coastwide average by computing approximate integral of intensity surface
    samples_overall <- apply(pred_df2, 2, FUN=function(x){
      as.numeric(by(x, INDICES = factor(df_tmp$year), FUN=mean))})

    #geom_mean_allregions <- mean(apply(samples_overall, 2, FUN=function(x){
    #  as.numeric(gm_mean(x))}))

    # reweight the samples of overall abundance by geometric mean
    #samples_overall <- samples_overall / geom_mean_allregions

    pred_df <- cbind(year=c(min_year:max_year), weight=1,
                     region='All',inlabru:::bru_summarise(samples_overall))

    # reweight the samples of regional abundance by regional geometric mean
    #pred_df2 <- pred_df2 / geom_mean_regions[df_tmp$region_INLA]

    samples_regional <- apply(pred_df2, 2, FUN=function(x){
      as.numeric(by(x, INDICES = list(factor(df_tmp$region_INLA), factor(df_tmp$year)), FUN=mean))})

    pred_df2 <- cbind(cprod(data.frame(region=levels(survey_boundaries$Region)),data.frame(year=c(min_year:max_year))),
                      inlabru:::bru_summarise(samples_regional) )

    pred_df <- rbind(pred_df, pred_df2)
    # need to remove the years 1996 and 1997 for WCVI as no data available
    # pred_df <-
    #   pred_df[!(pred_df$year %in% c(1996,1997) &
    #               pred_df$region %in% c('WCVI','All')),]

    if(plot)
    {
      print(
        ggplot(pred_df, aes(x=year ,y=mean, ymin=q0.025, ymax=q0.975)) +
          geom_point() + geom_errorbar() + ylab('Catch rate index') + facet_grid(~region) +
          ggtitle(paste0('Spatiotemporal Overdispersed ',ifelse(ICR_adjust,'ICR-Adjusted', ifelse(cprop<(1+1e-6),'Censored','')) ,' Poisson Index ',species),
                  subtitle = ifelse(cprop<(1+1e-6),paste0('censorship proportion ',cprop, ', censorship from data in upper ',max(0,(1-upper_bound_quantile)*100), '% of values removed'),''))
      )
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
                  pred_df_plot=pred_df_plot,
                  posterior_mean = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){inla.rmarginal(marginal=x, n=1000)}),
                  posterior_sample = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){rpois(lambda=inla.rmarginal(marginal=x, n=1000), n=1000)})))#, pred_poisson=pred_df_poisson))
    }
    if(keep == T)
    {
      return(list(mod=mod,
                  pred_spatiotemp=pred_df,
                  pred_df_plot=pred_df_plot,
                  posterior_mean = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){inla.rmarginal(marginal=x, n=1000)}),
                  posterior_sample = sapply(mod$marginals.fitted.values[grepl(rownames(mod$summary.fitted.values), pattern = 'AP')],FUN = function(x){rpois(lambda=inla.rmarginal(marginal=x, n=1000), n=1000)})))#, pred_poisson=pred_df_poisson))
    }
  }

}

spatio_temporal_plot_fun <- function(
  pred_df_plot, year, variable, hires_COAST, plot_figure=T, return_figure=F
)
{
  plot <-
  ggplot() + gg(pred_df_plot[which(pred_df_plot$year==year), c(variable)]) + gg(hires_COAST)
  if(plot_figure)
  {
    print(plot)
  }
  if(return_figure)
  {
    return(plot)
  }
}
