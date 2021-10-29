comp_factor_fun <- function(prop_hook, n_hooks)
{
  prop <- 1-prop_hook
  # if all hooks saturated - map to 1 hook
  prop[which(prop == 0)] <- 1 / n_hooks[which(prop == 0)]
  return(-log(prop)/(1-prop))
}

HT.prior = "expression:
  sigma = exp(-theta/2);
  nu = 7;
  log_dens = 0 - 0.5 * log(nu * pi) + (lgamma((nu+1)/2) - lgamma((nu)/2));
  log_dens = log_dens - 0.5 * (nu + 1) * log(1 + ((sigma * sigma)/nu));
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"

plotting_fun <-
  function(mod_lists,
           spatiotemporal = F,
           covs=NULL,
           cov_sp_list=NULL,
           times = 1,
           measure = 'mean',
           layout_mat = matrix(1, 1, 1),
           plot = T,
           n.samp=100,
           cov_plot=F,
           polys=NULL,
           fast=T){

    if(is.null(mod_lists))
    {
      return(NULL)
    }

    plot_pixels <- pixels(mesh, 150,150,mask=polys)

    plot_pixels_all <- cprod(plot_pixels,
                             data.frame(year_INLA = times))

    A_proj_plot <- inla.spde.make.A(
      mesh=mesh, loc=plot_pixels_all@coords
    )
    A_proj_st_plot <- inla.spde.make.A(
      mesh=mesh, loc=plot_pixels_all@coords,
      group=plot_pixels_all$year_INLA, group.mesh = mesh_time
    )

    cov_formula <- NULL
      if(!is.null(covs))
      {
        # build the covariate terms for inlabru
        for(i in covs)
        {
          cov_formula <- paste0(cov_formula,' + ',i,'*args$plot_pixels_all@data[,','"',i,'"',
                                ']')

          plot_pixels_all@data[,paste0(i)] <-
            over(plot_pixels_all, cov_sp_list[[i]])
        }
      }

    j <- 1
    plot_list <- vector('list', length=length(mod_lists))
    range_return <- vector()
    for (mod in mod_lists)
    {
      #plot_pixels_all2 <- plot_pixels_all
      if(is.null(mod))
      {
        plot_list[[j]] <- NULL
        next
      }
      if (spatiotemporal)
      {
        # for inlabru - FUTURE IMPLEMENTATION
        # RE_formula <- '~ Intercept_species +
        #     field_species + year_trend + stfield_species'
        RE_formula <- 'function(args){return(
        Intercept +
        as.numeric(args$A_proj_plot %*% as.numeric(field_species_idx)) +
        year_INLA[args$plot_pixels_all$year_INLA] +
        as.numeric(args$A_proj_st_plot %*% as.numeric(stfield_species_idx))'

          if(cov_plot)
          {
            RE_formula <- 'function(args){return(
        Intercept +'
          }

        if(fast)
        {
          all_pred1_samp <-
          inla.posterior.sample(n=n.samp,result=mod)
        }
        if(!fast)
        {
          all_pred1_samp <-
          inla.posterior.sample(n=n.samp,result=mod,seed=seed)
        }
        eval(parse(text=paste0('f <- ',
            RE_formula,
            cov_formula,')}')))

        all_pred1 <-
          inla.posterior.sample.eval(
            fun=f,
            samples = all_pred1_samp,
            args=list(
              plot_pixels_all=plot_pixels_all,
              A_proj_plot=A_proj_plot,
              A_proj_st_plot=A_proj_st_plot
          ))
        all_pred1 <- inlabru:::bru_summarise(all_pred1)

        plot_pixels_all2 <- plot_pixels_all
        plot_pixels_all2@data <- cbind(plot_pixels_all2@data,all_pred1)

        # for inlabru - FUTURE IMPLEMENTATION
        # all_pred1 <- predict(
        #   mod, seed=seed,
        #   plot_pixels_all,
        #   formula(paste0(
        #     RE_formula,
        #     cov_formula
        #   )),
        #   n.samples=n.samp
        # )
      }
      if (!spatiotemporal)
      {
        # for inlabru - FUTURE IMPLEMENTATION
        # RE_formula <- '~ Intercept_species +
        #     field_species + year_trend'

          RE_formula <- 'function(args){return(
        Intercept +
        as.numeric(args$A_proj_plot %*% as.numeric(field_species_idx)) +
        year_INLA[args$plot_pixels_all$year_INLA]'

          if(cov_plot)
          {
            RE_formula <- 'function(args){return(
        Intercept +'
          }

        all_pred1_samp <-
          inla.posterior.sample(n=n.samp,result=mod,seed=seed)

        eval(parse(text=paste0('f <- ',
            RE_formula,
            cov_formula,')}')))

        all_pred1 <-
          inla.posterior.sample.eval(
            fun=f,
            samples = all_pred1_samp,
            args=list(
              plot_pixels_all=plot_pixels_all,
              A_proj_plot=A_proj_plot
          )
          )
        all_pred1 <- inlabru:::bru_summarise(all_pred1)

        plot_pixels_all2 <- plot_pixels_all
        plot_pixels_all2@data <- cbind(plot_pixels_all2@data,all_pred1)

        # for inlabru - FUTURE IMPLEMENTATION
        # all_pred1 <- predict(mod, seed=seed,
        #                      plot_pixels_all,
        #                      formula(paste0(
        #     RE_formula,
        #     cov_formula
        #   )),
        #   n.samples=n.samp
        #   )
      }

      # for inlabru - FUTURE IMPLEMENTATION
      # range <- range(all_pred1$mean)
      # range_return <- range(range_return,range)
      # for (i in times)
      # {
      #   tmp <- all_pred1[which(all_pred1$year_INLA == i), measure]
      #   plot_list[[j]][[i]] <-
      #     ggplot() + gg(tmp) + gg(Coast_hires) + gg(survey_boundaries) + colsc(range)
      # }
      # if (plot)
      # {
      #   multiplot(plotlist = plot_list[[j]], layout = layout_mat)
      # }
      range <- range(plot_pixels_all2$mean)
      range_return <- range(range_return,range)
      for (i in times)
      {
        tmp <- plot_pixels_all2[which(plot_pixels_all2$year_INLA == i), measure]
        plot_list[[j]][[i]] <-
          ggplot() + gg(tmp) + gg(Coast_hires) + gg(survey_boundaries) + colsc(range)
      }
      if (plot)
      {
        multiplot(plotlist = plot_list[[j]], layout = layout_mat)
      }
      j <- j + 1
    }
    print(paste('range of plots across models is',range_return[1],range_return[2]))
    return(plot_list)
  }

cov_plotting_fun <-
  function(mod_lists,
           covs=NULL,
           plot_df=NULL,
           layout_mat = matrix(1, 1, 1),
           plot = T,
           log_scale=T,
           n.samp=100,
           fast=T){

    if(is.null(mod_lists))
    {
      return(NULL)
    }
    j <- 1
    plot_list <- vector('list', length=length(mod_lists))
    range_return <- vector()
    for (mod in mod_lists)
    {
      plot_df2 <- plot_df
      k <- 1
      cov_formula <- NULL
      if(!is.null(covs))
      {
        # build the covariate terms
        for(i in covs)
        {
          #inlabru - FUTURE IMPLEMENTATION
          #cov_formula <- paste0(cov_formula,' + ',i,'_eval(',i,')')
          if(k==1)
          {
            cov_formula <- paste0(i,'*args$plot_df[,"',i,'"]')
          }
           if(k!=1)
          {
          cov_formula <- paste0(cov_formula,' + ',i,'*args$plot_df[,"',i,'"]')
           }
          k <- k+1
        }
        if(log_scale==T)
        {
          #inlabru - FUTURE IMPLEMENTATION
          #cov_formula <- paste0('~ ',cov_formula)
          cov_formula <- paste0('function(args){return( ',cov_formula,')}')
        }
        if(log_scale==F)
        {
          #inlabru - FUTURE IMPLEMENTATION
          #cov_formula <- paste0('~ exp(',cov_formula,')')
          cov_formula <- paste0('function(args){return(exp( ',cov_formula,'))}')
        }
      }

      #inlabru   - FUTURE IMPLEMENTATION
      # all_pred1 <- predict(
      #     mod, seed=seed,
      #     plot_df,
      #     formula(cov_formula))

      if(fast)
      {
        all_pred1_samp <-
          inla.posterior.sample(n=n.samp,result=mod)
      }
      if(!fast)
      {
        all_pred1_samp <-
          inla.posterior.sample(n=n.samp,result=mod,seed=seed)
      }

        eval(parse(text=paste0('f <- ',cov_formula)))

        all_pred1 <-
          inla.posterior.sample.eval(
            fun=f,
            samples = all_pred1_samp,
            args=list(plot_df=plot_df)
          )
        all_pred1 <- inlabru:::bru_summarise(all_pred1)

      names(plot_df2)[1] <- 'x'

      plot_df2 <- cbind(plot_df2,all_pred1)

      range <- range(all_pred1[,c('q0.025','q0.975')])
      range_return <- range(range_return,range)

      plot_list[[j]] <-
        ggplot(plot_df2) +
  geom_line(aes(x, mean)) +
  geom_ribbon(aes(x,
                  ymin = q0.025,
                  ymax = q0.975),
              alpha = 0.2) +
  geom_ribbon(aes(x,
                  ymin = mean - 1 * sd,
                  ymax = mean + 1 * sd),
              alpha = 0.2) +
        xlab(covs[1])

      j <- j + 1
    }

    for(i in 1:length(plot_list))
    {
      plot_list[[i]] <- plot_list[[i]] + ylim(c(range_return[1],range_return[2]))
    }

    if (plot)
      {
        multiplot(plotlist = plot_list, layout = layout_mat)
      }

    return(plot_list)
  }

plotting_fun_diff <-
  function(mod_lists,
           spatiotemporal = F,
           covs=NULL,
           times = 1,
           measure = 'mean',
           plot = T,
           n.samp=100,
           cov_plot=F,
           log_scale=F,
           sensible_scale=T,
           polys=NULL){
    plot_pixels <- pixels(mesh, 150,150,mask=polys)
    plot_pixels_all <- cprod(plot_pixels,
                             data.frame(year_INLA = times))

    j <- 1
    all_pred1 <- vector('list', length=length(mod_lists))
    plot_list_return <- vector('list',length=length(times))
    names(plot_list_return) <- as.character(times + 2002)
    range_return <- matrix(NA, nrow=length(times), ncol=2)
    for (mod in mod_lists)
    {
      cov_formula <- NULL
      if(!is.null(covs))
      {
        # build the covariate terms for inlabru
        for(i in covs)
        {
          cov_formula <- paste0(cov_formula,' + ',i)
        }
        if(log_scale==F)
        {
          cov_formula <- paste0(cov_formula,')')
        }
      }
      if (spatiotemporal)
      {
        if(log_scale)
        {
          RE_formula <- '~ Intercept_species +
            field_species + year_trend + stfield_species'

          if(cov_plot)
          {
            RE_formula <- '~ '
          }
        }
        if(log_scale==F)
        {
          RE_formula <- '~ exp(Intercept_species +
            field_species + year_trend + stfield_species'

          if(cov_plot)
          {
            RE_formula <- '~ exp('
          }
        }


        # produce a matrix of Monte Carlo samples
        all_pred1[[j]] <- generate(
          mod,
          plot_pixels_all,
          formula(paste0(
            RE_formula,
            cov_formula
          )),
          n.samples=n.samp
        )
      }
      if (!spatiotemporal)
      {
        if(log_scale)
        {
          RE_formula <- '~ Intercept_species +
            field_species + year_trend'

          if(cov_plot)
          {
            RE_formula <- '~ '
          }
        }
        if(log_scale==F)
        {
          RE_formula <- '~ exp(Intercept_species +
            field_species + year_trend'

          if(cov_plot)
          {
            RE_formula <- '~ exp('
          }
        }

        all_pred1[[j]] <- generate(mod,
                             plot_pixels_all,
                             formula(paste0(
            RE_formula,
            cov_formula
          )),
          n.samples=n.samp
          )
      }

      j <- j + 1
    }
    i_3 <- 1
   for (i in times)
    {
        plot_list <- vector('list', length=(length(mod_lists)-1))

    for(i_1 in 1:(length(mod_lists)-1))
    {
      plot_list[[i_1]] <- vector('list', length = length(mod_lists)-i_1)

      for(i_2 in ((i_1+1):length(mod_lists)))
      {

        # compute all pairwise-differences across the years
         tmp <- all_pred1[[i_2]][which(plot_pixels_all$year_INLA == i), ] -
          all_pred1[[i_1]][which(plot_pixels_all$year_INLA == i), ]

        # Compute summary statistic of interest
         tmp <- apply(tmp, 1, measure)

        range <- range(tmp)
        range_return[i_3,] <- range(range_return[i_3,],range, na.rm=T)

        if(sensible_scale)
        {
          # threshold colour at 1% - 99% range of values
          range <- quantile(tmp, probs=c(0.01,0.99))
          range_return[i_3,] <- range(range_return[i_3,],range, na.rm=T)
        }

        tmp2 <- plot_pixels
        tmp2$measure <- tmp

        plot_list[[i_1]][[(i_2-1)]] <-
          ggplot() + gg(tmp2) + gg(Coast_hires) +
          gg(survey_boundaries) +colsc(range) +
          ggtitle(paste0('A plot of model ',i_2,' minus ',i_1, ' for year ',2002+i),
                  subtitle = paste0('Shown is the ',measure,' of the differences'))

      }
      }
      if (plot)
      {
     plot_list_combined <- do.call(c, plot_list)

     multiplot(plotlist = plot_list_combined,
               layout = matrix(1:((length(mod_lists)-1)^2),
                               (length(mod_lists)-1),(length(mod_lists)-1),
                               byrow = T))

      }

        plot_list_return[[i_3]] <- plot_list

   i_3 <- i_3 + 1
      }


    print(paste('range of plots across times is',range_return[,1],range_return[,2]))
    return(list(plot_list = plot_list_return, plot_ranges=range_return))
  }

index_generator_Alldata <- function(mod, covs, cov_sp_list, years, polys, n.samp=100, plot=T, species=NULL, spatiotemporal=F, fast=T)
{
  if(is.null(mod))
  {
    return(list(plot_list=NULL, index_regions=NULL, index_subregions=NULL))#,
                #rel_index_regions=NULL, rel_index_subregions=NULL))
  }
  n_regions <- length(polys)
  plot_list <- vector('list',n_regions)
  rel_plot_list <- vector('list',n_regions)

  # First create an index over the entire set of polygons
  # step 1 - generate spatialpixels
  ind_pixels <- pixels(mesh, 150,150, mask=polys)

  ind_pixels <- cprod(ind_pixels, data.frame(year_INLA=years))

  A_proj_plot <- inla.spde.make.A(
    mesh=mesh, loc=ind_pixels@coords
  )
  A_proj_st_plot <- inla.spde.make.A(
    mesh=mesh, loc=ind_pixels@coords,
    group=ind_pixels$year_INLA, group.mesh = mesh_time
  )
  A_proj_t_plot <- inla.spde.make.A(
    mesh=inla.mesh.1d(1:length(mod$summary.random$yearindex$mean)), loc=ind_pixels$year_INLA
  )

  # loop over the years and predict the average value
  # Remember - exponential scale is on count scale

  cov_formula <- NULL
  if(!is.null(covs))
  {
    # build the covariate terms
    for(i in covs)
    {
      cov_formula <- paste0(cov_formula,' + ',i,'*args$ind_pixels@data[,','"',i,'"',
                            ']')

      ind_pixels@data[,paste0(i)] <-
        over(ind_pixels, cov_sp_list[[i]])
    }
  }
  if (spatiotemporal)
  {
    # for inlabru - FUTURE IMPLEMENTATION
    # RE_formula <- '~ Intercept_species +
    #     field_species + year_trend + stfield_species'
    RE_formula <- 'function(args){return(exp(
        Intercept +
        as.numeric(args$A_proj_plot %*% as.numeric(field_species_idx)) +
        as.numeric(args$A_proj_t_plot %*% as.numeric(yearindex)) +
        as.numeric(args$A_proj_st_plot %*% as.numeric(stfield_species_idx))'
    #year_INLA[args$ind_pixels$year_INLA] +

    if(fast)
    {
      samples <-
        inla.posterior.sample(n=n.samp,result=mod)
    }
    if(!fast)
    {
      samples <-
        inla.posterior.sample(n=n.samp,result=mod,seed=seed)
    }

    eval(parse(text=paste0('f <- ',
                           RE_formula,
                           cov_formula,'))}')))

    samples <-
      inla.posterior.sample.eval(
        fun=f,
        samples = samples,
        args=list(
          A_proj_plot=A_proj_plot,
          A_proj_st_plot=A_proj_st_plot,
          A_proj_t_plot=A_proj_t_plot,
          ind_pixels=ind_pixels
        )
      )
  }
  if(!spatiotemporal)
  {
    # for inlabru - FUTURE IMPLEMENTATION
    # RE_formula <- '~ Intercept_species +
    #     field_species + year_trend + stfield_species'
    RE_formula <- 'function(args){return(exp(
        Intercept +
        as.numeric(args$A_proj_plot %*% as.numeric(field_species_idx)) +
        as.numeric(args$A_proj_t_plot %*% as.numeric(yearindex))'
        #year_INLA[args$ind_pixels$year_INLA]'

    if(fast)
    {
      samples <-
        inla.posterior.sample(n=n.samp,result=mod)
    }
    if(!fast)
    {
      samples <-
        inla.posterior.sample(n=n.samp,result=mod,seed=seed)
    }

    eval(parse(text=paste0('f <- ',
                           RE_formula,
                           cov_formula,'))}')))

    samples <-
      inla.posterior.sample.eval(
        fun=f,
        samples = samples,
        args=list(
          A_proj_plot=A_proj_plot,
          A_proj_t_plot=A_proj_t_plot,
          ind_pixels=ind_pixels
        )
      )
  }
  #browser()

  # compute geometric mean across all regions
  #gmmean_allregions <- apply(samples, 2, gm_mean)
  gmmean_allregions <- apply(
    as.matrix(do.call('cbind', by(samples, INDICES = factor(ind_pixels$year_INLA),
                                  FUN=function(x){
                                    aaply(as.matrix(x),2, mean,.drop = F)
                                  }, simplify = T))
    ), MARGIN = 1, FUN = gm_mean
  )
  #samples_subregions <- samples
  # scale the samples by the geometric mean
  #samples <- t(apply(samples, 1, FUN=function(x){x/gmmean_allregions}))

  # inlabru - FUTURE IMPLEMENTATION
  # cov_formula <- NULL
  #     if(!is.null(covs))
  #     {
  #       # build the covariate terms for inlabru
  #       for(i in covs)
  #       {
  #         cov_formula <- paste0(cov_formula,' + ',i)
  #       }
  #     }
  #
  # if(spatiotemporal==F)
  # {
  #   samples <-
  #   generate(mod, ind_pixels,
  #              formula(paste0(
  #              '~ exp(Intercept_species + field_species + year_trend',
  #              cov_formula,')')
  #              ),
  #            n.samples = n.samp)
  # }
  # if(spatiotemporal==T)
  # {
  #   samples <-
  #   generate(mod, ind_pixels,
  #            formula(paste0(
  #              '~ exp(Intercept_species +
  #              field_species + year_trend +
  #              stfield_species ',
  #              cov_formula,')')
  #              ),
  #            n.samples = n.samp)
  # }
  #

  for(i in years)
  {
    if(i == 1)
    {
      index_regions <-
        inlabru:::bru_summarise(
          t(aaply(samples[ind_pixels$year_INLA==i,],
                  2, mean,.drop = F))/
            gmmean_allregions)
      index_regions$year = i+1995
      index_regions$region = 'All'
    }
    if(i > 1)
    {
      tmp <- inlabru:::bru_summarise(
        t(aaply(samples[ind_pixels$year_INLA==i,],
                2, mean,.drop = F))/
          gmmean_allregions)
      tmp$year=i+1995
      tmp$region <- 'All'
      index_regions <- rbind(index_regions, tmp)
    }
  }

  index_allregions <- index_regions
  index_subregions <- NULL
  #rel_index_subregions <- NULL
  # If the polys object contain multiple boundaries, then compute indices for each
  regions_names=c('HS','QCS','WCHG','WCVI')
  if(n_regions > 1)
  {
    index_subregions <- vector('list',n_regions)
    #rel_index_subregions <- vector('list',n_regions)
    for(j in 1:n_regions)
    {
      #ind_pixels <- pixels(mesh, 150,150, mask=polys[j,])
      # find the indices corresponding to the jth subregion
      ind_j <- !is.na(over(ind_pixels, polys[j,]))

      g_mean_avg_subregion <- apply(
        as.matrix(do.call('cbind', by(samples[ind_j,], INDICES = factor(ind_pixels$year_INLA[ind_j]),
                                      FUN=function(x){
                                        aaply(as.matrix(x),2, mean,.drop = F)
                                      }, simplify = T))
        ), MARGIN = 1, FUN = gm_mean
      )

      # loop over the years and predict the average value
      # Remember - exponential scale is on count scale
      for(i in years)
      {

        if(i == 1)
        {
          index_subregions[[j]] <-
            inlabru:::bru_summarise(
              t(aaply(samples[ind_pixels$year_INLA==i & ind_j,],
                      2, mean,.drop = F))/
                g_mean_avg_subregion)
          index_subregions[[j]]$year = i+1995
          index_subregions[[j]]$region <- regions_names[j]
        }
        if(i > 1)
        {
          tmp <- inlabru:::bru_summarise(
            t(aaply(samples[ind_pixels$year_INLA==i & ind_j,],
                    2, mean,.drop = F))/
              g_mean_avg_subregion)
          tmp$year=i+1995
          tmp$region <- regions_names[j]
          index_subregions[[j]] <- rbind(index_subregions[[j]], tmp)
        }
      }

      index_allregions <- rbind(index_allregions, index_subregions[[j]])
      if(plot)
      {
        #  - Create non-relative indices
         plot_list[[j+1]] <- ggplot(index_subregions[[j]], aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
           geom_line() + geom_ribbon(alpha=0.4) + ylab(paste0('Index ',polys[j,]@polygons[[1]]@ID)) + ggtitle(paste0('Spatiotemporal index for species ',species))

        # rel_plot_list[[j+1]] <- ggplot(rel_index_subregions[[j]], aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
        #   geom_line() + geom_ribbon(alpha=0.4) + ylab(paste0('Relative index  ',polys[j,]@polygons[[1]]@ID)) + ggtitle(paste0('Spatiotemporal index for species ',species))
      }
    }

  }

  if(plot)
  {
     plot_list[[1]] <- ggplot(index_regions, aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
       geom_line() + geom_ribbon(alpha=0.4) + ylab('Index all survey boundaries') + ggtitle(paste0('Spatiotemporal index for species ',species))

    # rel_plot_list[[1]] <- ggplot(rel_index_regions, aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
    #   geom_line() + geom_ribbon(alpha=0.4) + ylab('Relative index all survey boundaries') + ggtitle(paste0('Spatiotemporal index for species ',species))
  }
  # return(list(plot_list=rel_plot_list,
  #             index_regions=rel_index_regions, index_subregions=rel_index_subregions))
   return(list(plot_list=plot_list, index_regions=index_allregions)) #rel_plot_list=rel_plot_list,
               #index_regions=index_regions, index_subregions=index_subregions))#, rel_index_regions=rel_index_regions, rel_index_subregions=rel_index_subregions))
}

multispecies_index_generator <- function(mod, covs, years, polys, n.samp=100, plot=T, species=NULL, species_ind=NULL, spatiotemporal=F, correction=T)
{
  n_regions <- length(polys)
  plot_list <- vector('list',n_regions)
  rel_plot_list <- vector('list',n_regions)

  # First create an index over the entire set of polygons
  # step 1 - generate spatialpixels
  ind_pixels <- pixels(mesh, 150,150, mask=polys)
  ind_pixels <- cprod(ind_pixels,
                      data.frame(year_INLA=years, species=species_ind))

  # loop over the years and predict the average value
  # Remember - exponential scale is on count scale

  cov_formula <- NULL
      if(!is.null(covs))
      {
        # build the covariate terms for inlabru
        for(i in covs)
        {
          cov_formula <- paste0(cov_formula,' + ',paste0(species,i))
        }
      }

  if(spatiotemporal==F)
  {
    # We subtract 1 to match mistake
    if(correction)
    {
      samples <-
    generate(mod, ind_pixels,
               formula(paste0(paste0('~ exp(-1 +',
                 paste0(c('Intercept_','field_'),
                 species, collapse = ' + '),
                ' + year_trend',
                cov_formula,collapse = " + "),')')),
             n.samples = n.samp)
    }
    if(!correction)
    {
      samples <-
    generate(mod, ind_pixels,
               formula(paste0(paste0('~ exp(',
                 paste0(c('Intercept_','field_'),
                 species, collapse = ' + '),
                ' + year_trend',
                cov_formula,collapse = " + "),')')),
             n.samples = n.samp)
    }
  }
  if(spatiotemporal==T)
  {
    if(correction)
    {
      samples <-
    generate(mod, ind_pixels,
               formula(paste0(paste0('~ exp(-1 +',
                 paste0(c('Intercept_','field_','stfield_'),
                 species, collapse = ' + '),
                ' + year_trend',
                cov_formula,collapse = " + "),')')),
             n.samples = n.samp)
    }
    if(!correction)
    {
      samples <-
    generate(mod, ind_pixels,
               formula(paste0(paste0('~ exp(',
                 paste0(c('Intercept_','field_','stfield_'),
                 species, collapse = ' + '),
                ' + year_trend',
                cov_formula,collapse = " + "),')')),
             n.samples = n.samp)
    }
  }

  for(i in years)
  {
    #ind_pixels$year_INLA <- i

    if(i == 1)
    {
      rel_index_regions <-
        inlabru:::bru_summarise(
      t(aaply(samples[ind_pixels$year_INLA %in% c(3,i),],
            2, .fun=function(x,year_ind){
              mean(x[year_ind==i])/mean(x[year_ind==3])
              },.drop = F,
            year_ind=ind_pixels$year_INLA[ind_pixels$year_INLA %in% c(3,i)])))
      rel_index_regions$year = i+2002

      # Generate non relative indices
    #   index_regions <-
    # inlabru:::bru_summarise(
    #   t(aaply(samples[ind_pixels$year_INLA==i,],
    #         2, mean,.drop = F)))
    #   index_regions$year = i+2002
    }
    if(i > 1)
    {
      rel_tmp <-
        inlabru:::bru_summarise(
      t(aaply(samples[ind_pixels$year_INLA %in% c(3,i),],
            2, .fun=function(x,year_ind){
              mean(x[year_ind==i])/mean(x[year_ind==3])
              },.drop = F,
            year_ind=ind_pixels$year_INLA[ind_pixels$year_INLA %in% c(3,i)])))
      rel_tmp$year = i+2002
      rel_index_regions <- rbind(rel_index_regions, rel_tmp)

       # Generate non relative indices
      # tmp <- inlabru:::bru_summarise(
      # t(aaply(samples[ind_pixels$year_INLA==i,],
      #       2, mean,.drop = F)))
      # tmp$year=i+2002
      # index_regions <- rbind(index_regions, tmp)
    }

  }

  index_subregions <- NULL
  # If the polys object contain multiple boundaries, then compute indices for each
  if(n_regions > 1)
  {
    index_subregions <- vector('list',n_regions)
    rel_index_subregions <- vector('list',n_regions)
    for(j in 1:n_regions)
    {
      #ind_pixels <- pixels(mesh, 150,150, mask=polys[j,])
      # find the indices corresponding to the jth subregion
      ind_j <- !is.na(over(ind_pixels, polys[j,]))

    # loop over the years and predict the average value
  # Remember - exponential scale is on count scale
  for(i in years)
  {
    if(i == 1)
    {
      rel_index_subregions[[j]] <-
        inlabru:::bru_summarise(
      t(aaply(samples[ind_pixels$year_INLA %in% c(3,i) & ind_j,],
            2, .fun=function(x,year_ind){
              mean(x[year_ind==i])/mean(x[year_ind==3])
              },.drop = F,
            year_ind=ind_pixels$year_INLA[ind_pixels$year_INLA %in% c(3,i) & ind_j])))
      rel_index_subregions[[j]]$year = i+2002

      # Generate non relative indices
    #   index_subregions[[j]] <-
    # inlabru:::bru_summarise(
    #   t(aaply(samples[ind_pixels$year_INLA==i & ind_j,],
    #         2, mean,.drop = F)))
    #   index_subregions[[j]]$year = i+2002
    }
    if(i > 1)
    {
      rel_tmp <-
        inlabru:::bru_summarise(
      t(aaply(samples[ind_pixels$year_INLA %in% c(3,i) & ind_j,],
            2, .fun=function(x,year_ind){
              mean(x[year_ind==i])/mean(x[year_ind==3])
              },.drop = F,
            year_ind=ind_pixels$year_INLA[ind_pixels$year_INLA %in% c(3,i) & ind_j])))
      rel_tmp$year = i+2002
      rel_index_subregions[[j]] <- rbind(rel_index_subregions[[j]], rel_tmp)
  }
      if(plot)
  {
    rel_plot_list[[j+1]] <- ggplot(rel_index_subregions[[j]], aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
      geom_line() + geom_ribbon(alpha=0.4) + ylab(paste0('Relative index  ',polys[j,]@polygons[[1]]@ID)) + ggtitle(paste0('Spatiotemporal index for species ',species),subtitle = 'Using multispecies model')
  }
  }

    }
  }
  if(plot)
  {
     # Generate non relative indices
    # plot_list[[1]] <- ggplot(index_regions, aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
    #   geom_line() + geom_ribbon(alpha=0.4) + ylab('Index all survey boundaries') + ggtitle(paste0('Spatiotemporal index for species ',species))

    rel_plot_list[[1]] <- ggplot(rel_index_regions, aes(x=year, y=mean, ymin=q0.025, ymax=q0.975)) +
      geom_line() + geom_ribbon(alpha=0.4) + ylab('Relative index all survey boundaries') + ggtitle(paste0('Spatiotemporal index for species ',species),subtitle = 'Using multispecies model')
  }
  return(list(plot_list=rel_plot_list,
              index_regions=rel_index_regions, index_subregions=rel_index_subregions))
}

index_comparison <- function(indices, subregion_indices, mod_names, species, polys)
{
  plot_list <- vector('list',length(indices[[1]]))

  for(i in 1:(1+length(subregion_indices[[1]])))
  {
    if(i==1)
    {
      indices[[1]]$mod <- mod_names[1]
    indices[[2]]$mod <- mod_names[2]

    combined_dataframe <- rbind(indices[[1]], indices[[2]])

    plot_list[[i]] <- ggplot(combined_dataframe, aes(x=year, y=mean, ymin=q0.025, ymax=q0.975, group=mod, color=mod, fill=mod)) +
      geom_line() + geom_ribbon(alpha=0.4) + ylab(paste0('Index across entire survey region')) + ggtitle(paste0('Spatiotemporal index for species ',species), subtitle = 'Using multispecies and single species models') + guides(color=NULL)
    }
     if(i>1)
    {
    subregion_indices[[1]][[i-1]]$mod <- mod_names[1]
    subregion_indices[[2]][[i-1]]$mod <- mod_names[2]

    combined_dataframe <- rbind(subregion_indices[[1]][[i-1]],
                                subregion_indices[[2]][[i-1]])

    plot_list[[i]] <- ggplot(combined_dataframe, aes(x=year, y=mean, ymin=q0.025, ymax=q0.975, group=mod, color=mod, fill=mod)) +
      geom_line() + geom_ribbon(alpha=0.4) + ylab(paste0('Index across survey boundary ',polys[i-1,]@polygons[[1]]@ID)) + ggtitle(paste0('Spatiotemporal index for species ',species), subtitle = 'Using multispecies and single species models') + guides(color=NULL)
    }

  }
  return(plot_list)
}

# Geometric Mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

index_comparison_quantitative <- function(true_ind, true_name='obs_mean', pred_ind, pred_pointest_name='mean', pred_LCL_name='q0.025', pred_UCL_name='q0.975')
{

  # First for the entire region
  ind_true <- na.omit(true_ind)
  # normalize
  #ind_true[,true_name] <- ind_true[,true_name]/max(ind_true[,true_name])

  ind_pred <- na.omit(pred_ind) # pre-normalized

  # subset ind_pred by the years seen in ind_true
  ind_pred <- ind_pred[which(ind_pred$year %in% unique(ind_true$year)),]

  results <- data.frame(bias=rep(0,length(unique(ind_pred$model))),
                        RMSE=rep(0,length(unique(ind_pred$model))),
                        coverage=rep(0,length(unique(ind_pred$model))),
                        model=unique(ind_pred$model),
                        #pearson_cor=rep(0,length(unique(ind_pred$model))),
                        spearman_cor=rep(0,length(unique(ind_pred$model))))
  for(i in 1:length(unique(ind_pred$model)))
  {
    #browser()
    ind_pred_sub <- ind_pred[ind_pred$model==unique(ind_pred$model)[i],]
    ind_true_sub <- ind_true[which(ind_true$year %in% unique(ind_pred_sub$year)),]
    # sort into matching order by year
    ind_true_sub <- ind_true_sub[pmatch(ind_pred_sub$year,ind_true_sub$year),]
    # compute bias, MSE, and coverage
    results$bias[i] <-
      mean(as.matrix(ind_pred_sub[,pred_pointest_name] - ind_true_sub[,true_name]))
    results$RMSE[i] <-
      sqrt(mean((as.matrix(ind_pred_sub[,pred_pointest_name] - ind_true_sub[,true_name]))^2))
    results$coverage[i] <-
      mean(ind_pred_sub[,pred_LCL_name] <= ind_true_sub[,true_name] &
             ind_pred_sub[,pred_UCL_name] >= ind_true_sub[,true_name])
    #results$pearson_cor[i] <- cor(as.matrix(ind_pred_sub[,pred_pointest_name]),
    #as.matrix(ind_true_sub[,true_name]), method='pearson')
    results$spearman_cor[i] <- cor(as.matrix(ind_pred_sub[,pred_pointest_name]),
                                   as.matrix(ind_true_sub[,true_name]), method='spearman')
  }
  return(results)
}

colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}
