#' @export
bootstrap_index_fun <- function(data, species, ICR_adjust=F, R=1000, return=F, ncpus=1, type='perc', subregion=NULL, plot=T)
{

  if(ICR_adjust)
  {
    rel_boot_species_ind <-
    boot(data@data[complete.cases(data@data$prop_removed),],
         statistic = function(x,ind){
           ind_est <-
             as.numeric(by(x[ind,],x$year[ind],
                           FUN=function(z){
                             mean((as.numeric(as.matrix(z[,paste0('N_it_',species)])[,1])/z$effSkateIPHC)*
                                    comp_factor_fun(z$prop_removed,z$obsHooksPerSet), na.rm=T)
                           }))
           ind_est <- ind_est / gm_mean(ind_est)
         },
         R=R, strata=data@data[complete.cases(data@data$prop_removed),]$year,
         parallel = 'multicore', ncpus=ncpus)
  }
  if(!ICR_adjust)
  {
    rel_boot_species_ind <-
    boot(data@data,
         statistic = function(x,ind){
           ind_est <-
             as.numeric(by(x[ind,],x$year[ind],
                           FUN=function(z){
                             mean(as.numeric(as.matrix(z[,paste0('N_it_',species)])[,1])/z$effSkateIPHC, na.rm=T)
                           }))
           ind_est <- ind_est / gm_mean(ind_est)
         },
         R=R, strata=data$year,
         parallel = 'multicore', ncpus=ncpus)
  }

  rel_boot_species_CI <- as.data.frame(matrix(0 ,length(rel_boot_species_ind$t0), 7))
  colnames(rel_boot_species_CI) <- c('year','region', 'mean', 'sd', 'q0.025', 'median', 'q0.975')

  if(is.null(subregion))
  {
    rel_boot_species_CI[,2] <- 'All'
  }
  if(!is.null(subregion))
  {
    rel_boot_species_CI[,2] <- subregion
  }
  rel_boot_species_CI[,1] <- sort(unique(data@data$year),decreasing = F)

  for(i in 1:length(rel_boot_species_ind$t0))
  {
    #rel_boot_species_CI[i,1] <- rel_boot_species_ind$t0[i]
    rel_boot_species_CI[i,3] <- median(rel_boot_species_ind$t[,i])
    rel_boot_species_CI[i,4] <- sd(rel_boot_species_ind$t[,i])
    rel_boot_species_CI[i,6] <- median(rel_boot_species_ind$t[,i])

    if(type=='bca')
    {

      rel_boot_species_CI[i,c(5,7)] <-
        boot.ci(rel_boot_species_ind, type='bca',
                index = i,
                t0=rel_boot_species_ind$t0[i],
                t=rel_boot_species_ind$t[,i])$bca[1,c(4,5)]
    }
    if(type=='perc')
    {

      rel_boot_species_CI[i,c(5,7)] <-
          boot.ci(rel_boot_species_ind, type='perc',
                index = i,
                t0=rel_boot_species_ind$t0[i],
                t=rel_boot_species_ind$t[,i])$perc[1,c(4,5)]
    }
  }

  rel_boot_species_CI <- as.data.frame(rel_boot_species_CI)
  #rel_boot_species_CI$year <- sort(unique(data@data$year),decreasing = F)#c(1996:2011,2013:2019)
  #rel_boot_species_CI[rel_boot_species_CI$year==2005,1:4] <- 1

  # if(plot)
  # {
  #   print(ggplot(rel_boot_species_CI, aes(x=year, y=mean, ymax=q0.975, ymin=q0.025)) +
  #           geom_point() + geom_errorbar() + ylab('Relative catch rate index') + ggtitle(paste0('Relative index ',species,ifelse(is.null(subregion),' across the entire region',paste0(' in subregion ',subregion) ))))
  # }

  if(return)
  {
    return(rel_boot_species_CI)
  }

}

bootstrap_index_fun_run <- function(
  species_vec,
  survey_boundaries,
  ICR_adjust=F,
  data,
  R=1000,
  return=F,
  ncpus=1,
  type='perc',
  plot=T
)
{
  species_boot <- vector('list', length(species_vec))
  #species_boot_subregions <- vector('list', length(species_vec))
  j <- 1
  for(i in species_vec)
  {
    species_boot[[j]] <- bootstrap_index_fun(data=data[which(!is.na(data$region_INLA)),], species=i, R=R, ICR_adjust=ICR_adjust, ncpus = ncpus, type=type, return=return, plot=plot)
    #species_boot[[j]] <- traditional_index_fun(data=data[which(!is.na(data$region_INLA)),], species=i, R=R, ncpus = ncpus, type=type, return=return, plot=plot)
    #species_boot_subregions[[j]] <-  vector('list', length(survey_boundaries))
    for(k in 1:length(survey_boundaries))
    {
      species_boot[[j]] <- rbind(
        species_boot[[j]],
        bootstrap_index_fun(data=data[which(data$region_INLA==k),],
                                species=i, ICR_adjust=ICR_adjust, R=R, ncpus = ncpus, type=type, return=return, plot=plot, subregion = survey_boundaries[k,]@polygons[[1]]@ID)
        )

        # species_boot_subregions[[j]][[k]] <-
        #   traditional_index_fun(data=data[which(data$region_INLA==k),],
        #                         species=i, R=R, ncpus = ncpus, type=type, return=return, plot=plot, subregion = survey_boundaries[k,]@polygons[[1]]@ID)

    }

    if(plot)
    {
      print(ggplot(species_boot[[j]], aes(x=year, y=mean, ymax=q0.975, ymin=q0.025)) + facet_grid(~region) +
              geom_point() + geom_errorbar() + ylab('Relative catch rate index') + ggtitle(paste0('Bootstrapped ',ifelse(ICR_adjust, 'ICR-adjusted ','') ,'Relative index ',species_vec[j]))
      )
    }

    j <- j + 1
  }

  return(list(
    species_boot=species_boot
    #species_boot_subregions=species_boot_subregions
  ))

}

