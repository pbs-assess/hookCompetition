# The actual workhorse of the bootstrap_index_fun_run
bootstrap_index_fun <- function(data, species, ICR_adjust=F, R=1000, return=F, ncpus=1, type='perc', subregion=NULL, plot=T, preserve_inter_regional_differences = F)
{

  if(ICR_adjust)
  {
    if(preserve_inter_regional_differences)
    {
      rel_boot_species_ind <-
        boot::boot(data[complete.cases(data$prop_removed),],
             statistic = function(x,ind){
               ind_est <-
                 as.numeric(by(x[ind,],x$year[ind],
                               FUN=function(z){
                                 mean((as.numeric(as.matrix(z[,paste0('N_it_',species)])[,1])/z$effSkateIPHC)*
                                        comp_factor_fun(z$prop_removed,z$obsHooksPerSet), na.rm=T)
                               }))
             },
             R=R, strata=data[complete.cases(data$prop_removed),]$year,
             parallel = 'multicore', ncpus=ncpus)
    }
    if(!preserve_inter_regional_differences)
    {
      rel_boot_species_ind <-
        boot::boot(data[complete.cases(data$prop_removed),],
             statistic = function(x,ind){
               ind_est <-
                 as.numeric(by(x[ind,],x$year[ind],
                               FUN=function(z){
                                 mean((as.numeric(as.matrix(z[,paste0('N_it_',species)])[,1])/z$effSkateIPHC)*
                                        comp_factor_fun(z$prop_removed,z$obsHooksPerSet), na.rm=T)
                               }))
               ind_est <- ind_est / gm_mean(ind_est)
             },
             R=R, strata=data[complete.cases(data$prop_removed),]$year,
             parallel = 'multicore', ncpus=ncpus)
    }
  }
  if(!ICR_adjust)
  {
    if(preserve_inter_regional_differences)
    {
      rel_boot_species_ind <-
        boot::boot(data,
             statistic = function(x,ind){
               ind_est <-
                 as.numeric(by(x[ind,],x$year[ind],
                               FUN=function(z){
                                 mean(as.numeric(as.matrix(z[,paste0('N_it_',species)])[,1])/z$effSkateIPHC, na.rm=T)
                               }))
             },
             R=R, strata=data$year,
             parallel = 'multicore', ncpus=ncpus)
    }
    if(!preserve_inter_regional_differences)
    {
      rel_boot_species_ind <-
        boot::boot(data,
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
  rel_boot_species_CI[,1] <- sort(unique(data$year),decreasing = F)

  for(i in 1:length(rel_boot_species_ind$t0))
  {
    #rel_boot_species_CI[i,1] <- rel_boot_species_ind$t0[i]
    rel_boot_species_CI[i,3] <- median(rel_boot_species_ind$t[,i])
    rel_boot_species_CI[i,4] <- sd(rel_boot_species_ind$t[,i])
    rel_boot_species_CI[i,6] <- median(rel_boot_species_ind$t[,i])

    if(type=='bca')
    {

      rel_boot_species_CI[i,c(5,7)] <-
        boot::boot.ci(rel_boot_species_ind, type='bca',
                index = i,
                t0=rel_boot_species_ind$t0[i],
                t=rel_boot_species_ind$t[,i])$bca[1,c(4,5)]
    }
    if(type=='perc')
    {

      rel_boot_species_CI[i,c(5,7)] <-
        boot::boot.ci(rel_boot_species_ind, type='perc',
                index = i,
                t0=rel_boot_species_ind$t0[i],
                t=rel_boot_species_ind$t[,i])$perc[1,c(4,5)]
    }
  }

  rel_boot_species_CI <- as.data.frame(rel_boot_species_CI)
  #rel_boot_species_CI$year <- sort(unique(data$year),decreasing = F)#c(1996:2011,2013:2019)
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

