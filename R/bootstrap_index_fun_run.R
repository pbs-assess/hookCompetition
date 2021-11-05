#' Computes Regional and Coastwide Bootstrapped Relative Abundance Indices
#' See the vignette for details
#'
#' @param data a sf points object containing the IPHC data
#' @param species a character of the species (linked to the variables `N_it` in the dataframe)
#' @param survey_boundaries a sf polygons object containing the survey boundary definitions
#' @param R specifies the number of bootstrap samples desired
#' @param ncpus specifies the number of cores to use in parallel
#' @param type specifies the type of bootstrap intervals to use (default percentile)
#' @param return when TRUE it returns the indices
#' @param plot when TRUE it returns the plot.
#' @param ICR_adjust if FALSE/TRUE it computes a CPUE-based/instantaneous catch rate-adjusted index
#' @param preserve_inter_regional_differences if TRUE estimated inter-regional differences in mean abundance are shown at a cost of higher variance. Does not affect coastwide index.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import stats
bootstrap_index_fun_run <- function(
  species,
  survey_boundaries,
  ICR_adjust=F,
  data,
  R=1000,
  return=F,
  ncpus=1,
  type='perc',
  plot=T,
  preserve_inter_regional_differences = F
)
{
  species_boot <- bootstrap_index_fun(data=data[which(!is.na(data$region_INLA)),], species=species, R=R, ICR_adjust=ICR_adjust, ncpus = ncpus, type=type, return=return, plot=plot)

  for(k in 1:dim(survey_boundaries)[1])
  {
    species_boot <- rbind(
      species_boot,
      bootstrap_index_fun(data=data[which(data$region_INLA==k),],
                          species=species, ICR_adjust=ICR_adjust, R=R, ncpus = ncpus, type=type, return=return, plot=plot, subregion = as.character(survey_boundaries$Region)[k], preserve_inter_regional_differences = preserve_inter_regional_differences)
    )
  }

  if(plot)
  {
    ind_plot <-
      ggplot2::ggplot(species_boot, ggplot2::aes(x=.data$year, y=.data$mean, ymax=.data$q0.975, ymin=.data$q0.025)) + ggplot2::facet_grid(~.data$region) +
      ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::ylab('Relative catch rate index') + ggplot2::ggtitle(paste0('Bootstrapped ',ifelse(ICR_adjust, 'ICR-adjusted ','') ,'Relative index ',species))
    print(
      ind_plot
    )
  }

  return(list(
    species_boot=species_boot,
    ind_plot = ind_plot,
    preserve_inter_regional_differences = preserve_inter_regional_differences
  ))

}
