#' Fit spatio-temporal indices using the GMRF-SPDE approximation of lindgren et al 2011
#' See the vignette for details
#'
#' @param pred_df_plot the sf object containing spatiotemporal predictions returned by `spatiotemp_censored_index_fun`
#' @param year integer specifying which year to plot
#' @param variable character specifying which variable to plot
#' @param hires_COAST a sf polygons object specifying the coastline used for plotting
#' @param plot_figure logical. Do you want to print the plot?
#' @param return_figure logical. Do you want to return the plot?
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import stats
spatio_temporal_plot_fun <- function(
  pred_df_plot, year, variable, hires_COAST, plot_figure=T, return_figure=F
)
{
  plot <-
    ggplot2::ggplot() +
    inlabru::gg(sf::as_Spatial(pred_df_plot)[which(pred_df_plot$year==year), c(variable)]) +
    inlabru::gg(sf::as_Spatial(hires_COAST)) +
    colsc(as.numeric(pred_df_plot[which(pred_df_plot$year==year), c(variable)]))
  if(plot_figure)
  {
    print(plot)
  }
  if(return_figure)
  {
    return(plot)
  }
}
