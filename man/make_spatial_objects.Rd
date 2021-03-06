% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_spatial_objects.R
\name{make_spatial_objects}
\alias{make_spatial_objects}
\title{load the R objects required for fitting spatio-temporal models.
These include: a delauney triangulation mesh for computing the GMRF
approximation to the SPDE (Lindgren et al 2011) within r-INLA (Rue et al 2005);
an INLA spde object with user-specified pc priors specified on the range and
standard deviation parameters.}
\usage{
make_spatial_objects(
  data,
  survey_boundaries = NULL,
  mesh_coast = NULL,
  hires_coast = NULL,
  prior_range = 20,
  prior_range_prob = 0.01,
  prior_sigma = 2,
  prior_sigma_prob = 0.01,
  nx = 150,
  ny = 150
)
}
\arguments{
\item{data}{a sf points object containing the IPHC data}

\item{survey_boundaries}{a sf polygons object containing the survey boundary definitions}

\item{mesh_coast}{a simplified (smoothed) sf polygons object over which the triangulation mesh will be created. If NULL, a simplified shapefile of BC Coastline will be returned.}

\item{hires_coast}{a high-resolution sf polygons object used for plotting. If NULL, a default BC coastline object will be returned.}

\item{prior_range}{Choose a value in kilometers such that the species density at two locations separated by distances greater than this will be `approximately' independent (typically a lower bound is chosen).}

\item{prior_range_prob}{What is you prior probability that the spatial range will be smaller than this value (typically a small probability)?}

\item{prior_sigma}{Choose a value for the marginal standard deviation of the random field. Typically an upper bound is used. Remember our linear predictor is on the log scale, so think about constraining exp(+- 2*SD) in terms of multiplicative regional differences.}

\item{prior_sigma_prob}{What is your prior probability that the true marginal standard deviation exceeds this?}

\item{nx}{How many columns of pixels to place over the survey_boundaries for plotting? Default 150, higher means smoother plots, but longer run times.}

\item{ny}{How many rows of pixels to place over the survey_boundaries for plotting? Default 150, higher means smoother plots, but longer run times.}
}
\description{
load the R objects required for fitting spatio-temporal models.
These include: a delauney triangulation mesh for computing the GMRF
approximation to the SPDE (Lindgren et al 2011) within r-INLA (Rue et al 2005);
an INLA spde object with user-specified pc priors specified on the range and
standard deviation parameters.
}
