#' load the R objects required for fitting spatio-temporal models.
#' These include: a delauney triangulation mesh for computing the GMRF
#' approximation to the SPDE (Lindgren et al 2011) within r-INLA (Rue et al 2005);
#' an INLA spde object with user-specified pc priors specified on the range and
#' standard deviation parameters.
#'
#' @param data a sf points object containing the IPHC data
#' @param survey_boundaries a sf polygons object containing the survey boundary definitions
#' @param mesh_coast a simplified (smoothed) sf polygons object over which the triangulation mesh will be created. If NULL, a simplified shapefile of BC Coastline will be returned.
#' @param hires_coast a high-resolution sf polygons object used for plotting. If NULL, a default BC coastline object will be returned.
#' @param prior_range Choose a value in kilometers such that the species density at two locations separated by distances greater than this will be `approximately' independent (typically a lower bound is chosen).
#' @param prior_range_prob What is you prior probability that the spatial range will be smaller than this value (typically a small probability)?
#' @param prior_sigma Choose a value for the marginal standard deviation of the random field. Typically an upper bound is used. Remember our linear predictor is on the log scale, so think about constraining exp(+- 2*SD) in terms of multiplicative regional differences.
#' @param prior_sigma_prob What is your prior probability that the true marginal standard deviation exceeds this?
#' @param nx How many columns of pixels to place over the survey_boundaries for plotting? Default 150, higher means smoother plots, but longer run times.
#' @param ny How many rows of pixels to place over the survey_boundaries for plotting? Default 150, higher means smoother plots, but longer run times.
#' @export
#' @importFrom magrittr %>%
#' @import stats
make_spatial_objects <- function(data, survey_boundaries=NULL, mesh_coast=NULL, hires_coast=NULL, prior_range=20, prior_range_prob=0.01, prior_sigma=2, prior_sigma_prob=0.01, nx = 150, ny = 150)
{
  # prior_range What value of range do we want to specify (in kilometres) for Gaussian random fields
  # prior_range_prob What is the prior probability that the range will be less than the specified range?
  # prior_sigma What value of standard deviation to we want to specify for Gaussian random fields?
  # prior_sigma_prob What is the prior probability that the SD will be greater than the specified SD?

  KM_CRS <- sf::st_crs(data)

  # Load the simplified shapefile of the coastline and plot
  if(is.null(mesh_coast))
  {
    Coast_simp <- sf::st_as_sf(hookCompetition::Coast_simp)
    # Coast_simp <- sp::spTransform(Coast_simp,KM_CRS)
    # Coast_simp <- sp::SpatialPolygons(
    #   Srl=lapply(Coast_simp@polygons[[1]]@Polygons,
    #              FUN = function(x){
    #                sp::Polygons(srl=list(x),ID=sample.int(1e9,size=1))
    #              }),
    #   proj4string = KM_CRS)
    #
    # Coast_simp <- Coast_simp[which(rgeos::gArea(Coast_simp, byid=T) > 100),]
    mesh_coast <- Coast_simp
  }

  if(is.null(hires_coast))
  {
    # Load the hi-resolution coastline
    # Coast_hires <- hookCompetition::Coast_hires
    # Coast_hires <- sp::spTransform(Coast_hires, KM_CRS)
    # Coast_hires <- rgeos::gSimplify(Coast_hires, 0.4)
    Coast_hires <- sf::st_as_sf(hookCompetition::Coastline)
  }

  #Now we define the computational mesh used to approximate the Gaussian random field. Next, we identify the indices of the triangles that lie on land. These indices will be used for implementing the barrier model.

  # define mesh
  mesh <- INLA::inla.mesh.2d(loc=sf::st_coordinates(data),
                       interior = sf::as_Spatial(mesh_coast[1,]),
                       max.edge = c(80),
                       min.angle = 21,
                       cutoff = c(20),
                       loc.domain = cbind(c(400,1080),c(1020,250)),
                       crs = KM_CRS)
  mesh$n #the number of triangle vertices - aim to keep < 300 for computation speed

  # Define the barrier model around polygons 4 and 9 (Vancouver Island and Haida Gwaii)
  # which triangle vertices are on land?
  # tl = length(mesh$graph$tv[,1])
  # # - the number of triangles in the mesh
  # posTri = matrix(0, tl, 2)
  # for (t in 1:tl){
  #   temp = mesh$loc[mesh$graph$tv[t, ], ]
  #   posTri[t,] = colMeans(temp)[c(1,2)]
  # }
  # posTri = sp::SpatialPoints(posTri, proj4string = KM_CRS)
  #
  # triangle_ind <- unlist(sp::over(mesh_coast[-1], posTri, returnList=T))
  triangle_ind <- NULL

  #Next, we define the barrier spde model using the function `inla.barrier.pcmatern()`. This defines an R-generic model that can now be handles in *inlabru* through the use of the `bru_mapper()` function. Then, for the purposes of exploratory analysis, we define a final count variable 'Other'. This will store the counts of caught fish other species.

  # Define the r-generic model
  # spde_mod <- INLA::inla.barrier.pcmatern(
  #   mesh = mesh,
  #   barrier.triangles = triangle_ind,
  #   prior.range = c(prior_range, prior_range_prob),
  #   prior.sigma = c(prior_sigma, prior_sigma_prob),
  #   range.fraction = 0.1)

  spde_mod <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(prior_range, prior_range_prob),
    prior.sigma = c(prior_sigma, prior_sigma_prob),
  )
  predict_pixels <- sf::st_as_sf(pixels_sf(mesh, mask=survey_boundaries, nx=nx, ny=ny))
  # extract centroids of the pixels
  # predict_points <-
  #   sp::SpatialPoints(
  #     coords = predict_pixels@coords,
  #     proj4string = predict_pixels@proj4string
  #   )

  predict_pixels$region_INLA <-
  sapply(sf::st_intersects(predict_pixels, survey_boundaries),
         function(z){ifelse(length(z)==0,NA_integer_,z[1])})

  # predict_pixels$region_INLA <-
  #   as.numeric(rgeos::gWithin(predict_points,
  #                      survey_boundaries,
  #                      byid=T, returnDense=F))
  predict_pixels$region_area <-
    as.numeric(sf::st_area(survey_boundaries, byid = T))[
      predict_pixels$region_INLA
    ]

  predict_pixels$region <-
    levels(survey_boundaries$Region)[predict_pixels$region_INLA]

  # predict_pixels$x <- predict_pixels@coords[,1]
  # predict_pixels$y <- predict_pixels@coords[,2]

return(list(
  spde_mod=spde_mod,
  mesh=mesh,
  triangle_ind=triangle_ind,
  Coast_hires=Coast_hires,
  Coast_simp=mesh_coast,
  predict_pixels=predict_pixels
))

}
