#' @export
make_spatial_objects <- function(data, survey_boundaries=NULL, mesh_coast=NULL, hires_coast=NULL, prior_range=20, prior_range_prob=0.01, prior_sigma=2, prior_sigma_prob=0.01)
{
  # prior_range What value of range do we want to specify (in kilometres) for Gaussian random fields
  # prior_range_prob What is the prior probability that the range will be less than the specified range?
  # prior_sigma What value of standard deviation to we want to specify for Gaussian random fields?
  # prior_sigma_prob What is the prior probability that the SD will be greater than the specified SD?

  KM_CRS <- data@proj4string

  # Load the simplified shapefile of the coastline and plot
  if(is.null(mesh_coast))
  {
    Coast_simp <- hookcompetition::Coast_simp
    Coast_simp <- spTransform(Coast_simp,KM_CRS)
    Coast_simp <- SpatialPolygons(
      Srl=lapply(Coast_simp@polygons[[1]]@Polygons,
                 FUN = function(x){
                   Polygons(srl=list(x),ID=sample.int(1e9,size=1))
                 }),
      proj4string = KM_CRS)

    Coast_simp <- Coast_simp[which(gArea(Coast_simp, byid=T) > 100),]
    mesh_coast <- Coast_simp
  }

  if(is.null(hires_coast))
  {
    # Load the hi-resolution coastline
    Coast_hires <- hookcompetition::Coast_hires
    Coast_hires <- spTransform(Coast_hires, KM_CRS)
    Coast_hires <- gSimplify(Coast_hires, 0.4)
    hires_coast <- Coast_hires
  }

  #Now we define the computational mesh used to approximate the Gaussian random field. Next, we identify the indices of the triangles that lie on land. These indices will be used for implementing the barrier model.

  # define mesh
  mesh <- inla.mesh.2d(loc=data@coords,
                       interior = mesh_coast[1],
                       max.edge = c(80),
                       min.angle = 21,
                       cutoff = c(20),
                       loc.domain = cbind(c(400,1080),c(1020,250)),
                       crs = KM_CRS)
  mesh$n #the number of triangle vertices - aim to keep < 300 for computation speed

  # Define the barrier model around polygons 4 and 9 (Vancouver Island and Haida Gwaii)
  # which triangle vertices are on land?
  tl = length(mesh$graph$tv[,1])
  # - the number of triangles in the mesh
  posTri = matrix(0, tl, 2)
  for (t in 1:tl){
    temp = mesh$loc[mesh$graph$tv[t, ], ]
    posTri[t,] = colMeans(temp)[c(1,2)]
  }
  posTri = SpatialPoints(posTri, proj4string = KM_CRS)

  triangle_ind <- unlist(over(mesh_coast[-1], posTri, returnList=T))

  #Next, we define the barrier spde model using the function `inla.barrier.pcmatern()`. This defines an R-generic model that can now be handles in *inlabru* through the use of the `bru_mapper()` function. Then, for the purposes of exploratory analysis, we define a final count variable 'Other'. This will store the counts of caught fish other species.

  # Define the r-generic model
  spde_mod <- inla.barrier.pcmatern(
    mesh = mesh,
    barrier.triangles = triangle_ind,
    prior.range = c(prior_range, prior_range_prob),
    prior.sigma = c(prior_sigma, prior_sigma_prob),
    range.fraction = 0.1)

  predict_pixels <- pixels(mesh, mask=survey_boundaries)
  predict_pixels$region_INLA <-
    as.numeric(gWithin(as(predict_pixels, 'SpatialPoints'),
                       survey_boundaries,
                       byid=T, returnDense=F))
  predict_pixels$region_area <-
    as.numeric(gArea(survey_boundaries, byid = T))[
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
  Coast_hires=hires_coast,
  Coast_simp=mesh_coast,
  predict_pixels=predict_pixels
))

}
