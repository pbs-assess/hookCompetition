# This computes the ICR-based hook-competition adjustment factor
comp_factor_fun <- function(prop_hook, n_hooks)
{
  prop <- 1-prop_hook
  # if all hooks saturated - map to 1 hook
  prop[which(prop == 0)] <- 1 / n_hooks[which(prop == 0)]
  return(-log(prop)/(1-prop))
}

# Geometric Mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Nice colour palette for map plots
colsc <- function(...) {
  ggplot2::scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}

# copied from inlabru author Finn Lindgren
bru_summarise <- function (data, x = NULL, cbind.only = FALSE)
{
  if (is.list(data)) {
    data <- do.call(cbind, data)
  }
  if (cbind.only) {
    smy <- data.frame(data)
    colnames(smy) <- paste0("sample.", 1:ncol(smy))
  }
  else {
    smy <- data.frame(apply(data, MARGIN = 1, mean, na.rm = TRUE),
                      apply(data, MARGIN = 1, sd, na.rm = TRUE), t(apply(data,
                                                                         MARGIN = 1, quantile, prob = c(0.025, 0.5, 0.975),
                                                                         na.rm = TRUE)), apply(data, MARGIN = 1, min,
                                                                                               na.rm = TRUE), apply(data, MARGIN = 1, max, na.rm = TRUE))
    colnames(smy) <- c("mean", "sd", "q0.025", "median",
                       "q0.975", "smin", "smax")
    smy$cv <- smy$sd/smy$mean
    smy$var <- smy$sd^2
  }
  if (!is.null(x)) {
    smy <- expand_to_dataframe(x, smy)
  }
  return(smy)
}

# copied from inlabru author Finn Lindgren
expand_to_dataframe <- function (x, data = NULL)
{
  if (is.null(data)) {
    data <- data.frame(matrix(nrow = NROW(x), ncol = 0))
  }
  only_x <- setdiff(names(x), names(data))
  if (length(only_x) < length(names(x))) {
    x <- x[!(names(x) %in% names(data))]
  }
  if (inherits(x, "SpatialPixels") && !inherits(x, "SpatialPixelsDataFrame")) {
    result <- sp::SpatialPixelsDataFrame(x, data = data)
  }
  else if (inherits(x, "SpatialGrid") && !inherits(x, "SpatialGridDataFrame")) {
    result <- sp::SpatialGridDataFrame(x, data = data)
  }
  else if (inherits(x, "SpatialLines") && !inherits(x, "SpatialLinesDataFrame")) {
    result <- sp::SpatialLinesDataFrame(x, data = data)
  }
  else if (inherits(x, "SpatialPolygons") && !inherits(x, "SpatialPolygonsDataFrame")) {
    result <- sp::SpatialPolygonsDataFrame(x, data = data)
  }
  else if (inherits(x, "SpatialPoints") && !inherits(x, "SpatialPointsDataFrame")) {
    result <- sp::SpatialPointsDataFrame(x, data = data)
  }
  else if (inherits(x, "Spatial")) {
    result <- sp::cbind.Spatial(x, data)
  }
  else {
    result <- cbind(x, data)
  }
  result
}
