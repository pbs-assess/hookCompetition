#' SpatialPolygonsDataFrame object containing 4 survey boundaries of BC (HS, QCS, WCHG, WCVI)
#'
#' SpatialPolygonsDataFrame
#'
#' @format SpatialPolygonsDataFrame contains:
#' \describe{
#'   \item{Region}{Factor variable with levels HS, QCS, WCHG, and WCVI}
#'   \item{spatio_temporal_CPUE_trajectory}{data.frame containing 10 posterior samples of regional yelloweye relative abundance from the spatio-temporal CPUE-based model}
#'   \item{spatio_temporal_CPUE_inddf}{data.frame containing posterior summaries of regional yelloweye relative abundance from the spatio-temporal CPUE-based model}
#'   \item{Other_Objects}{The same as above but for ICR-based and censored models}
#' }
"survey_boundaries"
