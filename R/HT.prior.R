#' The half t7 prior useful for the standard deviation of random effects
#'
#' @export
HT.prior = function(){
return(c("expression:
  sigma = exp(-theta/2);
  nu = 7;
  log_dens = 0 - 0.5 * log(nu * pi) + (lgamma((nu+1)/2) - lgamma((nu)/2));
  log_dens = log_dens - 0.5 * (nu + 1) * log(1 + ((sigma * sigma)/nu));
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"))
}
