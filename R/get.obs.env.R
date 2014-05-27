#'
#' @title Generates the observed environmental exposure data
#' @description Adds a set level of error to binary or quantitative data (the true data)
#' to obtain data with a larger variance (the observed data). The level of error is determined by
#' the misclassification rates in binary exposure and by the set level of variance in the quantitative
#' exposure.
#' @param env.data a vector of environmental measures that represents the true data.
#' @param env.model distribution of the exposure: binary=0 , normal=1 or uniform=2.
#' @param env.prev prevalence of the environmental exposure.
#' @param env.error misclassification rates: 1-sensitivity and 1-specificity.
#' @param env.reliability reliability of the assessment of quantitative exposure.
#' @return a dataframe with two coloumns/variables:
#' \code{true.environment} the error free exposure data (true data).
#' \code{observed.environment} the true esposure data with some added error (observed data).
#' @keywords internal
#' @author Amadou Gaye
#' 
get.obs.env <- function(env.data=NULL, env.model=0, env.sd=1, env.prev=0.2, env.error=c(0.1,0.1), env.reliability=0.9){
  
  if(env.model==0){
    obs.env <-  misclassify(binary.vector=env.data, error.1.0=env.error[1], error.0.1=env.error[2])
  }else{
    var.error <- (env.sd^2/env.reliability)-env.sd^2
    numsubs <- length(env.data)
    if(env.model==1){
      obs.env <- rnorm(numsubs, env.data, sqrt(var.error))
    }else{
      env.model.error <- rnorm(numsubs, 0, sqrt(var.error))
      obs.env <- env.data + env.model.error  
    }
  }
  
  return(obs.env)
}