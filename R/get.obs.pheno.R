#' 
#' @title Generates observed outcome data
#' @description Adds a set level of error to error free binary or quantitative data (the true phenotype data) 
#' to obtain data with a larger variance (the observed phenotype data).
#' @param seed 
#' @param phenotype outcome status.
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1
#' @param pheno.sd standard deviation of the outcome in the study population
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1 or uniform=2.
#' @param pheno.sd standard deviation of the outcome in the study the population
#' @param pheno.error misclassification rates: 1-sensitivity and 1-specificity
#' @param pheno.reliability reliability of the assessment for a quantitative outcome.
#' @return A dataframe containing:
#' \code{true.phenotype} the error free outcome data (true data).
#' \code{observed.phenotype} the true outcome data with some added error (observed data).
#' @keywords internal
#' @author Gaye A.
#'
get.obs.pheno <- function (phenotype=NULL, pheno.model=0, pheno.sd=1, pheno.error=c(0.05,0.05), pheno.reliability=0.9){ 
  
  if(is.null(phenotype)){
    cat("\n\n ALERT!\n")
    cat(" No phenotype data found.\n")
    cat(" Check the argument 'phenotype'\n")
    stop(" End of process!\n\n", call.=FALSE)
  } 
  if(is.null(pheno.model)){
    cat("\n\n ALERT!\n")
    cat(" No outcome  model provided\n")
    cat(" Check the argument 'pheno.model'\n")
    stop(" End of process!\n\n", call.=FALSE)
  } 
  
 true.phenotype <- phenotype
  
  # GET THE OBSERVED OUTCOME DATA
  if(pheno.model==0){ # IF THE OUTCOME IS BINARY
    
    observed.phenotype <- misclassify(binary.vector=phenotype, error.1.0=pheno.error[1], error.0.1=pheno.error[2])
    
  }else{ # IF THE OUTCOME IS CONTINUOUS NORMAL
    # USE THE RELIABITLITY OF PHENOTYPE ASSESSMENT TO COMPUTE THE VARIANCE OF MEASURED PHENOTYPES.
    # RELIABITY = (VAR.true.data/VAR.obs.data) AND VAR.obs.data = (VAR.true.data + VAR.measurement)
    # IT FOLLOWS THAT VAR.true.data + VAR.of.estimate = VAR.true.data/RELIABILITY AND THEN:
    # VAR.measurement = (VAR.true.data/RELIABILITY) - VAR.true.data
    var.m <- (pheno.sd^2/pheno.reliability)-(pheno.sd^2)
    
    # GENERATE THE NORMALLY DISTRIBUTED ERROR USING THE ABOVE COMPUTED VARIANCE
    num.obs <- length(phenotype)
    observed.phenotype <- rnorm(num.obs,phenotype,sqrt(var.m))
    var.m <- (pheno.sd^2/pheno.reliability)-(pheno.sd^2) 
    
    # ADD THE ERROR TO ORIGINAL PHENOTYPES TO GENERATE THE OBSERVED PHENOTYPE DATA
    num.obs <- length(true.phenotype)
    observed.phenotype <- rnorm(num.obs,true.phenotype,sqrt(var.m))
  } 
  
  # RETURN THE TRUE AND OBSERVED PHENOTYPE DATA AS A DATAFRAME
  df <- observed.phenotype
  
}
