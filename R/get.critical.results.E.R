#' 
#' @title Provides a summary of the main results
#' @description Gets the number of cases and controls or subjects and the 
#' empirical and theoretical power under each model and prints a summary on 
#' the screen.
#' @param scenario scenario number.
#' @param pheno.model outcome type; binary=0 and quantitative=1.
#' @param env.model model of the exposure: binary=0, quantitative-normal=1, quantitative-uniform=2.
#' @param sample.sizes.required number of cases and controls or number of subjects 
#' required to achieve the desired power.
#' @param empirical.power estimated empirical power.
#' @param modelled.power calculated theoretical power.
#' @param mean.beta mean beta value of each of the determinants.
#' @return a list that contains the sample size required to achieve the desired power and the
#' estimated modelled and theoretical power.
#' @keywords internal
#' @author Gaye A.
#'
get.critical.results.E <- function(scenario=1, pheno.model=0,env.model=0,sample.sizes.required=NULL,
                                   empirical.power=0.8,modelled.power=0.8,mean.beta=NULL){

  if(is.null(sample.sizes.required)){
    cat("\n\n ALERT!\n")
    cat(" The sample size required may have not been computed.\n")
    cat(" Check the output of the function 'samplsize.calc'.\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  if(is.null(mean.beta)){
    cat("\n\n ALERT!\n")
    cat(" The argument 'mean.beta' is empty.\n")
    cat(" This argument should average beta value'.\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  if(env.model==0){
    e.model <- "binary"
  }else{
    if(env.model==1){
      e.model <- "normal distribution"
    }else{
      e.model <- "uniform distribution"
    }
  }
  if(pheno.model == 0){
    numcases <- sample.sizes.required[[1]]
    numcontrols <- sample.sizes.required[[2]]
  }else{
    numsubjects <- sample.sizes.required[[1]]
  }

  if(pheno.model==0){
    # estimated ORs
    estimated.OR <- exp(mean.beta)
    
    cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
    cat("\nModels\n")
    cat("------\n")
    cat("  Outcome: binary \n")
    cat("  Environmental determinant: ",e.model)
    
    cat("\n\nNumber of cases required\n")
    cat("------------------------\n")
    cat(" ", numcases)
    
    cat("\n\nNumber of controls required\n")
    cat("---------------------------\n")
    cat(" ", numcontrols)
    
    cat("\n\nEmpirical power\n")
    cat("---------------\n")
    cat(" ",empirical.power)
    
    cat("\n\nModel power\n")
    cat("-----------\n")
    cat(" ",round(modelled.power,2))
    
    cat("\n\nEstimated ORs\n")
    cat("-----------\n")
    cat(" ",round(estimated.OR,2))
    
    
    cat("\n\n---- END OF SUMMARY ----\n")
    
    crit.res <- c(e.model,numcases,numcontrols,round(empirical.power,2),round(modelled.power,2),
                  round(estimated.OR,2))
    return(list(environment.model=crit.res[1],number.of.cases.required=crit.res[2],
                number.of.controls.required=crit.res[3],empirical.power=crit.res[4], 
                modelled.power=crit.res[5],estimated.OR=crit.res[6]))
    
  }else{
    # estimated ORs
    estimated.OR <- 'NA'
    
    cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
    cat("\nModels\n")
    cat("------\n")
    cat(" Outcome: quantitative \n")
    cat(" Environmental determinant: ",e.model)
    
    cat("\n\nNumber of subjects required\n")
    cat("------------------------\n")
    cat(" ",numsubjects)
    
    cat("\n\nEmpirical power\n")
    cat("---------------\n")
    cat(" ",empirical.power)
    
    cat("\n\nModel power\n")
    cat("-----------\n")
    cat(" ",round(modelled.power,2))
    
    cat("\n\nEstimated ORs\n")
    cat("-----------\n")
    cat(" ",estimated.OR)
    
    cat("\n\n---- END OF SUMMARY ----\n")
    
    crit.res <- c(e.model,numsubjects,round(empirical.power,2),round(modelled.power,2),estimated.OR)
    return(list(environment.model=crit.res[1],number.of.subjects.required=crit.res[2],
                empirical.power=crit.res[3], modelled.power=crit.res[4], estimated.OR=crit.res[5]))
  }
}

