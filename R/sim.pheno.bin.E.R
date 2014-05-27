#' 
#' @title Generates phenotype statuses
#' @description Generates affected and non-affected subjects.
#' @param num.obs number of observations to generate per iteration.
#' @param disease.prev prevalence of the binary outcome.
#' @param environment environmental exposure data.
#' @param subject.effect.data subject effect data, reflects the heterogenity in baseline disease risk.
#' @param env.OR odds ratios related to the exposure.
#' @return a vector binary vector, the phenotype data.
#' @keywords internal
#' @author Gaye A.
#'
sim.pheno.bin.E <- function(num.obs=10000, disease.prev=0.1, environment=NULL, subject.effect.data=NULL, env.OR=1.5){
  
  # IF GENOTYPE AND SUBJECT EFFECT DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
  if(is.null(environment)){
    cat("\n\n ALERT!\n")
    cat(" No environment data found.\n")
    cat(" Check the argument 'environment'\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  if(is.null(subject.effect.data)){
    cat("\n\n ALERT!\n")
    cat(" No baseline effect data found.\n")
    cat(" Check the argument 'subject.effect.data'\n")
    stop(" End of process!\n\n", call.=FALSE)
  }

  numobs <- num.obs
  pheno.prev <- disease.prev
  env.data <- environment
  s.efkt.data <- subject.effect.data
  e.OR <- env.OR
 
  # GET THE ALPHA AND BETA VALUES
  alpha <- log(pheno.prev/(1-pheno.prev))
  beta <- log(e.OR)
  
  # GENERATE THE LINEAR PREDICTOR
  lp <- alpha + beta*env.data + s.efkt.data

  # GET 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
  mu <- exp(lp)/(1 + exp(lp))

  # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRAME
  phenotype <- rbinom(numobs,1,mu)
  return(phenotype)
}

