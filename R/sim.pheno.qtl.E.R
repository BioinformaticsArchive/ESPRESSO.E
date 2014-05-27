#' 
#' @title Generates quantitative outcome data
#' @description The function uses the effects data of the determinants to construct a linear predictor(LP). 
#' The outcome is a normally distributed variable generated with a mean equal to the LP and a standard 
#' deviation of 1.
#' @param num.obs number of subjects to generate.
#' @param pheno.mean mean of the outcome variable in the study population
#' @param pheno.sd standard deviation of the outcome in the study population
#' @param environment environmental exposure data.
#' @param env.efkt effect size of the environmental exposure.
#' @return a numeric vector, the phenotype data
#' @keywords internal
#' @author Gaye A.
#'  
sim.pheno.qtl.E <- function(num.subjects=10000, pheno.mean=0, pheno.sd=1, environment=NULL, env.efkt=0.25){
  
   # IF GENOTYPE DATA ARE NOT SUPPLIED STOP AND ISSUE AN ALERT
   if(is.null(environment)){
			  cat("\n\n ALERT!\n")
			  cat(" No environmental exposure data found.\n")
			  cat(" Check the argument 'environment'\n")
			  stop(" End of process!\n\n", call.=FALSE)
	 }

   # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   alpha <- pheno.mean
   beta <-	env.efkt

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (beta*environment)

   # GENERATE THE TRUE PHENOTYPE DATA TO RETURN
   phenotype <- rnorm(num.subjects,lp,pheno.sd)
   return(phenotype)
}

