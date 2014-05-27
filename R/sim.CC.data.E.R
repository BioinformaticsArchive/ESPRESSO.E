#' 
#' @title Generates case and controls
#' @description Generates affected and non-affected subjects until the set sample size is achieved.
#' @param num.obs Number of observations to generate per iteration until the specified number 
#' of cases and controls is achieved.
#' @param numcases Number of cases to generate.
#' @param numcontrols Number of controls to generate.
#' @param allowed.sample.size Maximum number of observations allowed i.e. the total size of the population 
#' to sample from.
#' @param disease.prev Prevalence of the binary outcome.
#' @param env.model Model of the environmental exposure.
#' @param env.prev Prevelance of the environmental determinants.
#' @param env.mean Statistical mean under quantitative-normal model.
#' @param env.sd Standard deviation under quantitative-normal model.
#' @param env.low.lim Lower limit under quantitative-uniform model.
#' @param env.up.lim Upper limit under quantitative-uniform model.
#' @param env.OR Odds ratios of the environmental determinants.
#' @param baseline.OR Baseline odds ratio for subject on 95 percent population centile versus 5 percentile. 
#' This parameter reflects the heterogeneity in disease risk arising from determinates that have not been 
#' measured or have not been included in the model.
#' @param pheno.error misclassification rates: 1-sensitivity and 1-specificity
#' @return a list which holds a matrix, \code{data}, that contains the phenotype and exposure statuses
#' and an integer, \code{allowed.sample.size.exceeded}, which tells if the maximum population size has been 
#' exceeded.
#' @keywords internal
#' @author Gaye A.
#'
sim.CC.data.E <- function(num.obs=20000, numcases=2000, numcontrols=8000, allowed.sample.size=20000000, disease.prev=0.1,
env.model=0, env.prev=0.1, env.mean=0, env.sd=1, env.low.lim=0, env.up.lim=1, env.OR=1.5, baseline.OR=12.36, pheno.error=c(0.1,0.1)){
  
      # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
      complete <- 0
      complete.absolute <- 0
      cases.complete <- 0
      controls.complete <- 0
      block <- 0
      sample.size.excess <- 0
      
      # SET UP A MATRIX TO STORE THE GENERATED DATA
      sim.matrix <- matrix(numeric(0), ncol=2)
      
      # SET LOOP COUNTER
      numloops <- 0
      
      # LOOP UNTIL THE SET NUMBER OF CASES AND OR CONTROLS IS ACHIEVED OR THE 
      # THE SET POPULATION SIZE TO SAMPLE FROM IS REACHED
      while(complete==0 && sample.size.excess==0){
        
        # GENERATE THE TRUE EXPOSURE DATA
        if(env.model==0){
          env.U <- rbinom(num.obs, 1, env.prev)
          e.mean <- env.prev
          env.U <- env.U-e.mean
        }else{
          if(env.model==1){
            env.U <- rnorm(num.obs, env.mean, env.sd)
          }else{
            env.U <- rnorm(num.obs, env.low.lim, env.up.lim)
          }
        }
        env.data <- env.U
        
        # GENERATE SUBJECT EFFECT DATA THAT REFLECTS BASELINESET RISK: 
        # NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR WITH APPROPRIATE 
        # VARIANCE ON SCALE OF LOG-ODDS
        
        # CONVERT BASELINE ODDS RATIO FROM 5th TO 95th PERCENTILES INTO THE
        # CORRESPONDING VARIANCE FOR A NORMALLY DISTRIBUTED RANDOM EFFECT 
        baseline.variance <- (log(baseline.OR)/(2*qnorm(0.95)))^2
        
        # CREATE NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR
        # WITH APPROPRIATE VARIANCE ON SCALE OF LOG-ODDS
        subject.effect <- rnorm(num.obs,0,sqrt(baseline.variance))
        
        # GENERATE THE TRUE OUTCOME DATA
        # GET THE ALPHA AND BETA VALUES
        alpha <- log(disease.prev/(1-disease.prev))
        beta <- log(env.OR)
        
        # GENERATE THE LINEAR PREDICTOR
        lp <- alpha + beta*env.data + subject.effect
        
        # GET 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
        mu <- exp(lp)/(1 + exp(lp))
        
        # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRASETME
        pheno.data <- rbinom(num.obs,1,mu)
        
        # GENERATE THE OBSERVED OUTCOME DATA FROM WHICH WE SELECT CASES AND CONTROLS
        obs.phenotype <- misclassify(binary.vector=pheno.data, error.1.0=pheno.error[1], error.0.1=pheno.error[2])
        
        # STORE THE TRUE OUTCOME, GENETIC AND ALLELE DATA IN AN OUTPUT MATRIX 
        # WHERE EACH ROW HOLDS THE RECORDS OF ONE INDIVUDAL
        sim.matrix.temp <- cbind(obs.phenotype, env.data)
        
        # UPDATE THE MATRIX THAT HOLDS ALL THE DATA GENERATED SO FAR, AFTER EACH LOOP
        sim.matrix <- rbind(sim.matrix, sim.matrix.temp)
        
        # SELECT OUT CASES
        indxcases <- which(sim.matrix[,1]==1)
        sim.matrix.cases <- sim.matrix[indxcases,]
        
        # SELECT OUT CONTROLS
        indxcontrols <- which(sim.matrix[,1]==0)
        sim.matrix.controls <- sim.matrix[indxcontrols,]
        
        # COUNT THE NUMBER OF CASES AND CONTROLS THAT HAS BEEN GENERATED
        cases.simulated <- dim(sim.matrix.cases)[1]
        controls.simulated <- dim(sim.matrix.controls)[1]
        
        # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
        # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
        if(cases.simulated >= numcases)
        {
          sim.matrix.cases <- sim.matrix.cases[1:numcases,]
          cases.complete <- 1
        }
        
        # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
        # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
        if(controls.simulated>=numcontrols)
        {
          sim.matrix.controls <- sim.matrix.controls[1:numcontrols,]
          controls.complete <- 1
        }
        
        # HAVE WE NOW GENERATED THE SET NUMBER OF CASES AND CONTROLS?
        complete <- cases.complete*controls.complete  	
        
        # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
        block <- block + num.obs
        if(block >= allowed.sample.size) {sample.size.excess <- 1}
        
        # INCREMENT LOOP COUNTER
        numloops <- numloops + 1
      }
      
      # STACK FINAL DATA MATRIX WITH CASES FIRST
      sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
      totalnumrows <- dim(sim.matrix)[1]
      sim.matrix <- cbind(1:totalnumrows, sim.matrix)
      
      # NAME THE COLUMNS OF THE MATRIX AND RETURN IT AS A DATAFRAMEDATAFRAME
      colnames(sim.matrix) <- c("id", "phenotype", "environment")
      mm <- list(data=data.frame(sim.matrix), allowed.sample.size.exceeded=sample.size.excess) 
   return(mm)
}

