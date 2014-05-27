#'
#' @title Runs a full ESPRESSO analysis
#' @description This function calls the functions required to run a full ESPRESSO analysis
#' where the model consists of an outcome (binary or continuous) determined by a binary or
#' quantitative environmental determinant.
#' @details The function calls all the functions required to generate the error free data, add 
#' the error to obtain the observed data, run a glm analysis and calculate the sample size required and
#' the empirical and theoretical power. The functions called to carry the various tasks are internal.
#' @param simulation.params general parameters for the scenario(s) to analyse
#' @param pheno.params paramaters for the outcome variables
#' @param env.params parameters for the environmental determinant
#' @param scenarios2run the indices of the scenarios one wishes to analyse if there are more than
#' one scenario on the input tables.
#' @return a summary table that contains both the input parameters and the results of the analysis
#' @export
#' @author Gaye A.
#' @examples {
#'
#' # load the table that hold the input parameters; each of the table 
#' # hold parameters for 4 scenarios:
#' # scenario 1: a binary outcome determined by a binary exposure
#' # scenario 2: a binary outcome determined by a continuous exposure
#' # scenario 3: a quantitative outcome determined by a binary exposure
#' # scenario 4: a quantitative outcome determined by a continuous exposure
#' data(simulation.params)
#' data(pheno.params)
#' data(env.params)
#'
#' # run the function for the first two scenarios, two binomial models
#' run.espresso.E(simulation.params, pheno.params, env.params, scenarios2run=c(1,2))
#'
#' # run the function for the last two scenarios, two gaussian models
#' run.espresso.E(simulation.params, pheno.params, env.params, scenarios2run=c(3,4))
#' }
#'
run.espresso.E <- function(simulation.params=NULL, pheno.params=NULL, env.params=NULL, scenarios2run=1){
  
  # IF AN INPUT FILE IS NOT SUPPLIED LOAD THE DEFAULT TABLES WARNING
  if(is.null(simulation.params)){
    cat("\n WARNING!\n")
    cat(" No simulation parameters supplied\n")
    cat(" The default simulation parameters will be used\n")
    simulation.params <- data(simulation.params)
  }
  
  if(is.null(pheno.params)){
    cat("\n WARNING!\n")
    cat(" No outcome parameters supplied\n")
    cat(" The default outcome parameters will be used\n")
    pheno.params <- data(pheno.params)
  }
  
  if(is.null(env.params)){
    cat("\n WARNING!\n")
    cat(" No environmental parameters supplied\n")
    cat(" The default environmental parameters will be used\n")
    env.params <- data(env.params)
  }
  
  # MERGE INPUT FILES TO MAKE ONE TABLE OF PARAMETERS
  s.temp1 <- merge(simulation.params, pheno.params)
  s.parameters <- merge(s.temp1, env.params)
  
  
  #----------LOAD SET UP UP INITIAL PARAMETERS------------#
  
  # PRINT TRACER CODE EVERY Nth ITERATION
  # THIS ENSURES THAT YOU CAN SEE IF THE PROGRAM GRINDS TO A HALT FOR SOME REASON (IT SHOULDN'T)
  trace.interval <- 10
  
  
  # CREATE UP TO 20M SUBJECTS IN BLOCKS OF 20K UNTIL REQUIRED NUMBER OF
  # CASES AND CONTROLS IS ACHIEVED. IN GENERAL THE ONLY PROBLEM IN ACHIEVING THE
  # REQUIRED NUMBER OF CASES WILL OCCUR IF THE DISEASE PREVALENCE IS VERY LOW
  max.pop.size <- 20000000
  numobs <- 20000
  
  
  # DECLARE MATRIX THAT STORE THE RESULTS FOR EACH SCENARIO (ONE PER SCENARIO PER ROW)
  output.file <- "output.csv"
  output.matrix <- matrix(numeric(0), ncol=36)
  column.names <- c(colnames(s.parameters), "exceeded.sample.size?","numcases.required", "numcontrols.required", 
                    "numsubjects.required", "empirical.power", "modelled.power","estimated.OR")
  
  write(t(column.names),output.file,dim(output.matrix)[2],append=TRUE,sep=";")
  
  
  
  #-----------LOOP THROUGH THE SCENARIOS - DEALS WITH ONE SCENARIO AT A TIME-------------
  
  for(j in c(scenarios2run))
  {
    
    # SIMULATION PARAMETERS
    scenario <- s.parameters$scenario.id[j]  			     
    seed <- s.parameters$seed.val[j]					          
    nsims <- s.parameters$numsims[j]					            
    ncases <- s.parameters$numcases[j]					          
    ncontrols <- s.parameters$numcontrols[j]		
    nsubjects <- s.parameters$numsubjects[j]				     
    baseline.odds <- s.parameters$RR.5.95[j]					           
    pvalue <- s.parameters$p.val[j]						                
    tpower <- s.parameters$power[j]
    
    # OUTCOME PARAMETERS
    pheno.mod <- s.parameters$pheno.model[j]
    pheno.mean <- s.parameters$pheno.mean[j]
    pheno.sd <- s.parameters$pheno.sd[j]
    pheno.prev <- s.parameters$disease.prev[j]
    pheno.err <- c(1-s.parameters$pheno.sensitivity[j],1-s.parameters$pheno.specificity[j])
    pheno.rel <- s.parameters$pheno.reliability[j]    
    
    # ENVIRONMENTAL DETERMINANTS PARAMETERS
    e.mod<- s.parameters$env.model[j]
    e.prev <-  s.parameters$env.prevalence[j]    
    e.OR <- s.parameters$env.OR[j]
    e.efkt <- s.parameters$env.efkt[j]
    e.mean <- s.parameters$env.mean[j]
    e.sd <- s.parameters$env.sd[j]
    e.low.lim <- s.parameters$env.low.lim[j]
    e.up.lim <- s.parameters$env.up.lim[j]
    e.error <- c(1-s.parameters$env.sensitivity[j],1-s.parameters$env.specificity[j])
    e.rel <- s.parameters$env.reliability[j]
    
    # VECTORS TO HOLD BETA, SE AND Z VALUES AFTER EACH RUN OF THE SIMULATION
    beta.values <- rep(NA,nsims)
    se.values <- rep(NA,nsims)
    z.values<-rep(NA,nsims)
    
    
    # TRACER TO DETECT EXCEEDING MAX ALLOWABLE SAMPLE SIZE
    sample.size.excess <- 0
    
    # RANDOM NUMBER GENERATOR STARTS WITH SEED SET AS SPECIFIED 
    set.seed(seed)
    
    # GENERATE AND ANALYSE DATASETS ONE AT A TIME 
    for(s in 1:nsims)            # s from 1 to total number of simulations
    {
      
      #-----------------------------GENERATE "TRUE" EXPOSURE DATA AND "OBSERVED" OUTCOME DATA-----------------------------#
      
      if(pheno.mod == 0){ # UNDER BINARY OUTCOME MODEL
        # GENERATE CASES AND CONTROLS UNTILL THE REQUIRED NUMBER OF CASES, CONTROLS IS ACHIEVED 
        sim.data <- sim.CC.data.E(num.obs=numobs,numcases=ncases,numcontrols=ncontrols,allowed.sample.size=max.pop.size,disease.prev=pheno.prev,env.model=e.mod, 
                                  env.prev=e.prev,env.mean=e.mean,env.sd=e.sd,env.low.lim=e.low.lim,env.up.lim=e.up.lim,env.OR=e.OR,baseline.OR=baseline.odds,pheno.error=pheno.err)
        
        t.data <- sim.data$data
        
      }else{ # UNDER QUANTITATIVE OUTCOME MODEL
        # GENERATE THE SPECIFIED NUMBER OF SUBJECTS
        t.data <- sim.QTL.data.E(num.obs=nsubjects,env.model=e.mod,ph.mean=pheno.mean,ph.sd=pheno.sd,env.efkt=e.efkt,env.prev=e.prev,env.mean=e.mean,env.sd=e.sd,env.low.lim=e.low.lim,
                                 env.up.lim=e.up.lim)
      }
      
      #------------SIMULATE ERRORS AND ADD THEM TO THE TRUE COVARIATES DATA TO OBTAIN OBSERVED COVARIATES DATA-----------#
      
      # ADD APPROPRIATE ERRORS TO PRODUCE OBSERVED EXPOSURE MESEAUREMENTS 
      obs.env <- get.obs.env(env.data=t.data$environment, env.model=e.mod, env.sd=e.sd, env.prev=e.prev, env.error=e.error, env.reliability=e.rel)
      t.data$environment <- obs.env
      o.data <- t.data
      
      #--------------------------DATA ANALYSIS ----------------------------#
      glm.estimates <- glm.analysis.E(pheno.model=pheno.mod, observed.data=o.data)
      
      beta.values[s] <- glm.estimates[[1]]
      se.values[s] <- glm.estimates[[2]]
      z.values[s] <- glm.estimates[[3]]
      
      # PRINT TRACER AFTER EVERY Nth DATASET CREATED
      if(s %% trace.interval ==0)cat("\n",s,"of",nsims,"runs completed in scenario",scenario)
      
    }
    cat("\n\n")
    
    #------------------------ SUMMARISE RESULTS ACROSS ALL SIMULATIONS---------------------------#
    
    # SUMMARISE PRIMARY PARAMETER ESTIMATES
    # COEFFICIENTS ON LOG-ODDS SCALE
    m.beta <- mean(beta.values, na.rm=T)
    m.se <- sqrt(mean(se.values^2, na.rm=T))
    m.model.z <- m.beta/m.se
    
    
    #---------------------------POWER AND SAMPLE SIZE CALCULATIONS----------------------#
    
    # CALCULATE THE SAMPLE SIZE REQUIRED UNDER EACH MODEL
    sample.sizes.needed <- samplsize.calc(numcases=ncases, numcontrols=ncontrols, num.subjects=nsubjects, 
                                          pheno.model=pheno.mod, pval=pvalue, power=tpower, mean.model.z=m.model.z)
    
    # CALCULATE EMPIRICAL POWER AND THE MODELLED POWER 
    # THE EMPIRICAL POWER IS SIMPLY THE PROPORTION OF SIMULATIONS IN WHICH
    # THE Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC
    # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
    zvals <- z.values
    power <- power.calc(pval=pvalue, z.values=zvals, mean.model.z=m.model.z)
    
    
    #------------------MAKE FINAL A TABLE THAT HOLDS BOTH INPUT PARAMETERS AND OUTPUT RESULTS---------------#
    
    critical.res <- get.critical.results.E(scenario=j,pheno.model=pheno.mod,env.model=e.mod,sample.sizes.required=sample.sizes.needed,
                                           empirical.power=power$empirical,modelled.power=power$modelled,mean.beta=m.beta)
    
    #  WHEN OUTCOME IS BINARY INFORM IF RECORD EXCEEDED MAXIMUM SAMPLE SIZE
    if(pheno.mod==0){
      sample.size.excess <- sim.data$allowed.sample.size.exceeded
      if(sample.size.excess==1)
      {
        excess <- "yes"
        cat("\nTO GENERATE THE NUMBER OF CASES SPECIFIED AT OUTSET\n")
        cat("THE SIMULATION EXCEEDED THE MAXIMUM POPULATION SIZE OF ", max.pop.size,"\n")
      }else{
        excess <- "no"
      }
    }
    
    inparams <- s.parameters[j,]
    
    if(pheno.mod==0){
      mod <- "binary"
      if(e.mod==0){
        inparams [c(6,16,20,21,22,23,24,27)] <- "NA"
        inputs <- inparams
        outputs <- c(excess, critical.res[[2]], critical.res[[3]], "NA", critical.res[[4]], 
                     critical.res[[5]], critical.res[[6]])
      }else{
        if(e.mod==1){
          inparams [c(6,16,20,23,24,25,26)] <- "NA"
          inputs <- inparams
          outputs <- c(excess, critical.res[[2]], critical.res[[3]], "NA", critical.res[[4]], 
                       critical.res[[5]], critical.res[[6]])
        }else{
          inparams [c(6,16,20,23,24,25,26)] <- "NA"
          inputs <- inparams
          outputs <- c(excess, critical.res[[2]], critical.res[[3]], "NA", critical.res[[4]], 
                       critical.res[[5]], critical.res[[6]])
        }
      }
    }else{
      mod <- "quantitative"
      if(e.mod==0){
        inparams [c(4,5,11,14,15,19,21,22,23,24,27)] <- "NA"
        inputs <- inparams
        outputs <- c("NA", "NA", "NA", critical.res[[2]], critical.res[[3]], 
                     critical.res[[4]], critical.res[[5]])
      }else{
        if(e.mod==1){
          inparams <- s.parameters[j,]
          inparams [c(4,5,11,14,15,19,23,24,25,26)] <- "NA"
          inputs <- inparams
          outputs <- c("NA", "NA", "NA", critical.res[[2]], critical.res[[3]], 
                       critical.res[[4]], critical.res[[5]])
        }else{
          inparams [c(4,5,11,14,15,19,21,22,25,26)] <- "NA"
          inputs <- inparams
          outputs <- c("NA", "NA", "NA", critical.res[[2]], critical.res[[3]], 
                       critical.res[[4]], critical.res[[5]])
        }
      }
    }
    
    jth.row <- as.character(c(inputs,outputs))
    write(t(jth.row),output.file,dim(output.matrix)[2],append=TRUE,sep=";")
  }
}
