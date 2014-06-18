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
sim.CC.data.E <- function (num.obs = 20000, numcases = 2000, numcontrols = 8000, 
    allowed.sample.size = 2e+07, disease.prev = 0.1, env.model = 0, 
    env.prev = 0.1, env.mean = 0, env.sd = 1, env.low.lim = 0, 
    env.up.lim = 1, env.OR = 1.5, baseline.OR = 12.36, pheno.error = c(0.1, 
        0.1)) 
{
    complete <- 0
    complete.absolute <- 0
    cases.complete <- 0
    controls.complete <- 0
    block <- 0
    sample.size.excess <- 0
    sim.matrix <- matrix(numeric(0), ncol = 2)
    numloops <- 0
    while (complete == 0 && sample.size.excess == 0) {
        if (env.model == 0) {
            env.U <- rbinom(num.obs, 1, env.prev)
            e.mean <- env.prev
            env.U <- env.U - e.mean
        }
        else {
            if (env.model == 1) {
                env.U <- rnorm(num.obs, env.mean, env.sd)
            }
            else {
                env.U <- rnorm(num.obs, env.low.lim, env.up.lim)
            }
        }
        env.data <- env.U
        baseline.variance <- (log(baseline.OR)/(2 * qnorm(0.95)))^2
        subject.effect <- rnorm(num.obs, 0, sqrt(baseline.variance))
        alpha <- log(disease.prev/(1 - disease.prev))
        beta <- log(env.OR)
        lp <- alpha + beta * env.data + subject.effect
        mu <- exp(lp)/(1 + exp(lp))
        pheno.data <- rbinom(num.obs, 1, mu)
        obs.phenotype <- misclassify(binary.vector = pheno.data, 
            error.1.0 = pheno.error[1], error.0.1 = pheno.error[2])
        sim.matrix.temp <- cbind(obs.phenotype, env.data)
        sim.matrix <- rbind(sim.matrix, sim.matrix.temp)
        indxcases <- which(sim.matrix[, 1] == 1)
        sim.matrix.cases <- sim.matrix[indxcases, ]
        indxcontrols <- which(sim.matrix[, 1] == 0)
        sim.matrix.controls <- sim.matrix[indxcontrols, ]
        cases.simulated <- dim(sim.matrix.cases)[1]
        controls.simulated <- dim(sim.matrix.controls)[1]
        if (cases.simulated >= numcases) {
            sim.matrix.cases <- sim.matrix.cases[1:numcases, 
                ]
            cases.complete <- 1
        }
        if (controls.simulated >= numcontrols) {
            sim.matrix.controls <- sim.matrix.controls[1:numcontrols, 
                ]
            controls.complete <- 1
        }
        complete <- cases.complete * controls.complete
        block <- block + num.obs
        if (block >= allowed.sample.size) {
            sample.size.excess <- 1
        }
        numloops <- numloops + 1
    }
    sim.matrix <- rbind(sim.matrix.cases, sim.matrix.controls)
    totalnumrows <- dim(sim.matrix)[1]
    sim.matrix <- cbind(1:totalnumrows, sim.matrix)
    colnames(sim.matrix) <- c("id", "phenotype", "environment")
    mm <- list(data = data.frame(sim.matrix), allowed.sample.size.exceeded = sample.size.excess)
    return(mm)
}

