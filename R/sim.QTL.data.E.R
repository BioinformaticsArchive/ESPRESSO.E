#' 
#' @title Generates subjects for a continuous outcome
#' @description Generates the specified number of subjects for a quantitative outcome.
#' @param num.obs number of subjects to simulate.
#' @param ph.mean mean of the outcome variable in the study population
#' @param ph.sd standard deviation of the outcome in the study population
#' @param env.model model of the environmental exposure.
#' @param env.efkt effects of the environment determinants.
#' @param env.prev prevalence of the environmental exposure.
#' @param env.mean statistical mean for a normally distributed exposure.
#' @param env.sd standard deviation for a normally distributed exposure.
#' @param env.low.lim lower limit for a uniformly distributed exposure.
#' @param env.up.lim upper limit for a uniformly distributed exposure.
#' @param pheno.reliability reliability of the assessment for a quantitative outcome.
#' @return a matrix that holds the outcome (\code{phenotype}) and exposure (\code{environment}) data.
#' @keywords internal
#' @author Gaye, A.
#'
sim.QTL.data.E <- function (num.obs = 500, env.model = 0, ph.mean = 0, ph.sd = 1, 
    env.efkt = 0.5, env.prev = 0.1, env.mean = 0, env.sd = 1, 
    env.low.lim = 0, env.up.lim = 1, pheno.reliability = 0.9) 
{
    numobs <- num.obs
    e.mod <- env.model
    e.efkt <- env.efkt
    e.prev <- env.prev
    e.mean <- env.mean
    e.sd <- env.sd
    e.lowlim <- env.low.lim
    e.uplim <- env.up.lim
    pheno.rel <- pheno.reliability
    env.data <- sim.env.data(num.obs = numobs, env.model = e.mod, 
        env.prev = e.prev, env.mean = e.mean, env.sd = e.sd, 
        env.low.lim = e.lowlim, env.up.lim = e.uplim)
    pheno.data <- sim.pheno.qtl.E(num.subjects = num.obs, pheno.mean = ph.mean, 
        pheno.sd = ph.sd, environment = env.data, env.efkt = e.efkt)
    obs.phenotype <- get.obs.pheno(phenotype = pheno.data, pheno.model = 1, 
        pheno.sd = ph.sd, pheno.reliability = pheno.rel)
    phenotype <- obs.phenotype
    sim.matrix <- cbind(phenotype, env.data)
    totalnumrows <- dim(sim.matrix)[1]
    sim.matrix <- cbind(1:totalnumrows, sim.matrix)
    colnames(sim.matrix) <- c("id", "phenotype", "environment")
    mm <- data.frame(sim.matrix)
}
