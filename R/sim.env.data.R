#' 
#' @title Simulates cases and controls
#' @description Generates data for a binary, quantitative-normal, or quantitative-uniform environmental determinant.
#' @param num.obs number of observations to simulate.
#' @param env.model model of the exposure: binary=0, quantitative-normal=1, quantitative-uniform=2.
#' @param env.prev prevalence of the environmental exposure.
#' @param env.mean statisitical man under quantitative-normal model.
#' @param env.sd standard deviation under quantitative-normal model.
#' @param env.low.lim lower limit under quantitative-uniform model.
#' @param env.up.lim upper limit under quantitative-uniform model.
#' @return a vector of continuous or binary values.
#' @keywords internal
#' @author Gaye A.
#' 
sim.env.data <- function (num.obs = 10000, env.model = 0, env.prev = 0.1, env.mean = 0, 
    env.sd = 1, env.low.lim = 0, env.up.lim = 1) 
{
    numobs <- num.obs
    e.mod <- env.model
    e.prev <- env.prev
    e.sd <- env.sd
    e.low.lim <- env.low.lim
    e.up.lim <- env.up.lim
    if (e.mod == 0) {
        env.U <- rbinom(numobs, 1, e.prev)
        e.mean <- e.prev
        env.U <- env.U - e.mean
    }
    if (e.mod == 1) {
        e.mean <- env.mean
        env.U <- rnorm(numobs, e.mean, e.sd)
        env.U <- env.U - mean(env.U)
    }
    if (e.mod == 2) {
        if (env.low.lim >= env.up.lim) {
            stop("\n\nALERT!\n Uniform Distribution: The upper limit must be greater than the lower limit\n\n")
        }
        else {
            env.U <- runif(numobs, e.low.lim, e.up.lim)
        }
    }
    return(env.U)
}
