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
get.obs.pheno <- function (phenotype = NULL, pheno.model = 0, pheno.sd = 1, pheno.error = c(0.05, 
    0.05), pheno.reliability = 0.9) 
{
    if (is.null(phenotype)) {
        cat("\n\n ALERT!\n")
        cat(" No phenotype data found.\n")
        cat(" Check the argument 'phenotype'\n")
        stop(" End of process!\n\n", call. = FALSE)
    }
    if (is.null(pheno.model)) {
        cat("\n\n ALERT!\n")
        cat(" No outcome  model provided\n")
        cat(" Check the argument 'pheno.model'\n")
        stop(" End of process!\n\n", call. = FALSE)
    }
    true.phenotype <- phenotype
    if (pheno.model == 0) {
        observed.phenotype <- misclassify(binary.vector = phenotype, 
            error.1.0 = pheno.error[1], error.0.1 = pheno.error[2])
    }
    else {
        var.m <- (pheno.sd^2/pheno.reliability) - (pheno.sd^2)
        num.obs <- length(phenotype)
        observed.phenotype <- rnorm(num.obs, phenotype, sqrt(var.m))
        var.m <- (pheno.sd^2/pheno.reliability) - (pheno.sd^2)
        num.obs <- length(true.phenotype)
        observed.phenotype <- rnorm(num.obs, true.phenotype, 
            sqrt(var.m))
    }
    df <- observed.phenotype
}
