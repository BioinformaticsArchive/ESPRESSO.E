\name{run.espresso.E}
\alias{run.espresso.E}
\title{Runs a full ESPRESSO analysis}
\usage{
run.espresso.E(simulation.params = NULL, pheno.params = NULL,
  env.params = NULL, scenarios2run = 1)
}
\arguments{
  \item{simulation.params}{general parameters for the
  scenario(s) to analyse}

  \item{pheno.params}{paramaters for the outcome variables}

  \item{env.params}{parameters for the environmental
  determinant}

  \item{scenarios2run}{the indices of the scenarios one
  wishes to analyse if there are more than one scenario on
  the input tables.}
}
\value{
a summary table that contains both the input parameters and
the results of the analysis
}
\description{
This function calls the functions required to run a full
ESPRESSO analysis where the model consists of an outcome
(binary or continuous) determined by a binary or
quantitative environmental determinant.
}
\details{


The function calls all the functions required to generate
the error free data, add the error to obtain the observed
data, run a glm analysis and calculate the sample size
required and the empirical and theoretical power. The
functions called to carry the various tasks are internal.
}
\examples{
{

# load the table that hold the input parameters; each of the table
hold parameters for 4 scenarios:
# scenario 1: a binary outcome determined by a binary exposure
# scenario 2: a binary outcome determined by a continuous exposure
# scenario 3: a quantitative outcome determined by a binary exposure
# scenario 4: a quantitative outcome determined by a continuous exposure
data(simulation.params)
data(pheno.params)
data(env.params)

# run the function for the first two scenarios, two binomial models
run.espresso.E(simulation.params, pheno.params, env.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.E(simulation.params, pheno.params, env.params, scenarios2run=c(3,4))
}
}
\author{
Gaye A.
}

