\name{sim.CC.data.E}
\alias{sim.CC.data.E}
\title{Generates case and controls}
\usage{
sim.CC.data.E(num.obs = 20000, numcases = 2000, numcontrols = 8000,
  allowed.sample.size = 2e+07, disease.prev = 0.1, env.model = 0,
  env.prev = 0.1, env.mean = 0, env.sd = 1, env.low.lim = 0,
  env.up.lim = 1, env.OR = 1.5, baseline.OR = 12.36,
  pheno.error = c(0.1, 0.1))
}
\arguments{
  \item{num.obs}{Number of observations to generate per
  iteration until the specified number of cases and
  controls is achieved.}

  \item{numcases}{Number of cases to generate.}

  \item{numcontrols}{Number of controls to generate.}

  \item{allowed.sample.size}{Maximum number of observations
  allowed i.e. the total size of the population to sample
  from.}

  \item{disease.prev}{Prevalence of the binary outcome.}

  \item{env.model}{Model of the environmental exposure.}

  \item{env.prev}{Prevelance of the environmental
  determinants.}

  \item{env.mean}{Statistical mean under
  quantitative-normal model.}

  \item{env.sd}{Standard deviation under
  quantitative-normal model.}

  \item{env.low.lim}{Lower limit under quantitative-uniform
  model.}

  \item{env.up.lim}{Upper limit under quantitative-uniform
  model.}

  \item{env.OR}{Odds ratios of the environmental
  determinants.}

  \item{baseline.OR}{Baseline odds ratio for subject on 95
  percent population centile versus 5 percentile. This
  parameter reflects the heterogeneity in disease risk
  arising from determinates that have not been measured or
  have not been included in the model.}

  \item{pheno.error}{misclassification rates: 1-sensitivity
  and 1-specificity}
}
\value{
a list which holds a matrix, \code{data}, that contains the
phenotype and exposure statuses and an integer,
\code{allowed.sample.size.exceeded}, which tells if the
maximum population size has been exceeded.
}
\description{
Generates affected and non-affected subjects until the set
sample size is achieved.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

