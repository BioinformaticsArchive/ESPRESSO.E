\name{sim.QTL.data.E}
\alias{sim.QTL.data.E}
\title{Generates subjects for a continuous outcome}
\usage{
sim.QTL.data.E(num.obs = 500, env.model = 0, ph.mean = 0, ph.sd = 1,
  env.efkt = 0.5, env.prev = 0.1, env.mean = 0, env.sd = 1,
  env.low.lim = 0, env.up.lim = 1, pheno.reliability = 0.9)
}
\arguments{
  \item{num.obs}{number of subjects to simulate.}

  \item{ph.mean}{mean of the outcome variable in the study
  population}

  \item{ph.sd}{standard deviation of the outcome in the
  study population}

  \item{env.model}{model of the environmental exposure.}

  \item{env.efkt}{effects of the environment determinants.}

  \item{env.prev}{prevalence of the environmental
  exposure.}

  \item{env.mean}{statistical mean for a normally
  distributed exposure.}

  \item{env.sd}{standard deviation for a normally
  distributed exposure.}

  \item{env.low.lim}{lower limit for a uniformly
  distributed exposure.}

  \item{env.up.lim}{upper limit for a uniformly distributed
  exposure.}

  \item{pheno.reliability}{reliability of the assessment
  for a quantitative outcome.}
}
\value{
a matrix that holds the outcome (\code{phenotype}) and
exposure (\code{environment}) data.
}
\description{
Generates the specified number of subjects for a
quantitative outcome.
}
\details{

}
\author{
Gaye, A.
}
\keyword{internal}

