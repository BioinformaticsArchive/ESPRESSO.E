\name{sim.pheno.qtl.E}
\alias{sim.pheno.qtl.E}
\title{Generates quantitative outcome data}
\usage{
sim.pheno.qtl.E(num.subjects = 10000, pheno.mean = 0, pheno.sd = 1,
  environment = NULL, env.efkt = 0.25)
}
\arguments{
  \item{num.obs}{number of subjects to generate.}

  \item{pheno.mean}{mean of the outcome variable in the
  study population}

  \item{pheno.sd}{standard deviation of the outcome in the
  study population}

  \item{environment}{environmental exposure data.}

  \item{env.efkt}{effect size of the environmental
  exposure.}
}
\value{
a numeric vector, the phenotype data
}
\description{
The function uses the effects data of the determinants to
construct a linear predictor(LP). The outcome is a normally
distributed variable generated with a mean equal to the LP
and a standard deviation of 1.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

