\name{sim.pheno.bin.E}
\alias{sim.pheno.bin.E}
\title{Generates phenotype statuses}
\usage{
sim.pheno.bin.E(num.obs = 10000, disease.prev = 0.1, environment = NULL,
  subject.effect.data = NULL, env.OR = 1.5)
}
\arguments{
  \item{num.obs}{number of observations to generate per
  iteration.}

  \item{disease.prev}{prevalence of the binary outcome.}

  \item{environment}{environmental exposure data.}

  \item{subject.effect.data}{subject effect data, reflects
  the heterogenity in baseline disease risk.}

  \item{env.OR}{odds ratios related to the exposure.}
}
\value{
a vector binary vector, the phenotype data.
}
\description{
Generates affected and non-affected subjects.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

