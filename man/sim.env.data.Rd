\name{sim.env.data}
\alias{sim.env.data}
\title{Simulates cases and controls}
\usage{
sim.env.data(num.obs = 10000, env.model = 0, env.prev = 0.1,
  env.mean = 0, env.sd = 1, env.low.lim = 0, env.up.lim = 1)
}
\arguments{
  \item{num.obs}{number of observations to simulate.}

  \item{env.model}{model of the exposure: binary=0,
  quantitative-normal=1, quantitative-uniform=2.}

  \item{env.prev}{prevalence of the environmental
  exposure.}

  \item{env.mean}{statisitical man under
  quantitative-normal model.}

  \item{env.sd}{standard deviation under
  quantitative-normal model.}

  \item{env.low.lim}{lower limit under quantitative-uniform
  model.}

  \item{env.up.lim}{upper limit under quantitative-uniform
  model.}
}
\value{
a vector of continuous or binary values.
}
\description{
Generates data for a binary, quantitative-normal, or
quantitative-uniform environmental determinant.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

