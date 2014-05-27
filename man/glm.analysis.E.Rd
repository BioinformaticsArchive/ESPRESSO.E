\name{glm.analysis.E}
\alias{glm.analysis.E}
\title{Fits a Generalized Linear Model}
\usage{
glm.analysis.E(pheno.model = 0, observed.data = NULL)
}
\arguments{
  \item{pheno.model}{outcome type: binary=0 and normal=1}

  \item{observed.data}{a dataframe that contains the
  observerd outcome and covariate data}
}
\value{
a vector that holds the beta, standard-error and
z-statistic of each of the covariates
}
\description{
This function fits a conventional unconditional logistic
regression model on the observed outcome and exposure data.
In a Typical ESPRESSO the function is called at each
iteration and the estimates stored.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

