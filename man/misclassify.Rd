\name{misclassify}
\alias{misclassify}
\title{Adds some misclassification error to binary data}
\usage{
misclassify(binary.vector, error.1.0 = 0.05, error.0.1 = 0.05)
}
\arguments{
  \item{binary.vector}{a vector of binary values.}

  \item{error.1.0}{1 to 0 misclassification rate.}

  \item{error.0.1}{0 to 1 misclassification rate.}
}
\value{
a binary vector
}
\description{
Introduces some misclassification in binary data. The
number of values misclassified is determined by the
misclassification rates.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

