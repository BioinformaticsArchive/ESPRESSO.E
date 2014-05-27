\name{simulation.params}
\alias{simulation.params}
\docType{data}
\title{Main parameters of the simulations}
\description{A table with 9 columns: scenario.id, seed.val, numsims, numcases, numcontrols, numsubjects, RR.5.95, p.val, power}
\usage{data(simulation.params)}
\format{
  A table with 9 columns and 1 row for each scenario
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{seed.val}}{seed value}
    \item{\code{numsims}}{number of runs for each simulation}
    \item{\code{numcases}}{number of cases under binary outcome}
    \item{\code{numcontrols}}{number of controls under binary outcome}
    \item{\code{numsubjects}}{number of subjects under continuous outcome}
    \item{\code{RR.5.95}}{the baseline odds ratio for subjects on 95 percent population centile versus 5 percentile. This parameter reflects the heterogeneity in disease risk arising from determinants that have not been measured or have not been included in the model. If this parameter is set to 10, the implication is that a high risk subject (someone at the upper 95 percent centile of population risk) is, all else being equal, at 10 times the odds of developing disease compared to someone else who is at low risk (individual at the lower 5 percent centile of population risk)}
    \item{\code{p.val}}{cut-off p-value defining statistical significance}
    \item{\code{power}}{desired power}
  }
}


\examples{data(simulation.params)}
