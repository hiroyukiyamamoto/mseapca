\name{ora_bino}
\alias{ora_bino}
\title{Over-representation analysis with binomial resampling adjustment}
\description{
Performs ORA while adjusting for undetected metabolites by binomial resampling.
}
\usage{
ora_bino(SIG, DET, M, method = "naive", probs = c(0.025, 0.975), nsim = 1000, lambda = 5)
}
\arguments{
  \item{SIG}{Character vector of significant metabolites.}
  \item{DET}{Character vector of detected metabolites.}
  \item{M}{Named list of metabolite sets.}
  \item{method}{\code{"naive"}, \code{"weighted"}, or \code{"shrink"}.}
  \item{probs}{Quantiles for the confidence interval (default 95\%).}
  \item{nsim}{Number of binomial simulations (default 1000).}
  \item{lambda}{Shrinkage parameter used when \code{method = "shrink"} (default: 5).}
}
\value{
A list containing one matrix.
Rows = metabolite sets; columns = lower, median, and upper p-values.
}
\author{Hiroyuki Yamamoto}
\examples{
# Example: Metabolome data
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
DET <- fasting_mseapca$DET
M <- fasting_mseapca$pathway

# Perform ORA using detected metabolites only
B <- ora_bino(SIG, DET, M, method="naive", nsim = 10)

B$`Range of p-values for each pathway`
}
\keyword{ htest }
