\name{msea_ora_range}
\alias{msea_ora_range}
\title{Wrapper function for Over-Representation Analysis with p-value range estimation}
\description{
This function performs over-representation analysis (ORA) to assess metabolite set enrichment while considering uncertainty due to undetected metabolites. It wraps different methods for estimating a p-value range, including full enumeration and binomial resampling.
}
\usage{
msea_ora_range(SIG, DET = NULL, M,
               option = "ora_full",
               probs = c(0.025, 0.975),
               nsim = 1000,
               lambda = 5)
}
\arguments{
  \item{SIG}{A character vector of statistically significant metabolites.}
  \item{DET}{A character vector of all detected metabolites. Required for all methods except \code{ora_full}.}
  \item{M}{A named list of metabolite sets, each containing a character vector of metabolites.}
  \item{option}{Method to use for estimating the p-value range. One of \code{"ora_full"}, \code{"bino_naive"}, \code{"bino_weighted"}, or \code{"bino_shrink"}.}
  \item{probs}{Numeric vector of quantile probabilities for binomial simulation (e.g., \code{c(0.025, 0.975)}). Ignored for \code{ora_full}.}
  \item{nsim}{Number of simulations for binomial-based estimation. Ignored for \code{ora_full}.}
  \item{lambda}{Shrinkage parameter for \code{"bino_shrink"} option.}
}
\value{
A list containing a matrix with p-value range results for each metabolite set. Columns include lower, median, and upper p-values.
}
\details{
This wrapper function allows switching between multiple ORA implementations that estimate the uncertainty due to undetected metabolites. The \code{ora_full} method uses exhaustive enumeration of all possible detection patterns, while the other methods use binomial resampling with different estimation strategies (naive, weighted, or shrinkage-based).
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example: Metabolome data
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
DET <- fasting_mseapca$DET
M <- fasting_mseapca$pathway

# Perform ORA using detected metabolites only
B <- msea_ora_range(SIG, DET, M, option = "bino_naive", nsim=10)

B$`Range of p-values for each pathway`
}
\keyword{ htest }
