\name{ora_bino}
\alias{ora_bino}
\title{Over-representation analysis with binomial resampling adjustment}
\description{
Performs ORA while adjusting for undetected metabolites by binomial resampling.
}
\usage{
ora_bino(SIG, DET, M, method = "naive",
         probs = c(0.025, 0.975), nsim = 1000, lambda = 5)
}
\arguments{
  \item{SIG}{Character vector of significant metabolites.}
  \item{DET}{Character vector of detected metabolites (background).}
  \item{M}{Named list of metabolite sets.}
  \item{method}{\code{"naive"}, \code{"weighted"}, or \code{"shrink"}.}
  \item{probs}{Quantiles for the confidence interval (default 95\%).}
  \item{nsim}{Number of binomial simulations (default 1000).}
  \item{lambda}{Shrinkage constant used when \code{method = "shrink"}.}
}
\value{
A list containing one matrix.
Rows = metabolite sets; columns = lower, median, and upper p-values.
}
\author{Hiroyuki Yamamoto}
\examples{
M   <- list(set1 = c("A","B","C"), set2 = c("C","D","E"))
DET <- c("A","B","C","D")
SIG <- c("A","C")

# Binomial-adjusted ORA (shrink method)
ora_bino(SIG, DET, M, method = "shrink", nsim = 200)
}
\keyword{MSEA}
\keyword{ORA}
\keyword{simulation}
