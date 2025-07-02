\name{ora_full}
\alias{ora_full}
\title{Over-representation analysis with full enumeration of undetected metabolite patterns}
\description{
This function performs over-representation analysis (ORA) by enumerating all possible patterns of significant and non-significant assignments among undetected metabolites for each metabolite set. It returns the minimum, median, and maximum p-values from Fisher’s exact tests across these patterns, thereby estimating the full uncertainty range due to undetected metabolites.
}
\usage{
ora_full(SIG, DET, M)
}
\arguments{
  \item{SIG}{A character vector of statistically significant metabolite IDs.}
  \item{DET}{A character vector of all detected metabolite IDs (the background).}
  \item{M}{A named list of metabolite sets, where each element is a character vector of metabolite IDs.}
}
\value{
A list containing:
  \item{\code{Range of p-values}}{
    A matrix with rows corresponding to metabolite sets and three columns:
    \code{lower p-value}, \code{p-value(median)}, and \code{upper p-value}.
  }
}
\details{
For each metabolite set, the number of undetected metabolites is calculated.
The function then considers all possible numbers of significant metabolites (from 0 to the total number of undetected ones) among those undetected.
For each case, a 2x2 contingency table is constructed and Fisher’s exact test is applied.
The resulting p-values are aggregated to report the minimum, median, and maximum values.
}
\author{
Hiroyuki Yamamoto
}
\references{
Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.\\
Global functional profiling of gene expression. \emph{Genomics}. 2003;81(2):98–104.
}
\examples{
## Example: small set
M <- list(
  set1 = c("A", "B", "C"),
  set2 = c("C", "D", "E", "F"),
  set3 = c("G", "H")
)
DET <- c("A", "B", "C", "D", "E")  # Detected
SIG <- c("A", "C")                # Significant

ora_full(SIG, DET, M)$`Range of p-values`
}
\keyword{MSEA}
\keyword{ORA}
\keyword{enrichment}
