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
  \item{SIG}{A character vector of statistically significant metabolites.}
  \item{DET}{A character vector of all detected metabolites (the background).}
  \item{M}{A named list of metabolite sets, where each element is a character vector of metabolites.}
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
Yamamoto H, Fujimori T, Sato H, Ishikawa G, Kami K, Ohashi Y,
Statistical hypothesis testing of factor loading in principal component analysis and its application to metabolite set enrichment analysis.
BMC Bioinformatics, (2014) 15(1):51.

Yamamoto H.
Probabilistic Over-Representation Analysis for Metabolite Set Enrichment Analysis Considering Undetected Metabolites", Jxiv, (2024).
}
\examples{
# Example: Metabolome data
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
DET <- fasting_mseapca$DET
M <- fasting_mseapca$pathway

# Perform ORA using detected metabolites only
B <- ora_full(SIG, DET, M)
B$`Range of p-values`
}
\keyword{ htest }
