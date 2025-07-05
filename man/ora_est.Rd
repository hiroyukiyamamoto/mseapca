\name{ora_est}
\alias{ora_est}
\title{ORA adjusting for undetected metabolites}
\description{
This function performs metabolite set enrichment analysis using over-representation analysis (ORA), incorporating point estimates to adjust for potentially undetected metabolites.
It supports three estimation methods: naive, weighted, and shrinkage-based adjustment.
}
\usage{
ora_est(SIG, DET, M, method = "naive", lambda = 5)
}
\arguments{
  \item{SIG}{
Character vector of significant metabolites
}
  \item{DET}{
Character vector of detected metabolites (background set)
}
  \item{M}{
A named list, where each element is a metabolite set (e.g., pathway) containing character vectors of metabolites
}
  \item{method}{
A character string specifying the estimation method to use. One of \code{"naive"}, \code{"weighted"}, or \code{"shrink"}. Default is \code{"naive"}.
}
  \item{lambda}{
A numeric value used in the \code{"shrink"} method as a shrinkage parameter. Default is \code{5}.
}
}
\value{
A list containing:
\itemize{
  \item \code{Result of MSEA (ORA with adjustment)}: A matrix with raw p-values and adjusted q-values (BH correction) for each metabolite set.
  \item \code{significant metabolites}: A list of significant metabolites overlapping with each metabolite set.
  \item \code{Contingency tables}: A list of 2×2 contingency tables used in Fisher's exact tests.
}
}
\details{
This function estimates the impact of undetected metabolites on enrichment results. It builds upon the ORA results from detected metabolites, then adjusts the contingency tables by estimating how many undetected metabolites might be significant, based on a specified method.
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
B <- ora_est(SIG, DET, M)
B$`Result of MSEA (ORA with adjustment)`
}
\keyword{ htest }
