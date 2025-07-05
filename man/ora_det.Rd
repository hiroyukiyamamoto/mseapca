\name{ora_det}
\alias{ora_det}
\title{ORA using detected metabolites}
\description{
This function performs metabolite set enrichment analysis using over-representation analysis (ORA) under the assumption that only detected metabolites are used as the background.
A one-sided Fisher's exact test is applied to each metabolite set.
}
\usage{
ora_det(SIG, DET, M)
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
}
\value{
A list containing:
\itemize{
  \item \code{Result of MSEA(ORA)}: A matrix with raw p-values and adjusted q-values (BH correction) for each metabolite set.
  \item \code{significant metabolites}: A list of significant metabolites overlapping with each metabolite set.
  \item \code{Contingency tables}: A list of 2×2 contingency tables used in Fisher's exact tests.
}
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
B <- ora_det(SIG, DET, M)
B$`Result of MSEA(ORA)`
}
\keyword{ htest }
