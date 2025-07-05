\name{ora_all}
\alias{ora_all}
\title{ORA using all metabolites}
\description{
This function performs over-representation analysis (ORA) using all metabolites present in the given metabolite set list as the background, without specifying a reference metabolome. This corresponds to the behavior of MetaboAnalyst when no reference metabolome is uploaded.
}
\usage{
ora_all(SIG, M)
}
\arguments{
  \item{SIG}{Character vector of significant metabolites}
  \item{M}{Named list of metabolite sets}
}
\value{
A list with:
  \item{\code{Result of MSEA(ORA)}}{ Matrix of p-values and q-values}
  \item{\code{significant metabolites}}{ List of significant metabolites per set}
  \item{\code{Contingency tables}}{ A list of 2×2 contingency tables used in Fisher's exact tests.}
}
\references{
Yamamoto H, Fujimori T, Sato H, Ishikawa G, Kami K, Ohashi Y,
Statistical hypothesis testing of factor loading in principal component analysis and its application to metabolite set enrichment analysis.
BMC Bioinformatics, (2014) 15(1):51.

Yamamoto H.
Probabilistic Over-Representation Analysis for Metabolite Set Enrichment Analysis Considering Undetected Metabolites", Jxiv, (2024).
}
\author{Hiroyuki Yamamoto}
\examples{
# Example: Metabolome data
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
M <- fasting_mseapca$pathway

# Perform ORA using detected metabolites only
B <- ora_all(SIG, M)
B$`Result of MSEA(ORA)`
}
\keyword{ htest }
