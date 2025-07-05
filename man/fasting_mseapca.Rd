\name{fasting_mseapca}
\alias{fasting_mseapca}
\docType{data}
\title{Metabolome data from a fasting mouse study}
\description{
Example data derived from a fasting mouse liver metabolomics study,
originally published by Yamamoto et al. and distributed in the \pkg{loadings} package as \code{fasting}.
This version is reformatted and annotated to reproduce the analysis setup used in MetaboAnalyst.
It is intended for demonstrating ORA and MSEA functions in the \pkg{mseapca} package.
}
\format{
A named \code{list} with the following components:

\describe{
  \item{\code{SIG}}{Character vector of significant metabolites (e.g., \eqn{p < 0.05}).}
  \item{\code{DET}}{Character vector of all detected metabolites.}
  \item{\code{pathway}}{Named list of metabolite sets (pathways). Each element is a character vector of metabolite IDs belonging to that pathway.}
}
}
\usage{
data(fasting_mseapca)
}
\details{
The dataset is intended for vignette examples:

\itemize{
  \item \code{SIG} — Significant metabolites (\eqn{p < 0.05}) obtained from
        statistical analysis of a fasting vs control comparison.
  \item \code{DET} — Background list of all detected metabolites.
  \item \code{pathway} — Pathway definitions from MetaboAnalyst, reused for mouse data.
}

The object can be used directly with functions such as
\code{\link{msea_ora}}, \code{\link{ora_det}}, and \code{\link{ora_bino}}.
}
\source{
Derived from the \code{fasting} dataset in the \pkg{loadings} package.

Original reference:
Yamamoto H., Fujimori T., Sato H., Ishikawa G., Kami K., Ohashi Y. (2014).
"Statistical hypothesis testing of factor loading in principal component analysis and its application to metabolite set enrichment analysis".
\emph{BMC Bioinformatics}, 15(1):51.
}
\examples{
## Load data
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
DET <- fasting_mseapca$DET
M   <- fasting_mseapca$pathway

## Simple ORA on detected metabolites
res <- ora_det(SIG, DET, M)
head(res$`Result of MSEA(ORA)`)
}
\keyword{ datasets }
