\name{setlabel}
\alias{setlabel}
\title{Binary labeling of metabolite IDs across metabolite sets}
\description{
Assigns binary labels (1 or 0) to indicate whether each metabolite ID is included in each metabolite set. This function is used for over-representation analysis (ORA) in metabolomics or other omics data.
}
\usage{
setlabel(M_ID, M, option = NULL)
}
\arguments{
  \item{M_ID}{A character vector of metabolite IDs (e.g., detected or significant metabolites).}
  \item{M}{A named list of metabolite sets, where each element is a character vector of metabolite IDs.}
  \item{option}{Optional argument. If \code{"anno"}, input IDs are assumed to be annotation strings (e.g., "C00031,C00022") and are split by semicolon after replacing commas. If \code{NULL} (default), IDs are used as-is.}
}
\value{
A binary matrix with rows representing input metabolite IDs and columns representing metabolite sets. Each cell contains 1 if the metabolite ID is found in the set, and 0 otherwise.
}
\details{
When \code{option = "anno"}, each element of \code{M_ID} may contain multiple IDs separated by commas (e.g., from an annotation field). These are split and matched against the metabolite sets.

This function is used internally in various ORA methods (e.g., \code{ora_det}, \code{ora_all}, \code{ora_est}) to compute contingency tables for enrichment analysis.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example usage
M <- list(
  set1 = c("A", "B", "C"),
  set2 = c("C", "D"),
  set3 = c("E", "F")
)

M_ID <- c("A", "C", "E", "G")

setlabel(M_ID, M)
}
\keyword{MSEA}
\keyword{labeling}
\keyword{ORA}
