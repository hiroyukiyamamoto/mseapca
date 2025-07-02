\name{ora_det}
\alias{ora_det}
\title{Over-representation analysis using detected metabolites only}
\description{
This function performs metabolite set enrichment analysis using over-representation analysis (ORA) under the assumption that only detected metabolites are used as the background.
A one-sided Fisher's exact test is applied to each metabolite set.
}
\usage{
ora_det(SIG, DET, M)
}
\arguments{
  \item{SIG}{
Character vector of metabolite IDs considered statistically significant.
}
  \item{DET}{
Character vector of all detected metabolite IDs (background set).
}
  \item{M}{
A named list, where each element is a metabolite set (e.g., pathway) containing character vectors of metabolite IDs.
}
}
\value{
A list containing:
  \itemize{
    \item{\code{Result of MSEA(ORA)}: A matrix with raw p-values and adjusted q-values (BH correction) for each metabolite set.}
    \item{\code{significant metabolites}: A list of significant metabolites overlapping with each metabolite set.}
    \item{\code{TAB}: A list of 2x2 contingency tables used for Fisher's exact tests.}
  }
}
\author{
Hiroyuki Yamamoto
}
\references{
Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.\\
Global functional profiling of gene expression.\\
\emph{Genomics}. 2003 Feb;81(2):98-104.
}
\examples{
# Example: Simple ORA with dummy data
data(fasting)
data(pathway)

# Detected and significant metabolites
metabolites <- colnames(fasting$X)
pca <- prcomp(fasting$X, scale=TRUE)
pca <- pca_loading(pca)
SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]
DET <- metabolites
M <- pathway$fasting

# Perform ORA using detected metabolites only
B <- ora_det(SIG, DET, M)
B$`Result of MSEA(ORA)`
}
\keyword{msea}
