\name{ora_est}
\alias{ora_est}
\title{Over-representation analysis adjusting for undetected metabolites}
\description{
This function performs metabolite set enrichment analysis using over-representation analysis (ORA), incorporating point estimates to adjust for potentially undetected metabolites.
It supports three estimation methods: naive, weighted, and shrinkage-based adjustment.
}
\usage{
ora_est(SIG, DET, M, method = "naive", lambda = 5)
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
    \item{\code{Result of MSEA (ORA with adjustment)}: A matrix with raw p-values and adjusted q-values (BH correction) for each metabolite set.}
  }
}
\details{
This function estimates the impact of undetected metabolites on enrichment results. It builds upon the ORA results from detected metabolites, then adjusts the contingency tables by estimating how many undetected metabolites might be significant, based on a specified method.
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
# Example using metabolomics data
data(fasting)
data(pathway)

# Compute PCA loadings
pca <- prcomp(fasting$X, scale=TRUE)
pca <- pca_loading(pca)

# Detected and significant metabolites
metabolites <- colnames(fasting$X)
SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]
DET <- metabolites
M <- pathway$fasting

# Run adjusted ORA (e.g., shrinkage method)
B <- ora_est(SIG, DET, M, method = "shrink", lambda = 5)
B$`Result of MSEA (ORA with adjustment)`
}
\keyword{msea}
