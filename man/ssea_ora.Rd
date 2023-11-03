\name{ssea_ora}
\alias{ssea_ora}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single sample enrichment analysis by over representation analysis}
\description{
This function performs single sample enrichment analysis (SSEA) by over representation analysis (ORA). SSEA performs MSEA by ORA between detected and not detected metabolites in each sample."
}
\usage{
ssea_ora(det_list, det_all, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{det_list}{
  metabolite names of detected metabolites
  }
  \item{det_all}{
  metabolite names of all metabolites
  }
  \item{M}{
  list of metabolite set and metabolite names
  }
}
\details{
The threshold for determining whether a metabolite is detected or not is typically set by the signal-to-noise (S/N) ratio. If the S/N ratio is unavailable, one might consider using the signal intensity or peak area for each metabolite as an alternative. In such cases, all values below the threshold can be set to 0.
}
\value{
A matrix where each row represents a sample and each column represents a set of metabolites.
}
\author{
Hiroyuki Yamamoto
}
\references{
Yamamoto H., Single sample enrichment analysisfor mass spectrometry-based omics data, Jxiv.(2023)
}
\examples{
\dontrun{
data(fasting)
data(pathway)

det_list <- pathway$data$fasting
M <- pathway$fasting
det_all <- unique(c(colnames(fasting$X), as.character(unlist(M)))) 

# SSEA
Z <- ssea_ora(det_list, det_all, M)

## PCA for SSEA score
pca <- prcomp(Z, scale=TRUE)
pca <- pca_loading(pca)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }
