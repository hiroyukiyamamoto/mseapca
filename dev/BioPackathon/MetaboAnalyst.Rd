\name{MetaboAnalyst}
\alias{MetaboAnalyst}
\docType{data}
\title{Example metabolomics data extracted from MetaboAnalyst (fasting study)}
\description{
A list containing example data for demonstrating the ORA and MSEA functions
implemented in the \pkg{mseapca} package.
The data originate from the fasting study included in the MetaboAnalyst
web platform.
}
\format{
A named \code{list} with one top-level element, \code{$fasting},
which itself is a list with the following components:

\describe{
  \item{\code{$data}}{A list with two character vectors:
    \describe{
      \item{\code{SIG}}{Metabolite IDs judged significant in the fasting study.}
      \item{\code{DET}}{All metabolite IDs detected in the experiment.}
    }}
  \item{\code{$pathway}}{A named list of metabolite sets (pathways).
    Each element is a character vector of metabolite IDs belonging to that pathway.
    The names of the list elements correspond to pathway names.}
}
}
\usage{
data(MetaboAnalyst)
}
\details{
The dataset is intended for vignette examples:

\itemize{
  \item \code{SIG} — Significant metabolites (\eqn{p < 0.05}) obtained from
        statistical analysis of a fasting vs.~control comparison.
  \item \code{DET} — Background list of all detected metabolites.
  \item \code{pathway} — Human pathway definitions curated in MetaboAnalyst.
}

The object can be used directly with functions such as
\code{\link{msea_ora}}, \code{\link{ora_det}}, and \code{\link{ora_bino}}.
}
\source{
Downloaded from the MetaboAnalyst web server
(\url{https://www.metaboanalyst.ca/}).
}
\examples{
## Load data
data(MetaboAnalyst)

SIG <- MetaboAnalyst$fasting$data$SIG
DET <- MetaboAnalyst$fasting$data$DET
M   <- MetaboAnalyst$fasting$pathway

## Simple ORA on detected metabolites
res <- ora_det(SIG, DET, M)
head(res$`Result of MSEA(ORA)`)
}
\keyword{datasets}
\keyword{metabolomics}
