\name{msea_ora}
\alias{msea_ora}
\title{Wrapper for several over-representation analysis (ORA) flavours}
\description{
\code{msea_ora} is a convenience wrapper that calls one of five ORA
implementations, allowing the user to choose how undetected metabolites are
handled.

\itemize{
  \item \code{det}  – ORA using \emph{detected metabolites only} as background
        (\code{\link{ora_det}}).
  \item \code{all}  – ORA using \emph{all metabolites appearing in the
        pathway list} as background (\code{\link{ora_all}}).
  \item \code{est_naive}, \code{est_weighted}, \code{est_shrink} –
        ORA adjusted for undetected metabolites by the
        “naive”, “weighted”, or “shrinkage” estimator
        (all via \code{\link{ora_est}}).
}
}
\usage{
msea_ora(SIG, DET = NULL, M,
         option = "det", lambda = NULL)
}
\arguments{
  \item{SIG}{Character vector of significant metabolite IDs.}
  \item{DET}{Character vector of detected metabolite IDs (background). Required for every option except \code{"all"}.}
  \item{M}{Named list of metabolite sets; each element is a character vector of metabolite IDs.}
  \item{option}{One of \code{"det"}, \code{"all"},
                \code{"est_naive"}, \code{"est_weighted"}, or \code{"est_shrink"}.}
  \item{lambda}{Shrinkage constant used only when \code{option = "est_shrink"}.}
}
\value{
Whatever object is returned by the chosen back-end function:
\describe{
  \item{\code{det} / \code{all}}{A list with p-values, q-values, tables, etc.\ from
        \code{\link{ora_det}} or \code{\link{ora_all}}.}
  \item{\code{est_*}}{A list produced by \code{\link{ora_est}}.}
}
}
\seealso{
\code{\link{ora_det}}, \code{\link{ora_all}}, \code{\link{ora_est}}.
}
\author{Hiroyuki Yamamoto}
\examples{
## Toy data ----------------------------------------------------
M   <- list(set1 = c("A","B","C"),
            set2 = c("C","D","E"))
DET <- c("A","B","C","D")
SIG <- c("A","C")

## Plain ORA with detected metabolites
msea_ora(SIG, DET, M, option = "det")

## ORA with shrinkage adjustment
msea_ora(SIG, DET, M, option = "est_shrink", lambda = 5)
}
\keyword{MSEA}
\keyword{ORA}
\keyword{enrichment}
