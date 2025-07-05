\name{setlabel}
\alias{setlabel}
\title{Generate binary label matrix of metabolite set}
\description{
This function generates binary label matrix of metabolites and metabolite sets.
This is mainly used to be called by other functions, and used
to count the number of metabolites in a specific metabolite set.
}
\usage{
setlabel(MET, M)
}
\arguments{
  \item{MET}{A character vector of metabolites (e.g., detected or significant metabolites).}
  \item{M}{A named list of metabolite sets, where each element is a character vector of metabolites.}
}
\value{
binary label matrix of metabolites in metabolite sets
}
\details{
This function is used internally in various ORA methods (e.g., \code{ora_det}, \code{ora_all}, \code{ora_est}) to compute contingency tables for enrichment analysis.
}
\author{
Hiroyuki Yamamoto
}
\examples{
# Example
data(fasting)
data(pathway)

MET <- colnames(fasting$X) # detected metabolites
M <- pathway$fasting # metabolite set list

# Fix for multiple annotations
MET[MET == "UDP-glucose ; UDP-galactose"] <- "UDP-glucose"
MET[MET == "trans-Glutaconic acid ; Itaconic acid"] <- "Itaconic acid"
MET[MET == "Isonicotinamide ; Nicotinamide"] <- "Nicotinamide"
MET[MET == "Isobutyric acid ; Butyric acid"] <- "Isobutyric acid"
MET[MET == "GDP-mannose ; GDP-galactose"] <- "GDP-mannose"
MET[MET == "ADP-glucose ; GDP-fucose"] <- "ADP-glucose"
MET[MET == "1-Methylhistidine ; 3-Methylhistidine"] <- "3-Methylhistidine"

L <- setlabel(MET, M)

# Example 2
data(fasting_mseapca)

MET <- fasting_mseapca$DET
M <- fasting_mseapca$pathway

L <- setlabel(MET, M)
}
\keyword{utilities}
