\name{setlabel}
\alias{setlabel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate binary label matrix of metabolite set}
\description{
This function generates binary label matrix of metabolite names and metabolite sets.
This is mainly used to be called by other functions, and used
to count the number of metabolites in a specific metabolite set.
}
\usage{
setlabel(M_ID, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M_ID}{
detected metabolites
}
  \item{M}{
list of metabolite set and metabolite names
}
}
\details{
If single peak has multiple metabolite IDs in M_ID, split by "," or ";".
}
\value{
binary label matrix of metabolite names in metabolite sets
}
\author{
Hiroyuki Yamamoto
}
\examples{
data(fasting)
data(pathway)

M_ID <- colnames(fasting$X) # detected metabolites
M <- pathway$fasting # metabolite set list

L <- setlabel(M_ID,M)	# binary label matrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
