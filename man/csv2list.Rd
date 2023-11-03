\name{csv2list}
\alias{csv2list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert metabolite set / csv to list}
\description{
This function converts your own metabolite set (csv file to list).
}
\usage{
csv2list(filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
file path of metabolite set (csv file)
}
}
\details{
The first row of csv file are "metabolite set name" and "metabolite IDs" as header.
The first column must be metabolite IDs and second column must be metabolite set name.
}
\value{
list of metabolite set name and metabolite IDs
}
\author{
Hiroyuki Yamamoto
}
\examples{
\dontrun{
	# ---------------------------
	#  Convert csv file to list
	# ---------------------------
	filepath <- "C:/pathway.csv"	# filepath of csv file
	N <- csv2list(filepath)	# convert csv file to list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
