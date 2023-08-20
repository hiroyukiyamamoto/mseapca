\name{read_pathway}
\alias{read_pathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read metabolite set file (*.xml)}
\description{
This function generates metabolite set list from metabolite set file (XML). 
This is mainly used to be called by other functions.
}
\usage{
read_pathway(fullpath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullpath}{
	file path of metabolite set (XML)
}
}
\value{
	list of metabolite set name and metabolite IDs.
}
\author{
Hiroyuki Yamamoto
}
\examples{
	# filename <- "C:/R/pathway.xml"	# load metabolite set file
	# M <- read_pathway(filename)		# Convert XML to metabolite set (list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }
