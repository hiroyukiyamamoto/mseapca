\name{list2xml}
\alias{list2xml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Save compound set as XML file}
\description{
This function save compound set of list format as XML file.
}
\usage{
list2xml(filepath, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
  filepath of XML file to save
}
  \item{M}{
  list fomat of compound set and compound names
  }
}
\details{
This function is used to store a compound set. Saved xml file can be read using the read_pathway function.
}
\value{
filepath of saved XML file
}
\author{
Hiroyuki Yamamoto
}
\examples{
\dontrun{
	data(pathway)
	M <- pathway$fasting
	xml_file <- "pathway_fasting.xml"
	N <- list2xml(xml_file, M)
	# XML::saveXML(N,filepath)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
