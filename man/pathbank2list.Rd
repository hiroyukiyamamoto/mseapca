\name{pathbank2list}
\alias{pathbank2list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate metabolite set list from PathBank database}
\description{
This function generates metabolite set list of PathBank database by referencing the AHPathbankDbs Bioconductor package.
}
\usage{
pathbank2list(tbl_pathbank, subject, id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbl_pathbank}{
  tibble from AHPathbankDbs
  }
  \item{subject}{
  Pathway subject (Metabolic, Disease, etc.) in tibble
  }
  \item{id}{
  database ID (HMDB ID, Uniprot ID, etc.) used for analysis
  }
}
\details{
AHPathbankDbs needs to be installed separately.
}
\value{
list of metabolite or protein set
}
\author{
Hiroyuki Yamamoto
}
\examples{
\dontrun{
## PathBank
library(AnnotationHub)

ah <- AnnotationHub()
qr <- query(ah, c("pathbank", "Homo sapiens"))

#tbl_pathbank <- qr[[1]] # metabolomics
tbl_pathbank <- qr[[2]] # proteomics

ids <- names(tbl_pathbank)[-c(1:4)]
id <- ids[1] # Uniprot ID

subs <- unique(tbl_pathbank$`Pathway Subject`)
subject <- subs[6] # Protein

M <- pathbank2list(tbl_pathbank, subject, id)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
