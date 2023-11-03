\name{msea_ora}
\alias{msea_ora}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MSEA by over representation analysis}
\description{
This function performs metabolite set enrichment analysis by over representation analysis (ORA).
Statistical hypothesis test of cross tabulation is performed by one-sided Fisher's exact test.
}
\usage{
msea_ora(SIG, ALL, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SIG}{
Metabolite names of significant metabolites
}
  \item{ALL}{
Metabolite names of all detected metabolites
}
  \item{M}{
list of metabolite set name and metabolite name
}
}
\value{
list of p-value and q-value for metabolite set and selected (significant) metabolite IDs for each metabolite set
}
\references{
Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.
Global functional profiling of gene expression.
Genomics. 2003 Feb;81(2):98-104.
}
\author{
Hiroyuki Yamamoto
}
\examples{
## Example1 : Metabolome data
data(fasting)
data(pathway)

# pca and pca loading
pca <- prcomp(fasting$X, scale=TRUE)
pca <- pca_loading(pca)

# all detected metabolites
metabolites <- colnames(fasting$X)

# statistically significant negatively correlated metabolites in PC1 loading
SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]
ALL <- metabolites #all detected metabolites

# metabolite set list
M <- pathway$fasting

# MSEA by over representation analysis
B <- msea_ora(SIG, ALL, M)
B$`Result of MSEA(ORA)`

## Example2 : Proteome data
data(covid19)
data(pathway)

X <- covid19$X$proteomics
Y <- covid19$Y
D <- covid19$D
tau <- covid19$tau

protein_name <- colnames(X)

# pls-rog and pls-rog loading
plsrog <- pls_rog(X,Y,D)
plsrog <- plsrog_loading(plsrog)

# statistically significant proteins
index_prot <- which(plsrog$loading$R[,1]>0 & plsrog$loading$p.value[,1]<0.05)
sig_prot <- protein_name[index_prot]

# detected proteins
protein_name <- colnames(X)

# protein set list
M <- pathway$covid19$proteomics

# MSEA by over representation analysis
B <- msea_ora(sig_prot, protein_name, M)
B$`Result of MSEA(ORA)`

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
