\name{ora_all}
\alias{ora_all}
\title{ORA using all metabolites in the pathway}
\description{
This function performs over-representation analysis (ORA) using all metabolites present in the given metabolite set list as the background, without specifying a reference metabolome. This corresponds to the behavior of MetaboAnalyst when no reference metabolome is uploaded.
}
\usage{
ora_all(SIG, M)
}
\arguments{
  \item{SIG}{Character vector of significant metabolites}
  \item{M}{Named list of metabolite sets}
}
\value{
A list with:
  \item{\code{Result of MSEA(ORA)}}{Matrix of p-values and q-values}
  \item{\code{significant metabolites}}{List of significant IDs per set}
  \item{\code{Contingency tables}}{A list of 2×2 contingency tables used in Fisher's exact tests.}
}
\references{
Draghici S, Khatri P, Martins RP, Ostermeier GC, Krawetz SA.
Global functional profiling of gene expression.
Genomics. 2003 Feb;81(2):98-104.

Yamamoto H, Fujimori T, Sato H, Ishikawa G, Kami K, Ohashi Y,
Statistical hypothesis testing of factor loading in principal component analysis and its application to metabolite set enrichment analysis.
BMC Bioinformatics, (2014) 15(1):51.

Yamamoto H.
Probabilistic Over-Representation Analysis for Metabolite Set Enrichment Analysis Considering Undetected Metabolites", Jxiv, (2024).

}
\author{Hiroyuki Yamamoto}
\examples{
# Example1 : Metabolome data
data(fasting)
data(pathway)

# pca and pca loading
pca <- prcomp(fasting$X, scale=TRUE)
pca <- pca_loading(pca)

# all detected metabolites
metabolites <- colnames(fasting$X)

# statistically significant negatively correlated metabolites in PC1 loading
SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]

# metabolite set list
M <- pathway$fasting

# MSEA by over representation analysis
B <- ora_all(SIG, M)
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

# protein set list
M <- pathway$covid19$proteomics

# MSEA by over representation analysis
B <- ora_all(sig_prot, M)
B$`Result of MSEA(ORA)`

# Example3: Metabolome data
data(MetaboAnalyst)

SIG <- MetaboAnalyst$fasting$data$SIG
M <- MetaboAnalyst$fasting$pathway

# Perform ORA using detected metabolites only
B <- ora_all(SIG, M)
B$`Result of MSEA(ORA)`

}
\keyword{ htest }
