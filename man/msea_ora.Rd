\name{msea_ora}
\alias{msea_ora}
\title{Wrapper function for over-representation analysis (ORA)}
\description{
\code{msea_ora} is a wrapper function that calls ORA
implementations, allowing the user to choose how undetected metabolites are
handled.

\itemize{
  \item \code{"det"} – ORA using \emph{detected metabolites only} as background
        (\code{\link{ora_det}}).
  \item \code{"all"} – ORA using \emph{all metabolites appearing in the
        pathway list} as background (\code{\link{ora_all}}).
  \item \code{"est_naive"}, \code{"est_weighted"}, \code{"est_shrink"} –
        ORA adjusted for undetected metabolites by the
        "naive", "weighted", or "shrinkage" estimator
        (all via \code{\link{ora_est}}).
}
}
\usage{
msea_ora(SIG, DET, M, option = "det", lambda=NULL)
}
\arguments{
  \item{SIG}{Character vector of significant metabolites}
  \item{DET}{Character vector of detected metabolites. Required for all options except \code{"all"}.}
  \item{M}{A named list, where each element is a metabolite set (e.g., pathway) containing character vectors of metabolites.}
  \item{option}{One of \code{"det"}, \code{"all"},
                \code{"est_naive"}, \code{"est_weighted"}, or \code{"est_shrink"}.}
  \item{lambda}{Shrinkage parameter used only when \code{option = "est_shrink"}.}
}
\value{
A list with:
  \item{\code{Result of MSEA(ORA)}}{ Matrix of p-values and q-values}
  \item{\code{significant metabolites}}{ List of significant metabolites per set}
  \item{\code{Contingency tables}}{ A list of 2×2 contingency tables used in Fisher's exact tests.}
}
\seealso{
\code{\link{ora_det}}, \code{\link{ora_all}}, \code{\link{ora_est}}
}
\references{
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
DET <- metabolites

# Fix for multiple annotations
DET[DET == "UDP-glucose ; UDP-galactose"] <- "UDP-glucose"
DET[DET == "Isonicotinamide ; Nicotinamide"] <- "Nicotinamide"
DET[DET == "1-Methylhistidine ; 3-Methylhistidine"] <- "3-Methylhistidine"

# metabolite set list
M <- pathway$fasting

# MSEA by over representation analysis
B <- msea_ora(SIG, DET, M)
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
B <- msea_ora(sig_prot, protein_name, M)
B$`Result of MSEA(ORA)`

## Example3: Metabolome data
data(fasting_mseapca)

SIG <- fasting_mseapca$SIG
DET <- fasting_mseapca$DET
M <- fasting_mseapca$pathway

# Perform ORA using detected metabolites only
B <- msea_ora(SIG, DET, M)
B$`Result of MSEA(ORA)`

}
\keyword{ htest }
