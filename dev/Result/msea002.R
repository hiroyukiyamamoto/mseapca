library(loadings)

data(fasting)
X <- fasting$X
pca <- prcomp(X, scale=TRUE)
pca <- pca_loading(pca)
R1 <- pca$loading$R[,1] # PC loading

# ------------------------------------------

# 検出された代謝物と有意代謝物の情報
library(mseapca)

metabolites <- colnames(fasting$X)

ALL <- unique(unlist(pathway$fasting))
M <- pathway$fasting

index <- which(!ALL %in% names(R1))
## 検出できない物質に対して、ランダムサンプリング
R2 <- sample(R1,length(index), replace = TRUE)
names(R2) <- ALL[index]

R <- c(R1,R2)

### SubramanianのGSEAを計算
source("C:/Users/yamamoto/Documents/R/msea/msea_sub_core.R")

ALLB <- NULL
for(i in 1:10){
  print(i)
  B <- msea_sub_core(R, M)
  ALLB[[i]] <- B[[1]]
}

BB <- NULL
for(i in 1:10){
  BB <- cbind(BB,ALLB[[i]][,3])
}

B0 <- msea_sub_core(R1, M)
