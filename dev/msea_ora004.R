rm(list=ls(all=TRUE))

source("C:/Users/hyama/Documents/R/msea/msea_ora.R")

library(mseapca)

## Example1 : Metabolome data
data(fasting)
data(pathway)

# pca and pca loading
pca <- prcomp(fasting$X, scale=TRUE)
pca <- pca_loading(pca)

# all detected metabolites
metabolites <- colnames(fasting$X)

DET0 <- metabolites

# statistically significant negatively correlated metabolites in PC1 loading
SIG <- metabolites[pca$loading$R[,1] < 0 & pca$loading$p.value[,1] < 0.05]
SIG0 <- SIG

# Annotation
file_annotation <- "C:/Users/hyama/Documents/R/msea/MetaboAnalyst/fasting_annotation.csv"
anno_fasting <- read.csv(file_annotation)

### 対象物質の選択(有意)
index_SIG <- (anno_fasting$Query %in% SIG) & !is.na(anno_fasting$Match) & anno_fasting$HMDB!=""
SIG <- anno_fasting$Match[index_SIG]

DSIG <- anno_fasting$Query[anno_fasting$Query %in% SIG0 & !index_SIG] # 3物質
SIG <- c(SIG,DSIG) # 89物質

### 対象物質の選択(検出物質)
DET <- metabolites

index_DET <- (anno_fasting$Query %in% DET) & !is.na(anno_fasting$Match) & anno_fasting$HMDB!=""
DET <- anno_fasting$Match[index_DET] # 検出された物質の名前(HMDB)
DIF <- anno_fasting$Query[!index_DET] # 検出された物質の中で、データベースに名前がない物質
DET <- c(DET,DIF) # 89物質

# ----- 以下、MetaboAnalystR ----- #

library(MetaboAnalystR)

mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
pathway_hmdb <- current.msetlib$member # HMDBの名前

ALL0 <- unique(unlist(pathway_hmdb)) # 対象物質の情報
ALL <- unique(c(ALL0, DET))

# metabolite set list
M <- pathway_hmdb

# setlabel関数を修正しないといけなさそう

setlabel <- function (M_ID, M) 
{
  L <- NaN
  for (i in 1:length(M)) {
    m <- as.character(unlist(M[i]))
    m <- unique(m)
    l <- NaN
    for (j in 1:length(M_ID)) {
      b <- M_ID[j]
      c <- b %in% m
      l[j] <- 0
      if (sum(c) >= 1) {
        l[j] <- 1
      }
    }
    L <- cbind(L, l)
  }
  L <- L[, -1]
  colnames(L) <- names(M)
  return(L)
}

M <- append(M,list(metabolites[!index_DET]))
names(M)[length(M)] <- "not match metabolites" 

# ORAの計算
#B0 <- msea_ora(SIG, ALL, M) # MetaboAnalyst

B <- msea_ora(SIG, DET, M)
M <- M[names(M) %in% rownames(B$`Result of MSEA(ORA)`) ]

# Data
# SIG, DET, ALL, M

M <- M[-length(M)]
save(SIG, DET, ALL, M, file="C:/R/ora_fasting.RData")

### ORA
source("C:/Users/hyama/Documents/R/msea/msea_ora_point.R")
B1 <- msea_ora_point(SIG, DET, ALL, M)

source("C:/Users/hyama/Documents/R/msea/msea_ora_all_combination.R")
B2 <- msea_ora_all_combination(SIG, DET, ALL, M)
#B3 <- msea_ora_all_combination(SIG, DET, ALL, M[1:3], "range")

source("C:/Users/hyama/Documents/R/msea/msea_ora_binomial_ci.R")
B4 <- msea_ora_binomial_ci(SIG, DET, ALL, M)

source("C:/Users/hyama/Documents/R/msea/msea_ora_binomial_ci.R")
B5 <- msea_ora_binomial_ci(SIG, DET, ALL, M, c(0.2,0.8))

#source("C:/Users/hyama/Documents/R/msea/msea_ora_binomial_ci2.R")
#B5 <- msea_ora_binomial_ci2(SIG, DET, ALL, M)

#source("C:/Users/hyama/Documents/R/msea/msea_ora_beta_ci.R")
#B5 <- msea_ora_beta_ci(SIG, DET, ALL, M)




