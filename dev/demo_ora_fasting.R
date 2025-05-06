rm(list=ls(all=TRUE))

load(file="C:/R/ora_fasting.RData")
source("C:/Users/hyama/Documents/R/msea/msea_ora.R")
source("C:/Users/hyama/Documents/R/msea/setlabel.R")

### ORA
source("C:/Users/hyama/Documents/R/msea/msea_ora_point.R")
B1 <- msea_ora_point(SIG, DET, ALL, M)

source("C:/Users/hyama/Documents/R/msea/msea_ora_all_combination.R")
B2 <- msea_ora_all_combination(SIG, DET, ALL, M)
#B3 <- msea_ora_all_combination(SIG, DET, ALL, M[1:3], "range")

source("C:/Users/hyama/Documents/R/msea/msea_ora_binomial_ci.R")
B4 <- msea_ora_binomial_ci(SIG, DET, ALL, M)

source("C:/Users/hyama/Documents/R/msea/msea_ora_beta_ci.R")
B5 <- msea_ora_beta_ci(SIG, DET, ALL, M)