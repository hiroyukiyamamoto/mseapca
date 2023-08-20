ssea_ora <- function(detected, M){

  ## SSEA
  ALL <- detected # detected metabolites
  umet <- unique(as.character(unlist(M))) # all metabolites in pathways

  # ORA for each samples
  Z <- NULL
  for(i in 1:nrow(X)){
    index <- which(X[i,]>0) # threshold = 0 (detected metabolites)

    B <- msea_ora(ALL[index], umet, M)
    z <- -log(B$`Result of MSEA(ORA)`[,1])

    Z <- rbind(Z,z)
  }

  return(Z)
}
