ssea_ora <- function (det_list, det_all, M) 
{
  Z <- NULL
  for (i in 1:length(det_list)) {
    B <- msea_ora(det_list[[i]], det_all, M)
    z <- -log(B$`Result of MSEA(ORA)`[, 1])
    Z <- rbind(Z, z)
  }
  return(Z)
}

