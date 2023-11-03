setlabel <- function (M_ID, M) 
{
  L <- NaN
  for (i in 1:length(M)) {
    m <- as.character(unlist(M[i]))
    m <- unique(m)
    l <- NaN
    for (j in 1:length(M_ID)) {
      a <- chartr(",", ";", M_ID[j])
      b <- unlist(strsplit(a, ";"))
      b <- unique(b)
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
