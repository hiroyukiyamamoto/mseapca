setlabel <- function(MET, M) {
  n <- length(MET)
  p <- length(M)

  L <- matrix(0, nrow = n, ncol = p)
  colnames(L) <- names(M)

  for (i in seq_len(p)) {
    m <- unique(as.character(unlist(M[[i]])))

    for (j in seq_len(n)) {
      b <- MET[j]
      if (b %in% m) {
        L[j, i] <- 1
      }
    }
  }
  return(L)
}
