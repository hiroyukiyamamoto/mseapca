setlabel <- function(M_ID, M, option = NULL) {
  n <- length(M_ID)
  p <- length(M)

  L <- matrix(0, nrow = n, ncol = p)
  colnames(L) <- names(M)

  for (i in seq_len(p)) {
    m <- unique(as.character(unlist(M[[i]])))

    for (j in seq_len(n)) {
      if (is.null(option)){
        b <- M_ID[j]
      }
      else if (option == "anno") {
        a <- chartr(",", ";", M_ID[j])
        b <- unique(unlist(strsplit(a, ";")))
      } else {
        stop("Unknown option: ", option)
      }

      if (any(b %in% m)) {
        L[j, i] <- 1
      }
    }
  }

  return(L)
}
