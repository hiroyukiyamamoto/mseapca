msea_ora <- function (SIG, ALL, M) 
{
  ALL <- as.character(as.matrix(ALL))
  SIG <- as.character(as.matrix(SIG))
  num_all <- length(ALL)
  num_sig <- length(SIG)
  Lall0 <- setlabel(ALL, M)
  Lall <- Lall0[, colSums(Lall0) != 0]
  if (ncol(Lall) < 2) {
    stop("more than two metabolite set are necessary")
    return()
  }
  Lsig <- setlabel(SIG, M)
  Lsig <- Lsig[, colSums(Lall0) != 0]
  l <- colSums(Lall0) != 0
  P <- NaN
  TAB <- NULL
  for (i in 1:sum(l)) {
    a1 <- sum(Lsig[, i])
    a2 <- sum(Lall[, i]) - sum(Lsig[, i])
    a3 <- length(SIG) - a1
    a4 <- (length(ALL) - length(SIG)) - a2
    tab <- t(matrix(c(a1, a2, a3, a4), 2))
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- resfish$p.value
    TAB[i] <- list(tab)
  }
  Q <- p.adjust(P, method = "BH")
  LES <- NaN
  for (i in 1:ncol(Lsig)) {
    les <- SIG[Lsig[, i] == 1]
    LES[i] <- list(les)
  }
  names(LES) <- colnames(Lsig)
  PQ <- cbind(P, Q)
  rownames(PQ) <- colnames(Lsig)
  colnames(PQ) <- c("p.value", "q.value")
  RES <- list(PQ, LES,TAB)
  names(RES) <- c("Result of MSEA(ORA)", "significant metabolites", "TAB")
  return(RES)
}