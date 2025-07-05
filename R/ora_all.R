# ora_all: over-representation analysis under the same conditions as MetaboAnalyst,
# when no reference metabolome is specified.
ora_all <- function (SIG, M)
{
  ALL <- unique(as.character(unlist(M)))
  SIG <- as.character(as.matrix(SIG))
  SIG <- SIG[SIG %in% ALL]
  num_all <- length(ALL)
  num_sig <- length(SIG)
  Lall <- setlabel(ALL, M)
  if (ncol(Lall) < 2) {
    stop("more than two metabolite set are necessary")
    return()
  }
  Lsig <- setlabel(SIG, M)
  P <- NULL
  TAB <- NULL
  for (i in 1:length(M)) {
    a1 <- sum(Lsig[, i])
    a2 <- sum(Lall[, i]) - sum(Lsig[, i])
    a3 <- length(SIG) - a1
    a4 <- (length(ALL) - length(SIG)) - a2
    tab <- t(matrix(c(a1, a2, a3, a4), 2))
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- resfish$p.value
    TAB[i] <- list(tab)
  }
  names(TAB) <- colnames(Lsig)
  Q <- p.adjust(P, method = "BH")
  LES <- NULL
  for (i in 1:ncol(Lsig)) {
    les <- SIG[Lsig[, i] == 1]
    LES[i] <- list(les)
  }
  names(LES) <- colnames(Lsig)
  PQ <- cbind(P, Q)
  rownames(PQ) <- colnames(Lsig)
  colnames(PQ) <- c("p.value", "q.value")
  RES <- list(PQ, LES, TAB)
  names(RES) <- c("Result of MSEA(ORA)", "significant metabolites", "Contingency tables")
  return(RES)
}
