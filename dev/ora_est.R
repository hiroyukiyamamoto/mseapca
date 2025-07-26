ora_est <- function (SIG, DET, M, method = "naive", lambda = 5) 
{
  ALL0 <- unique(unlist(M))
  ALL <- unique(c(ALL0, DET))
  L1 <- setlabel(ALL, M)
  L2 <- setlabel(DET, M)
  Lsig <- setlabel(SIG, M)
  B <- ora_det(SIG, DET, M)
  P <- NULL
  P_range <- NULL
  p <- length(SIG)/length(DET)
  TAB <- NULL
  for (i in 1:length(M)) {
    l1 <- colSums(L1)[i]
    l2 <- colSums(L2)[i]
    l3 <- colSums(Lsig)[i]
    n <- l1 - l2
    if (method == "naive") {
      r <- if (l2 > 0) 
        l3/l2
      else p
    }
    else if (method == "weighted") {
      r <- (l3 + n * p)/l1
    }
    else if (method == "shrink") {
      r <- if (l2 + lambda > 0) 
        (l3 + lambda * p)/(l2 + lambda)
      else p
    }
    a <- round(B$`Contingency tables`[[i]][1, 1] + n * r)
    a <- min(a, l1)
    #b <- round(B$`Contingency tables`[[i]][1, 2] + n * (1 - r))
    b <- l1 - a
    c <- max(0, round(length(ALL) * p - a))
    d <- max(0, round(length(ALL) * (1 - p) - b))
    tab <- t(matrix(c(a, b, c, d), nrow = 2))
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
  names(RES) <- c("Result of MSEA (ORA with adjustment)", "significant metabolites", 
                  "Contingency tables")
  return(RES)
}