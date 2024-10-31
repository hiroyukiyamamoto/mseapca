library(mseapca)

### Subramanian QEA
msea_sub2 <- function (R, M, maxiter = 100) 
{

  S_name <- names(M)
  M_ID  <- names(R)
  L <- setlabel(M_ID, M)
  L <- L[, colSums(L) != 0]
  if (ncol(L) < 2) {
    stop("more than two metabolite set are necessary")
    return()
  }
  ESMAT <- NaN

  # ---
  ro <- R
  d <- sort.list(ro, decreasing = TRUE)
  r <- sort(ro, decreasing = TRUE)
  Lo <- L[d, ]
  Mo <- M_ID[d]
  p <- 1
  NR <- (abs(r)^(p)) %*% Lo
  ESo <- NaN
  Z <- NaN
  for (S in 1:length(NR)) {
    NH <- colSums(Lo)
    N <- nrow(Lo)
    Phit <- NaN
    Pmiss <- NaN
    for (i in 1:length(r)) {
      Phit[i] <- (abs(r[1:i])^(p)) %*% Lo[1:i, S]/NR[S]
      Pmiss[i] <- sum(Lo[1:i, S] == 0) * (1/(N - NH[S]))
    }
    z <- Phit - Pmiss
    Z <- cbind(Z, z)
    maxdev <- max(abs(Phit - Pmiss))
    maxdev_index <- which.max(abs(Phit - Pmiss))
    ESo[S] = sign(Phit[maxdev_index] - Pmiss[maxdev_index]) * maxdev
  }
  
  # ---
  
  Z <- Z[, -1]
  for (iter in 1:maxiter) {
    Lr <- NaN
    for (i in 1:ncol(L)) {
      Lr <- cbind(Lr, L[sample(d), i])
    }
    Lr <- Lr[, -1]
    NR <- (abs(r)^(p)) %*% Lr
    ESr <- NaN
    for (S in 1:length(NR)) {
      NH <- colSums(Lr)
      N <- nrow(Lr)
      Phit <- NaN
      Pmiss <- NaN
      for (i in 1:length(r)) {
        Phit[i] <- (abs(r[1:i])^(p)) %*% Lr[1:i, S]/NR[S]
        Pmiss[i] <- sum(Lr[1:i, S] == 0) * (1/(N - NH[S]))
      }
      maxdev <- max(abs(Phit - Pmiss))
      maxdev_index <- which.max(abs(Phit - Pmiss))
      ESr[S] <- sign(Phit[maxdev_index] - Pmiss[maxdev_index]) * 
        maxdev
    }
    ESMAT <- cbind(ESMAT, ESr)
  }
  ESMAT <- ESMAT[, -1]
  P <- NaN
  for (i in 1:nrow(ESMAT)) {
    esmat <- ESMAT[i, ]
    if (ESo[i] >= 0) {
      p_esmat <- esmat[esmat >= 0]
      P[i] <- sum(p_esmat >= ESo[i])/length(p_esmat)
    }
    if (ESo[i] < 0) {
      n_esmat <- esmat[esmat <= 0]
      P[i] <- sum(n_esmat <= ESo[i])/length(n_esmat)
    }
  }
  NESMAT <- NaN
  nesmat <- NaN
  NES <- NaN
  for (i in 1:nrow(ESMAT)) {
    esmat <- ESMAT[i, ]
    p_esmat <- esmat[esmat >= 0]
    n_esmat <- esmat[esmat <= 0]
    nesmat[esmat >= 0] <- p_esmat/mean(p_esmat)
    nesmat[esmat < 0] <- n_esmat/abs(mean(n_esmat))
    NESMAT <- cbind(NESMAT, nesmat)
    if (ESo[i] >= 0) {
      NES[i] <- ESo[i]/mean(p_esmat)
    }
    if (ESo[i] < 0) {
      NES[i] <- ESo[i]/abs(mean(n_esmat))
    }
  }
  NESMAT <- NESMAT[, -1]
  Q1 <- NaN
  Q2 <- NaN
  K <- length(NES)
  for (l in 1:length(NES)) {
    if (ESo[l] >= 0) {
      Q1[l] <- sum(sum((NES[l] <= NESMAT)))/sum(sum((NESMAT >= 0)))
    }
    if (ESo[l] < 0) {
      Q1[l] <- sum(sum((NES[l] >= NESMAT)))/sum(sum((NESMAT <= 0)))
    }
    if (NES[l] >= 0) {
      Q2[l] <- sum(NES[l] <= NES)/sum(NES >= 0)
    }
    if (NES[l] < 0) {
      Q2[l] <- sum(NES[l] >= NES)/sum(NES <= 0)
    }
  }
  Q <- Q1/Q2
  Q[Q > 1] <- 1
  PQ <- cbind(ESo, NES, P, Q)
  colnames(PQ) <- c("ES","normalized enrichment score", "p-value", 
                    "q-value")
  rownames(PQ) <- colnames(Lo)
  M <- NaN
  for (i in 1:ncol(Z)) {
    if (ESo[i] > 0) {
      m <- max(Z[, i])
      n <- which.max(Z[, i])
    }
    if (ESo[i] < 0) {
      m <- min(Z[, i])
      n <- which.min(Z[, i])
    }
    M[i] <- n
  }
  LES <- NaN
  for (i in 1:ncol(Lo)) {
    if (ESo[i] > 0) {
      m <- Mo[1:M[i]]
      les <- m[Lo[1:M[i], i] == 1]
    }
    if (ESo[i] < 0) {
      m <- Mo[-(1:M[i] - 1)]
      les <- m[Lo[M[i]:nrow(Lo), i] == 1]
    }
    LES[i] <- list(les)
  }
  names(LES) <- colnames(Lo)
  g <- list(PQ, LES)
  names(g) <- c("Result of MSEA", "leading edge subset")
  return(g)
}

