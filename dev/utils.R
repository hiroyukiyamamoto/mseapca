## ===== utils.R =====
.get_counts <- function(SIG, DET, M){
  N_vec <- vapply(M, function(s) length(unique(s)), integer(1))
  n_vec <- vapply(M, function(s) sum(unique(DET) %in% s), integer(1))
  k_vec <- vapply(M, function(s) sum(unique(SIG) %in% s), integer(1))
  lenALL <- length(unique(c(unique(unlist(M)), unique(DET))))
  list(N=N_vec, n=n_vec, k=k_vec, lenALL=lenALL)
}

.compute_m <- function(k_vec, n_vec, type=c("global","loo"), eps=1e-9){
  type <- match.arg(type)
  if (type=="global"){
    m <- sum(k_vec)/pmax(sum(n_vec), 1L)
    rep(pmin(pmax(m, eps), 1-eps), length(k_vec))
  } else {
    Kall <- sum(k_vec); Nall <- sum(n_vec)
    m_vec <- (Kall - k_vec) / pmax(Nall - n_vec, eps)
    pmin(pmax(m_vec, eps), 1-eps)
  }
}

fisher_from_phat <- function(p_hat, k, n, m, M, DET, N_vec){
  lenALL <- length(unique(c(unique(unlist(M)), unique(DET))))
  vapply(seq_along(M), function(i){
    n_miss <- N_vec[i]-n[i]; r <- p_hat[i]
    a <- round(k[i] + n_miss*r)
    b <- round(n[i]-k[i] + n_miss*(1-r))
    c <- max(0, round(lenALL*m[i]     - a))
    d <- max(0, round(lenALL*(1-m[i]) - b))
    tab <- t(matrix(c(a,b,c,d), nrow=2))
    fisher.test(tab, alternative="greater")$p.value
  }, numeric(1))
}
