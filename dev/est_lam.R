est_lam <- function(k, n, m0 = NULL, method = "fct2", Kmin = 3) {
  if (method == "fct2") {
    j <- (n >= 2) & (k >= Kmin)
    t <- mean((k[j]*(k[j]-1))/(n[j]*(n[j]-1)), na.rm=TRUE)
    lam <- (m0 - t)/(t - m0^2)
  }
  
  if(method == "var") {
    j <- (n > 0) & (k >= Kmin)
    sbar <- mean(((k[j]-n[j]*m0)^2)/(n[j]*m0*(1-m0)), na.rm=TRUE)
    nbar <- mean(n[j], na.rm=TRUE)
    lam <- (nbar - sbar)/(sbar - 1)
  }
  
  return(lam)
}

### demo
## フィッシャー正確確率検定
#m0 <- length(SIG)/length(DET)
#lam <- est_lam(k0, n0, m0)
#B <- ora_est(SIG, DET, M, method="shrink", lambda=lam)
#sum(B$`Result of MSEA (ORA with adjustment)`[,1]<0.05)

## ベイズ事後確率
#a_post <- lam * m0 + k0
#b_post <- lam * (1 - m0) + (n0 - k0)
#post_enrich <- 1 - pbeta(m0, a_post, b_post) ## 3)