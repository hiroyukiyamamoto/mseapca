fit_cs <- function(k,n,N,rho=1e-3){
  m <- length(SIG)/length(DET) # 事前分布の割合p0
  
  i <- is.finite(k)&is.finite(n)&is.finite(N)&(n>=2)
  
  y <- (k[i]*(k[i]-1))/(n[i]*(n[i]-1)) # 観測の2次モーメント
  U <- N[i]-n[i]
  
  f <- function(p){
    c <- exp(p[1]); s <- exp(p[2])
    lam <- s + c*U # lambda
    v <- mean((y - (m*(m*lam+1))/(lam+1))^2, na.rm=TRUE) + rho*mean(lam,na.rm=TRUE)
  }
  optim(c(0,0), f, method="L-BFGS-B", control=list(maxit=1000))
}

res  <- fit_cs(k,n,N, rho=1e-10)
cHat <- exp(res$par[1]); sHat <- exp(res$par[2])
lam  <- sHat + cHat*pmax(N-n,0)
