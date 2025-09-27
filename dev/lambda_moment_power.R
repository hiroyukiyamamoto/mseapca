## ===== method_moment_power.R =====
compute_lambda_moment_power <- function(k, n, N, m, ...){
  U <- pmax(N-n,0); eps <- 1e-9
  v_model <- function(lam) n*m*(1-m)*((lam+n)/pmax(lam+1, eps))
  par0  <- c(log(1e-6), log(1), 2.0, log(0.5))  # (lmin, c, alpha, s)
  obj <- function(par){
    lmin <- exp(par[1]); c <- exp(par[2]); alpha <- par[3]; s <- exp(par[4])
    lam <- pmax(lmin + c*(pmax(U+s,1e-3))^alpha, eps)
    v_th <- v_model(lam); resid <- (k - n*m)^2 - v_th
    w <- (1/pmax(v_th,eps))^0.5; w <- pmin(w,1e3)
    sum((sqrt(w)*resid)^2)
  }
  fit <- optim(par0, obj, method="L-BFGS-B",
               lower=c(log(1e-6), log(1e-6), 0.5, log(1e-3)),
               upper=c(log(1e3),  log(1e6),  8,   log(5)))
  lmin <- exp(fit$par[1]); c <- exp(fit$par[2]); alpha <- fit$par[3]; s <- exp(fit$par[4])
  pmax(lmin + c*(pmax(U+s,1e-3))^alpha, eps)
}
