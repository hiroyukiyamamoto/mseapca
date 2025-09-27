## ===== transforms.R =====
build_lambda_from_weighted <- function(
    U, target_mean=mean(U), target_var_frac=0.10,
    model=c("affine","power","scale"), floor=1e-9,
    power_alpha_range=c(0.5,3.5), power_s=NULL
){
  model <- match.arg(model); U <- pmax(U,0)
  vU <- stats::var(U); mU <- mean(U)
  if (!is.finite(vU) || vU<=0) stop("U の分散が 0")
  target_var <- target_var_frac * vU
  if (model=="scale"){
    a <- target_mean/pmax(mU, .Machine$double.eps)
    lambda <- pmax(a*U, floor)
    return(list(lambda=lambda, params=list(a=a)))
  }
  if (model=="affine"){
    b <- sqrt(pmax(target_var,0)/vU); a <- target_mean - b*mU
    lambda <- pmax(a + b*U, floor)
    return(list(lambda=lambda, params=list(a=a,b=b)))
  }
  if (is.null(power_s)){
    qs <- stats::quantile(U, probs=0.05, na.rm=TRUE)
    power_s <- max(as.numeric(qs)/5, 1e-3)
  }
  obj <- function(alpha){
    X <- (U+power_s)^alpha; c <- target_mean/mean(X); lam <- c*X
    abs(stats::var(lam) - target_var)
  }
  opt <- stats::optimize(obj, interval=power_alpha_range)
  alpha <- opt$minimum; X <- (U+power_s)^alpha; c <- target_mean/mean(X)
  lambda <- pmax(c*X, floor)
  list(lambda=lambda, params=list(c=c, alpha=alpha, s=power_s))
}

## 0–1 変換（必要なら）
lambda_transform01 <- function(raw_lambda, N, n,
                               transform=c("none","ratio","minmax","logistic","lambda_over_c"),
                               logistic_alpha=0.2, lambda_over_c=10, minmax_eps=1e-12){
  transform <- match.arg(transform); U <- pmax(N-n,0)
  if (transform=="none") return(raw_lambda)
  if (transform=="ratio") return(ifelse(N>0, U/N, 0))
  if (transform=="minmax"){
    rmin <- min(raw_lambda, na.rm=TRUE); rmax <- max(raw_lambda, na.rm=TRUE)
    rng <- max(rmax-rmin, minmax_eps); return((raw_lambda-rmin)/rng)
  }
  if (transform=="logistic"){
    z <- pmax(pmin(logistic_alpha*raw_lambda,60),-60); return(1/(1+exp(-z)))
  }
  if (transform=="lambda_over_c"){
    cst <- max(lambda_over_c, minmax_eps); return(raw_lambda/(raw_lambda+cst))
  }
}
