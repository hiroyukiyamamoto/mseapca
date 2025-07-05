# msea_ora_range : wrapper function of over-representation analysis using one of the following options:
#          "ora_full","bino_naive","bino_weighted", "bino_shrink"
msea_ora_range <- function(SIG, DET = NULL, M, option = "ora_full", probs = c(0.025, 0.975), nsim = 1000, lambda = 5) {
  if (option == "ora_full") {
    B <- ora_full(SIG, DET, M)
  } else if (option == "bino_naive") {
    B <- ora_bino(SIG, DET, M, method="naive", probs = probs, nsim = nsim)
  } else if (option == "bino_weighted") {
    B <- ora_bino(SIG, DET, M, method="weighted", probs = probs, nsim = nsim)
  } else if (option == "bino_shrink") {
    B <- ora_bino(SIG, DET, M, method="shrink", probs = probs, nsim = nsim, lambda = lambda)
  } else {
    stop("Invalid option. Use 'ora_full', 'bino_naive', 'bino_weighted', or 'bino_shrink'.")
  }
  return(B)
}
