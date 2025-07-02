# msea_ora : wrapper function of over-representation analysis using the following options:
#          "det","all","est_naive", "est_weighted", or "est_shrink"
msea_ora <- function(SIG, DET = NULL, M, option = "det", lambda = NULL) {
  if (option == "det") {
    B <- ora_det(SIG, DET, M)
  } else if (option == "all") {
    B <- ora_all(SIG, M)
  } else if (option == "est_naive") {
    B <- ora_est(SIG, DET, M, method = "naive")
  } else if (option == "est_weighted") {
    B <- ora_est(SIG, DET, M, method = "weighted")
  } else if (option == "est_shrink") {
    B <- ora_est(SIG, DET, M, method = "shrink", lambda = lambda)
  } else {
    stop("Invalid option. Use 'det', 'all', 'est_naive', 'est_weighted', or 'est_shrink'.")
  }
  return(B)
}
