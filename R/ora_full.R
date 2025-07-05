# ora_full : ORA that enumerates every possible hit/miss pattern among
#            undetected metabolites and returns the minimum, median,
#            and maximum p-values for each pathway.
ora_full <- function(SIG, DET, M) {

  ALL <- unique(as.character(unlist(M)))

  # Step 1: Label assignment for each metabolite group
  L1 <- setlabel(ALL, M) # All metabolites
  L2 <- setlabel(DET, M) # Detected metabolites
  L3 <- setlabel(SIG, M) # Significant metabolites

  # Step 2: Perform base ORA on detected metabolites
  B <- ora_det(SIG, DET, M)

  # Initialize vectors to store results
  P <- NULL
  P_range <- NULL

  # Step 3: Loop through each pathway
  for (i in 1:length(M)) {
    # Counts for each pathway
    l1 <- colSums(L1)[i]  # Total metabolites
    l2 <- colSums(L2)[i]  # Detected metabolites
    l3 <- colSums(L3)[i]  # Significant metabolites

    # Proportion of significant among detected, and estimate for undetected
    r <- l3 / l2
    n <- l1 - l2  # Count of undetected metabolites

    # Construct adjusted 2x2 table
    a <- round(B$`Contingency tables`[[i]][1,1] + n * r)
    b <- round(B$`Contingency tables`[[i]][1,2] + n * (1 - r))

    # Baseline significance proportion from detected metabolites
    p <- length(SIG) / length(DET)

    # Count of non-significant detected and total substances
    c <- max(0, round(length(ALL) * p - a))
    d <- max(0, round(length(ALL) * (1 - p) - b))

    # Conditional branching based on option
    # Calculate p-value range
    possible_values <- 0:n  # Full range of undetected metabolites
    p_values <- sapply(possible_values, function(x) {
      # Construct 2x2 table for all patterns
      a_var <- l3 + x
      b_var <- l2 - l3 + (n - x)
      c_var <- max(0, round(length(ALL) * p - a_var))
      d_var <- max(0, round(length(ALL) * (1 - p) - b_var))

      tab_var <- matrix(c(a_var, b_var, c_var, d_var), nrow = 2)

    # Fisher's exact test for each pattern if valid table
      if (all(tab_var >= 0) && all(is.finite(tab_var))) {
        fisher.test(tab_var, alternative = "greater")$p.value
      } else {
        NA  # Invalid table entries
      }
    })

    # Obtain the range of p-values (minimum and maximum)
    p_values <- na.omit(p_values)  # Remove NAs from invalid tables
    p_min <- min(p_values)
    p_max <- max(p_values)

    # Store the range result
    P_range <- rbind(P_range, c(p_min, median(p_values), p_max))
  }

  # Output the range of lower and upper p-values
  rownames(P_range) <- names(M)
  colnames(P_range) <- c("lower p-value", "p-value(median)", "upper p-value")

  result <- list(P_range)
  names(result) <- c("Range of p-values")

  return(result)
}
