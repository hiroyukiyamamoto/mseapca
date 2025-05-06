# msea_ora_all_combination: Adjusted Over-Representation Analysis (ORA) for Metabolite Set Enrichment Analysis (MSEA)
# This function performs an adjusted ORA that accounts for undetected metabolites in each pathway.
# It builds upon the base ORA function and applies a correction based on estimated significance of undetected metabolites.

msea_ora_all_point <- function(SIG, DET, ALL, M) {
  
  # Step 1: Label assignment for each metabolite group
  L1 <- setlabel(ALL, M) # All metabolites
  L2 <- setlabel(DET, M) # Detected metabolites
  L3 <- setlabel(SIG, M) # Significant metabolites
  
  # Step 2: Perform base ORA on detected metabolites
  B <- msea_ora(SIG, DET, M)
  
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
    a <- round(B$TAB[[i]][1,1] + n * r)
    b <- round(B$TAB[[i]][1,2] + n * (1 - r))
    
    # Baseline significance proportion from detected metabolites
    p <- length(SIG) / length(DET)
    
    # Count of non-significant detected and total substances
    c <- max(0, round(length(ALL) * p - a))
    d <- max(0, round(length(ALL) * (1 - p) - b))
    
    # Perform Fisher's test in the default case
    tab <- matrix(c(a, b, c, d), nrow = 2)
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- resfish$p.value

  }
  
  # Adjust p-values for multiple testing (default)
  Q <- p.adjust(P, method = "BH")
  PQ <- cbind(P, Q)
  rownames(PQ) <- names(M)
  colnames(PQ) <- c("p.value", "q.value")
  result <- list(PQ)
  names(result) <- c("Result of MSEA (ORA with adjustment)")

  return(result)
}