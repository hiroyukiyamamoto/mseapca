# msea_ora2: Adjusted Over-Representation Analysis (ORA) for Metabolite Set Enrichment Analysis (MSEA)
# This function performs an adjusted ORA that accounts for undetected metabolites in each pathway.
# It builds upon the base ORA function and applies a correction based on estimated significance of undetected metabolites.

msea_ora2 <- function (SIG, DET, ALL, M) {
  
  # Step 1: Label assignment for each metabolite group
  # L1: All metabolites (in ALL)
  # L2: Detected metabolites (in DET)
  # L3: Significant metabolites (in SIG)
  L1 <- setlabel(ALL, M)
  L2 <- setlabel(DET, M)
  L3 <- setlabel(SIG, M)
  
  # Step 2: Perform base ORA on detected metabolites
  # B includes basic ORA results used for adjustment
  B <- msea_ora(SIG, DET, M)
  
  # Initialize vector to store p-values for each pathway
  P <- NULL
  
  # Step 3: Loop through each pathway and calculate adjusted counts
  for (i in 1:length(M)) {
    
    # Total, detected, and significant counts for each pathway
    l1 <- colSums(L1)[i]  # Total metabolites
    l2 <- colSums(L2)[i]  # Detected metabolites
    l3 <- colSums(L3)[i]  # Significant metabolites
    
    # Proportion of significant among detected and estimate for undetected
    r <- l3 / l2
    n <- l1 - l2  # Undetected metabolites count
    
    # Adjusted contingency table counts
    a <- round(B$TAB[[i]][1,1] + n * r)     # Adjusted significant count
    b <- round(B$TAB[[i]][1,2] + n * (1 - r))  # Adjusted non-significant count
    
    # Calculate baseline significance proportion from detected metabolites
    p <- length(SIG) / length(DET)
    
    # Counts for non-significant detected and total substances
    c <- round(length(ALL) * p - a)
    d <- round(length(ALL) * (1 - p) - b)
    
    # Fisher's exact test for the adjusted contingency table
    tab <- t(matrix(c(a, b, c, d), 2))
    resfish <- fisher.test(tab, alternative = "greater")
    print(resfish)  # Display test results for each pathway
    
    # Store p-value
    P[i] <- resfish$p.value
  }
  
  # Step 4: Adjust p-values for multiple testing using Benjamini-Hochberg correction
  Q <- p.adjust(P, method = "BH")
  PQ <- cbind(P, Q)
  
  # Set column and row names for clarity
  colnames(PQ) <- c("p.value", "q.value")
  rownames(PQ) <- names(M)
  
  # Final results as a named list
  RES <- list(PQ)
  names(RES) <- c("Result of MSEA (ORA with adjustment)")
  
  return(RES)
}
