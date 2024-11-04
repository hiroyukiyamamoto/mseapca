# msea_ora_binomial_ci: Adjusted Over-Representation Analysis (ORA) for Metabolite Set Enrichment Analysis (MSEA)
# This function performs an adjusted ORA that accounts for undetected metabolites in each pathway.
# It uses a binomial distribution to simulate potential undetected significant metabolites and computes
# the range of p-values based on a 95% confidence interval. The results provide both the baseline p-value
# and the range of possible p-values based on the binomial resampling approach.

msea_ora_binomial_ci <- function(SIG, DET, ALL, M, num_simulations = 1000) {
  
  # Set a random seed for reproducibility
  set.seed(1)
  
  # Perform base ORA on detected metabolites
  B <- msea_ora(SIG, DET, M)
  
  # Initialize vectors to store results
  P <- NULL
  P_range <- NULL
  
  # Perform calculations for each pathway
  for (i in 1:length(M)) {
    
    # Get metabolite counts for each pathway
    l1 <- sum(ALL %in% M[[i]])  # Total metabolites in the pathway
    l2 <- sum(DET %in% M[[i]])   # Detected metabolites in the pathway
    l3 <- sum(SIG %in% M[[i]])   # Significant metabolites in the pathway
    
    # Directly use the proportion of significant metabolites among detected metabolites
    r <- ifelse(l2 > 0, l3 / l2, 0)  # Proportion of significant metabolites
    n <- l1 - l2  # Number of undetected metabolites in the pathway
    
    # Initial 2x2 table based on undetected metabolites
    a <- B$TAB[[i]][1,1]
    b <- B$TAB[[i]][1,2]
    
    # Overall proportion of significant metabolites
    p <- length(SIG) / length(DET)
    
    # Overall count of non-significant metabolites
    c <- round(length(ALL) * p - a)
    d <- round(length(ALL) * (1 - p) - b)
    
    # Resampling based on binomial distribution and p-value range calculation
    simulated_p_values <- numeric(num_simulations)
    
    for (j in 1:num_simulations) {
      # Resample significant metabolites among undetected ones using a binomial distribution
      sampled_significant <- rbinom(1, size = n, prob = r)
      
      # Reconstruct the 2x2 table based on resampling
      a_var <- a + sampled_significant
      b_var <- b + (n - sampled_significant)
      c_var <- round(length(ALL) * p - a_var)
      d_var <- round(length(ALL) * (1 - p) - b_var)
      
      tab_var <- matrix(c(a_var, b_var, c_var, d_var), nrow = 2)
      
      # Calculate p-value using Fisher's exact test
      simulated_p_values[j] <- fisher.test(tab_var, alternative = "greater")$p.value
    }
    
    # Obtain the range of p-values from simulations
    p_min <- quantile(simulated_p_values, probs = 0.025)
    p_max <- quantile(simulated_p_values, probs = 0.975)
    
    # Calculate the default p-value using Fisher's exact test
    tab <- matrix(c(a, b, c, d), nrow = 2)
    resfish <- fisher.test(tab, alternative = "greater")
    P[i] <- mean(simulated_p_values)
    
    # Store the range of p-values
    P_range <- rbind(P_range, c(p_min, P[i], p_max))
  }
  
  # Adjust p-values for multiple testing
  Q <- p.adjust(P, method = "BH")
  PQ <- cbind(P, Q)
  rownames(PQ) <- names(M)
  colnames(PQ) <- c("p.value", "q.value")
  
  # Set row and column names for p-value range output
  rownames(P_range) <- names(M)
  colnames(P_range) <- c("lower p-value", "p-value(mean)","upper p-value")
  
  # Display results
  list("Result of MSEA (ORA with adjustment)" = PQ, 
       "Range of p-values for each pathway (95% confidence interval)" = P_range)
}
