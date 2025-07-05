# msea_ora_beta_ci: Adjusted Over-Representation Analysis (ORA) using Beta distribution
# This function performs an adjusted ORA that accounts for undetected metabolites in each pathway.
# It uses a Beta distribution to estimate the proportion of significant metabolites among undetected ones
# and computes the range of p-values based on the 95% confidence interval for the proportion.

msea_ora_beta_ci <- function(SIG, DET, ALL, M, alpha_prior = 1, beta_prior = 1, num_simulations = 1000) {
  
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
    
    # Set Beta distribution parameters
    alpha <- alpha_prior + l3  # Posterior alpha
    beta <- beta_prior + l2 - l3  # Posterior beta
    
    # Number of undetected metabolites
    n <- l1 - l2
    
    # Simulate the proportion of significant metabolites using the Beta distribution
    simulated_proportions <- rbeta(num_simulations, alpha, beta)
    
    # Resample and compute p-values based on the simulated proportions
    simulated_p_values <- numeric(num_simulations)
    for (j in 1:num_simulations) {
      sampled_significant <- round(n * simulated_proportions[j])

      # Reconstruct 2x2 table based on resampling
      a_var <- B$TAB[[i]][1, 1] + sampled_significant
      b_var <- B$TAB[[i]][1, 2] + (n - sampled_significant)
      
      # Calculate overall non-significant counts
      p <- length(SIG) / length(DET)
      c_var <- round(length(ALL) * p - a_var)
      d_var <- round(length(ALL) * (1 - p) - b_var)
      
      tab_var <- matrix(c(a_var, b_var, c_var, d_var), nrow = 2)
      
      # Calculate p-value using Fisher's exact test
      simulated_p_values[j] <- fisher.test(tab_var, alternative = "greater")$p.value
    }
    
    # Obtain the 95% confidence interval for the p-values
    p_min <- quantile(simulated_p_values, probs = 0.025)
    p_median <- median(simulated_p_values)
    p_max <- quantile(simulated_p_values, probs = 0.975)
    
    # Calculate the default p-value using Fisher's exact test
    tab <- matrix(c(B$TAB[[i]][1, 1], B$TAB[[i]][1, 2], round(length(ALL) * p - B$TAB[[i]][1, 1]), round(length(ALL) * (1 - p) - B$TAB[[i]][1, 2])), nrow = 2)
    resfish <- fisher.test(tab, alternative = "greater")
    #P[i] <- mean(simulated_p_values)

    # Store the range of p-values
    P_range <- rbind(P_range, c(p_min, p_median, p_max))
  }
  
  # Set row and column names for p-value range output
  rownames(P_range) <- names(M)
  colnames(P_range) <- c("lower p-value", "p-value(median)","upper p-value")
  
  # Display results
  list("Range of p-values for each pathway (95% CI using Beta distribution)" = P_range)
  
}