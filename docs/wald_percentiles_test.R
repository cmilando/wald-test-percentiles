Wald_percentiles_test <- function(x1, x2, p, n.B, 
                                  show.plot = T, 
                                  signif.digits = 2) {
  #' Arguments:
  #'    @x1   = numeric vector, data from group 1
  #'    @x2   = numeric vector, data from group 2
  #'    @p    = numeric vector, percentiles to investigate
  #'    @n.B  = single integer, number of bootstrap replicates
  
  # -- Step 0 --
  # library requirements and global variables and checks
  require(ggplot2)
  require(boot)
  require(MASS)
  require(Matrix)
  require(ggrepel)
  
  n.p <- length(p)
  len.x1 <- length(x1)
  len.x2 <- length(x2)
  
  # error checking
  stopifnot(is.double(x1))
  stopifnot(is.double(x2))
  
  stopifnot(all(p >= 0))
  stopifnot(all(p <= 1))
  stopifnot(n.p <= 5)
  
  stopifnot(is.numeric(n.B))
  stopifnot(n.B %% 1 == 0)
  stopifnot(length(n.B) == 1)
  if(n.B < 500) warning('n.B should be > 500')
  
  # -- Step 1 & 2 --
  # percentiles profile, q, for each population
  get_quantile <- function(x, p = p) quantile(x, probs = p, names = F)
  q1 <- get_quantile(x1, p)
  q2 <- get_quantile(x2, p)
  q <- c(q1, q2)
  
  # plot
  if(show.plot) {
    
    p_df1 <- data.frame(p, x = q1, fct = 'x1', 
                        z = paste0(p*100, "th: ", signif(q1, signif.digits)))
    
    p_df2 <- data.frame(p, x = q2, fct = 'x2', 
                        z = paste0(p*100, "th: ", signif(q2, signif.digits)))
    
    p_df  <- rbind(p_df1, p_df2)
    
    data_df <- rbind(
      data.frame(x = x1, fct = 'x1'),
      data.frame(x = x2, fct = 'x2')
    )
    
    pp <- ggplot() + xlab("X") + ylab("Density") +
      geom_density(aes(x = x, fill = fct), data = data_df, alpha = 0.5) +
      scale_fill_manual(values = c('x1' = "pink", 'x2' = 'lightblue')) +
      geom_point(data = p_df, aes(x = x, y = 0, color = fct), show.legend = F) + 
      geom_label_repel(data = p_df, aes(x = x, y = 0, label = z, color = fct), 
                       point.padding = 1, direction = 'y', nudge_y = 0.1, 
                       show.legend = F, size = 3) + 
      scale_color_manual(values = c('x1' = 'red', 'x2' = 'blue'))
    
    print(pp)
  }
  
  # -- Step 3 --
  # function to get B bootstrap sample estimates of Q
  get_bootstrap_Q <- function(x) {
    
    # update from previous method and use tsboot function
    x.ts <- as.ts(x) # convert to time-series object
    block.len <- round(length(x)^(1/3))
    
    set.seed(1)
    Q.out <- tsboot(tseries = x.ts, 
                    statistic = get_quantile, 
                    sim = "geom", 
                    R = n.B, l = block.len, 
                    n.sim = length(x.ts), orig.t = F,
                    p = p) 
    
    return(Q.out$t)
    
  }
  
  # tilda - qi, from step 3
  q.B.x1 <- get_bootstrap_Q(x = x1)
  q.B.x2 <- get_bootstrap_Q(x = x2)
  
  # -- Step 4 -- 
  # Compute the p x p covariance matrix Vi for each population
  ones.B <- matrix(rep(1, n.B), ncol = 1)
  internal.mat <- diag(n.B) - n.B ^ -1 * (ones.B %*% t(ones.B))
  
  V.x1 <- (n.B - 1) ^ -1 * ( t(q.B.x1) %*% internal.mat %*% q.B.x1 )
  V.x2 <- (n.B - 1) ^ -1 * ( t(q.B.x2) %*% internal.mat %*% q.B.x2 )
  
  # then combine into the Kp x Kp block matrix V
  V <- as.matrix(bdiag(V.x1, V.x2))
  
  # -- Step 5 --
  # Contrast matrix A = [Ip | -Ip]
  # A is dimension c x Kp, where c is the degrees of freedom
  # and 1 <= c <= (Kp - 1)
  A <- cbind(diag(n.p), -1 * diag(n.p))
  deg.f <- dim(A)[1]
  
  # Compute the Wald statistic, and compare to chi-squared
  # W > qchisq means differences between percentile profiles are
  # statistically significant
  W <- t(q) %*% t(A) %*% ginv( A %*% V %*% t(A)) %*% A %*% q
  chi <- qchisq(0.95, df = deg.f)
  note <- ifelse(W > chi, 
                 "* Differences btwn percentile profiles are statistically significant", 
                 "Percentile profiles are statistically similar")
  cat("W > qchisq? ", W > chi, note, "\n")
  
  # -- Step 6 --
  # Compute confidence intervals
  if(W > chi) {
    z <- 1.96 # critical value of the standard normal distribution for CI 95%
    for (i in 1:deg.f){
      est <- A[i, ] %*% q
      mod <- z * (A[i, ] %*% V %*% t(A)[, i]) ^ 0.5
      lower <- est - mod
      higher <- est + mod
      note <- ifelse(sign(lower) == sign(higher), "* Significant difference", "")
      cat("CI95 of diff at p=", signif(p[i], signif.digits), 
          ": [", signif(est - mod, signif.digits), ", ",
          signif(est + mod, signif.digits), "] ", note, "\n")
    }
  }
  
}