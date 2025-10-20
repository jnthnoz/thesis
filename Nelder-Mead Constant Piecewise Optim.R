# ==============================================================================
# Piecewise Constant Change-Point Model
# ==============================================================================
rm(list = ls(all = TRUE))
graphics.off()
set.seed(123)

# Needed functions
count_less_than_t <- function(T, t) {sum(T < t, na.rm = TRUE)}
count_less_equal_than_t <- function(T, t) {as.integer(T <= t)}
count_more_than_t <- function(T, t) {as.integer(T > t)}

# ==============================================================================
# ==============================================================================
# Piecewise Constant Change-Point Model
# ==============================================================================
# ==============================================================================
summary_list <- list()
lambda0 <- .2
lambda1 <- .3
lambda2 <- .4
lambda3 <- .5
tau0    <-  0
tau1    <-  1
tau2    <-  2
tau3    <-  3

# ==============================================================================
#   Survival times data through the inverse cumulative distribution function   #
# ==============================================================================
Y = rexp(400, rate = 1)
T = c()

for (i in 1:length(Y)) 
{
  if (Y[i] >= 0 && Y[i] < lambda0*tau1)
  {T[i] = Y[i]/lambda0}
  else if (Y[i] >= lambda0*tau1 && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1))
  {T[i] = tau1 + (Y[i]-lambda0*tau1)/lambda1}
  else if (Y[i] >= lambda0*tau1 + lambda1*(tau2-tau1) && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1) + lambda2*(tau3-tau2))
  {T[i] = tau2 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1))/lambda2}
  else {
    T[i] = tau3 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1)-lambda2*(tau3-tau2))/lambda3
  }
}

fmsfundata <- structure(list(data = T), class = 'optimbase.functionargs')

# Objective function for likelihood
negLL <- function(x, fmsfundata){
  tau1 <- x[1]
  tau2 <- x[2]
  tau3 <- x[3]
  
  out <- -sum((log(lambda0*exp(-lambda0*fmsfundata$data))) *
                ifelse(fmsfundata$data < x[1], 1, 0) +
                
                (log(lambda1*exp(-lambda0*x[1] - 
                                   lambda1*(fmsfundata$data-x[1])))) *
                ifelse(x[1] <= fmsfundata$data & fmsfundata$data < x[2], 1, 0) + 
                
                (log(lambda2*exp(-lambda0*x[1] - 
                                   lambda1*(x[2]-x[1]) - 
                                   lambda2*(fmsfundata$data-x[2])))) *
                ifelse(x[2] <= fmsfundata$data & fmsfundata$data < x[3], 1, 0) + 
                
                (log(lambda3*exp(-lambda0*x[1] - 
                                   lambda1*(x[2]-x[1]) - 
                                   lambda2*(x[3]-x[2]) - 
                                   lambda3*(fmsfundata$data-x[3])))) *
                ifelse(x[3] <= fmsfundata$data, 1, 0))
  
  return(out)
}

# Number of repetitions
n_runs <- 5000

# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix <- matrix(NA, nrow = n_runs, ncol = 3)
hazard_matrix <- matrix(NA, nrow = n_runs, ncol = 4)

# Run optimization 5000 times
for (j in 1:n_runs) {
  Y = rexp(400, rate = 1)
  T = c()
  
  for (i in 1:length(Y)) 
  {
    if (Y[i] >= 0 && Y[i] < lambda0*tau1)
    {T[i] = Y[i]/lambda0}
    else if (Y[i] >= lambda0*tau1 && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1))
    {T[i] = tau1 + (Y[i]-lambda0*tau1)/lambda1}
    else if (Y[i] >= lambda0*tau1 + lambda1*(tau2-tau1) && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1) + lambda2*(tau3-tau2))
    {T[i] = tau2 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1))/lambda2}
    else {
      T[i] = tau3 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1)-lambda2*(tau3-tau2))/lambda3
    }
  }
  
  fmsfundata <- structure(list(data = T), class = 'optimbase.functionargs')
  
  # Initial guess for tau1, tau2, tau3
  x0 <- runif(3, 0, 4)
  
  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,                          # Starting values for tau1, tau2, tau3
    fn = negLL,                        # The negative log-likelihood function
    fmsfundata = fmsfundata,           # Data passed to the function
    method = "Nelder-Mead",            # Optimization method set to Nelder-Mead
    control = list(maxit = 50000)      # Control parameters (max iterations)
  )
  
  # Store the results in the result_matrix (tau1, tau2, tau3)
  result_matrix[j, ] <- optim_result$par
  
  # Estimate the lambda values
  hazard_matrix[j,1] <-  count_less_than_t(T, optim_result$par[1])/sum(pmin(T,optim_result$par[1]))
  hazard_matrix[j,2] <- (count_less_than_t(T, optim_result$par[2]) - count_less_than_t(T, optim_result$par[1])) / (pmin(T,optim_result$par[2])-optim_result$par[1])%*%count_more_than_t(T, optim_result$par[1])
  hazard_matrix[j,3] <- (count_less_than_t(T, optim_result$par[3]) - count_less_than_t(T, optim_result$par[2])) / (pmin(T,optim_result$par[3])-optim_result$par[2])%*%count_more_than_t(T, optim_result$par[2])
  hazard_matrix[j,4] <- (length(T) - count_less_than_t(T, optim_result$par[3])) / ((T-optim_result$par[3])%*%count_more_than_t(T, optim_result$par[3]))
}

# Calculate the mean for each column (tau1, tau2, tau3)
colMeans(hazard_matrix)
colMeans(result_matrix)

# Calculate the standard deviation for each column (tau1, tau2, tau3)
apply(hazard_matrix, 2, sd)
apply(result_matrix, 2, sd)
# ==============================================================================















# ==============================================================================
# Adding 15% censorship
# ==============================================================================
set.seed(123)

# Number of repetitions
n_runs <- 5000

# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix_15 <- matrix(NA, nrow = n_runs, ncol = 3)
hazard_matrix_15 <- matrix(NA, nrow = n_runs, ncol = 4)

num <- c()
den <- c()

# Run optimization 5000 times
for (j in 1:n_runs) {
  Y = rexp(400, rate = 1)
  T = c()
  
  for (i in 1:length(Y)) 
  {
    if (Y[i] >= 0 && Y[i] < lambda0*tau1)
    {T[i] = Y[i]/lambda0}
    else if (Y[i] >= lambda0*tau1 && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1))
    {T[i] = tau1 + (Y[i]-lambda0*tau1)/lambda1}
    else if (Y[i] >= lambda0*tau1 + lambda1*(tau2-tau1) && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1) + lambda2*(tau3-tau2))
    {T[i] = tau2 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1))/lambda2}
    else {
      T[i] = tau3 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1)-lambda2*(tau3-tau2))/lambda3
    }
  }
  
  Z = rexp(400, rate = 1)
  C = c()
  
  for (i in 1:length(Z)) 
  {
    if (Z[i] >= 0 && Z[i] < lambda0*tau1)
    {C[i] = runif(1, tau0, tau1 + 10)}
    else if (Z[i] >= lambda0*tau1 && Z[i] < lambda0*tau1 + lambda1*(tau2-tau1))
    {C[i] = runif(1, tau1, tau2 + 10)}
    else if (Z[i] >= lambda0*tau1 + lambda1*(tau2-tau1) && Z[i] < lambda0*tau1 + lambda1*(tau2-tau1) + lambda2*(tau3-tau2))
    {C[i] = runif(1, tau2, tau3 + 10)}
    else {
      C[i] = runif(1, tau3, max(T))
    }
  }
  
  X = pmin(T,C)

  fmsfundata_15 <- structure(list(data = pmin(T,C)),
                             class = 'optimbase.functionargs')
  
  
  # Initial guess for tau1, tau2, tau3
  x0 <- runif(3, 0, 4)
  
  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,                          # Starting values for tau1, tau2, tau3
    fn = negLL,                        # The negative log-likelihood function
    fmsfundata = fmsfundata_15,        # Data passed to the function
    method = "Nelder-Mead",            # Optimization method set to Nelder-Mead
    control = list(maxit = 50000)      # Control parameters (max iterations)
  )
  
  # Store the results in the result_matrix (tau1, tau2, tau3)
  result_matrix_15[j, ] <- optim_result$par
  
  # Estimate the lambda values
  hazard_matrix_15[j,1] <-  count_less_equal_than_t(T, optim_result$par[1])%*%count_less_equal_than_t(T, C) / sum(pmin(X,optim_result$par[1]))
  hazard_matrix_15[j,2] <- (count_less_equal_than_t(T, optim_result$par[2])%*%count_less_equal_than_t(T, C) - count_less_equal_than_t(T, optim_result$par[1])%*%count_less_equal_than_t(T, C)) / (pmin(X,optim_result$par[2])-optim_result$par[1])%*%count_more_than_t(X, optim_result$par[1])
  hazard_matrix_15[j,3] <- (count_less_equal_than_t(T, optim_result$par[3])%*%count_less_equal_than_t(T, C) - count_less_equal_than_t(T, optim_result$par[2])%*%count_less_equal_than_t(T, C)) / (pmin(X,optim_result$par[3])-optim_result$par[2])%*%count_more_than_t(X, optim_result$par[2])
  
  num[j] <- (length(T) - sum(X==C) - count_less_equal_than_t(T, optim_result$par[3])%*%count_less_equal_than_t(T, C))
  den[j] <- ((X-optim_result$par[3])%*%count_more_than_t(X, optim_result$par[3]))
  hazard_matrix_15[j,4] <- if (den[j] == 0) lambda3 else num[j] / den[j]
}

# Calculate the mean for each column (tau1, tau2, tau3)
colMeans(hazard_matrix_15)
colMeans(result_matrix_15)

# Calculate the standard deviation for each column (tau1, tau2, tau3)
apply(hazard_matrix_15, 2, sd)
apply(result_matrix_15, 2, sd)


















# ==============================================================================
# Adding 30% censorship
# ==============================================================================
# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix_30 <- matrix(NA, nrow = n_runs, ncol = 3)
hazard_matrix_30 <- matrix(NA, nrow = n_runs, ncol = 4)

num <- c()
den <- c()

# Run optimization 5000 times
for (j in 1:n_runs) {
  Y = rexp(400, rate = 1)
  T = c()
  
  for (i in 1:length(Y)) 
  {
    if (Y[i] >= 0 && Y[i] < lambda0*tau1)
    {T[i] = Y[i]/lambda0}
    else if (Y[i] >= lambda0*tau1 && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1))
    {T[i] = tau1 + (Y[i]-lambda0*tau1)/lambda1}
    else if (Y[i] >= lambda0*tau1 + lambda1*(tau2-tau1) && Y[i] < lambda0*tau1 + lambda1*(tau2-tau1) + lambda2*(tau3-tau2))
    {T[i] = tau2 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1))/lambda2}
    else {
      T[i] = tau3 + (Y[i]-lambda0*tau1-lambda1*(tau2-tau1)-lambda2*(tau3-tau2))/lambda3
    }
  }
  
  Z = rexp(400, rate = 1)
  C = c()
  
  for (i in 1:length(Z)) 
  {
    if (Z[i] >= 0 && Z[i] < lambda0*tau1)
    {C[i] = runif(1, tau0, tau1 + 3)}
    else if (Z[i] >= lambda0*tau1 && Z[i] < lambda0*tau1 + lambda1*(tau2-tau1))
    {C[i] = runif(1, tau1, tau2 + 3)}
    else if (Z[i] >= lambda0*tau1 + lambda1*(tau2-tau1) && Z[i] < lambda0*tau1 + lambda1*(tau2-tau1) + lambda2*(tau3-tau2))
    {C[i] = runif(1, tau2, tau3 + 3)}
    else {
      C[i] = runif(1, tau3, max(T))
    }
  }
  
  X = pmin(T,C)

  fmsfundata_30 <- structure(list(data = pmin(T,C)),
                             class = 'optimbase.functionargs')
  
  
  # Initial guess for tau1, tau2, tau3
  x0 <- runif(3, 0, 4)
  
  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,                          # Starting values for tau1, tau2, tau3
    fn = negLL,                        # The negative log-likelihood function
    fmsfundata = fmsfundata_30,        # Data passed to the function
    method = "Nelder-Mead",            # Optimization method set to Nelder-Mead
    control = list(maxit = 50000)      # Control parameters (max iterations)
  )
  
  # Store the results in the result_matrix (tau1, tau2, tau3)
  result_matrix_30[j, ] <- optim_result$par
  
  # Estimate the lambda values
  hazard_matrix_30[j,1] <-  count_less_equal_than_t(T, optim_result$par[1])%*%count_less_equal_than_t(T, C) / sum(pmin(X,optim_result$par[1]))
  hazard_matrix_30[j,2] <- (count_less_equal_than_t(T, optim_result$par[2])%*%count_less_equal_than_t(T, C) - count_less_equal_than_t(T, optim_result$par[1])%*%count_less_equal_than_t(T, C)) / (pmin(X,optim_result$par[2])-optim_result$par[1])%*%count_more_than_t(X, optim_result$par[1])
  hazard_matrix_30[j,3] <- (count_less_equal_than_t(T, optim_result$par[3])%*%count_less_equal_than_t(T, C) - count_less_equal_than_t(T, optim_result$par[2])%*%count_less_equal_than_t(T, C)) / (pmin(X,optim_result$par[3])-optim_result$par[2])%*%count_more_than_t(X, optim_result$par[2])
  
  num[j] <- (length(T) - sum(X==C) - count_less_equal_than_t(T, optim_result$par[3])%*%count_less_equal_than_t(T, C))
  den[j] <- ((X-optim_result$par[3])%*%count_more_than_t(X, optim_result$par[3]))
  hazard_matrix_30[j,4] <- if (den[j] == 0) lambda3 else num[j] / den[j]
}

# Calculate the mean for each column (tau1, tau2, tau3)
colMeans(hazard_matrix_30)
colMeans(result_matrix_30)

# Calculate the standard deviation for each column (tau1, tau2, tau3)
apply(hazard_matrix_30, 2, sd)
apply(result_matrix_30, 2, sd)