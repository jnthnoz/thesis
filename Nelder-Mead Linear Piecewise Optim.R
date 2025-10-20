rm(list = ls(all = TRUE))
graphics.off()
pacman::p_load(purrr)
set.seed(123)

# Number of repetitions
n_runs <- 5000

# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix <- matrix(NA, nrow = n_runs, ncol = 8)

lower=c(0.00,0.10,0.20,0.30,0.7,1.7,0.30,0.50)
upper=c(0.20,0.35,0.45,0.50,1.3,2.3,0.70,0.90)
# ==============================================================================
# ==============================================================================
# Piecewise Linear Change-Point Model
# ==============================================================================
# ==============================================================================
lambda  <- .1
lambda0 <- .2
lambda1 <- .3
lambda2 <- .4
tau0    <-  0
tau1    <-  1
tau2    <-  2
beta1   <- .5
beta2   <- .7

negLL <- function(x, fmsfundata, Z) {
  lambda  <- x[1]
  lambda0 <- x[2]
  lambda1 <- x[3]
  lambda2 <- x[4]
  tau1    <- x[5]
  tau2    <- x[6]
  beta1   <- x[7]
  beta2   <- x[8]
  
  out <- - sum(
    ((x[1] + x[2] * fmsfundata$data + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*fmsfundata$data)-1))*
      ifelse (fmsfundata$data < x[5], 1, 0) +
      
      ((x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3])*fmsfundata$data)-exp((x[2]+x[3])*x[5])))*
      ifelse(fmsfundata$data >= x[5] & fmsfundata$data < x[6], 1, 0) +
      
      ((x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + x[4] * (fmsfundata$data - x[6]) + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3])*x[6])-exp((x[2]+x[3])*x[5])) - (1/(x[2]+x[3]+x[4]))*exp(x[1] - x[3]*x[5] - x[4]*x[6] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3]+x[4])*fmsfundata$data)-exp((x[2]+x[3]+x[4])*x[6])))*
      ifelse(fmsfundata$data >= x[6], 1, 0))
  return(out)
}


# Run optimization n_runs-times
for (j in 1:n_runs) {
  Y <- rexp(400, rate = 1)
  T <- c()
  Z <- c(rbinom(1, size = 1, p = 0.7), runif(1, min = 0, max = 1))
  
  A <- (1 / lambda0) * (exp(lambda + lambda0 * tau1 + (Z[1] * beta1 + Z[2] * beta2)) - exp(lambda + (Z[1] * beta1 + Z[2] * beta2)))
  B <- (1 / (lambda0 + lambda1)) * (exp(lambda - lambda1 * tau1 + lambda0 * tau2 + lambda1 * tau2 + (Z[1] * beta1 + Z[2] * beta2)) - exp(lambda + lambda0 * tau1 + (Z[1] * beta1 + Z[2] * beta2)))
  
  for (i in 1:length(Y)) {
    if (Y[i] >= 0 && Y[i] < A) {
      T[i] <- (1 / lambda0) * (log(lambda0 * Y[i] + exp(lambda + (Z[1] * beta1 + Z[2] * beta2))) - lambda - (Z[1] * beta1 + Z[2] * beta2))
    } else if (Y[i] >= A && Y[i] < A + B) {
      T[i] <- (1 / (lambda0 + lambda1)) * log((lambda0 + lambda1) * (Y[i] - A) + exp(lambda + lambda0 * Y[i] + (Z[1] * beta1 + Z[2] * beta2))) + (1 / (lambda0 + lambda1)) * (-lambda + lambda1 * tau1 - (Z[1] * beta1 + Z[2] * beta2))
    } else {
      T[i] <- (1 / (lambda0 + lambda1 + lambda2)) * (log((lambda0 + lambda1 + lambda2) * (Y[i] - A - B) + exp(lambda - lambda1 * tau1 + (lambda0 + lambda1) * tau2 + (Z[1] * beta1 + Z[2] * beta2)))) + (1 / (lambda0 + lambda1 + lambda2)) * (-lambda + lambda1 * tau1 * lambda2 * tau2 - (Z[1] * beta1 + Z[2] * beta2))
    }
  }
  
  # For integration, define each step separately to facilitate R calculations.
  fmsfundata <- structure(list(data = T), class = 'optimbase.functionargs')
  
  # Initial guess for lambda, lambda0, lambda1, lambda2, tau1, tau2, beta1 and beta2.
  x0 <- c(.5,.5,.5,.5,.5,1,.5,.5)
  
  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,
    fn = negLL,
    fmsfundata = fmsfundata,
    Z = Z,
    method = "Nelder-Mead",
    control = list(maxit = 50000)
  )
  
  # Putting in place constraint bound arguments.
  for (i in 1:length(result_matrix[j,]))
  {
    if (optim_result$par[i] <= lower[i]) { result_matrix[j, i] <- lower[i] }
    else if (optim_result$par[i] >= upper[i]) { result_matrix[j, i] <- upper[i] }
    else {result_matrix[j, i] <- optim_result$par[i]}
  }
  
}

colMeans(result_matrix)
apply(result_matrix, 2, sd)
# 
# > colMeans(result_matrix)
# [1] 0.1331248 0.1797567 0.3013423 0.4767679 1.0793983 1.7547334 0.4797054 0.6540960
# > apply(result_matrix, 2, sd)
# [1] 0.09222838 0.10947940 0.11871972 0.06229786 0.22395743 0.14654383 0.18838705 0.18780053





















# ==============================================================================
# Adding 15% censorship
# ==============================================================================
# Number of repetitions
n_runs <- 5000

# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix_15 <- matrix(NA, nrow = n_runs, ncol = 8)

lower=c(0.00,0.10,0.20,0.30,0.7,1.7,0.30,0.50)
upper=c(0.20,0.35,0.45,0.50,1.3,2.3,0.70,0.90)
# ==============================================================================
# ==============================================================================
# Piecewise Linear Change-Point Model
# ==============================================================================
# ==============================================================================
lambda  <- .1
lambda0 <- .2
lambda1 <- .3
lambda2 <- .4
tau0    <-  0
tau1    <-  1
tau2    <-  2
beta1   <- .5
beta2   <- .7

negLL <- function(x, fmsfundata, d, Z) {
  lambda  <- x[1]
  lambda0 <- x[2]
  lambda1 <- x[3]
  lambda2 <- x[4]
  tau1    <- x[5]
  tau2    <- x[6]
  beta1   <- x[7]
  beta2   <- x[8]
  
  out <- - sum(
    (d*(x[1] + x[2] * fmsfundata$data + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*fmsfundata$data)-1))*
      ifelse (fmsfundata$data < x[5], 1, 0) +
      
      (d*(x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3])*fmsfundata$data)-exp((x[2]+x[3])*x[5])))*
      ifelse(fmsfundata$data >= x[5] & fmsfundata$data < x[6], 1, 0) +
      
      (d*(x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + x[4] * (fmsfundata$data - x[6]) + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3])*x[6])-exp((x[2]+x[3])*x[5])) - (1/(x[2]+x[3]+x[4]))*exp(x[1] - x[3]*x[5] - x[4]*x[6] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3]+x[4])*fmsfundata$data)-exp((x[2]+x[3]+x[4])*x[6])))*
      ifelse(fmsfundata$data >= x[6], 1, 0))
  return(out)
}


# Run optimization n_runs-times
for (j in 1:n_runs) {
  Y <- rexp(400, rate = 1)
  W = rexp(400, rate = 0.15)

  T <- c()
  C <-
  Z <- c(rbinom(1, size = 1, p = 0.7), runif(1, min = 0, max = 1))
  
  A <- (1 / lambda0) * (exp(lambda + lambda0 * tau1 + (Z[1] * beta1 + Z[2] * beta2)) - exp(lambda + (Z[1] * beta1 + Z[2] * beta2)))
  B <- (1 / (lambda0 + lambda1)) * (exp(lambda - lambda1 * tau1 + lambda0 * tau2 + lambda1 * tau2 + (Z[1] * beta1 + Z[2] * beta2)) - exp(lambda + lambda0 * tau1 + (Z[1] * beta1 + Z[2] * beta2)))
  
  for (i in 1:length(Y)) {
    if (Y[i] >= 0 && Y[i] < A) {
      T[i] <- (1 / lambda0) * (log(lambda0 * Y[i] + exp(lambda + (Z[1] * beta1 + Z[2] * beta2))) - lambda - (Z[1] * beta1 + Z[2] * beta2))
    } else if (Y[i] >= A && Y[i] < A + B) {
      T[i] <- (1 / (lambda0 + lambda1)) * log((lambda0 + lambda1) * (Y[i] - A) + exp(lambda + lambda0 * Y[i] + (Z[1] * beta1 + Z[2] * beta2))) + (1 / (lambda0 + lambda1)) * (-lambda + lambda1 * tau1 - (Z[1] * beta1 + Z[2] * beta2))
    } else {
      T[i] <- (1 / (lambda0 + lambda1 + lambda2)) * (log((lambda0 + lambda1 + lambda2) * (Y[i] - A - B) + exp(lambda - lambda1 * tau1 + (lambda0 + lambda1) * tau2 + (Z[1] * beta1 + Z[2] * beta2)))) + (1 / (lambda0 + lambda1 + lambda2)) * (-lambda + lambda1 * tau1 * lambda2 * tau2 - (Z[1] * beta1 + Z[2] * beta2))
    }
  }
  
  for (i in 1:length(W)) {
    if (W[i] >= 0 && W[i] < A) {
      C[i] <- (1 / lambda0) * (log(lambda0 * W[i] + exp(lambda + (Z[1] * beta1 + Z[2] * beta2))) - lambda - (Z[1] * beta1 + Z[2] * beta2))
    } else if (W[i] >= A && W[i] < A + B) {
      C[i] <- (1 / (lambda0 + lambda1)) * log((lambda0 + lambda1) * (W[i] - A) + exp(lambda + lambda0 * W[i] + (Z[1] * beta1 + Z[2] * beta2))) + (1 / (lambda0 + lambda1)) * (-lambda + lambda1 * tau1 - (Z[1] * beta1 + Z[2] * beta2))
    } else {
      C[i] <- (1 / (lambda0 + lambda1 + lambda2)) * (log((lambda0 + lambda1 + lambda2) * (W[i] - A - B) + exp(lambda - lambda1 * tau1 + (lambda0 + lambda1) * tau2 + (Z[1] * beta1 + Z[2] * beta2)))) + (1 / (lambda0 + lambda1 + lambda2)) * (-lambda + lambda1 * tau1 * lambda2 * tau2 - (Z[1] * beta1 + Z[2] * beta2))
    }
  }
  
  X = pmin(T,C)
  d = ifelse(pmin(T,C)==T,1,0)
  
  # For integration, define each step separately to facilitate R calculations.
  fmsfundata <- structure(list(data = pmin(T,C)), class = 'optimbase.functionargs')
  
  # Initial guess for lambda, lambda0, lambda1, lambda2, tau1, tau2, beta1 and beta2.
  x0 <- c(.5,.5,.5,.5,.5,1,.5,.5)
  
  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,
    fn = negLL,
    fmsfundata = fmsfundata,
    d = d,
    Z = Z,
    method = "Nelder-Mead",
    control = list(maxit = 50000)
  )
  
  # Putting in place constraint bound arguments.
  for (i in 1:length(result_matrix_15[j,]))
  {
    if (optim_result$par[i] <= lower[i]) { result_matrix_15[j, i] <- lower[i] }
    else if (optim_result$par[i] >= upper[i]) { result_matrix_15[j, i] <- upper[i] }
    else {result_matrix_15[j, i] <- optim_result$par[i]}
  }
  
}

colMeans(result_matrix_15)
# [1] 0.1237555 0.1831214 0.2983130 0.4585877 1.0401764 1.7989002 0.4837043 0.6613592
apply(result_matrix_15, 2, sd)
# [1] 0.09485384 0.11040445 0.11858852 0.07942841 0.24370931 0.20586139 0.19018502 0.18964134



















# ==============================================================================
# Adding 30% censorship
# ==============================================================================
# Number of repetitions
n_runs <- 5000

# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix_30 <- matrix(NA, nrow = n_runs, ncol = 8)

lower=c(0.00,0.10,0.20,0.30,0.7,1.7,0.30,0.50)
upper=c(0.20,0.35,0.45,0.50,1.3,2.3,0.70,0.90)
# ==============================================================================
# ==============================================================================
# Piecewise Linear Change-Point Model
# ==============================================================================
# ==============================================================================
lambda  <- .1
lambda0 <- .2
lambda1 <- .3
lambda2 <- .4
tau0    <-  0
tau1    <-  1
tau2    <-  2
beta1   <- .5
beta2   <- .7

negLL <- function(x, fmsfundata, d, Z) {
  lambda  <- x[1]
  lambda0 <- x[2]
  lambda1 <- x[3]
  lambda2 <- x[4]
  tau1    <- x[5]
  tau2    <- x[6]
  beta1   <- x[7]
  beta2   <- x[8]
  
  out <- - sum(
    (d*(x[1] + x[2] * fmsfundata$data + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*fmsfundata$data)-1))*
      ifelse (fmsfundata$data < x[5], 1, 0) +
      
      (d*(x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3])*fmsfundata$data)-exp((x[2]+x[3])*x[5])))*
      ifelse(fmsfundata$data >= x[5] & fmsfundata$data < x[6], 1, 0) +
      
      (d*(x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + x[4] * (fmsfundata$data - x[6]) + (Z[1] * x[7] + Z[2] * x[8])) - (1/x[2])*exp(x[1] + (Z[1] * x[7] + Z[2] * x[8]))*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3])*x[6])-exp((x[2]+x[3])*x[5])) - (1/(x[2]+x[3]+x[4]))*exp(x[1] - x[3]*x[5] - x[4]*x[6] + (Z[1] * x[7] + Z[2] * x[8]))*(exp((x[2]+x[3]+x[4])*fmsfundata$data)-exp((x[2]+x[3]+x[4])*x[6])))*
      ifelse(fmsfundata$data >= x[6], 1, 0))
  return(out)
}


# Run optimization n_runs-times
for (j in 1:n_runs) {
  Y <- rexp(400, rate = 1)
  W = rexp(400, rate = 0.4)
  
  T <- c()
  C <- c()
  Z <- c(rbinom(1, size = 1, p = 0.7), runif(1, min = 0, max = 1))
  
  A <- (1 / lambda0) * (exp(lambda + lambda0 * tau1 + (Z[1] * beta1 + Z[2] * beta2)) - exp(lambda + (Z[1] * beta1 + Z[2] * beta2)))
  B <- (1 / (lambda0 + lambda1)) * (exp(lambda - lambda1 * tau1 + lambda0 * tau2 + lambda1 * tau2 + (Z[1] * beta1 + Z[2] * beta2)) - exp(lambda + lambda0 * tau1 + (Z[1] * beta1 + Z[2] * beta2)))
  
  for (i in 1:length(Y)) {
    if (Y[i] >= 0 && Y[i] < A) {
      T[i] <- (1 / lambda0) * (log(lambda0 * Y[i] + exp(lambda + (Z[1] * beta1 + Z[2] * beta2))) - lambda - (Z[1] * beta1 + Z[2] * beta2))
    } else if (Y[i] >= A && Y[i] < A + B) {
      T[i] <- (1 / (lambda0 + lambda1)) * log((lambda0 + lambda1) * (Y[i] - A) + exp(lambda + lambda0 * Y[i] + (Z[1] * beta1 + Z[2] * beta2))) + (1 / (lambda0 + lambda1)) * (-lambda + lambda1 * tau1 - (Z[1] * beta1 + Z[2] * beta2))
    } else {
      T[i] <- (1 / (lambda0 + lambda1 + lambda2)) * (log((lambda0 + lambda1 + lambda2) * (Y[i] - A - B) + exp(lambda - lambda1 * tau1 + (lambda0 + lambda1) * tau2 + (Z[1] * beta1 + Z[2] * beta2)))) + (1 / (lambda0 + lambda1 + lambda2)) * (-lambda + lambda1 * tau1 * lambda2 * tau2 - (Z[1] * beta1 + Z[2] * beta2))
    }
  }
  
  for (i in 1:length(W)) {
    if (W[i] >= 0 && W[i] < A) {
      C[i] <- (1 / lambda0) * (log(lambda0 * W[i] + exp(lambda + (Z[1] * beta1 + Z[2] * beta2))) - lambda - (Z[1] * beta1 + Z[2] * beta2))
    } else if (W[i] >= A && W[i] < A + B) {
      C[i] <- (1 / (lambda0 + lambda1)) * log((lambda0 + lambda1) * (W[i] - A) + exp(lambda + lambda0 * W[i] + (Z[1] * beta1 + Z[2] * beta2))) + (1 / (lambda0 + lambda1)) * (-lambda + lambda1 * tau1 - (Z[1] * beta1 + Z[2] * beta2))
    } else {
      C[i] <- (1 / (lambda0 + lambda1 + lambda2)) * (log((lambda0 + lambda1 + lambda2) * (W[i] - A - B) + exp(lambda - lambda1 * tau1 + (lambda0 + lambda1) * tau2 + (Z[1] * beta1 + Z[2] * beta2)))) + (1 / (lambda0 + lambda1 + lambda2)) * (-lambda + lambda1 * tau1 * lambda2 * tau2 - (Z[1] * beta1 + Z[2] * beta2))
    }
  }
  
  X = pmin(T,C)
  d = ifelse(pmin(T,C)==T,1,0)
  
  # For integration, define each step separately to facilitate R calculations.
  fmsfundata <- structure(list(data = pmin(T,C)), class = 'optimbase.functionargs')
  
  # Initial guess for lambda, lambda0, lambda1, lambda2, tau1, tau2, beta1 and beta2.
  x0 <- c(.5,.5,.5,.5,.5,1,.5,.5)
  
  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,
    fn = negLL,
    fmsfundata = fmsfundata,
    d = d,
    Z = Z,
    method = "Nelder-Mead",
    control = list(maxit = 50000)
  )
  
  # Putting in place constraint bound arguments.
  for (i in 1:length(result_matrix_30[j,]))
  {
    if (optim_result$par[i] <= lower[i]) { result_matrix_30[j, i] <- lower[i] }
    else if (optim_result$par[i] >= upper[i]) { result_matrix_30[j, i] <- upper[i] }
    else {result_matrix_30[j, i] <- optim_result$par[i]}
  }
  
}

colMeans(result_matrix_30)
# [1] 0.1172180 0.1807136 0.3063051 0.4437772 1.0452985 1.8332076 0.4908861 0.6646095
apply(result_matrix_30, 2, sd)
# [1] 0.09612666 0.11001374 0.12093692 0.08849528 0.25501462 0.22822982 0.19183822 0.19049030