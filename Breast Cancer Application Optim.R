rm(list = ls(all = TRUE))
graphics.off()
# install.packages("pacman")
pacman::p_load(neldermead, data.table, ggplot2, tidyverse, readxl, survival, muhaz, numDeriv, dplyr, strucchange)

# Leemos los datos.
data <- read_xlsx(file.choose())

# Nos quedamos solo con los que son malignos de primera vez y los que se conoce su fecha de seguimiento.
data1 <- data %>% filter(First_malignant == 'Yes',
                         Survival_months != 'Unknown')

# Convertir Survival_months a numérico y filtrar valores no numéricos
data1 <- data1 %>% 
            mutate(Survival_months = as.numeric(Survival_months)) %>%
            filter(!is.na(Survival_months))

# Marcamos como censura todos los cánceres que no son de pecho. Tiempo de supervivencia en años.
data2 <- data1 %>% 
            mutate(Status = case_when(Cause_of_Death == 'Dead (attributable to this cancer dx)' & Site == 'Breast' ~ 1, TRUE ~ 0), Survival_years = Survival_months/12)

# Cantidad de censurados:
sum(data2$Status==0)
# 463078

# Porporción de no censurados:
sum(data2$Status==1)/length(data2$Status)
# 0.1940793

# Continuous plot follow-up time against hazard rate until follow-up time 30:
# mod <- muhaz(data2$Survival_years, data2$Status, min.time = min(data2$Survival_years), max.time = 30)
plot(mod, xlim=c(0,30))

# Hazard function estimation with Nelson-Aalen estimator:
na <- basehaz(coxph(Surv(Survival_years, Status) ~ 1, data = data2), centered = FALSE)
na_hazard <- diff(na$hazard) / diff(na$time)
time_points <- na$time[-1]  # Remover el primer punto de tiempo
plot(time_points, na_hazard, type = "s", xlim = c(0,30), xlab = "Years", ylim = c(0,0.04), ylab = "Hazard", main = "Hazard Function (Nelson-Aalen Estimator)")

# Hazard function estimation with Kaplan-Meier estimator:
km <- survfit(Surv(Survival_years, Status) ~ 1, data = data2)
km_hazard <- -diff(log(km$surv)) / diff(km$time)
time_points <- km$time[-1]  # Remover el primer punto de tiempo
plot(time_points, km_hazard, type = "s", xlim = c(0,30), xlab = "Years", ylim = c(0,0.04), ylab = "Hazard", main = "Hazard Function (Kaplan-Meier Estimator)")

#===============================================================================
# Continuous hazard function against Nelson-Aalen (red) and Kaplan-Meier (blue) estimators:
#===============================================================================
plot(mod,
     main = "Hazard Functions Comparison",
     xlab = "Years",
     ylab = "Hazard",
     xlim = c(0,30),
     ylim = c(0,0.04),
     col = "black")
lines(time_points, na_hazard, type = "s", col = "red")
lines(time_points, km_hazard, type = "s", col = "blue")
legend("topright", 
       legend = c("Limit-Product", "Nelson-Aalen", "Kaplan-Meier"), 
       col = c("black", "red", "blue"), 
       lty = 1,
       cex = 0.8)


#===============================================================================
# How many change points should there be? And where are they?
#===============================================================================
library(strucchange)

haz <- mod$haz.est
t   <- mod$est.grid

# Allow up to 6 change points:
bp <- breakpoints(haz ~ t, breaks = 5)

# Summary for change points 0 to 5:
summary(bp)

# Graphing how many breaking points minimize the BIC:
plot(bp)
title(ylab = "BIC")

# Selected manually according to the graph:
opt_breaks <- 2

# Which positions give the optimal model?
cat("Indices of change points:", bp$breakpoints[1:opt_breaks], "\n")
cat("Times of change points:", t[bp$breakpoints[1:opt_breaks]], "\n")


#===============================================================================
# Checking the hazard function estimation for the data
#===============================================================================
# Number of repetitions
n_runs <- 100

# Create an empty matrix to store results (500 rows, 3 columns for tau1, tau2, tau3)
result_matrix <- matrix(NA, nrow = n_runs, ncol = 6)

lower <- c(-3.3,-0.3, 0.0,-0.1, 3.9, 11.4)
upper <- c(-3.1,-0.1, 0.1, 0.1, 4.5, 12.0)
# -3.22 -0.125 0.021 0.065 4.2 11.7
# ==============================================================================
# ==============================================================================
# Piecewise Linear Change-Point Model
# ==============================================================================
# ==============================================================================
lambda  <- -3
lambda0 <- -0.1
lambda1 <- +0.1
lambda2 <- +0.1
tau0    <- 0
tau1    <- 4.2
tau2    <- 11.7

negLL <- function(x, fmsfundata, d) {
  lambda  <- x[1]
  lambda0 <- x[2]
  lambda1 <- x[3]
  lambda2 <- x[4]
  tau1    <- x[5]
  tau2    <- x[6]

  out <- - sum(
    (d*(x[1] + x[2] * fmsfundata$data) - (1/x[2])*exp(x[1])*(exp(x[2]*fmsfundata$data)-1))*
      ifelse (fmsfundata$data < x[5], 1, 0) +
      
      (d*(x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5])) - (1/x[2])*exp(x[1])*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5])*(exp((x[2]+x[3])*fmsfundata$data)-exp((x[2]+x[3])*x[5])))*
      ifelse(fmsfundata$data >= x[5] & fmsfundata$data < x[6], 1, 0) +
      
      (d*(x[1] + x[2] * fmsfundata$data + x[3] * (fmsfundata$data - x[5]) + x[4] * (fmsfundata$data - x[6])) - (1/x[2])*exp(x[1])*(exp(x[2]*x[5])-1) - (1/(x[2]+x[3]))*exp(x[1] - x[3]*x[5])*(exp((x[2]+x[3])*x[6])-exp((x[2]+x[3])*x[5])) - (1/(x[2]+x[3]+x[4]))*exp(x[1] - x[3]*x[5] - x[4]*x[6])*(exp((x[2]+x[3]+x[4])*fmsfundata$data)-exp((x[2]+x[3]+x[4])*x[6])))*
      ifelse(fmsfundata$data >= x[6], 1, 0))
  return(out)
}

# Run optimization n_runs-times
for (j in 1:n_runs) {
  set.seed(j)
  data3 <- sample_n(data2, 5000, replace = FALSE)
  
  T <- data3$Survival_years
  # cat("Mean survival this run:", mean(T), "\n")

  d <- data3$Status

  # For integration, define each step separately to facilitate R calculations.
  fmsfundata <- structure(list(data = T), class = 'optimbase.functionargs')

  # Initial guess for lambda, lambda0, lambda1, lambda2, tau1, tau2, beta1 and beta2.
  x0 <- runif(8, lower, upper)

  # Use optim to minimize the negative log-likelihood
  optim_result <- optim(
    par = x0,
    fn = negLL,
    fmsfundata = fmsfundata,
    d = d,
    method = "Nelder-Mead",
    control = list(maxit = 50000)
  )
  
  # cat("Run", j, "mean survival =", mean(T), "log-lik =", optim_result$value, "\n")
  
  # Putting in place constraint bound arguments.
  for (i in 1:length(result_matrix[j,]))
  {
    if (optim_result$par[i] <= lower[i]) { result_matrix[j, i] <- lower[i] }
    else if (optim_result$par[i] >= upper[i]) { result_matrix[j, i] <- upper[i] }
    else {result_matrix[j, i] <- optim_result$par[i]}
  }

}
colMeans(result_matrix)
# apply(result_matrix, 2, sd)
# -3.22 -0.125 0.021 0.065 4.2 11.7


# Your coefficients
lambda  <- colMeans(result_matrix)[1]
lambda0 <- colMeans(result_matrix)[2]
lambda1 <- colMeans(result_matrix)[3]
lambda2 <- colMeans(result_matrix)[4]
tau1    <- colMeans(result_matrix)[5]
tau2    <- colMeans(result_matrix)[6]


# Your existing plot
plot(mod, xlim = c(0, 30), ylim = c(0,0.03), main="Estimated piecewise linear multiple change-point hazard function")
specific_points <- c(colMeans(result_matrix)[5],colMeans(result_matrix)[6])
abline(v = specific_points, col = "red", lty = 2, lwd = 1)
text(x = c(colMeans(result_matrix)[5],colMeans(result_matrix)[6]), y = max(haz) * 0.10, labels = c("4.2", "11.7"), pos = 4, col = "red")
grid()

hazard_fn <- function(t) {
  ifelse(t < tau1,
         exp(lambda + lambda0 * t),
         ifelse(t < tau2,
                exp(lambda + lambda0 * t + lambda1 * pmax(t - tau1, 0) ),
                exp(lambda + lambda0 * t + lambda1 * pmax(t - tau1, 0) + lambda2 * pmax(t - tau2, 0) )
         )
  )
}

# Add the piecewise hazard line
curve(hazard_fn(x), from = 0, to = 30, col = "blue", lwd = 2, add = TRUE)