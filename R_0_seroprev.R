sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)



## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- seq(from = 0.05, to = 0.5, by = 0.05)
## input a value for recovery rate, gamma
gamma <- c(1/5, 1/10, 1/15, 1/20, 1/25, 1/30, 1/35, 1/40)

paras <- expand.grid(beta, gamma)
beta_op <- paras$Var1
gamma_op <- paras$Var2
R_0_op <- beta_op/gamma_op
filter <- R_0_op > 1.5
beta_vals <- beta_op[filter]
gamma_vals <- gamma_op[filter]
R_0 <- beta_vals/gamma_vals

## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0

## input values for the age dependent death rate

mu_1m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 100000

## input the time period that you wish to run the model for
time_period <- 3600
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

adult_seropositivity <- vector(length = length(R_0))
juvenile_seropositivity <- vector(length = length(R_0))
total_seropositivity <- vector(length = length(R_0))
max_adult_seropositivity <- vector(length = length(R_0))
max_juvenile_seropositivity <- vector(length = length(R_0))
max_total_seropositivity <- vector(length = length(R_0))
min_adult_seropositivity <- vector(length = length(R_0))
min_juvenile_seropositivity <- vector(length = length(R_0))
min_total_seropositivity <- vector(length = length(R_0))

for(i in 1:(length(R_0))){
message(i)  
  ## include any user-defined parameters as arguments here
  x <- sir_model(alpha = alpha, beta = beta_vals[i], gamma = gamma_vals[i], Ab_susc = Ab_susc, 
                 mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
                 importation_rate = importation_rate)
  
  ## run the model
  x_res <- as.data.frame(replicate(100, x$run(1:3600)[,c(1,117:119,124,127)]))
 
  adult_seropositivity[i] <- mean(rowSums(x_res[time_window, seq(from = 2, to = 2+(99*6), by = 6)]))/100
  juvenile_seropositivity[i] <- mean(rowSums(x_res[time_window, seq(from = 3, to = 3+(99*6), by = 6)]))/100
  total_seropositivity[i] <- mean(rowSums(x_res[time_window, seq(from = 4, to = 4+(99*6), by = 6)]))/100
  
  max_adult_seropositivity[i] <- max(rowSums(x_res[time_window, seq(from = 2, to = 2+(99*6), by = 6)]))/100
  max_juvenile_seropositivity[i] <- max(rowSums(x_res[time_window, seq(from = 3, to = 3+(99*6), by = 6)]))/100
  max_total_seropositivity[i] <- max(rowSums(x_res[time_window, seq(from = 4, to = 4+(99*6), by = 6)]))/100
 
  min_adult_seropositivity[i] <- min(rowSums(x_res[time_window, seq(from = 2, to = 2+(99*6), by = 6)]))/100
  min_juvenile_seropositivity[i] <- min(rowSums(x_res[time_window, seq(from = 3, to = 3+(99*6), by = 6)]))/100
  min_total_seropositivity[i] <- min(rowSums(x_res[time_window, seq(from = 4, to = 4+(99*6), by = 6)]))/100  
}
  
plot(x = R_0, y = total_seropositivity, pch = 19, lwd = 2, ylim = c(20, 90), 
     ylab = "mean seropositivity", xlab = "R0")
