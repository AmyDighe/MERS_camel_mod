sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 0.3

## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0

## input values for the age dependent death rate

mu_1m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 10000000

## input the time period that you wish to run the model for
time_period <- 36000 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0

start_time <- Sys.time()
## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate)

x_res <- as.data.frame(replicate(100, x$run(0:36000)[,]))
end_time <- Sys.time()

time_elapsed <- end_time - start_time

print(time_elapsed)


#for just 11 key variable
###########################################

N_0 <- 10000000

## input the time period that you wish to run the model for
time_period <- 36000 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

#start_time <- Sys.time()
## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = 2, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate)

x_res <- as.data.frame(replicate(10, x$run(0:720)[,]))
end_time <- Sys.time()

time_elapsed <- end_time - start_time

print(time_elapsed)

###############################################
sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 0.3

## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0

## input values for the age dependent death rate

mu_1m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 10000000

## input the time period that you wish to run the model for
time_period <- 3600 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

start_time <- Sys.time()
## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate)

x_res <- as.data.frame(replicate(100, x$run(0:3600)[,94:97]))
sir_col <- c("aquamarine2", "#cc0044", "#8c8cd9", "black", "pink", "gold", "dimgrey")
library(scales)
matplot(0:3600, x_res, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Discrete SIR model - stochastic - 100 realisations", type = "l", col = rep(alpha(sir_col[1:4], 0.1), 100),
        lwd = 3, lty = 1, xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 3600, by = 360), labels = seq(from = 0, to = 10, by = 1))
legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))
text("T", x = 6000000, y = 4000, cex = 1.2)

