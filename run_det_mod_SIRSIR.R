sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 0.3

## input a value for recovery rate, gamma

infectious_period <- 14 ## duration of infection in days
gamma <- 1/infectious_period

## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0

## input values for the age dependent death rate

mu_1m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 10000

## input the time period that you wish to run the model for
time_period <- 10000 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate)

## run the model
mod_run <- x$run(t)
out <- as.data.frame(mod_run)
