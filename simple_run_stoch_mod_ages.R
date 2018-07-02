# simple run script #

sir_model <- odin::odin("stoch_mod_extra_ages.R", verbose = FALSE, skip_cache = TRUE)


## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 3/14

## input a value for recovery rate, gamma

infectious_period <- 14 ## duration of infection in days
gamma <- 1/infectious_period


## input a value for R --> S2, sigma
sigma <- 1/180


## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0
mAb_susc <- 0

## input values for the age dependent death rate

mu_6m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 100000

## input the time period that you wish to run the model for
time_period <- 4700 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01
imp_t <- 100000000000000

## input a strength of seasonality of births (1 = max, 0 = no seasonality)

delta = 1

## run model

x <- sir_model(alpha = alpha, beta = beta[i], gamma = gamma, sigma = sigma, Ab_susc = Ab_susc[k], mAb_susc = mAb_susc[k], 
               mu_6m = mu_6m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate, imp_t = imp_t, delta = delta)

out <- as.data.frame(x$run(t))

## run multiple iterations of the model

x_res <- as.data.frame(replicate(100, x$run(0:5800)[,c(165:191,205:206)]))

