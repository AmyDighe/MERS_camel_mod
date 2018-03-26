## A simple stochastic discrete-time model with ageing

N_age <- 6 ## number of age categories between 0 and 6 months
N_age2 <- 7 ## number of age categories above 6 months

## the proportion of individuals in each age category
## 6m-9m, 9m-12m, 12m-15m, 15m-18m, 18m-21, 21-24m, 2yr+
den_age2[1:N_age]<- 0.1
den_age2[N_age2] <- 0.4

## rates of transition

alpha <- user(0) # birth rate, user-defined, default = 0
beta <- user(0.3) # infection rate, user-defined, default =  0.3
rate_infection <- beta * sum(I[1:N_age2]) / N
mu <- user(0) # death rate, user-defined, default = 0
gamma <- user(0.05) # recovery rate, user-defined, default = 0.1

## converting rates to probabilities

p_alpha <- 1 - exp(-alpha)
p_infection <- 1 - exp(-rate_infection)
p_mu <- 1 - exp(-mu)
p_gamma <- 1 - exp(-gamma)

p_S <- p_infection + p_mu # probability of leaving S (with the exception of through ageing)
p_I <- p_gamma + p_mu # probability of leaving I (with the exception of through ageing)

## birth process: any new individual will enter 'M[1]'

birth_rate <- N_0 * alpha # took '*' function out of rpois 
new_births <- rpois(birth_rate) 

# outflows (due to infection, recovery or death - ageing is dealt with seperately)
# removed the '1' for the number of random draws argument incase it is not supported...?

outflow_M[1:N_age] <- rbinom(M[i], prob = p_mu)
outflow_S[1:(N_age2)] <- rbinom(S[i], prob = p_S)
outflow_I[1:(N_age2)] <- rbinom(I[i], prob = p_I)
outflow_R[1:(N_age2)] <- rbinom(R[i], prob = p_mu)

# number of individuals leaving S through infection and I through recovery
# new_infections[1:(N_age2)] <- rmultinom(1, outflow_S[i], prob = c(rate_infection, mu))[1]
# new_recoveries[1:(N_age2)] <- rmultinom(1, outflow_I[i], prob = c(gamma, mu))[1]

# # reassingning arrays for use in stochastic functions compatible with odin
# p__S[1] <- p_infection
# p__S[2] <- p_mu
# dim(p__S) <- 2
# 
# p__I[1] <- p_gamma
# p__I[2] <- p_mu
# dim(p__I) <- 2
# 
# rmn_S[1:N_age2, ] <- rmultinom(outflow_S[i], p__S)
# new_infections[1:N_age2] <- rmn_S[i,1]
# 
# rmn_I[1:N_age2, ] <- rmultinom(outflow_I[i], p__I)
# new_recoveries[1:N_age2] <- rmn_I[i,1]

# messing without arrays to check if it works
foo <- 10
moo[1] <- mu
moo[2] <- gamma
dim(moo) <- 2

rmn_S[] <- rmultinom(foo, moo)
new_infections[1] <- rmn_S[1]

rmn_I[] <- rmultinom(foo, moo)
new_recoveries[1] <- rmn_I[1]

dim(rmn_S) <- 2
dim(rmn_I) <- 2

# number of individuals leaving each compartment through ageing

aged_M[1:N_age] <- if(tt %% 30 == 0) M[i] - outflow_M[i] else 0

aged_S[1:(N_age2 - 1)] <- if(tt %% 91 == 0) S[i] - outflow_S[i] else 0

aged_I[1:(N_age2 - 1)] <- if(tt %% 91 == 0) I[i] - outflow_I[i] else 0

aged_R[1:(N_age2 - 1)] <- if(tt %% 91 == 0) R[i] - outflow_R[i] else 0


## equations for transitions between compartments
## M = individuals protected by maternal antibodies (<6 months of age)
## S = susceptible individuals
## I = infectious individuals
## R = recovered individuals

update(M[1]) <- M[1] - outflow_M[1] - aged_M[1] + new_births
update(M[2:N_age]) <- M[i] - outflow_M[i] - aged_M[i] + aged_M[i-1]

update(S[1]) <- S[1] - outflow_S[1] - aged_S[1] + aged_M[N_age]
update(S[2:(N_age2 - 1)]) <-S[i] - outflow_S[i] - aged_S[i] + aged_S[i - 1]
update(S[N_age2]) <- S[N_age2] - outflow_S[N_age2] + aged_S[N_age2 - 1]

update(I[1]) <- if(tt == ttt) I[1] - outflow_I[1] - aged_I[1] + new_infections[1] + imported_cases else I[1] - outflow_I[1] - aged_I[1] + new_infections[1] 
update(I[2:(N_age2 - 1)]) <- I[i] - outflow_I[i] - aged_I[i] + new_infections[i] + aged_I[i - 1]
update(I[N_age2]) <-  I[N_age2] - outflow_I[N_age2] + new_infections[N_age2] + aged_I[N_age2 - 1]

update(R[1]) <- R[1] - outflow_R[1] - aged_R[1] + new_recoveries[1]
update(R[2:(N_age2 - 1)]) <- R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i] + aged_R[i - 1]
update(R[N_age2]) <- R[N_age2] - outflow_R[N_age2] + new_recoveries[N_age2] + aged_R[(N_age2 - 1)]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size

N <- (sum(M[1:N_age]) + sum(S[1:(N_age2)]) + sum(I[1:(N_age2)]) + sum(R[1:(N_age2)]))

## initial states

initial(M[1:N_age]) <- M_ini # will be user-defined
initial(S[1:N_age2]) <- S_ini[i] # will be user-defined
initial(I[1:N_age2]) <- I_ini # will be user-defined
initial(R[1:N_age2]) <- 0
initial(tt) <- 1

## initial population size for use in birthrate

N_0 <- sum(S_ini[1:7]) + I_ini*(N_age2)

## importation of cases
ttt <- user(1500) # time of importation, user-defined, default day 1500
imported_cases <- user(0) # imported cases, user-defined, default = 0

M_ini <- user(0) # user-defined, default = 0
S_ini[1:7] <- 1000 * den_age2[i]
I_ini <- user(1) # user-defined, default = 1

## useful outputs

output(MM) <- sum(M[1:N_age]) # total number of individuals still protected by maternal immunity
output(SS) <- sum(S[1:N_age2]) # total number of susceptible individuals
output(II) <- sum(I[1:N_age2]) # total number of infectious individuals
output(RR) <- sum(R[1:N_age2]) # total number of recovered individuals (modelled to be immune)
output(N) <- N # total number of individuals

output(SJ) <- sum(S[1:(N_age2 - 1)]) # total number of susceptible juveniles
output(IJ) <- sum(I[1:(N_age2 - 1)]) # total number of infectious juveniles
output(RJ) <- sum(R[1:(N_age2 - 1)]) # total number of recovered juveniles (modelled to be immune)

output(SA) <- S[N_age2] # total number of susceptible adults
output(IA) <- I[N_age2] # total number of infectious adults
output(RA) <- R[N_age2] # total number of recovered adults (modelled to be immune)

## dim calls needed for arrays

dim(M) <- N_age
dim(S) <- N_age2
dim(I) <- N_age2
dim(R) <- N_age2
dim(S_ini) <- N_age2
dim(den_age2) <- N_age2
dim(aged_M) <- N_age
dim(new_infections) <- N_age2
dim(new_recoveries) <- N_age2
dim(aged_S) <- N_age2
dim(aged_I) <- N_age2
dim(aged_R) <- N_age2
dim(outflow_M) <- N_age
dim(outflow_S) <- N_age2
dim(outflow_I) <- N_age2
dim(outflow_R) <- N_age2
#dim(rmn_S)<- c(N_age2,2)
#dim(rmn_I)<- c(N_age2,2)