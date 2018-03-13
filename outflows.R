## A simple deterministic discrete-time model with ageing

N_age <- 6 ## number of age categories between 0 and 6 months
N_age2 <- 6 ## number of age categories between 6 months and 2 years

## the proportion of individuals in each age category
## 6m-9m, 9m-12m, 12m-15m, 15m-18m, 18m-21, 21-24m, 2yr+
den_age2[1:6]<- 0.1
den_age2[7] <- 0.4


# number of individuals leaving S through infection and I through recovery
new_infections[1:(N_age2 + 1)] <- beta * S[i] * sum(I[1:(N_age2 + 1)]) / N
new_recoveries[1:(N_age2 + 1)] <- gamma * I[i]

# number of individuals leaving each compartment through ageing
aged_M[1:N_age] <- if(tt %% 30 == 0) M[i] - mu * M[i] else 0

aged_S[1:N_age2] <- if(tt %% 91 == 0) S[i] - new_infections[i] - mu * S[i] else 0

aged_I[1:(N_age2)] <- if(tt %% 91 == 0) I[i] - new_recoveries[i] - mu * I[i] else 0

aged_R[1:(N_age2)] <- if(tt %% 91 == 0) R[i] - mu * R[i] else 0

# outflows 
outflow_M[1:N_age] <- aged_M[i] + mu * M[i]
outflow_S[1:N_age2] <- aged_S[i] + new_infections[i] + mu * S[i]
outflow_S[N_age2 + 1] <- new_infections[N_age2 + 1] + mu * S[N_age2 + 1]
outflow_I[1:N_age2] <- aged_I[i] + new_recoveries[i] + mu * I[i]
outflow_I[N_age2 + 1] <- new_recoveries[N_age2 + 1] + mu * I[N_age2 + 1]
outflow_R[1:N_age2] <- aged_R[i] + mu * R[i]
outflow_R[N_age2 + 1] <- mu * R[N_age2 + 1]

## equations for transitions between compartments
## M = individuals protected by maternal antibodies (<6 months of age)
## S = susceptible individuals
## I = infectious individuals
## R = recovered individuals

update(M[1]) <- M[1] - outflow_M[1] + alpha * N_0
update(M[2:N_age]) <- M[i] - outflow_M[i] + aged_M[i-1]

update(S[1]) <- S[1] - outflow_S[1] + aged_M[N_age]
update(S[2:N_age2]) <-S[i] - outflow_S[i] + aged_S[i - 1]
update(S[N_age2 + 1]) <- S[N_age2 + 1] - outflow_S[N_age2 + 1] + aged_S[N_age2]

update(I[1]) <- if(tt == ttt) I[1] - outflow_I[1] + new_infections[1] + imported_cases else I[1] - outflow_I[1] + new_infections[1] 
update(I[2:N_age2]) <- I[i] - outflow_I[i] + new_infections[i] + aged_I[i - 1]
update(I[N_age2 + 1]) <-  I[N_age2 + 1] - outflow_I[N_age2 + 1] + new_infections[N_age2 + 1]

update(R[1]) <- R[1] - outflow_R[1] + new_recoveries[1]
update(R[2:N_age2]) <- R[i] - outflow_R[i] + new_recoveries[i] + aged_R[i - 1]
update(R[N_age2 + 1]) <- R[N_age2 + 1] - outflow_R[N_age2 + 1] + new_recoveries[N_age2 + 1]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size

N <- (sum(M[1:N_age]) + sum(S[1:(N_age2+1)]) + sum(I[1:(N_age2+1)]) + sum(R[1:(N_age2+1)]))

## initial states

initial(M[1:N_age]) <- M_ini # will be user-defined
initial(S[1:(N_age2 + 1)]) <- S_ini[i] # will be user-defined
initial(I[1:(N_age2 + 1)]) <- I_ini # will be user-defined
initial(R[1:(N_age2 + 1)]) <- 0
initial(tt) <- 1

## initial population size for use in birthrate

N_0 <- sum(S_ini[1:7]) + I_ini*(N_age2 + 1)

M_ini <- user(0) # user-defined, default = 0
S_ini[1:7] <- 1000 * den_age2[i]
I_ini <- user(1) # user-defined, default = 1

## rates of transition

alpha <- user(0) # birth rate, user-defined, default = 0
beta <- user(0.3) # infection rate, user-defined, default =  0.3
mu <- user(0) # death rate, user-defined, default = 0
gamma <- user(0.05) # recovery rate, user-defined, default = 0.1

ttt <- user(1500) # time of importation of cases, user-defined, default day 1500
imported_cases <- user(0) # imported cases, user-defined, default = 0

## useful outputs

output(MM) <- sum(M[1:N_age]) # total number of individuals still protected by maternal immunity
output(SS) <- sum(S[1:(N_age2 + 1)]) # total number of susceptible individuals
output(II) <- sum(I[1:(N_age2 + 1)]) # total number of infectious individuals
output(RR) <- sum(R[1:(N_age2 + 1)]) # total number of recovered individuals (modelled to be immune)
output(N) <- N # total number of individuals

output(SJ) <- sum(S[1:N_age2]) # total number of susceptible juveniles
output(IJ) <- sum(I[1:N_age2]) # total number of infectious juveniles
output(RJ) <- sum(R[1:N_age2]) # total number of recovered juveniles (modelled to be immune)

output(SA) <- S[N_age2 + 1] # total number of susceptible adults
output(IA) <- I[N_age2 + 1] # total number of infectious adults
output(RA) <- R[N_age2 + 1] # total number of recovered adults (modelled to be immune)

## dim calls needed for arrays

dim(M) <- N_age
dim(S) <- 7
dim(I) <- 7
dim(R) <- 7
dim(S_ini) <- 7
dim(den_age2) <- 7
dim(aged_M) <- N_age
dim(new_infections) <- 7
dim(new_recoveries) <- 7
dim(aged_S) <- N_age2
dim(aged_I) <- N_age2
dim(aged_R) <- N_age2
dim(outflow_M) <- N_age
dim(outflow_S) <- 7
dim(outflow_I) <- 7
dim(outflow_R) <- 7
