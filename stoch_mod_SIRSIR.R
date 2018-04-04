## A simple stochastic discrete-time model with ageing

N_age <- 14 ## number of age categories


## the proportion of individuals in each age category
## 6m-9m, 9m-12m, 12m-15m, 15m-18m, 18m-21, 21-24m, 2yr+
den_age[1:13]<- 0.05 # COME BACK TO
den_age[14] <- 0.65

## rates of transition

alpha <- user(0) # birth rate, user-defined, default = 0

beta <- user(0.3) # infection rate, user-defined, default =  0.3
Ab_protec <- user(0) # proportion of susceptibility experienced if Abs present. default = 0
beta_2 <- Ab_protec * beta # reinfection rate (infection rate in the presence of Ab protection)
rate_infection <- beta * (sum(I[1:N_age]) / N)
rate_reinfection <- beta_2 * (sum(I[1:N_age]) / N)

mu <- user(0) # death rate, user-defined, default = 0
gamma <- user(0.05) # recovery rate, user-defined, default = 0.1

## converting rates to probabilities

#p_alpha <- 1 - exp(-alpha)
p_infection <- 1 - exp(-rate_infection)
p_reinfection <- 1 - exp(-rate_reinfection)
p_mu <- 1 - exp(-mu)
p_gamma <- 1 - exp(-gamma)


p_S <- p_infection + p_mu # probability of leaving S (with the exception of through ageing)
p_I <- p_gamma + p_mu # probability of leaving I (with the exception of through ageing)
p_R <- p_reinfection + p_mu # probability of leaving R (or S for those protected by maternal Abs)

## birth process: any new individual will enter 'M[1]'
pi <- 3.14159
birth_rate <- N_0 * alpha * (1 + cos(3 * cos(2 * pi * tt / 365)))
new_births <- rpois(birth_rate) 

# outflows (due to infection, recovery or death - ageing is dealt with seperately)
outflow_S[1:6] <- rbinom(S[i], prob = p_R) #same as R because protected by mAbs
outflow_I[1:6] <- rbinom(I[i], prob = p_I)
outflow_R[1:6] <- rbinom(R[i], prob = p_mu)
outflow_S[7:N_age] <- rbinom(S[i], prob = p_S)
outflow_I[7:N_age] <- rbinom(I[i], prob = p_I)
outflow_R[7:N_age] <- rbinom(R[i], prob = p_mu)

# number of individuals leaving S through infection and I through recovery
norm_p_infection <- p_infection/(p_infection + p_mu)
norm_p_reinfection <- p_reinfection/(p_reinfection + p_mu)
norm_p_gamma <- p_gamma/(p_gamma + p_mu)

new_infections[1:6] <- rbinom(outflow_S[i], prob = norm_p_reinfection)
new_infections[7:N_age] <- rbinom(outflow_S[i], prob = norm_p_infection)
new_recoveries[1:N_age] <- rbinom(outflow_I[i], prob = norm_p_gamma)
#new_reinfections[] <- rbinom(outflow_R[i], prob = norm_p_reinfection)

# number of individuals leaving each compartment through ageing

aged_S[1:12] <- if(tt %% 30 == 0) S[i] - outflow_S[i] else 0

aged_S[13] <- if(tt %% 360 == 0) S[13] - outflow_S[13] else 0 #13th comp has width 1 yr

aged_I[1:12] <- if(tt %% 30 == 0) I[i] - outflow_I[i] else 0

aged_I[13] <- if(tt %% 360 == 0) I[13] - outflow_I[13] else 0 #13th comp has width 1 yr

aged_R[1:12] <- if(tt %% 30 == 0) R[i] - outflow_R[i] else 0

aged_R[13] <- if(tt %% 360 == 0) R[13] - outflow_R[13] else 0 #13th comp has width 1 yr


## equations for transitions between compartments
## S[1:6] = individuals protected by maternal antibodies (<6 months of age)
## S[7:N_age] = fully susceptible individuals
## I = infectious individuals
## R = recovered individuals

update(S[1]) <- S[1] - outflow_S[1] - aged_S[1] + new_births
update(S[2:(N_age - 1)]) <- S[i] - outflow_S[i] - aged_S[i] + aged_S[i-1]
update(S[N_age]) <- S[N_age] - outflow_S[N_age] + aged_S[(N_age - 1)]

update(I[1]) <- I[1] - outflow_I[1] - aged_I[1] + new_infections[1] 
update(I[2:(N_age - 1)]) <- I[i] - outflow_I[i] - aged_I[i] + new_infections[i] + aged_I[i - 1]
update(I[N_age]) <-  if(tt == ttt) I[N_age] - outflow_I[N_age] + new_infections[N_age] + aged_I[(N_age - 1)] + imported_cases else I[N_age] - outflow_I[N_age] + new_infections[N_age] + aged_I[(N_age - 1)]

update(R[1]) <- R[1] - outflow_R[1] - aged_R[1] + new_recoveries[1]
update(R[2:(N_age - 1)]) <- R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i] + aged_R[i - 1]
update(R[N_age]) <- R[N_age] - outflow_R[N_age] + new_recoveries[N_age] + aged_R[(N_age - 1)]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size

N <- (sum(S[1:(N_age)]) + sum(I[1:(N_age)]) + sum(R[1:(N_age)]))

## initial states

initial(S[1:N_age]) <- S_ini[i] # will be user-defined
initial(I[1:N_age]) <- I_ini # will be user-defined
initial(R[1:N_age]) <- 0
initial(tt) <- 1

## initial population size for use in birthrate

N_0 <- sum(S_ini[N_age]) + I_ini * N_age # the adult population size (>2 yrs)

## importation of cases

ttt <- user(500) # time of importation, user-defined, default day 500
imported_cases <- user(0) # imported cases, user-defined, default = 0
S_ini[1:N_age] <- 2000 * den_age[i]
I_ini <- 1


## useful outputs

output(MM) <- sum(S[1:6]) + sum(I[1:6]) + sum(R[1:6]) # total number of individuals still protected by maternal immunity
output(SS) <- sum(S[1:N_age]) # total number of susceptible individuals
output(II) <- sum(I[1:N_age]) # total number of infectious individuals
output(RR) <- sum(R[1:N_age]) # total number of recovered individuals (modelled to be immune)
output(N) <- N # total number of individuals

output(SC) <- sum(S[7:13]) # total number of susceptible juveniles
output(IC) <- sum(I[7:13]) # total number of infectious juveniles
output(RC) <- sum(R[7:13]) # total number of recovered juveniles

output(SJ) <- sum(S[7:13]) # total number of susceptible juveniles w/0 mAbs
output(IJ) <- sum(I[7:13]) # total number of infectious juveniles w/0 mAbs
output(RJ) <- sum(R[7:13]) # total number of recovered juveniles w/0 mAbs (modelled to be immune)

output(SA) <- S[N_age] # total number of susceptible adults
output(IA) <- I[N_age] # total number of infectious adults
output(RA) <- R[N_age] # total number of recovered adults (modelled to be immune)

output(reinf) <- rate_reinfection
output(inf) <- rate_infection
output(births) <- new_births

## dim calls needed for arrays
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(S_ini) <- N_age
#dim(I_ini) <- N_age
dim(den_age) <- N_age
dim(new_infections) <- N_age
dim(new_recoveries) <- N_age
dim(aged_S) <- 13
dim(aged_I) <- 13
dim(aged_R) <- 13
dim(outflow_S) <- N_age
dim(outflow_I) <- N_age
dim(outflow_R) <- N_age