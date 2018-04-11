## A simple stochastic discrete-time model with ageing

N_age <- 14 ## number of age categories


## the proportion of individuals in each age category
## first year of life broken into age categories width 1 month, 
## then 2nd year of life width 1 yr, then 2yr+

# den_age[1:13]<- 0.05 # COME BACK TO
# den_age[14] <- 0.65

## rates of transition

#alpha <- user(0) # birth rate, user-defined, default = 0
#alpha_array[1:N_age] <- mu[i] * (S_ini[i]/sum(S_ini[1:N_age]))
#alpha = sum(alpha_array[1:N_age])
alpha <- user(0.0005) #(v rough estimate as use of above was a circular definition with S_ini)
beta <- user(0.3) # infection rate, user-defined, default =  0.3
Ab_protec <- user(0) # proportion of susceptibility experienced if Abs present. default = 0
beta_2 <- Ab_protec * beta # reinfection rate (infection rate in the presence of Ab protection)
rate_infection <- beta * (sum(I[1:N_age]) / N)
rate_infection_mAb <- beta_2 * (sum(I[1:N_age]) / N)
mu[1] <- 0.005 # death rate for 1st month of life
mu[2:(N_age - 1)] <- 0.001 # death rate for the rest of the first 2 yrs of life
mu[N_age] <- 0.0005 # death rate in adulthood (2 yrs +)
gamma <- user(0.05) # recovery rate, user-defined, default = 0.1
sigma <- 4/(8*365) ## waning immunity, 1/mean age of mothers.. shady data... # user-defined, default = 0

## converting rates to probabilities

p_alpha <- 1 - exp(-alpha)
p_infection <- 1 - exp(-rate_infection)
p_infection_mAb <- 1 - exp(-rate_infection_mAb)
p_mu[1:N_age] <- 1 - exp(-mu[i])
p_gamma <- 1 - exp(-gamma)
p_sigma <- 1 - exp(-sigma)

p_S[1:6] <- 1 - exp(- (rate_infection_mAb + mu[i])) # probability of leaving S for those <6m protected by maternal Abs
p_S[7:N_age] <- 1 - exp(- (rate_infection + mu[i])) # probability of leaving S (with the exception of through ageing)
p_I[1:N_age] <- 1 - exp(- (gamma + mu[i])) # probability of leaving I (with the exception of through ageing)
p_R[1:N_age] <- 1 - exp(- (sigma + mu[i])) # probability of leaving R (with the exception of through ageing)

## birth process: any new individual will enter 'M[1]'
pi <- 3.14159
birth_rate <- N_0 * p_alpha * (1 + cos(3 * cos(2 * pi * tt / 360)))
new_births <- rpois(birth_rate) 

# outflows (due to infection, recovery or death - ageing is dealt with seperately)
outflow_S[1:6] <- rbinom(S[i], prob = p_S[i])
outflow_I[1:6] <- rbinom(I[i], prob = p_I[i])
outflow_R[1:6] <- rbinom(R[i], prob = p_R[i])
outflow_S[7:N_age] <- rbinom(S[i], prob = p_S[i])
outflow_I[7:N_age] <- rbinom(I[i], prob = p_I[i])
outflow_R[7:N_age] <- rbinom(R[i], prob = p_R[i])

# number of new infections and recoveries

#normalising the probabilities 
norm_p_infection[1:6] <- p_infection_mAb/(p_infection_mAb + p_mu[i])
norm_p_infection[7:N_age] <- p_infection/(p_infection + p_mu[i])
norm_p_gamma[1:N_age] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[1:N_age] <- p_sigma/(p_sigma + p_mu[i])

new_infections[1:N_age] <- rbinom(outflow_S[i], prob = norm_p_infection[i])
new_recoveries[1:N_age] <- rbinom(outflow_I[i], prob = norm_p_gamma[i])
new_waned[1:N_age] <- rbinom(outflow_R[i], prob = norm_p_sigma[i])

# number of individuals leaving each compartment through ageing
# individuals which have entered age compartment 14 remain there until they die

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

update(S[1]) <- S[1] - outflow_S[1] - aged_S[1] + new_births + new_waned[1]
update(S[2:(N_age - 1)]) <- S[i] - outflow_S[i] - aged_S[i] + aged_S[i-1] + new_waned[i]
update(S[N_age]) <- S[N_age] - outflow_S[N_age] + aged_S[(N_age - 1)] + new_waned[i]

update(I[1]) <- I[1] - outflow_I[1] - aged_I[1] + new_infections[1]
update(I[2:(N_age - 1)]) <- I[i] - outflow_I[i] - aged_I[i] + new_infections[i] + aged_I[i - 1]
update(I[N_age]) <-  if(tt == ttt) I[N_age] - outflow_I[N_age] + new_infections[N_age] + aged_I[(N_age - 1)] + imported_cases else I[N_age] - outflow_I[N_age] + new_infections[N_age] + aged_I[(N_age - 1)]

update(R[1]) <- R[1] - outflow_R[1] - aged_R[1] + new_recoveries[1]
update(R[2:(N_age - 1)]) <- R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i] + aged_R[i - 1]
update(R[N_age]) <- R[N_age] - outflow_R[N_age] + new_recoveries[N_age] + aged_R[(N_age - 1)]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size

N <- (sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]))

## initial states

initial(S[1:N_age]) <- S_ini[i] # will be user-defined
initial(I[1:N_age]) <- I_ini # will be user-defined
initial(R[1:N_age]) <- 0
initial(tt) <- 1

## initial population size for use in birthrate

#N_0 <- sum(S_ini[N_age]) + I_ini * N_age # the adult population size (>2 yrs)
N_0 <- 1000
## importation of cases

ttt <- user(500) # time of importation, user-defined, default day 500
imported_cases <- user(0) # imported cases, user-defined, default = 0

## setting initial conditions using the equilibrium solution for age distribution
births_det[1:360] <- N_0 * p_alpha * (1 + cos(3 * cos(2 * pi * i / 360)))

## if we start the model with the equilibrium amount in each of the first month-wide compartments,
## and no camels in the 2nd year of life (they would have just moved into the adult compartment),
## then from here camels will start filling the yr 2 compartment every month and then eveyr year this will
## empty into the adult compartment. Birthrate will be set to balance summed death rate of this age distribution.

a_max <- 10 ## estimated max age of camels in years (as death rate of adults +2yrs is 1/5.5 years #7.5 years on av)
S_ini[1] <- 0
#S_ini[2:12] <- (sum(births_det[1:(360 - ((i - 2) * 30))]) - 
#               sum(births_det[1:(360 - ((i - 1) * 30))])) * exp(- (30 * sum(mu[1:(i - 1)])))
S_ini[2] <- (sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * mu[1]))
S_ini[3] <- (sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:2])))
S_ini[4] <- (sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:3])))
S_ini[5] <- (sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:4])))
S_ini[6] <- (sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:5])))
S_ini[7] <- (sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:6])))
S_ini[8] <- (sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:7])))
S_ini[9] <- (sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:8])))
S_ini[10] <- (sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:9])))
S_ini[11] <- (sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:10])))
S_ini[12] <- (sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:11])))
S_ini[13] <- 0
S_ini[N_age] <- a_max * (sum(births_det[1:360]) * exp(- ((30 * sum(mu[1:12])) + 360 * mu[13])))

I_ini <- 0


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

output(inf) <- rate_infection
output(births) <- new_births

output(deaths_S) <- sum(outflow_S[1:N_age]) - sum(new_infections[1:N_age])
output(death_per) <- (sum(outflow_S[1:N_age]) - sum(new_infections[1:N_age])) / sum(S[1:(N_age)])
output(death_S1) <- (outflow_S[1] - new_infections[1])

## dim calls needed for arrays
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(S_ini) <- N_age
#dim(I_ini) <- N_age
#dim(den_age) <- N_age
dim(new_infections) <- N_age
dim(new_recoveries) <- N_age
dim(new_waned) <- N_age
dim(aged_S) <- 13
dim(aged_I) <- 13
dim(aged_R) <- 13
dim(outflow_S) <- N_age
dim(outflow_I) <- N_age
dim(outflow_R) <- N_age
#dim(alpha_array) <- N_age
dim(births_det) <- 360

# because mu varies with age
dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_S) <- N_age
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(norm_p_infection) <- N_age
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age