## A simple stochastic discrete-time model with ageing
## In this model 1 year is approximated as 360 days

N_age <- 27 ## number of age categories


## rates of transition


alpha <- user(0.0005) # birth rate, user-defined, default = 0.0005 v rough estimate as use of above was a circular definition with S_ini

beta <- user(0.3) # infection rate, user-defined, default =  0.3
mAb_susc <- user(0) # proportion of susceptibility experienced if mAbs present. default = 0
Ab_susc <- user(0) # proportion of susceptibility experienced if previoulsy infected. default = 0
reduced_shed <- user(1) # proportion of shedding seen in reinfections default = 1 (no difference)
beta_mAb <- mAb_susc * beta # reinfection rate (infection rate in the presence of Ab protection)
beta_Ab <- Ab_susc * beta
beta_I2 <- reduced_shed * beta
beta_I2_mAb <- mAb_susc * reduced_shed * beta
beta_I2_Ab <- Ab_susc * reduced_shed * beta

rate_infection <- beta * (sum(I[1:N_age]) / N)  + beta_I2 * (sum(I2[1:N_age]) / N)
rate_infection_mAb <- beta_mAb * (sum(I[1:N_age]) / N) + beta_I2_mAb * (sum(I2[1:N_age]) / N)
rate_reinfection <- beta_Ab * (sum(I[1:N_age]) / N) + beta_I2_Ab * (sum(I2[1:N_age]) / N)

mu_6m <- user(0.005) # death rate for 1st 6 months of life, user-defined, default = 0.005
mu_2y <- user(0.001) # death rate for the rest of the 1st 2 yrs of life
mu_adult <- user(0.0005) # death rate in adulthood (2 yrs +)
mu[1:6] <- mu_6m
mu[7:24] <- mu_2y
mu[25:27] <- mu_adult

gamma <- user(0.05) # recovery rate, user-defined, default = 0.1
sigma <- user(0.0005) ## waning immunity, 1/mean age of mothers.. shady data... # user-defined, default = 0

## converting rates to probabilities

p_alpha <- 1 - exp(-alpha)
p_infection <- 1 - exp(-rate_infection)
p_infection_mAb <- 1 - exp(-rate_infection_mAb)
p_reinfection <- 1 - exp(-rate_reinfection)
p_mu[1:N_age] <- 1 - exp(-mu[i])
p_gamma <- 1 - exp(-gamma)
p_sigma <- 1 - exp(-sigma)

p_S[1:6] <- 1 - exp(- (rate_infection_mAb + mu[i])) # probability of leaving S for those <6m protected by maternal Abs
p_S[7:N_age] <- 1 - exp(- (rate_infection + mu[i])) # probability of leaving S (with the exception of through ageing)
p_I[1:N_age] <- 1 - exp(- (gamma + mu[i])) # probability of leaving I (with the exception of through ageing)
p_R[1:N_age] <- 1 - exp(- (sigma + mu[i])) # probability of leaving R (with the exception of through ageing)

p_S2[1:N_age] <- 1 - exp(- (rate_reinfection + mu[i]))
p_I2[1:N_age] <- 1 - exp(- (gamma + mu[i]))
p_R2[1:N_age] <- 1 - exp(- (sigma + mu[i]))

## birth process: any new individual will enter 'S[1]'
delta <- user() #modulates the seasonality of births
pi <- 3.14159
birth_rate <- N_0 * alpha * (1 + delta *(cos(3 * cos(pi * tt / 360))))
new_births <- rpois(birth_rate) 

## importation process
importation_rate <- user(0.01)
imported_cases <- rpois(importation_rate)

# outflows (due to infection, recovery or death - ageing is dealt with seperately)
outflow_S[1:N_age] <- rbinom(S[i], prob = p_S[i])
outflow_I[1:N_age] <- rbinom(I[i], prob = p_I[i])
outflow_R[1:N_age] <- rbinom(R[i], prob = p_R[i])
outflow_S2[1:N_age] <- rbinom(S2[i], prob = p_S2[i])
outflow_I2[1:N_age] <- rbinom(I2[i], prob = p_I2[i])
outflow_R2[1:N_age] <- rbinom(R2[i], prob = p_R2[i])

# number of new infections and recoveries

#normalising the probabilities 
norm_p_infection[1:6] <- p_infection_mAb/(p_infection_mAb + p_mu[i])
norm_p_infection[7:N_age] <- p_infection/(p_infection + p_mu[i])
norm_p_reinfection[1:N_age] <- p_reinfection/(p_reinfection + p_mu[i])
norm_p_gamma[1:N_age] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[1:N_age] <- p_sigma/(p_sigma + p_mu[i])

new_infections[1:N_age] <- rbinom(outflow_S[i], prob = norm_p_infection[i])
new_recoveries[1:N_age] <- rbinom(outflow_I[i], prob = norm_p_gamma[i])
new_waned[1:N_age] <- rbinom(outflow_R[i], prob = norm_p_sigma[i])
new_reinfections[1:N_age] <- rbinom(outflow_S2[i], prob = norm_p_reinfection[i])
new_recoveries_2[1:N_age] <- rbinom(outflow_I2[i], prob = norm_p_gamma[i]) 
new_waned_2[1:N_age] <- rbinom(outflow_R2[i], prob = norm_p_sigma[i])

# number of individuals leaving each compartment through ageing alone (but remaining in S or I or R)
# individuals which have entered age compartment 27 remain there until they die

aged_S[1:24] <- if(tt %% 30 == 0) S[i] - outflow_S[i] else 0

aged_S[25:26] <- if(tt %% 360 == 0) S[i] - outflow_S[i] else 0 #13th comp has width 1 yr

aged_I[1:24] <- if(tt %% 30 == 0) I[i] - outflow_I[i] else 0

aged_I[25:26] <- if(tt %% 360 == 0) I[i] - outflow_I[i] else 0 #13th comp has width 1 yr

aged_R[1:24] <- if(tt %% 30 == 0) R[i] - outflow_R[i] else 0

aged_R[25:26] <- if(tt %% 360 == 0) R[i] - outflow_R[i] else 0 #13th comp has width 1 yr

aged_S2[1:24] <- if(tt %% 30 == 0) S2[i] - outflow_S2[i] else 0

aged_S2[25:26] <- if(tt %% 360 == 0) S2[i] - outflow_S2[i] else 0 #13th comp has width 1 yr

aged_I2[1:24] <- if(tt %% 30 == 0) I2[i] - outflow_I2[i] else 0

aged_I2[25:26] <- if(tt %% 360 == 0) I2[i] - outflow_I2[i] else 0 #13th comp has width 1 yr

aged_R2[1:24] <- if(tt %% 30 == 0) R2[i] - outflow_R2[i] else 0

aged_R2[25:26] <- if(tt %% 360 == 0) R2[i] - outflow_R2[i] else 0 #13th comp has width 1 yr

## equations for transitions between compartments
## S[1:6] = individuals protected by maternal antibodies (<6 months of age)
## S[7:N_age] = fully susceptible individuals
## I = infectious individuals
## R = recovered individuals
imp_t <- user()

update(S[1]) <- S[1] - outflow_S[1] - aged_S[1] + new_births
update(S[2:(N_age - 1)]) <- S[i] - outflow_S[i] - aged_S[i] + aged_S[i-1]
update(S[N_age]) <- S[N_age] - outflow_S[N_age] + aged_S[(N_age - 1)]

update(I[1]) <-  if (tt %% 30 ==0) I[1] - outflow_I[1] - aged_I[1] else I[1] - outflow_I[1] - aged_I[1] + new_infections[1]
update(I[2:24]) <- if(tt %% 30 == 0) I[i] - outflow_I[i] - aged_I[i] + new_infections[i - 1] + aged_I[i - 1] else I[i] - outflow_I[i] - aged_I[i] + new_infections[i] + aged_I[i - 1]
update(I[25]) <- if(tt %% 30 == 0) I[25] - outflow_I[25] - aged_I[25] + new_infections[24] + aged_I[24] + imported_cases else if(tt == imp_t) 1 +  I[25] - outflow_I[25] - aged_I[25] + new_infections[25] + aged_I[24] + imported_cases else I[25] - outflow_I[25] - aged_I[25] + new_infections[25] + aged_I[24] + imported_cases
update(I[26]) <- if(tt %% 360 == 0) I[26] - outflow_I[26] - aged_I[26] + new_infections[25] + aged_I[25] else I[26] - outflow_I[26] - aged_I[26] + new_infections[26] + aged_I[25]
update(I[N_age]) <- if(tt %% 360 == 0) I[N_age] - outflow_I[N_age] + sum(new_infections[26:N_age]) + aged_I[(N_age - 1)] else I[N_age] - outflow_I[N_age] + new_infections[N_age] + aged_I[(N_age - 1)]

update(R[1]) <- if(tt %% 30 == 0) R[1] - outflow_R[1] - aged_R[1] else R[1] - outflow_R[1] - aged_R[1] + new_recoveries[1]
update(R[2:25]) <- if(tt %% 30 == 0) R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i - 1] + aged_R[i - 1] else R[i] - outflow_R[i] - aged_R[i] + new_recoveries[i] + aged_R[i - 1]
update(R[26]) <- if(tt %% 360 == 0) R[26] - outflow_R[26] - aged_R[26] + new_recoveries[25] + aged_R[25] else R[26] - outflow_R[26] - aged_R[26] + new_recoveries[26] + aged_R[25]
update(R[N_age]) <- if(tt %% 360 == 0) R[N_age] - outflow_R[N_age] + sum(new_recoveries[26:N_age]) + aged_R[(N_age - 1)] else R[N_age] - outflow_R[N_age] + new_recoveries[N_age] + aged_R[(N_age - 1)]

update(S2[1]) <- if(tt %% 30 == 0) S2[1] - outflow_S2[1] - aged_S2[1] else S2[1] - outflow_S2[1] - aged_S2[1] + new_waned[1] + new_waned_2[1]
update(S2[2:25]) <- if(tt %% 30 == 0) S2[i] - outflow_S2[i] - aged_S2[i] + aged_S2[i-1] + new_waned[i - 1] + new_waned_2[i - 1] else S2[i] - outflow_S2[i] - aged_S2[i] + aged_S2[i - 1] + new_waned[i] + new_waned_2[i]
update(S2[26]) <- if(tt %% 360 == 0)S2[i] - outflow_S2[26] - aged_S2[26] + aged_S2[25] + new_waned[26 - 1] + new_waned_2[25] else S2[26] - outflow_S2[26] - aged_S2[26] + aged_S2[25] + new_waned[26] + new_waned_2[26]
update(S2[N_age]) <- if(tt %% 360 == 0) S2[N_age] - outflow_S2[N_age] + aged_S2[(N_age - 1)] + sum(new_waned[26:N_age]) + sum(new_waned_2[26:N_age]) else S2[N_age] - outflow_S2[N_age] + aged_S2[(N_age - 1)] + new_waned[N_age] + new_waned_2[N_age]

update(I2[1]) <- if(tt %% 30 == 0) I2[1] - outflow_I2[1] - aged_I2[1] else I2[1] - outflow_I2[1] - aged_I2[1] + new_reinfections[1]
update(I2[2:25]) <- if(tt %% 30 == 0) I2[i] - outflow_I2[i] - aged_I2[i] + new_reinfections[i - 1] + aged_I2[i - 1] else I2[i] - outflow_I2[i] - aged_I2[i] + new_reinfections[i] + aged_I2[i - 1]
update(I2[26]) <- if(tt %% 360 == 0)I2[26] - outflow_I2[26] - aged_I2[26] + new_reinfections[25] + aged_I2[25] else I2[26] - outflow_I2[26] - aged_I2[26] + new_reinfections[26] + aged_I2[25]
update(I2[N_age]) <- if(tt %% 360 == 0) I2[N_age] - outflow_I2[N_age] + sum(new_reinfections[26:N_age]) + aged_I2[(N_age - 1)] else I2[N_age] - outflow_I2[N_age] + new_reinfections[N_age] + aged_I2[(N_age - 1)]

update(R2[1]) <- if(tt %% 30 == 0) R2[1] - outflow_R2[1] - aged_R2[1] else R2[1] - outflow_R2[1] - aged_R2[1] + new_recoveries_2[1]
update(R2[2:25]) <- if(tt %% 30 == 0) R2[i] - outflow_R2[i] - aged_R2[i] + new_recoveries_2[i - 1] + aged_R2[i - 1] else R2[i] - outflow_R2[i] - aged_R2[i] + new_recoveries_2[i] + aged_R2[i - 1]
update(R2[26]) <- if(tt %% 360 == 0) R2[26] - outflow_R2[26] - aged_R2[26] + new_recoveries_2[25] + aged_R2[25] else R2[26] - outflow_R2[26] - aged_R2[26] + new_recoveries_2[26] + aged_R2[25]
update(R2[N_age]) <- if(tt %% 360 == 0) R2[N_age] - outflow_R2[N_age] + sum(new_recoveries_2[26:N_age]) + aged_R2[(N_age - 1)] else R2[N_age] - outflow_R2[N_age] + new_recoveries_2[N_age] + aged_R2[(N_age - 1)]


update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work


## record total population size

N <- (sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))

update(seroprevalence[1:N_age]) <- (I[i] + R[i] + S2[i] + I2[i] + R2[i]) / (S[i] + I[i] + R[i] + S2[i] + I2[i] + R2[i])

## initial states

initial(S[1:N_age]) <- S_ini_p[i] # will be user-defined
initial(I[1:N_age]) <- I_ini[i] # will be user-defined
initial(R[1:N_age]) <- 0
initial(S2[1:N_age]) <- 0
initial(I2[1:N_age]) <- 0
initial(R2[1:N_age]) <- 0

initial(tt) <- 1
initial(seroprevalence[1:N_age]) <- 0

## initial population size for use in birthrate

N_0 <- user(1000) # user-defined, default 1000

## setting initial conditions using the equilibrium solution for age distribution
births_detr[1:360] <- 10000 * p_alpha * (1 + cos(3 * cos(pi * i / 360))) # change to fixed N_0 = 100 to avoid NaNs
births_det[1:360] <- rpois(births_detr[i]) 

## if we start the model with the equilibrium amount in each of the first month-wide compartments,
## and no camels in the 2nd year of life (they would have just moved into the adult compartment),
## then from here camels will start filling the yr 2 compartment every month and then eveyr year this will
## empty into the adult compartment. Birthrate will be set to balance summed death rate of this age distribution.

a_max <- 14 ## estimated max age of camels in years (as death rate of adults +2yrs is 1/5.5 years #7.5 years on av)

S_ini[1] <- 0
S_ini[2] <- round((sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * mu[1])), 0)
S_ini[3] <- round((sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:2]))),0)
S_ini[4] <- round((sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:3]))),0)
S_ini[5] <- round((sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:4]))),0)
S_ini[6] <- round((sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:5]))),0)
S_ini[7] <- round((sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:6]))),0)
S_ini[8] <- round((sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:7]))),0)
S_ini[9] <- round((sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:8]))),0)
S_ini[10] <- round((sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:9]))),0)
S_ini[11] <- round((sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:10]))),0)
S_ini[12] <- round((sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:11]))),0)
S_ini[13] <- round(sum(births_det[1:(360 - 330)]) * exp(- (30 * sum(mu[1:12]))),0)
S_ini[14] <- round((sum(births_det[1:360]) - sum(births_det[1:(360 - 30)])) * exp(- (30 * sum(mu[1:13]))), 0)
S_ini[15] <- round((sum(births_det[1:(360 - 30)]) - sum(births_det[1:(360 - 60)])) * exp(- (30 * sum(mu[1:14]))),0)
S_ini[16] <- round((sum(births_det[1:(360 - 60)]) - sum(births_det[1:(360 - 90)])) * exp(- (30 * sum(mu[1:15]))),0)
S_ini[17] <- round((sum(births_det[1:(360 - 90)]) - sum(births_det[1:(360 - 120)])) * exp(- (30 * sum(mu[1:16]))),0)
S_ini[18] <- round((sum(births_det[1:(360 - 120)]) - sum(births_det[1:(360 - 150)])) * exp(- (30 * sum(mu[1:17]))),0)
S_ini[19] <- round((sum(births_det[1:(360 - 150)]) - sum(births_det[1:(360 - 180)])) * exp(- (30 * sum(mu[1:18]))),0)
S_ini[20] <- round((sum(births_det[1:(360 - 180)]) - sum(births_det[1:(360 - 210)])) * exp(- (30 * sum(mu[1:19]))),0)
S_ini[21] <- round((sum(births_det[1:(360 - 210)]) - sum(births_det[1:(360 - 240)])) * exp(- (30 * sum(mu[1:20]))),0)
S_ini[22] <- round((sum(births_det[1:(360 - 240)]) - sum(births_det[1:(360 - 270)])) * exp(- (30 * sum(mu[1:21]))),0)
S_ini[23] <- round((sum(births_det[1:(360 - 270)]) - sum(births_det[1:(360 - 300)])) * exp(- (30 * sum(mu[1:22]))),0)
S_ini[24] <- round((sum(births_det[1:(360 - 300)]) - sum(births_det[1:(360 - 330)])) * exp(- (30 * sum(mu[1:23]))),0)
S_ini[25] <- 0
S_ini[26] <- round(sum(births_det[1:360]) *  exp(-(30*(sum(mu[1:25])) + 360*mu[25])), 0)
S_ini[N_age] <- round(a_max * (sum(births_det[1:360]) * exp(- ((30 * sum(mu[1:24])) + 360 * sum(mu[25:26])))),0)

I_ini[12] <- 0
I_ini[1:11] <- 0
I_ini[13:14] <- 0

S_ini_p[1:N_age] <- round(((S_ini[i] / sum(S_ini[1:N_age])) * N_0), 0)





## useful outputs

#########################################
## number of individuals in each state ##
#########################################

# total number of individuals still protected by maternal immunity
output(M) <- sum(S[1:6])+sum(I[1:6])+sum(R[1:6])+sum(S2[1:6])+sum(I2[1:6])+sum(R2[1:6])

output(S_1) <- sum(S[1:N_age]) # susceptible individuals never infected
output(I_1) <- sum(I[1:N_age]) # individuals infectious for the 1st time
output(R_1) <- sum(R[1:N_age]) # individuals recovered from a 1st infection
output(S_2) <- sum(S2[1:N_age]) # susceptible individuals whose immunity has waned
output(I_2) <- sum(I2[1:N_age]) # individuals infectious for the 2nd+ time
output(R_2) <- sum(R2[1:N_age]) # individuals recovered from 2nd+ infections

output(Stot) <- sum(S[1:N_age]) + sum(S2[1:N_age]) # total number of susceptible individuals
output(Itot) <- sum(I[1:N_age]) + sum(I2[1:N_age]) # total number of infectious individuals
output(Rtot) <- sum(R[1:N_age]) + sum(R2[1:N_age]) # total number of recovered individuals 

output(N) <- N # total number of individuals
output(N_J) <- sum(S[1:13]) + sum(S2[1:13]) + sum(I[1:13]) + sum(I2[1:13]) + sum(R[1:13]) + sum(R2[1:13])
output(N_A) <- sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age])

# output(SA) <- S[N_age] + S2[N_age] # total number of susceptible adults
# output(IA) <- I[N_age] + I2[N_age] # total number of infectious adults
# output(RA) <- R[N_age] + R2[N_age] # total number of recovered adults (modelled to be immune)

####################
## seroprevalence ##
####################

output(seropoz_A) <- 100 * (sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age]))/(sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]) + sum(R2[25:N_age]))
output(seropoz_J) <- 100 * (sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]) + sum(R2[1:24])) / (sum(S[1:24]) + sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]) + sum(R2[1:24]))
output(seropz_tot) <- 100 * (sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))/(sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]) + sum(R2[1:N_age]))


# output(outflowS2) <- outflow_S[2]
# output(outflowS3) <- outflow_S[3]
# output(outflowS4) <- outflow_S[4]
# output(outflowS5) <- outflow_S[5]
# output(outflowS6) <- outflow_S[6]
# output(outflowS7) <- outflow_S[7]
# output(outflowS8) <- outflow_S[8]
# output(outflowS9) <- outflow_S[9]
# output(outflowS10) <- outflow_S[10]
# output(outflowS11) <- outflow_S[11]
# output(outflowS12) <- outflow_S[12]
# output(outflowS13) <- outflow_S[13]

############################
## age at first infection ##
############################
output(reinf_1) <- new_reinfections[1]
output(reinf_2) <- new_reinfections[2]
output(incidence_new_inf) <- sum(new_infections[1:N_age])
output(incidence_indig_inf) <- sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age])
output(total_incidence) <- (sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age]))
###########
## other ##
###########
output(inf) <- rate_infection
output(birthrate) <- birth_rate
output(births) <- new_births
output(importations) <- imported_cases
output(inf_6m_1yr) <- (sum(I[7:13]) + sum(I2[7:13]))
## dim calls needed for arrays
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(S2) <- N_age
dim(I2) <- N_age
dim(R2) <- N_age
dim(S_ini) <- N_age
dim(S_ini_p) <- N_age
dim(I_ini) <- N_age
dim(new_infections) <- N_age
dim(new_recoveries) <- N_age
dim(new_waned) <- N_age
dim(new_reinfections) <- N_age
dim(new_recoveries_2) <- N_age
dim(new_waned_2) <- N_age
dim(aged_S) <- 26
dim(aged_I) <- 26
dim(aged_R) <- 26
dim(outflow_S) <- N_age
dim(outflow_I) <- N_age
dim(outflow_R) <- N_age
dim(aged_S2) <- 26
dim(aged_I2) <- 26
dim(aged_R2) <- 26
dim(outflow_S2) <- N_age
dim(outflow_I2) <- N_age
dim(outflow_R2) <- N_age
dim(births_det) <- 360
dim(births_detr) <- 360
dim(seroprevalence) <- N_age

# because mu varies with age
dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_S) <- N_age
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- N_age
dim(p_I2) <- N_age
dim(p_R2) <- N_age
dim(norm_p_infection) <- N_age
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_reinfection) <- N_age