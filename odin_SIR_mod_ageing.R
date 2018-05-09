## A simple deterministic discrete time model with ageing

N_age <- 6 ## number of age categories under 6 months


## equations for transitions between compartments
## M are individuals protected by maternal antibodies

update(M[1]) <- if((tt %% 30) == 0) M[1] - M[1] + alpha * N_0 - mu * M[1] else M[1] + alpha * N_0 - mu * M[1]
update(M[2:N_age]) <- if(tt %% 30 == 0) M[i] - M[i] + M[i-1] - mu * M[i] else M[i] - mu * M[i]
update(S) <- if((tt %% 30) == 0) S - beta * S * I / N - mu * S + M[N_age] else S - beta * S * I / N - mu * S
update(I) <- I + beta * S * I / N - gamma * I - mu * I
update(R) <- R + gamma * I - mu * R
update(tt) <- tt + 1

## record total population size

N <- (sum(M[1:N_age]) + S + I + R)


## initial states

initial(M[1:N_age]) <- M_ini # will be user-defined
initial(S) <- S_ini # will be user-defined
initial(I) <- I_ini # will be user-defined
initial(R) <- 0
initial(tt) <- 0

## initial population size for use in birthrate

N_0 <- S_ini + I_ini

M_ini <- user(0) # user-defined, default = 0
S_ini <- user(1000) # user-defined, default = 1000
I_ini <- user(1) # user-defined, default = 1
alpha <- user(0) # user-defined, default = 0
beta <- user(0.3) # user-defined, default =  0.3
mu <- user(0) # user-defined, default = 0
gamma <- user(0.05) # user-defined, default = 0.1

## useful outputs

output(N) <- N
output(mm) <- sum(M[1:N_age])

## dim calls needed for arrays

dim(M) <- N_age
