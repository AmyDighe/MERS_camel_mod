## A simple deterministic discrete time model 


## equations for transitions between compartments
## M are individuals protected by maternal antibodies

update(M) <- if ((tt %% 30) == 0) M - M + alpha * N_0 else M + alpha * N_0
update(S) <- if ((tt %% 30) == 0) S - beta * S * I / N - mu * S + M else S - beta * S * I / N - mu * S
update(I) <- I + beta * S * I / N - gamma * I - mu * I
update(R) <- R + gamma * I - mu * R
update(tt) <- tt + 1

## record total population size

N <- M + S + I + R


## initial states

initial(M) <- M_ini # will be user-defined
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
gamma <- user(0.1) # user-defined, default = 0.1

## useful outputs

output(N) <- N