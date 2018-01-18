#states
S = 90    #susceptible
I = 10    #infected
R = 0     #recovered
N = S + I + R    #total camel population

#probabilities of transitions on a given day
b = 0.002 #probability of a birth (per camel)
mu = 0.002 #probability of death
lamda = 0.01 #probability of infection S-->I
sigm = 0.001 #probability of recovery to an immune, non-infectious state I-->R
omega = 0.0002 #probability of immunity waning R-->S

#length of time in days for model to run
dy = 365
time = seq(1:(3*dy))

#defining outflows
OR <- c()
OS <- c()
OI <- c()

  for(i in 1:(length(time)-1)){
    
    OR[i] <- rbinom(n = 1, size = R[i], prob = omega + mu)
    OS[i] <- rbinom(n = 1, size = S[i], prob = lamda + mu)
    OI[i] <- rbinom(n = 1, size = I[i], prob = sigm + mu)
    
S[i+1] = S[i] + rpois(n = 1, lambda = b)*(S[i]+I[i]+R[i]) + rbinom(n = 1, size = OR[i], prob = omega/(omega+mu)) - OS[i]
I[i+1] = I[i] + rbinom(n = 1, size = OS[i], prob= lamda/(lamda+mu)) - OI[i]
R[i+1] = R[i] + rbinom(n = 1, size = OI[i], prob = sigm/(sigm+mu)) - OR[i]
N[i+1] = S[i+1] + I[i+1] + R[i +1]
  }
  out <- data.frame(time, S, I, R, N)
  matplot(x = out$time, y = out[c(2,3,4,5)], type = "l", lty = 1, lwd = 3, xlab = "time", ylab = "number of individuals")
  legend("topright", lwd = 3, lty = 1, col = c("black", "red", "green", "blue"), legend = c("S", "I", "R", "N"))

#####################################################################################################################  
## THOUGHTS ##
  
  # to look at basic behavious of the infection want to start with an equal birth and death rate ie. constant pop?
  # just remove stochasticity from these two parameters and have a look
  # DOESN'T WORK BECAUSE CREATES NUMBERS THAT ARE NOT INTEGERS MEANING THE STOCHASTIC BITS CAN'T WORK!!
  
  #states
  S = 90    #susceptible
  I = 10    #infected
  R = 0     #recovered
  N = S + I + R    #total camel population
  
  #probabilities of transitions on a given day
  b = 0.002 #probability of a birth (per camel)
  mu = 0.002 #probability of death
  lamda = 0.01 #probability of infection S-->I
  sigm = 0.001 #probability of recovery to an immune, non-infectious state I-->R
  omega = 0.0002 #probability of immunity waning R-->S
  
  #length of time in days for model to run
  dy = 365
  time = seq(1:(3*dy))
  
  #defining outflows
  infect <- c()
  rec <- c()
  wan <- c()
  
  for(i in 1:(length(time)-1)){
    
    infect[i] <- rbinom(n = 1, size = (S[i]-(mu*S[i])), prob = lamda)
    rec[i] <- rbinom(n = 1, size = (I[i]-(mu*I[i])), prob = sigm)
    wan[i] <- rbinom(n = 1, size = (R[i]-(mu*R[i])), prob = omega)
    
    S[i+1] = S[i] + b*(S[i]+I[i]+R[i]) -infect[i] + wan[i] - mu*S[i]
    I[i+1] = I[i] + infect[i] - rec[i] - mu*I[i]
    R[i+1] = R[i] + rec[i] - wan[i] - mu*R[i]
    N[i+1] = S[i+1] + I[i+1] + R[i +1]
  }
  out_constant_demog <- data.frame(time, S, I, R, N)
  matplot(x = out_constant_demog$time, y = out_constant_demog[c(2,3,4,5)], type = "l", lty = 1, lwd = 3)
  