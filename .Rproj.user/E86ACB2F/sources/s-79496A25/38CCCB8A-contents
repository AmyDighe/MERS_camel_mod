#states
S = 90
I = 10
R = 0
N = S + I + R

#probabilities of transitions on a given day
b = 0.002 #probability of a birth (per camel)
mu = 0.002 #probability of death
lamda = 0.01 #probability of infection S-->I
sigm = 0.001 #probability of recovery to an immune, non-infectious state I-->R
omega = 0.0002 #probability of immunity waning R-->S

#length of time in days for model to run
dy = 365
time = seq(1:(3*dy))

#outflows
#OR <- rbinom(n = 1, size = R, prob = omega + mu)
#OS <- rbinom(n = 1, size = S, prob = lamda + mu)
#OI <- rbinom(n = 1, size = I, prob = sigm + mu)
OR <- c()
OS <- c()
OI <- c()

DROM <- function(S, I, R, N, b, mu, lamda, sigm, omega, OR, OS, OI, time){
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
matplot(x = out$time, y = out[c(2,3,4,5)], type = "l", lty = 1, lwd = 3)
}

DROM(S = 90, I = 10, R = 0, N = 100, b = 0.002, mu = 0.002, lamda = 0.01, sigm = 0.001, omega = 0.0002, OR = OR, OS = OS, OI = OI, time = seq(1:3*365))

