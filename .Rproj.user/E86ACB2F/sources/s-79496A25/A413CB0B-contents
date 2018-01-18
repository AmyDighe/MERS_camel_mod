# variables have been renamed and subsetted to make it more user friendly
# the random sampling from a poisson distribution for camel births has been corrected

#states
S = 90    #susceptible
I = 10    #infected
R = 0     #recovered
N = S + I + R    #total camel population

#rates of transitions on a given day
b = 0.002 #probability of a birth (per camel)
mu = 0.002 #probability of death
lamda = 0.01 #probability of infection S-->I
sigm = 0.001 #probability of recovery to an immune, non-infectious state I-->R
omega = 0.0002 #probability of immunity waning R-->S


#rates --> probabilities
r2p <- function(r){
  p = 1 - exp(-r)
  print(p)
}

b_prob <- r2p(b)
mu_prob <- r2p(mu)
lamda_prob <- r2p(lamda)
sigm_prob <- r2p(sigm)
omega_prob <- r2p(omega)


#length of time in days for model to run
dy = 365
time = seq(1:(3*dy))

#defining outflows
Outflow_R <- c()
Outflow_S <- c()
Outflow_I <- c()

for(i in 1:(length(time)-1)){
  
  Outflow_R[i] <- rbinom(n = 1, size = R[i], prob = omega_prob + mu_prob)
  Outflow_S[i] <- rbinom(n = 1, size = S[i], prob = lamda_prob + mu_prob)
  Outflow_I[i] <- rbinom(n = 1, size = I[i], prob = sigm_prob + mu_prob)
  
  new_birth <- sum(rpois(n = (S[i]+I[i]+R[i]), lambda = b_prob)) 
  new_waned <- rbinom(n = 1, size = Outflow_R[i], prob = omega_prob/(omega_prob+mu_prob))  
  new_S <- new_birth + new_waned
  dead_S <- rbinom(n = 1, size = Outflow_S[i], prob = mu_prob/(lamda+mu))
  new_infectious <- rbinom(n = 1, size = Outflow_S[i], prob= lamda_prob/(lamda_prob+mu_prob))
  dead_I <- rbinom(n = 1, size = Outflow_I[i], prob = mu_prob/(sigm_prob+mu_prob))
  new_recovered <- rbinom(n = 1, size = Outflow_I[i], prob = sigm_prob/(sigm_prob+mu_prob))
  dead_R <- rbinom(n = 1, size = Outflow_R[i], prob = mu_prob/(omega_prob+mu_prob))
  
  S[i+1] = S[i] + new_birth + new_waned - Outflow_S[i]
  I[i+1] = I[i] + new_infectious - Outflow_I[i]
  R[i+1] = R[i] + new_recovered - Outflow_R[i]
  N[i+1] = S[i+1] + I[i+1] + R[i +1]
  
  cat(new_birth)
}

# by sampling again from the distribution, even if the size is the size of the outflow,
# the sum of the two types of transitions could exceed or not reach the outflow as they are
# not sampled simultaneously
# will multinomial fix this?

out <- data.frame(time, S, I, R, N)
matplot(x = out$time, y = out[c(2,3,4,5)], type = "l", lty = 1, lwd = 3, xlab = "time", ylab = "number of individuals")
legend("topright", lwd = 3, lty = 1, col = c("black", "red", "green", "blue"), legend = c("S", "I", "R", "N"))
