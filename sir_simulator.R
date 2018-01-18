
## This is a simple SIR simulator to play with
## Copyright: Amy Dighe 2018

## beta: individual rate of infection
## alpha: birth rate; defaults to 0
## mu: death rate; defaults to 0

sir_sim <- function(beta, omega, sigma = 1 / 20, alpha = 0, mu = 0, 
                    ini_S = 100, ini_I = 1, ini_R = 0
                    duration = 365) {

    ## This function checks that rates have acceptable values.
    check_rate <- function(r) {
        if (!is.finite(r)) stop("rate is not a finite number")
        if (r < 0) stop("rate is negative")
    }

    
    ## Make a few checks of inputs
    check_rate(beta)
    check_rate(alpha)
    check_rate(mu)
    check_rate(duration)
    
    
    ## This function changes a rate lambda into a probabilitie p, using:
    ## p = 1 - exp(-lambda)
    
    r2p <- function(r) {
        check_rate(r)
        p <- 1 - exp(-r)
        return(p)
    }

    alpha_prob <- r2p(alpha)
    mu_prob <- r2p(mu)
    beta_prob <- r2p(beta)
    sigma_prob <- r2p(sigma)
    omega_prob <- r2p(omega)

    stay_S_prob <- 1 - (mu_prob + beta_prob)
    stay_I_prob <- 1 - (mu_prob + sigma_prob)
    stay_R_prob <- 1 - (mu_prob + omega_prob)


    ## defining outflows etc.

    Outflow_R <- data.frame(waned = 0, dead = 0, remain = 0)
    Outflow_S <- data.frame(infected = 0, dead = 0, remain = 0)
    Outflow_I <- data.frame(recovered = 0, dead = 0,remain = 0)

    S <- I <- R <- N <- integer(duration)
    S[1] <- ini_S
    I[1] <- ini_I
    R[1] <- ini_R
    
    N[1]  <- S[1] + I[1] + R[1]

    time <- seq_len(duration)

    for (i in time[-1]) {

        ## handle new susceptibles (S)
        new_susceptibles <- rpois(1, S[i - 1] * alpha) # birth
        S[i] <- S[i - 1] + new_susceptibles
                    
    
        ## handle new infected (I)
        rate_infection <- beta * I[i - 1] / N[i - 1]
        
        outflow_S <- rbinom(S[i - 1], prob = r2p(rate_infection + mu))
        Outflow_R[i,] <-rmultinom(n = 1, size = R[i], prob = c(rate_infection, mu_prob, stay_S_prob))[,1]
        
        ## handle changes in recovered (R)
        
        
        
        ## 
        ## Outflow_S[i,] <- (rmultinom(n = 1, size = S[i], prob = c(beta_prob, mu_prob, stay_I_prob)))[,1]
        ## Outflow_I[i,] <- (rmultinom(n = 1, size = I[i], prob = c(sigma_prob, mu_prob, stay_R_prob)))[,1]
        
        ## new_birth <- sum(rpois(n = (S[i]+I[i]+R[i]), lambda = alpha_prob)) 
        ## new_waned <- Outflow_R$waned[i]  
        ## new_S <- new_birth + new_waned
        ## dead_S <- Outflow_S$dead[i]
        ## new_infectious <- Outflow_S$infected[i]
        ## dead_I <- Outflow_I$dead[i]
        ## new_recovered <- Outflow_I$recovered[i]
        ## dead_R <- Outflow_R$dead[i]
        
        ## S[i+1] = S[i] + new_birth + new_waned - new_infectious - dead_S
        ## I[i+1] = I[i] + new_infectious - new_recovered - dead_I
        ## R[i+1] = R[i] + new_recovered - new_waned - dead_R
        ## N[i+1] = S[i+1] + I[i+1] + R[i +1]
        
    }

    out <- data.frame(duration, S, I, R, N)
    matplot(x = out$duration, y = out[c(2,3,4,5)], type = "l", lty = 1, lwd = 3, xlab = "duration", ylab = "number of individuals")
    legend("topright", lwd = 3, lty = 1, col = c("black", "red", "green", "blue"), legend = c("S", "I", "R", "N"))


}




## example

sir_sim(S = 90, I = 9, R = 1, alpha = 0.002, mu = 0.002, beta = 0.01, sigma = 0.001, omega = 0.0002, time = 3*365)


