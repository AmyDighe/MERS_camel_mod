
## This is a simple SIR simulator to play with
## Copyright: Amy Dighe 2018

## beta: individual rate of infection
## alpha: birth rate; defaults to 0
## mu: death rate; defaults to 0

sir_sim <- function(beta, omega, sigma = 1 / 20, alpha = 0, mu = 0, 
                    ini_S = 100, ini_I = 1, ini_R = 0,
                    duration = 365) {

    ## This function checks that rates have acceptable values.
    check_rate <- function(r) {
        if (is.null(r)) stop("rate is NULL")
        if (!is.finite(r)) stop("rate is not a finite number:", r)
        if (r < 0) stop("rate is negative:", r)
    }

    
    ## Make a few checks of inputs
    check_rate(beta)
    check_rate(alpha)
    check_rate(mu)
    check_rate(duration)
    
    
    ## This function changes a rate lambda into a probability p, using:
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

    #stay_S_prob <- 1 - (mu_prob + beta_prob)
    #stay_I_prob <- 1 - (mu_prob + sigma_prob)
    #stay_R_prob <- 1 - (mu_prob + omega_prob)


    ## defining outflows etc.

    #outflow_R <- data.frame(waned = 0, dead = 0)
    #outflow_S <- data.frame(infected = 0, dead = 0)
    #outflow_I <- data.frame(recovered = 0, dead = 0)

    S <- I <- R <- N <- integer(duration)
    S[1] <- ini_S
    I[1] <- ini_I
    R[1] <- ini_R 
    N[1]  <- S[1] + I[1] + R[1]

    time <- seq_len(duration)

    
    for (i in time[-1]) {
        
        ## cat("iteration:", i)
        ## browser()
        
        ## transitions into and out of susceptible state (S) in a time step
        new_births <- rpois(1, S[i - 1] * alpha)
        
        #new_waned <- (rmultinom(1, outflow_R, prob = c(omega_prob, mu))[1])
        new_susceptibles <- new_births #+ new_waned
        
        ## transitions into and out of infectious state (I) in a time step
        rate_infection <- beta * I[i - 1] / N[i - 1]

        outflow_S <- rbinom(1, S[i - 1],
                            prob = r2p(rate_infection + mu))
                
        S[i] <- S[i - 1] + new_susceptibles - outflow_S
        
        new_infectious <- rmultinom(1, size = outflow_S,
                                    prob = c(rate_infection, mu))[1]
        
        outflow_I <- rbinom(1, I[i - 1], prob = r2p(sigma + mu))
        
        I[i] <- I[i - 1] + new_infectious - outflow_I
        
        ## transitions into and out of recovered state (R) in a time step
        new_recovered <- rmultinom(1, outflow_I, prob = c(sigma, mu))[1]
        
        outflow_R <- rbinom(1, R[i-1], prob = r2p(omega + mu))
        
        R[i] <- R[i - 1] + new_recovered - outflow_R  
   
        N[i]  <- S[i] + I[i] + R[i]

    }

    out <- data.frame(time, S, I, R, N)
    
}





## example

out <- sir_sim(ini_S = 90, ini_I = 9, ini_R = 1, alpha = 0.002, mu = 0.002,
               beta = 0.01, omega = 0.0002)



matplot(x = out$time, y = out[c(2,3,4,5)], type = "l", lty = 1, lwd = 3,
        xlab = "Time", ylab = "number of individuals")

legend("topright", lwd = 3, lty = 1,
       col = c("black", "red", "green", "blue"),
       legend = c("S", "I", "R", "N"))
