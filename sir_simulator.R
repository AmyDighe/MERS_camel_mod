
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

    # initialising 
    # N_0 is used to set the birth rate
    
    S <- I <- R <- N <- integer(duration)
    S[1] <- ini_S
    I[1] <- ini_I
    R[1] <- ini_R 
    N[1]<- S[1] + I[1] + R[1]
    N_0  <- ini_S + ini_I + ini_R
    
    ## For clarity we structure the code in terms of what leaves compartments,
    ## whenever this applies. When individuals can leave a compartment through
    ## different competing processes, the strategy is:

    ## i) draw number of individuals leaving using a binomial distribution

    ## ii) decide where they go using a multinomial distribution

    ## iii) update the respective compartments accordingly


    i <- 1
    while ((N[i] > 0) && (i < duration)) {
      

  
        ## increment time step
        i <- i + 1
        
        ## birth process: any new individual will enter 'S'
        
        new_births <- rpois(1, N_0 * alpha) # changed from N[i -1] to N_0
        
        
        ## individuals leaving 'S', either by becoming 'I' (new infections) or
        ## dying

        if (N[i - 1] == 0) {
            message("N is zero at iteration ", i-1)
        }
        rate_infection <- beta * I[i - 1] / N[i - 1]
        outflow_S <- rbinom(1, S[i - 1],
                            prob = r2p(rate_infection + mu))            
        S[i] <- S[i - 1] + new_births - outflow_S
        new_infectious <- rmultinom(1, size = outflow_S,
                                    prob = c(rate_infection, mu))[1]

        ## individuals leaving 'I', either by becoming 'R' (new recovered) or
        ## dying
        
        outflow_I <- rbinom(1, I[i - 1], prob = r2p(sigma + mu)) 
        I[i] <- I[i - 1] + new_infectious - outflow_I
        new_recovered <- rmultinom(1, outflow_I, prob = c(sigma, mu))[1]


        ## individuals leaving 'R', either by becoming 'S' (waning immunity) or
        ## dying
        
        outflow_R <- rbinom(1, R[i-1], prob = r2p(omega + mu))
        new_waned <- rmultinom(1, outflow_R, prob = c(omega, mu))[1]
        R[i] <- R[i - 1] + new_recovered - outflow_R  
        S[i] <- S[i] + new_waned

        ## update the current population size
        N[i]  <- S[i] + I[i] + R[i]
        
    } # end of the while loop
    


    ## put results into shape and return
    time <- seq_len(i)

    out <- data.frame(time,
                      S[time],
                      I[time],
                      R[time],
                      N[time],
                      new_B[time])

    class(out) <- c("sir_sim", "data.frame")
    return(out)
    
}



## plot method for sir_sim objects

plot.sir_sim <- function(x, ...) {
    sirn_colors <- c("#7575a3", "#b3003b", "#669999", "#8a8a5c")
    
    matplot(x = x$time, y = x[, -1], type = "l", lty = 1, lwd = 3,
            xlab = "Time", ylab = "number of individuals", col = sirn_colors)

    legend("topright", lwd = 3, lty = 1,
           col = sirn_colors, bg = "white",
           legend = c("S", "I", "R", "N"))

}






## example

out <- sir_sim(ini_S = 900, ini_I = 10, ini_R = 0, alpha = 0.002, mu = 0.002,
               beta = 0.1, omega = 0.002)

plot(out)



