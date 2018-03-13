
sir_modell <- odin::odin("outflows.R", verbose = FALSE)

alpha <- 0.001 ## input a value for birth rate (default = 0)
mu <- 0.001 ## input a value for death rate (default = 0)
ttt <- 2000 ## time for importation of cases
imported_cases <- 1 ## input a value for imported cases at time, ttt
time_period <- 2000 ## input the time period that you wish to run the model for
I_ini <- 0

t <- seq(0:time_period)
x <- sir_modell(mu = mu, alpha = alpha, imported_cases = imported_cases, ttt = ttt, I_ini = I_ini) ## include any updated paramters as arguments in the model function

mod_run <- x$run(t) ## run the model



###########
## PLOTS ##
###########

sir_col <- c("aquamarine2", "#8c8cd9", "#cc0044", "black", "pink", "gold", "dimgrey")

matplot(x = mod_run[, 1], y = mod_run[, 30:34], xlab = "Time (days)", ylab = "Number of individuals",
        main = "Discrete MSIR model - deterministic", type = "l", col = sir_col,
        lwd = 3, lty = 1, xlim = c(0, 2000))
legend("topright", lwd = 3, col = sir_col, legend = c("M", "S", "I", "R", "N"))

matplot(x = mod_run[, 1], y = mod_run[, c(3,4,5,6,9)], xlab = "Time (days)", ylab = "Number of individuals",
        main = "Discrete MSIR model - deterministic", type = "l", col = sir_col,
        lwd = 3, lty = 1, ylim = c(0,100), xlim = c(365, 2000))
legend("topright", lwd = 3, col = sir_col, legend = c("S1", "S2", "S3", "S4", "S7"))