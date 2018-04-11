sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

alpha <- 0.00076 ## input a value for birth rate (default = 0.0005)
#mu
ttt <- 1170 ## time for importation of cases
imported_cases <- 1 ## input a value for imported cases at time, ttt
time_period <- 10000 ## input the time period that you wish to run the model for


t <- seq(0:time_period)
x <- sir_model(alpha = alpha, imported_cases = imported_cases, ttt = ttt) ## include any updated paramters as arguments in the model function
mod_run <- x$run(t) ## run the model
out <- as.data.frame(mod_run)


###########
## PLOTS ##
###########

sir_col <- c("aquamarine2", "#8c8cd9", "#cc0044", "black", "pink", "gold", "dimgrey")

matplot(x = out[, 1], y = out[, 46:49], xlab = "Time (days)", ylab = "Number of individuals",
        main = "Discrete SIR model - stochastic", type = "l", col = sir_col,
        lwd = 3, lty = 1, xlim = c(0, 10000))
legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))