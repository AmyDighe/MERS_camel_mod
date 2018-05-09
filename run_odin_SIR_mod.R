
sir_model <- odin::odin("odin_SIR_mod.R", verbose = FALSE)
sir_col <- c("aquamarine2", "#8c8cd9", "#cc0044", "#999966", "black")

alpha <- 0.001 # input a value for birth rate (default = 0)
mu <- 0.001 # input a value for death rate (default = 0)

time_period <- 365 # input the time period that you wish to run the model for

t <- seq(0:time_period)

x <- sir_model(mu = mu, alpha = alpha)
#x$run(0:10)
mod_run <- x$run(t)

matplot(x = mod_run[, 1], y = mod_run[, c(2:5,7)], xlab = "Time", ylab = "Number of individuals",
        main = "Discrete MSIR model - deterministic", type = "l", col = sir_col,
        lwd = 3, lty = 1)
legend("topright", lwd = 3, col = sir_col, legend = c("M", "S", "I", "R", "N"))