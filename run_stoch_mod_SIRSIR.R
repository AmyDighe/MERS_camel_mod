sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for infectivity (default = 0.3)
beta <- 0.003

## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0

## input values for the age dependent death rate

mu_1m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 10000

## input the time period that you wish to run the model for
time_period <- 10000 
t <- seq(0:time_period)

## introduce infectious individuals
ttt <- 1171 # time of introduction 
imported_cases <- 11 ## number of imported cases at time, ttt, (default = 0)

## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult,
               imported_cases = imported_cases, ttt = ttt)

## run the model
mod_run <- x$run(t)
out <- as.data.frame(mod_run)


###########
## PLOTS ##
###########

## PLOT OF BASIC DYNAMICS ##
sir_col <- c("aquamarine2", "#cc0044", "#8c8cd9", "black", "pink", "gold", "dimgrey")

matplot(x = out[, 1], y = out[, 94:97], xlab = "Time (days)", ylab = "Number of individuals",
        main = "Discrete SIR model - stochastic", type = "l", col = sir_col,
        lwd = 3, lty = 1, xlim = c(0, 10000))
legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))

## PLOT OF AGE CATEGORIES ##
par(las = 1)
matplot(x = out[, 1], y = out[,3:16], col = (rainbow(n = 20))[1:14], type = 'l', lty = 1,
        xlim = c(0, 3600), xlab = "years", ylab = "number of camels", xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 3600, by = 360), labels = seq(from = 0, to = 10, by = 1))
legend(title = "Age", "right", col = (rainbow(n = 20))[1:14], lty = 1, cex = 0.6, bty= "n",
       legend = c("1m", "2m", "3m", "4m", "5m", "6m", "7m", "8m", "9m", "10m", "11m", "12m", "1-2yrs", "adults"))
arrows(x0 = 720, y0 = 4000, x1= 1080, y1 = 5600, lwd = 3, col = "darkgray")
text(x = 700, y = 3700, "adults")

## HISTOGRAM OF AGE DISTRIBUTION ##
age_dist <- vector(length = 13)
width <- c(rep(30, 12), 360)
for(i in 1:13){
  age_dist[i] <- sum(out[720:1079, (i + 2)])/width[i]
}  
age_dist <- round(age_dist, 0)
namess <- c("<1m","1<x>2m", "2<x>3m","3<x>4m","4<x>5m","5<x>6m","6<x>7m","7<x>8m","8<x>9m","9<x>10m","10<x>11m","11<x>12m", "1<x>2yrs")
par(las = 2)
library(scales)
barplot(height = age_dist, width = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12), ylim = c(0, 1800), ylab = "frequency density", axes = FALSE,
        cex.names = 1, space = 0)
#lines(x = c(0,24), y = c(age_dist[2], 800), lty = 2, lwd = 2, col = "darkgrey")
#lines(x = c(12, 24), y = c(age_dist[12], age_dist[12]), lty = 2, lwd = 2, col = "darkgrey")
#lines(x = c(12, 24), y = c(age_dist[12] - 0.5*(age_dist[12] - age_dist[13]), age_dist[12] - 0.5*(age_dist_12 - age_dist[13])), lty = 2, lwd = 2, col = "darkgrey")
#polygon(x = c(12, 12, 24, 24), y = c(0, age_dist[12] - 0.5*(age_dist[12] - age_dist[13]), age_dist[12] - 0.5*(age_dist[12] - age_dist[13]), 0), border = NA, col = alpha("lightpink", 0.3))
axis(side = 2, at = seq(from = 0, to = 2000, by = 200))
end_point = 0.5 + 11
text(c(seq(0.5,end_point,by=1), 18), par("usr")[3]-0.25, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = namess, cex=1)

## PLOT OF SEASONAL BIRTHS ##
par(las = 1, oma = c(0, 0, 0, 2))
library(scales)
tttt <- seq(1:360)
birthratezz <- alpha * (1 + cos(3 * cos(pi * tttt / 360)))
plot(x = out$tt, y = out$births, xlim = c(0, 360), xlab = "time of year", ylab = "births", 
     pch = 19, col = alpha("lightsteelblue3",0.4), xaxt = "n")
lines(10000 * birthratezz, col = "cornflowerblue", lwd = 3)
axis(side = 1, at = c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360), tick = FALSE,
     labels = c("July", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June"))
axis(side = 1, at = c(15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345, 375), tick = TRUE, labels = FALSE)
axis(side = 4, at = pretty(10000*birthratezz), labels = c(0, 0.0005, 0.0010, 0.0015, 0.0020), 
     col = "cornflowerblue", lwd = 3)
polygon(x = c(105, 105, 255, 255), y = c(-2, 32, 32, -2), col = alpha("lightpink", 0.3), border = NA)     
legend("topright", bty = "n", pch = c(NA, 19), lty = c(1, NA), lwd = c(3, 1), 
       col = c("cornflowerblue", "lightsteelblue3"), legend = c("birth rate", "births (stochastic)"))



## PLOT OF REPEATS ##
sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 
## input the time period that you wish to run the model for
time_period <- 10000 
t <- seq(0:time_period)
## introduce infectious individuals
ttt <- 250 # time of introduction 
imported_cases <- 1 ## number of imported cases at time, ttt, (default = 0)

## include any user-set parameters as arguments here:
x <- sir_model(alpha = alpha, imported_cases = imported_cases, ttt = ttt)

## run for multiple stochastic realisations and plot here:
x_res <- as.data.frame(replicate(100, x$run(0:3600)[, 94:97]))
matplot(0:3600, x_res, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Discrete SIR model - stochastic - 100 realisations", type = "l", col = rep(alpha(sir_col[1:4], 0.1), 100),
        lwd = 3, lty = 1, xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 3600, by = 360), labels = seq(from = 0, to = 10, by = 1))
legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))
text("T", x = 250, y = 4000, cex = 1.2)



