sir_model <- odin::odin("stoch_mod_SIRSIR.R", verbose = FALSE, skip_cache = TRUE)

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 0.3

## input a value for recovery rate, gamma

infectious_period <- 20 ## duration of infection in days
gamma <- 1/infectious_period

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
time_period <- 36000 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate)

## run the model
mod_run <- x$run(t)
out <- as.data.frame(mod_run)








###########
## PLOTS ##
###########

## PLOT OF BASIC DYNAMICS ##
sir_col <- c("aquamarine2", "#cc0044", "#8c8cd9", "black", "pink", "gold", "dimgrey")

matplot(x = out[, 1], y = out[, 108:111], xlab = "Time (days)", ylab = "Number of individuals",
        main = "Discrete SIR model - stochastic", type = "l", col = sir_col,
        lwd = 3, lty = 1, xlim = c(0, 36000))
legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))

## PLOT OF AGE CATEGORIES ## !!! FOR USE WITH AN IMPORTATION RATE OF 0 !!!

## set importation rate for introducing infectious individuals
importation_rate <- 0

## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0,
               importation_rate = importation_rate)

## run the model
mod_run <- x$run(t)
out <- as.data.frame(mod_run)

par(las = 1)
matplot(x = out[, 1], y = out[,3:16], col = (rainbow(n = 20))[1:14], type = 'l', lty = 1,
        xlim = c(0, 3600), xlab = "years", ylab = "number of camels", xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 3600, by = 360), labels = seq(from = 0, to = 10, by = 1))
legend(title = "Age", "right", col = (rainbow(n = 20))[1:14], lty = 1, cex = 0.6, bty= "n",
       legend = c("1m", "2m", "3m", "4m", "5m", "6m", "7m", "8m", "9m", "10m", "11m", "12m", "1-2yrs", "adults"))
arrows(x0 = 720, y0 = 4000, x1= 1080, y1 = 5600, lwd = 3, col = "darkgray")
text(x = 700, y = 3700, "adults")

## HISTOGRAM OF AGE DISTRIBUTION ## !!! FOR USE WITH AN IMPORTATION RATE OF 0 !!!
age_dist <- vector(length = 13)
width <- c(rep(30, 12), 360)

for(i in 1:13){
  age_dist[i] <- sum(out[720:1079, (i + 2)])/width[i]
}
namess <- c("<1m","1<x>2m", "2<x>3m","3<x>4m","4<x>5m","5<x>6m","6<x>7m","7<x>8m","8<x>9m","9<x>10m","10<x>11m","11<x>12m", "1<x>2yrs")
par(las = 2)
library(scales)
barplot(height = age_dist, width = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12), ylim = c(0, 1800), ylab = "frequency density", axes = FALSE,
        cex.names = 1, space = 0)
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

## PLOT OF SEROPREVALENCE WITH INCIDENCE CURVE ##
par(las = 1)
out$importations[out$importations == 0] <- NA

seroprevalence <- 100 * (out$Itot + out$Rtot + out$S_2) / (out$N)

plot(seroprevalence[0:3600], type = 'l', lwd = 3, col = "cadetblue", xaxt = 'n', xlab = "time (years)",
     ylab = "seroprevalence (%)", ylim = c(0,100))
axis(side = 1, at = c(seq(from = 0, to = 3600, by = 360)), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
lines(out$incidence_new_inf, lwd = 2, col = alpha("firebrick", 0.6))
points(out$importations, pch = 8, col = "black")
legend("topright", legend = c("seroprevalence", "incidence", "importations"), col = c("cadetblue", "firebrick", "black"), lty = c(1, 1, NA), pch = c(NA, NA, 8), lwd = c(3, 1, NA), bty= "n")

## PLOT OF SEROPREVALENCE WITH AGE ##

plot(out$seropoz_A[0:36000], type = 'l', lwd = 3, col = "mediumpurple1", xaxt = 'n', xlab = "time (years)",
     ylab = "seroprevalence (%)", ylim = c(0,100))
axis(side = 1, at = c(seq(from = 0, to = 36000, by = 3600)), labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"))
lines(out$seropoz_J, col = "lightpink", lwd = 3)
legend(x = 25000, y = 67, legend = c("seroprevalence in adults", "seroprevalence in calves"), col = c("mediumpurple1", "lightpink"), lty = c(1, 1), lwd = c(3, 3), bty= "n")

############################
## AGE AT FIRST INFECTION ##
############################
new_infections_tot <- sum(out[1080:2159, 87:100])
inf_prop <- vector(length = 14)
inf_prop_dens <- vector(length = 14)

namess <- c("<1m","1<x>2m", "2<x>3m","3<x>4m",
            "4<x>5m","5<x>6m","6<x>7m","7<x>8m",
            "8<x>9m","9<x>10m","10<x>11m","11<x>12m",
            "1<x>2yrs", "2yrs +")

for(i in 1:14){
  inf_prop[i] <- sum(out[1080:2159, (i + 86)]) / new_infections_tot
}
barplot(height = inf_prop)

inf_prop_dens[1:12] <- inf_prop[1:12]
inf_prop_dens[13] <- inf_prop[13] / 12
inf_prop_dens[14] <- inf_prop[14] / (12*5)

barplot(height = inf_prop_dens, width = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12, 60),
        las = 2, ylab = "proportion of infections", ylim = c(0, 0.1), main = "age at first infection", 
        names.arg = namess)

#####################
## PLOT OF REPEATS ##
#####################

## input a value for average daily birth rate (default = 0.0005)
alpha <- 0.00076 

## input a value for transmission rate (aka effective contact rate - 
## rate of contact between S and I and rate of transmission upon contact) (default = 0.3)
beta <- 0.3

## input a value between 0 and 1 for susceptibility experienced by individuals with maternal antibodies 
## 0 would mean mAbs afford complete protection from MERS (default)
## 1 would mean mAbs afford no protection at all
Ab_susc <- 0

## input values for the age dependent death rate

mu_1m <- 0.001 # death rate for 1st month of life
mu_2y <- 0.001 # death rate for the rest of the 1st 2 yrs of life
mu_adult <- 0.0005 # death rate in adulthood (2 yrs +)

## input an initial population size

N_0 <- 10000000

## input the time period that you wish to run the model for
time_period <- 36000 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult,
               importation_rate = importation_rate)

## run for multiple stochastic realisations and plot here:
x_res <- (replicate(100, x$run(0:36000)[,]))
matplot(0:3600, x_res, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Discrete SIR model - stochastic - 100 realisations", type = "l", col = rep(alpha(sir_col[1:4], 0.1), 100),
        lwd = 3, lty = 1, xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 3600, by = 360), labels = seq(from = 0, to = 10, by = 1))
legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))
text("T", x = 250, y = 4000, cex = 1.2)

x_res <- as.data.frame(replicate(100, x$run(0:36000)[,]))
########################################################################################################################
## run from here for check_model function ##############################################################################
########################################################################################################################

check_model <- function(n = 100, t = 0:3600, ...) {
  model <- sir_model(...)
  
  res <- as.data.frame(replicate(n, model$run(t)[, 94:97]))
  matplot(t, res, xlab = "Time", ylab = "Number of individuals",
          main = "SIR model", type = "l",
          col = rep(alpha(sir_col[1:4], 0.1), n),
          lwd = 3, lty = 1)
  legend("right", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"), bty = 'n')
}

check_model(importation_rate = 0)