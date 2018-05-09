par(mfrow = c(3,2))

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

N_0 <- c(100, 1000, 10000, 100000, 1000000, 10000000)

## input the time period that you wish to run the model for
time_period <- 36000 
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

for(i in 1:(length(N_0))){
  
## include any user-defined parameters as arguments here
x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, Ab_susc = Ab_susc, 
               mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0[i],
               importation_rate = importation_rate)

## run the model
mod_run <- x$run(t)
out <- as.data.frame(mod_run)

message(i)
#prop_imp_epi <- 

adult_seropositivity <- mean(out$seropoz_A[1080:2160])
juvenile_seropositivity <- mean(out$seropoz_J[1080:2160])
total_seropositivity <- mean(out$seropz_tot[1080:2160])

adult_RNA_pos <-100 * max(out$IA[1080:2160]/out$N_A[1080:2160])
juvenile_RNA_pos <-100 * max((out$Itot[1080:2160] - (out$IA[1080:2160] + out$importations[1080:2160]))/out$N_J[1080:2160])
total_RNA_pos <-100 * max(out$Itot[1080:2160]/ out$N[1080:2160])

print(adult_seropositivity)
print(adult_RNA_pos)

print(juvenile_seropositivity)
print(juvenile_RNA_pos)

new_infections_tot <- sum(out[1080:2159, 87:100])
age_at_1st_inf <- vector(length = 14) 
age_at_1st_inf_dens <- vector(length = 14) 
titlee <- c("N_0 = 100", "N_0 = 1000", "N_0 = 10,000", "N_0 = 100,000", "N_0 = 1,000,000", "N_0 = 10,000,000")
for(j in 1:14){
  age_at_1st_inf[j] <- sum(out[1080:2159, (j + 86)]) / new_infections_tot
}

age_at_1st_inf_dens[1:12] <- age_at_1st_inf[1:12]
age_at_1st_inf_dens[13] <- age_at_1st_inf[13] / 12
age_at_1st_inf_dens[14] <- age_at_1st_inf[14] / (12*5)

namess <- c("<1m","1<x>2m", "2<x>3m","3<x>4m",
            "4<x>5m","5<x>6m","6<x>7m","7<x>8m",
            "8<x>9m","9<x>10m","10<x>11m","11<x>12m",
            "1<x>2yrs", "2yrs +")  
barplot(height = age_at_1st_inf_dens, width = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12, 60),
        las = 2, ylab = "age at 1st infection (proportion)", ylim = c(0, 0.1), 
        names.arg = namess, main = titlee[i])
text(x = 70, y = 0.08, paste("total seroprevalence:", round(total_seropositivity, 1), "%"))
text(x = 70, y = 0.075, paste("total RNA positivity:     ", round(total_RNA_pos, 1), "%"))
text(x = 70, y = 0.06, paste("adult seroprevalence:", round(adult_seropositivity, 1), "%"))
text(x = 70, y = 0.055, paste("adult RNA positivity:     ", round(adult_RNA_pos, 1), "%"))
text(x = 70, y = 0.04, paste("<2yrs seroprevalence:", round(juvenile_seropositivity, 1), "%"))
text(x = 70, y = 0.035, paste("<2yrs RNA positivity:   ", round(juvenile_RNA_pos, 1), "%"))
}

