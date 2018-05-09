par(mgp = c(3, 0.5, 0), las = 2, mar = c(4, 3, 0.7, 0.5))
layout(matrix(c(1,3, 5, 2, 4, 6, 7, 9, 11, 8, 10, 12), byrow = TRUE, ncol = 3, nrow = 4))
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
time_period <- 3600
t <- seq(0:time_period)

## set importation rate for introducing infectious individuals
importation_rate <- 0.01

for(i in 1:(length(N_0))){
  
  ## include any user-defined parameters as arguments here
  x <- sir_model(alpha = alpha, beta = beta, gamma = gamma, Ab_susc = Ab_susc, 
                 mu_1m = mu_1m, mu_2y = mu_2y, mu_adult = mu_adult, N_0 = N_0[i],
                 importation_rate = importation_rate)
  
  ## run the model
  x_res <- as.data.frame(replicate(100, x$run(1:3600)[,c(1, 87:100, 108:113, 115, 117:119, 124, 127)]))
  
  message(i)
  #prop_imp_epi <- 
  time_window <- 1080:2160
  adult_seropositivity <- mean(rowSums(x_res[time_window, seq(from = 23, to = 23+(99*27), by = 27)]))/100
  juvenile_seropositivity <- mean(rowSums(x_res[time_window, seq(from = 24, to = 24+(99*27), by = 27)]))/100
  total_seropositivity <- mean(rowSums(x_res[time_window, seq(from = 25, to = 25+(99*27), by = 27)]))/100
  
  colMax <- function(data) sapply(data, max, na.rm = TRUE)
  
  IA_df <- x_res[time_window, seq(from = 22, to = 22+(99*27), by = 27)]
  NA_df <- x_res[time_window, seq(from = 21, to = 21+(99*27), by = 27)]
  IA_NA_df <- IA_df/NA_df
  adult_RNA_pos <- 100 * mean(colMax(IA_NA_df))
  
  Itot_df <- x_res[time_window, seq(from = 17, to = 17+(99*27), by = 27)]
  Importations_df <- x_res[time_window, seq(from = 26, to = 26+(99*27), by = 27)]
  NJ_df <- x_res[time_window, seq(from = 20, to = 20+(99*27), by = 27)]
  IJ_NJ_df <- (Itot_df - (IA_df + Importations_df)) / NJ_df
  juvenile_RNA_pos <-100 * mean(colMax(IJ_NJ_df))
  
  N_df <- x_res[time_window, seq(from = 19, to = 19+(99*27), by = 27)]
  
  total_RNA_pos <-100 * max(out$Itot[time_window]/ out$N[time_window])
  
  print(adult_seropositivity)
  print(adult_RNA_pos)
  
  print(juvenile_seropositivity)
  print(juvenile_RNA_pos)
  
  age_at_1st_inf <- matrix(dim(14,100)) 
  age_at_1st_inf_dens <- matrix(dim(14,100))   
  new_infections_tot <- vector(length = 100)
  for(k in 1:100){
  new_infections_tot[k] <- sum(x_res[time_window, (c(2:15) + ((k-1) * 27))])
  
  titlee <- c("N_0 = 100", "N_0 = 1000", "N_0 = 10,000", "N_0 = 100,000", "N_0 = 1,000,000", "N_0 = 10,000,000")
  
  for(j in 1:14){
    age_at_1st_inf[j] <- sum(x_res[1080:2159, (j + 1 + ((k-1)*27))]) / new_infections_tot[k]
  }
  
  age_at_1st_inf_dens[1:12,k] <- age_at_1st_inf[1:12]
  age_at_1st_inf_dens[13,k] <- age_at_1st_inf[13] / 12
  age_at_1st_inf_dens[14,k] <- age_at_1st_inf[14] / (12*5)
  
  }
 
  namess <- c("<1m","1<x>2m", "2<x>3m","3<x>4m",
              "4<x>5m","5<x>6m","6<x>7m","7<x>8m",
              "8<x>9m","9<x>10m","10<x>11m","11<x>12m",
              "1<x>2yrs", "2yrs +")  
  barplot(height = age_at_1st_inf_dens, width = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12, 60),
          las = 2, ylab = "age at 1st infection (proportion)", ylim = c(0, 0.1), 
          names.arg = namess, main = titlee[i])
  text(x = 70, y = 0.08, paste("total seroprevalence:", round(total_seropositivity, 1), "%"))
  text(x = 70, y = 0.075, paste("max RNA positivity:     ", round(total_RNA_pos, 1), "%"))
  text(x = 70, y = 0.06, paste("adult seroprevalence:", round(adult_seropositivity, 1), "%"))
  text(x = 70, y = 0.055, paste("max adult RNA positivity:     ", round(adult_RNA_pos, 1), "%"))
  text(x = 70, y = 0.04, paste("<2yrs seroprevalence:", round(juvenile_seropositivity, 1), "%"))
  text(x = 70, y = 0.035, paste("max <2yrs RNA positivity:   ", round(juvenile_RNA_pos, 1), "%"))
  
  library(scales)
  noo <- seq(from = 0, to = 99, by = 1)
  col_select_SIRN <- c(16 + noo*27, 17 + noo*27, 18 + noo*27, 19 + noo*27)
  col_select_SIRN_sort <- sort(col_select_SIRN, decreasing = FALSE)
  sir_col <- c("aquamarine2", "#cc0044", "#8c8cd9", "black", "pink", "gold", "dimgrey")
  
  matplot(x_res[time_window,col_select_SIRN_sort], xlab = "", ylab = "Number of individuals", 
          type = "l", col = rep(alpha(sir_col[1:4], 0.1), 100), lwd = 3, lty = 1, xaxt = "n")
  axis(side = 1, at = seq(from = 0, to = 3600, by = 360), labels = seq(from = 0, to = 10, by = 1))
  mtext(side = 1, line = 2, "Time (years)")
  legend("topright", lwd = 3, col = sir_col, legend = c("S", "I", "R", "N"))
}