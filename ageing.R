
## Adding ageing into a simple SIR simulator
## Copyright: Amy Dighe 2018

    ## set the age increment for ageing in the first 6 months of life
    ## individuals will move into the next age compartment 
    ## every time this increment comes around
    
      unit_of_age_1 <- 30 #days

    ## set the age increment for between 6 and 24 months
    ## individuals will move into the next age compartment 
    ## every time this increment comes around
    
      unit_of_age_2 <- 91 #days

      age_cats_pre_6m <- 6 #number of age categories dividing calves under 6m
      age_cats_post_6m <- 6 #number of age categories dividing calves between 6m-2yrs    


    ## birth process: any new individual will enter 'Rp[,1]' where they are 
    ## protected from MERS-CoV by maternal antibodies for the first 6 months of life, 
    
      new_births <- rpois(1, N_0 * alpha) # where N_0 is the total number of adult camels initially
    
    ## deaths per time step in the first 6 months of life (as a stochastic process)
    
      deaths <- rbinom(1, Rp[i, j], prob = r2p(mu))
    

    Rp[i, 1] <- if (t %% unit_of_age_1 = 0) 0 + new_births else Rp[i-1] + new_births - deaths

    Rp[i, 2:age_cats_pre_6m] <- if (t %% unit_of_age_1 = 0) 0 + Rp[i-1, j-1] - deaths else Rp[i-1, j] - deaths

    ## individuals leave 'Rp' either by losing maternal 
    ## antibody protection and becoming susceptible 'S' or by dying
    ## individuals leave 'S' by becoming infected 'I' or dying

      rate_infection <- beta * I[i - 1] / N[i - 1]
      outflow_S[j] <- rbinom(1, S[i - 1, j],
                    prob = r2p(rate_infection + mu))   

    S[i, 1] <- S[i - 1] - outflow_S[1] # removing new infections and deaths

    S[i, 1] <- if (t %% 91 = 0) S[i, 1] - S[i, 1] else S[i, 1] # removing calves that have aged

    S[i, 1] <- if (t %% 30 = 0) S[i, 1] + Rp[i - 1, age_cats_pre_6m] else S[i, 1] # add calves from Rp whose maternal immunity has waned

      new_infectious[j] <- rmultinom(1, size = outflow_S[j],
                            prob = c(rate_infection, mu))[1]

    S[i, 2:age_cats_post_6m] <- s[i - 1, j] - outflow_S[j]
    
    S[i, 2:age_cats_post_6m] <- if (t %% 91 = 0) S[i, j] - S[i, j] + S[i, j-1]

      outflow_I <- rbinom(1, I[i - 1], prob = r2p(sigma + mu)) 
    I[i] <- I[i - 1] + new_infectious - outflow_I
    new_recovered <- rmultinom(1, outflow_I, prob = c(sigma, mu))[1]
    
    
I[i, 1] <- I[i - 1, 1] + new_infectious[1] - outflow_I[1] # removing those whose immunity has waned and deaths

I[i, 1] <- if (t %% 91 = 0) I[i, 1] - I[i, 1] + S[i, age_cats_post_6m]
  
I[i, 2:age_cats_post_6m] <- I[i - 1, j] + new_infectious[j] - outflow_I[j]

I[i, 2:age_cats_post_6m] <- if (t %% 91 = 0)

I_A[i]
  
R[i, 1] <- 

R[i, 2:age_cats_post_6m] <- 

R_A[i]