
###############
## FUNCTIONS ##
###############



## converting rates (r) to probabilities (p)
## using p = 1 - exp(-r)

r2p <- function(r) {
  p <- 1 - exp(-r)
  return(p)
}