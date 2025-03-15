# To create a daily sequence of transition matrices
# considering daily thermal mortality and plant scarcity
# INPUT
## A: reference transition matrix (A)
## mort: vector with daily thermal mortalities
# OUTPUT
## list with an entry per each daily transition matrix to feed the population projection

change_juvs_mult <- function(A, mort, larvmort_robtest = 1) {
  
  # Dimensions of the transition matrix
  ## maximum juvenile age
  m <- length(which(A[1,] == 0))
  ## maximum adult age
  a <- ncol(A) - m
  
  # Sequence of transition matrices
  A_new <- list()
  
  
  for(j in seq_along(mort)) { #for each day of the projection
    
    A_new[[j]] <- A
    
    for (i in 2:(m+a)) { #then, add thermal mortality at juvenile parameters
      if(i <= m) {
        #surviving and resting as juveniles
        A_new[[j]][i, i-1] <- A[i, i-1]*(1-mort[j])*larvmort_robtest
        
      } 
      if(i == (m+1)) {
        A_new[[j]][i, 1:(m)] <- A[i, 1:(m)]*(1-mort[j])*larvmort_robtest
      }
    }
  }
  
  return(A_new)
}




