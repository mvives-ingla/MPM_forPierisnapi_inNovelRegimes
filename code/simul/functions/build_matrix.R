# A function to construct a daily-based matrix projection model
# and extract its main properties

# No restriction in minimum age for pupal eclosion

build_matrix <- function(juvsurvs, juvecls, adfecs, adsurvs) {
  
  a <- length(adfecs) # max adult longevity
  m <- length(juvsurvs) # max juvenile longevity
  
  A <- matrix(0, nrow = m+a, ncol = m+a) # transition matrix initialization
  
  for (i in 1:(m+a)) {
    if (i == 1) {
      A[i,(m+1):(m+a)] <- adfecs
      # fecundity vector (row 1, adult columns)
    } else if(i <= m) {
      A[i, i-1] <- juvsurvs[i-1]*(1-juvecls[i-1])
      # subdiagonal of juvenile survival and persist as a juvenile
      # survival*probability of not ecloding
    } else if(i == (m+1)) {
      A[i, 1:m] <- juvsurvs*juvecls
      # vector of pupal eclosion (row m+1, juvenile columns)
      # survival*probabilty of eclosion
    } else {
      if (length(adsurvs) == 1) {
        A[i, i-1] <- adsurvs
        # constant adult survival
      } else {
        A[i, i-1] <- adsurvs[i-m-1]
        # varying adult survival
      }
      #subdiagonal of adult survival
    }
  }
  
  ## Avoiding spontaneous generation of individuals and probabilities > 1
  col_prob <- which(colSums(A[-1,]) > 1)
  if(length(col_prob) > 0) {
    if(length(col_prob) > 1) {
      col_sums <- colSums(A[-1, col_prob])
    } else {
      col_sums <- sum(A[-1, col_prob])
    }
    for(i in col_prob) {
      if(i <= m) {
        A[i+1, i] <- A[i+1, i]/col_sums[which(col_prob == i)]
        # this is the same as (1- probability of eclosion)
        A[m+1, i] <- A[m+1, i]/col_sums[which(col_prob == i)]
        # this is the same as probability of eclosion
      } else {
        A[i+1, i] <- A[i+1, i]/col_sums[which(col_prob == i)]
      }
    }
  }
  
  
  colnames(A) <- rownames(A) <-  c(paste0("j0", 1:9),
                                   paste0("j", 10:m),
                                   paste0("a0", 1:9),
                                   paste0("a", 10:a))
  
  ir_test <- isPrimitive(A)+isIrreducible(A)+isErgodic(A)
  if(ir_test == 3) {
    return(A)
  } else {
      return("try-error") # should this be the case, simul_matrix() would take a new bootstrap sample
    } 

}
