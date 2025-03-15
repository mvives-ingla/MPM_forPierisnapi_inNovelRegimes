# functions to estimate thermal mortality depending on the event.
## group: type of thermal event


thermal_mortality <- function(group) {
  
  termor <- as.numeric()
  
  for(i in seq_along(group)){
    if(group[i] == "00") { #non-extreme event
      tx <- runif(1, 17, 35)
      termor[i] <- exp(-23.1705859770024+0.501511610449519*tx)
      } else if (group[i] == "10") { #non-LLHI extreme thermal event
        tx <- runif(1, 35, 40)
        termor[i] <- exp(-22.6557794795911+0.484215947292227*tx)
        } else if (group[i] == "11") { #LLHI extreme thermal event
          tx <- runif(1, 40, 45)
          termor[i] <- 1/(1+exp(-(-37.9697676066161+0.873108327565261*tx)))
        }
    }
  
  return(termor)

}









