# To simulate the dynamics of a population during time time-steps using the project function from popdemo package
# INPUT
## A: the transition matrix or list of matrices to use for the simulation
## vector: Stage distribution at time = 0
#### n to try all stage-specific distributions
#### a vector to try a particular distribution
## Aseq: the order of use of the transition matrices contained in A
# OUTPUT
## The array with the stage distributions for each time-step and initial distributions

project_pop <- function(A, vector = "n", time = 100, Aseq = "unif") {
  # projection (popdemo::project)
  transenv <- project(A,
                      time = time,
                      vector = vector,
                      Aseq = Aseq,
                      return.vec = T)
  vec_transenv <- vec(transenv) # Stage distributions
  
  return(vec_transenv)
}

# To calculate the net reproductive rate between the initial and the next projected broods
# INPUT
## vec_transenv: the array with stage distributions for each time step and initial distribution
## start_name: in case we did the projection with a single initial distribution, which was the name of the starting stage?
## Plot (T/F) showing how abundance of adults of age1 (a01) changes with time
## OUTPUT: a data.frame with the growth rates between the first and the second simulated broods

calculate_growth <- function(vec_transenv, start_name = NULL, plot = F) {
  
  if (length(dim(vec_transenv)) == 3) { # if tried different initial distribution, vec_transenv will be an array of dim = 3
    popintime <- vec_transenv %>% 
      array_branch(3) %>% # a list with a matrix per each initial distribution
      map(as.data.frame) %>% # a column for each stage and a row for each day of the projection
      map(rownames_to_column,
          var = "day") %>% 
      map(mutate,
          day = as.numeric(day),
          day = day-1) %>% 
      map(select,
          day,
          a01) %>% # taking only the number of adults of age 1 (a01)
      bind_rows(.id = "start") %>% 
      mutate(start = str_sub(start, start = -3, end = -1),
             start_stage = str_sub(start, start = 1, end = 1),
             start_age = str_sub(start, start = -2, end = -1)) 
  } else {
    popintime <- vec_transenv %>% # if we only started with an initial distribtuion vec_transenv is a matrix (cols: stages, rows: days)
      as.data.frame() %>% 
      rownames_to_column(var = "day") %>% 
      mutate(day = as.numeric(day),
             day = day-1) %>% 
      select(day, a01) %>%  
      mutate(start = start_name, # the stage we started from
             start_stage = str_sub(start, start = 1, end = 1),
             start_age = str_sub(start, start = -2, end = -1))
  }
  
  # considering the beginning of the next brood in the day where abundance of a01 begins to rise again after a general trend of decrease (3 or more days of decrease in a period of 4)
  popintime <- popintime %>% 
    group_by(start, start_stage, start_age) %>% 
    mutate(N_past = lag(a01, default = 1),
           ch = (a01-N_past)/N_past*100,
           ch = if_else(day < 10, NA_real_, ch),
           q95 = quantile(ch, probs = .95, na.rm = T),
           q5 = quantile(ch, probs = .05, na.rm = T),
           ch = if_else(ch == Inf, NA_real_, ch),
           ch = if_else(ch < q5 | ch > q95, NA_real_, ch),
           ch_past_3 = rollmean(ch, k = 3, align = "right", fill = 0, na.rm = T),
           ch_next_3 = rollmean(ch, k = 3, align = "left", fill = 0, na.rm = T),
           # incr = a01 - N_past,
           # incr_past = lag(incr),
           # sign = (incr > 0) & (incr > incr_past),
           # sign_past = lag(sign),
           sign = a01 > N_past, # T if pop is rising
           sign_past = lag(sign), # was the population rising in the preceding time-step?
           sign_last5 = rollsum(sign, k = 5, align = "right", fill = 1),
           sign_next5 = rollsum(sign, k = 5, align = "left", fill = 1),
           brood = ((sign & !sign_past) &
                      ((sign_last5 <= 2 & ch_next_3 > 0) |
                                           (sign_next5 >= 4 &
                                              -5 < ch_past_3 & ch_past_3 < 5 &
                                              0 < ch_next_3 & ch_next_3 < 5)) &
                      day > 10 |
                      day == 1), #preceded or followed by a clear trend of increase or decrease
           # (abs(ch_past) < 10 | ch_past == -100),
           brood = ifelse(is.na(brood), F, brood),
           # brood = ifelse(brood & day > 2 & day < 20, F, brood), #forcing the next brood starts beyond day 20
           brood = cumsum(brood)+1) 
  
  brood_border <- popintime %>% 
    filter(brood == 2) %>% 
    slice_max(day)
  
  if (plot) { #plotting how adults of age 1 changed with time
      timeseries <- popintime %>% 
      filter(start_stage == "j" | start  == "a01", day < 50, day > 0) %>% 
      ggplot(aes(x = day, y = a01)) +
      geom_line(aes(group = start), color = "grey") +
      geom_vline(data = brood_border, aes(xintercept = day), color = "blue") +
      facet_wrap(vars(start_stage), nrow = 2, scales = "free",
                 labeller = label_both)
    
    print(timeseries)
  }
  
  # estimating the growth rate between the first two generations of the projection
  suppressMessages(
    growth_rates <- popintime  %>% 
      filter(brood %in% 1:2) %>% 
      group_by(start, start_stage, start_age, brood) %>% 
      summarise(N = sum(a01)) %>% 
      mutate(brood = paste0("g", brood)) %>% 
      pivot_wider(id_cols = c(start, start_stage, start_age),
                  names_from = "brood",
                  values_from = "N") %>% 
      mutate(gr = g2/g1,
             gr = if_else(brood_border$day[1] < 50, gr, NA_real_)) # it's equivalent to the net reproductive rate (how many individuals of one age an individual of the same age gives birth during one generation period)
  ) 
  
  return(growth_rates)
  
}
