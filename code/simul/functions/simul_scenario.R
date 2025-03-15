## To simulate different scenarios with thermal mortality and plant scarcity
## where every day a thermal extreme can happen with prob_ext probability
## and with an intensity from a poisson distribution of lambda therm_mor
## therm_mor can also be a vector of lambda's,
## in which case each extreme day a lambda is previously chosen at random (sample, uniform dist.)
## Each simulation uses a bootstrapped transition matrix
## INPUT
## nsims: number of simulations
## ndays: number of days of the population projection
## broods: the simulated period
## data: list with juvenile, fecundity & adult data + plant scarcity data
## prob_35: probability of a thermal extreme event
## prob_40: probability of a thermal LLHI extreme event
## prob_drought: probability of extreme drought
## rob_test: survival after predatory pressure applied


## output: tibble with the estimated net reproductive rate of each simulation
## and the associated information of the simulation

simul_scenario <- function(nsims = 1000,
                           ndays = 91,
                           broods,
                           data,
                           prob_35 = 0,
                           prob_40 = 0,
                           prob_drought = 0,
                           rob_test = 1) {
  
    plant_exp <- data$plsc  # data of the plant experiment
  
    plant <- rbinom(nsims, 1, prob = prob_drought)
    
    sims <- expand_grid(sim = 1:nsims,
                        days = 1:ndays) %>%
      mutate(tx35 = rbinom(nrow(.), 1, prob = prob_35), # an extreme event occurs?
             tx40 = rbinom(nrow(.), 1, prob = prob_40/prob_35), #a very extreme and rare event occurs?
             tx40 = if_else(tx35 == 0, 0, as.double(tx40)),
             mort = thermal_mortality(group = paste0(tx35, tx40)),  #the intensity of the mortality event
             sim_with_plant = rep(plant, each = ndays)) %>%  # whether the simulation is with plant scarcity or not
      nest(sim_df = -c(sim, sim_with_plant)) %>%
      mutate(A = future_map(sim,
                            ~ simul_matrix(data_juv = data$juv,
                                           data_repr = data$repr,
                                           data_ad = data$ad,
                                           brood_period = broods)),
             # we create a transition matrix per each simulation and return it with the associated bootstrapped samples
             # future_map should be applied to an ungrouped df, as now!
             dims = map(A, pluck, "dims"), #dimensions of the transition matrix (m & a)
             A = future_map_if(.x = A,
                               .p = sim_with_plant == 1,
                               ~ plant_scarcity(A = .x,
                                                data_plant = plant_exp)), # a bootstrapped transition matrix if plant scarcity occurs
             A = map_if(.x = A, .p = sim_with_plant == 0, pluck, "A"), # in simulations where there's no plant scarcity, take just A
             A_seq = map2(sim_df, A,
                          ~ change_juvs_mult(A = .y,
                                             mort = .x$mort,
                                             larvmort_robtest = rob_test)), # creating a sequence of daily transition matrices depending on thermal mortality and plant scarcity
             growth_rate = future_map2(A_seq, dims,
                                       ~ project_pop(
                                         A = .x,
                                         Aseq = 1:length(.x),
                                         vector = c(rep(0, .y["m"]), 1, rep(0, .y["a"] - 1))),
                                       .options = furrr_options(chunk_size = 30)), # projection of the population with the given sequence of transition matrices
             A_seq = NULL,
             growth_rate = future_map(growth_rate,
                                      ~ calculate_growth(vec_transenv = ., start_name = "a01"))) %>% # calculation of the net reproductive rate of the population
       unnest(growth_rate)

  return(sims)
}
