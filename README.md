# MPM_forPierisnapi_inNovelRegimes
Data and code associated with the manuscript
Vives-Ingla, M., Capdevila, P., Clements, C. F., Stefanescu, C., & Carnicer, J. (2025). Novel regimes of extreme climatic events trigger negative population rates in a common insect. Global Change Biology.

<b>contact:</b>  
m.vives[at]creaf.uab.cat  
mariavivesingla[at]gmail.com

---
## Abstract  
The IPCC predicts that events at the extreme tail of the probability distribution will increase at a higher rate relative to less severe but still abnormal events. Such outlier events are of particular concern due to nonlinear physiological and demographic responses to climatic exposure, meaning that these events are expected to have disproportionate impacts on populations over the next decades (so called low-likelihood, high-impact events —LLHI). Because such events are historically rare, forecasting how biodiversity will respond requires mechanistic models that integrate the fundamental processes driving biological responses to our changing climate. Here we built a matrix population model (MPM) from long-term monitored populations of an insect model species in a Mediterranean area. The model simultaneously integrates the effects of extreme microclimatic heat exposure and drought-induced host-plant scarcity on early life stages, a key methodological step forward because these understudied life stages are usually very susceptible to climatic events. This model for the first time allowed us to forecast the demographic impacts that LLHI events will have on a well-known insect considering their whole life cycle. We found that juveniles were the life stage with the largest relative contribution to population dynamics. In line with field observations, simulated population rates in current climatic regimes were importantly determined by drought impacts, producing a regional mosaic of non-declining and declining populations. The simulations also indicated that, in future climate scenarios not meeting the Paris Agreement, LLHI heat extremes triggered regionally-widespread and severe declines in this currently abundant species. Our results suggest that LLHI events could thus emerge as a critical new —but overlooked— driver of the declines in insect populations, risking the crucial ecosystem functions they perform. We suggest that process-based and whole-cycle modelling approaches are a fundamental tool with which to understand the true impacts of climate change.


---
## __`data/`__

### __`clima/`__

#### raw data  
  - __`microclimate_data.csv`__: hourly temperatures recorded in the rearing microhabitats of the two types of <i>Pieris napi</i> populations.  
  - __`tdt_experiments.txt`__: larval survival time recorded in the static assays of <i>P. napi</i> heat tolerance.  
  - __`interactive_atlas/`__: macroclimatic data for the Mediterranean region obtained from the IPCC interactive atlas (http://interactive-atlas.ipcc.ch/, DOI: 10.1038/s41597-022-01739-y) organised in monthly folders (from March to September) that include: 
    - __`CMIP6 - *`__: Predicted number of days with the CMIP6 models in which maximum temperature exceed 35 or 40 ºC, with bias adjustment or not, for different Global Warming Levels.
    - __`CORDEX Europe - *`__: Predicted number of days with the CORDEX Europe models in which maximum temperature exceed 35 or 40 ºC, with bias adjustment or not, for different Global Warming Levels.
    - __`historical/*`__: Estimated number of days with the CMIP6 or CORDEX Europe models in which maximum temperature exceeded 35 or 40 ºC, with bias adjustment or not, for different historical periods.
#### processed data
  - __`future_mortalities.csv`__: predicted values of daily thermal mortality obtained with the TDT dynamic model from a random subset of microclimatic records.   
  - __`*_extreme_prob.csv`__: estimated values for current and future probabilities of heat microclimatic events.  

- __`more_info_on_data.txt`__: more information about the variables in each data file.

### __`mpm`__

#### raw data
  - __`growth__chamber_experiments.csv`__: observations from the growth chamber experiments monitoring the offspring of 11 <i>P. napi</i> females.  
  - __`mean_daily_fec_articles.csv`__: fecundity rates of <i>P. napi</i> females extracted from various published articles with the \textsc{metaDigitise r} package.
  - __`ad_survival_oleracea`__: adult survival rates of <i>P. oleracea</i> females extracted from Kerr et al. 2020 Popul. Ecol. 62(1) with the \textsc{metaDigitise r} package.

#### processed data
  - __`A_*.RData`__: transition matrices obtained from experimental and bibliographic data for spring (`g1g2`), summer (`g2g3`) and dry summer (`plant`) generations.
  - __`sims_output/`__: simulation outputs obtained can be found in: https://doi.org/10.5281/zenodo.15033204 


__`more_info_on_data.txt`__: more information about the variables in each data file.



## __`code/`__

### __`clima/`__

  - __`thermal-mortality-from-tmax.R`__: to obtain the model equations relating predicted thermal mortality and daily maximum temperature.
  - __`Thermal_landscape_functions_mod.R`__: to calculate the TDT curves from experimental data and apply the TDT dynamical model. Modified from Rezende et al. 2020 (https://datadryad.org/dataset/doi:10.5061/dryad.stqjq2c1r).
  - __`extreme-event-probabilities.R`__: to estimate current and future probabilities of heat extremes.

### __`simul/`__
  - __`future_sims_imp0.R`__: to project the populations in main conditions (main set of simulations).
  - __`calcul_elassens.R`__: to calculate elasticities and sensitivities of simulates matrices.
  - __`future_sims.R`__: to project the populations without imposing a minimum pupal age of eclosion.
  - __`future_sims_diffmort.R`__: to project the populations at varying predation rates.
  - __`future_sims_difflong.R`__: to project the populations at varying juvenile and adult longevities.
  - __`future_sims_imp0_macro.R`__: to project the populations without microclimatic thermal buffering effects.
  - __`validation.R`__: to project the populations over recorded microclimatic series.
  - __`functions/`__: functions used for the simulations
    - __`build_matrix_*.R`__: to build transition matrices given survival and eclosion probabilities and fecundities, imposing a mininum age of pupal eclosion (`imp0`) or not.
    - __`simul_matrix_*.R`__: to simulate a transition matrix by bootstraping from survival, eclosion and fecundity data, with fixed adult and juvenile longevity or not (`long`).
    - __`plant_scarcity_*.R`__: to change juvenile mortality and eclosion due to drought-induced plant scarcity effects, imposing a mininum age of pupal eclosion (`imp0`) or not.
    - __`thermal_mortality.R`__: to estimate thermal mortality depending on maximum temperature.
    - __`change_juvs_mult.R`__: to change the parameters of juvenile mortality in the transition matrix given an external mortality pulse.
    - __`project_pop_alt.R`__: to project the dynamics of a population given a transition matrix and scenario.
    - __`simul_scenario.R`__: to perform several simulations of population dynamics in the specified conditions.

  
  
  