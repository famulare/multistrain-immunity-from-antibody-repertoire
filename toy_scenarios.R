# toy scenarios

library(tidyverse)
source('immune_system_life_history_model.R')

######### Scenario #1: 70 years of life with poisson exposure to three pathogens ######
######### pathogens 1 and 2 are closely related (from the same serogroup)
######### pathogen 3 is unrelated

# configure
Duration = 70*12 # months
N_expected_antibodies_per_pathogen = 1e2 # this is really meant to be few epitopes per pathogen times many antibodies per epitope...
N_pathogens = 3

# antibody correlation matrix that defines serogroups: 
# - pathogens 1 and 2 are different strains of the same serogroup
# - pathogen 3 is a different serogroup 
antibody_correlation_matrix = diag(1, N_pathogens)
antibody_correlation_matrix[1,2] <- 0.8 -> antibody_correlation_matrix[2,1]
rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
antibody_correlation_matrix

# initialize!
set.seed(10)
pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)
exposures = initialize_poisson_exposure(exposure_rate=1/48,Duration,N_pathogens)

# run
person = immune_system_life_history(pathogens,exposures,Duration)

## plot some fun stuff!

# serum titers over time
gg_serum_titers(person)

# sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)

# within-serogroup comparisons

serum_plot = data.frame(year=person$serum_NAb$year,
                        relative_titer_2over1 = person$serum_NAb$pathogen_2/person$serum_NAb$pathogen_1) |>
  left_join(person$infections |> filter(pathogen!=3) |> mutate(pathogen = as.character(pathogen)))

ggplot(serum_plot) +
  geom_line(aes(x=year,y=relative_titer_2over1)) +
  geom_point(data = serum_plot |> filter(!is.na(infected)),aes(x=year,y=relative_titer_2over1,color=pathogen)) +
  theme_bw() +
  scale_y_continuous(trans='log2')



######### Scenario #2: let's do something like polio vaccination ######
######### pathogens 1 and 2 are closely related (from the same serogroup)
######### pathogen 3 is unrelated

# configure
Duration = 20*12 # months
N_expected_antibodies_per_pathogen = 1e2 # this is really meant to be few epitopes per pathogen times many antibodies per epitope...
N_pathogens = 3

# antibody correlation matrix that defines serogroups: 
# - pathogens 1 and 2 are different strains of the same serogroup
# - pathogen 3 is a different serogroup 
antibody_correlation_matrix = diag(1, N_pathogens)
antibody_correlation_matrix[1,2] <- 0.3 -> antibody_correlation_matrix[2,1]
antibody_correlation_matrix[1,3] <- 0.1 -> antibody_correlation_matrix[3,1]
antibody_correlation_matrix[2,3] <- 0.1 -> antibody_correlation_matrix[3,2]
rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
antibody_correlation_matrix

# initialize!
# tOPV-ish  
  set.seed(10)
  pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)

  exposures = data.frame(time_exposed = c(2,3,4, 5,6,7, 8,9,10, 18,19,20, 60,61,62),
                         pathogen_exposed = c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)) # don't want to think about co-infection
  
  # what should I use for gamma?
  # x=2^seq(0,9,by=0.1)
  # VE_approx_OPV = 10^(-2.3/(x^0.44))
  # lazy_VE = 1-x^(-.15)
  # plot(log2(x),VE_approx_OPV)
  # lines(log2(x),lazy_VE)
  
  # run
  person = immune_system_life_history(pathogens,exposures,Duration,gamma=0.15,shape=0.72, mean_decay_time=1)
  
  ## plot some fun stuff!
  
  # serum titers over time
  gg_serum_titers(person)
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)


# bOPV-ish
  # set.seed(10) # no cross reaction
  # set.seed(1) # big cross reaction
  set.seed(5) # moderate cross reaction
  
  pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)

  exposures = data.frame(time_exposed = c(2,3, 5,6, 8,9, 18,19, 60,61),
                         pathogen_exposed = c(1,3,1,3,1,3,1,3,1,3)) # don't want to think about co-infection
  
  # run
  person = immune_system_life_history(pathogens,exposures,Duration,gamma=0.15,shape=0.72, mean_decay_time=30/21)
  
  ## plot some fun stuff!
  
  # serum titers over time
  gg_serum_titers(person)
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)
  
  
# mOPV1-ish
  set.seed(5)
  pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)

  exposures = data.frame(time_exposed = c(2, 5, 8, 18, 60),
                         pathogen_exposed = c(1,1,1,1,1)) # don't want to think about co-infection
  
  # run
  person = immune_system_life_history(pathogens,exposures,Duration,gamma=0.15,shape=0.72, mean_decay_time=30/21)
  
  ## plot some fun stuff!
  
  # serum titers over time
  gg_serum_titers(person)
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)



######## Scenario #3: let's do something COVID-like ######
######### pathogens 1 and 2 wuhan/WA-1 and Delta-ish
######### pathogen 3 is omicron-ish
  
  # configure
  Duration = 26 # months
  N_expected_antibodies_per_pathogen = 1e2
  N_pathogens = 3
  
  # antibody correlation matrix that defines serogroups: 
  antibody_correlation_matrix = diag(1, N_pathogens)
  antibody_correlation_matrix[1,2] <- 0.9 -> antibody_correlation_matrix[2,1]
  antibody_correlation_matrix[1,3] <- 0.4 -> antibody_correlation_matrix[3,1]
  antibody_correlation_matrix[2,3] <- 0.4 -> antibody_correlation_matrix[3,2]
  rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
  colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
  antibody_correlation_matrix
  
  # original mrna-like
  set.seed(10)
  pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix,
                                  alpha=1)

  exposures = data.frame(time_exposed = c(11,12, 20),
                         pathogen_exposed = c(1,1,1)) 
  
  # run
  # response model loosely informed by https://www.nejm.org/doi/full/10.1056/NEJMc2119912
  person = immune_system_life_history(pathogens,exposures,Duration,
                                      gamma=0, # IM vaccine doesn't protect from itself
                                      mu = 8.3,sigma=2,#sigma=3.3,
                                      max_log2_NAb=30)
  
  # serum titers over time
  gg_serum_titers(person) + scale_y_continuous(trans='log10',limits=10^c(-0.1,5), breaks = 10^c(0:5) ) #
   
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)
  

  # cross-responses
  serum_plot = data.frame(year=person$serum_NAb$year,
                          relative_titer_delta_over_wuhan = person$serum_NAb$pathogen_2/person$serum_NAb$pathogen_1,
                          relative_titer_omicron_over_wuhan = person$serum_NAb$pathogen_3/person$serum_NAb$pathogen_1) |>
    pivot_longer(-year,names_to = "strain",
                 names_pattern = "^relative_titer_(.*)_over_wuhan$",
                 values_to = "relative_titer"
    ) |>
    mutate(reference = "wuhan") |>
    left_join(person$infections |> mutate(pathogen = as.character(pathogen)))
  
  ggplot(serum_plot) +
    geom_line(aes(x=year,y=1/relative_titer,group=strain,color=strain)) +
    geom_point(data = serum_plot |> filter(!is.na(infected)),aes(x=year,y=1/relative_titer,color=strain)) +
    theme_bw() +
    scale_y_continuous(trans='log2') +
    scale_color_discrete(name='vs. WA/1')+
    ylab('fold-difference') 

  
  
  # original mrna and then omicron monovalent booster
  set.seed(10)
  pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)

  exposures = data.frame(time_exposed = c(11,12, 20),
                         pathogen_exposed = c(1,1,3)) 
  
  # run
  # response model loosely informed by https://www.nejm.org/doi/full/10.1056/NEJMc2119912
  person = immune_system_life_history(pathogens,exposures,Duration,
                                      gamma=0, # IM vaccine doesn't protect from itself
                                      max_log2_NAb=20,mu = 9, sigma=2)
  
  # serum titers over time
  gg_serum_titers(person) + scale_y_continuous(trans='log10',limits=10^c(-0.1,5), breaks = 10^c(0:5))
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)
  
  # does reproduce the OAS-like finding "for free"
  
  
  
  # original mrna and then wuhan-omicron bivalent booster
  set.seed(10)
  pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)
  
  # need to play with "dose" and solve why booster isn't as strong as it should be...
  pathogens$immunogenicity = cbind(pathogens$immunogenicity, 1*(pathogens$immunogenicity[,1]+pathogens$immunogenicity[,3]))
  pathogens$sensitivity = cbind(pathogens$sensitivity, 1*(pathogens$sensitivity[,1]+pathogens$sensitivity[,3]))

  exposures = data.frame(time_exposed = c(11,12, 20),
                         pathogen_exposed = c(1,1,4)) 
  
  # run
  # response model loosely informed by https://www.nejm.org/doi/full/10.1056/NEJMc2119912
  person = immune_system_life_history(pathogens,exposures,Duration,
                                      gamma=0, # IM vaccine doesn't protect from itself
                                      max_log2_NAb=20,mu = 9, sigma=2)
  
  person_plot = person
  person_plot$serum_NAb = person_plot$serum_NAb[,-5]
  # serum titers over time
  gg_serum_titers(person_plot) + scale_y_continuous(trans='log10',limits=10^c(-0.1,5), breaks = 10^c(0:5))
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)



  