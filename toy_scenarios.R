# toy scenarios

library(tidyverse)
source('immune_system_life_history_model.R')

######### Scenario #1: 70 years of life with poisson exposure to three pathogens ######
######### pathogens 1 and 2 are closely related (from the same serogroup)
######### pathogen 3 is unrelated

# configure
Duration = 70*12 # months
N_expected_antibodies_per_pathogen = 1e3 # this is really meant to be few epitopes per pathogen times many antibodies per epitope...
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
immune_system = initialize_immune_system(pathogens,Duration,shape=1, mean_decay_time=30/21)
exposures = initialize_poisson_exposure(exposure_rate=1/60,Duration,N_pathogens)

# run
person = immune_system_life_history(immune_system,pathogens,exposures,Duration)

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
N_expected_antibodies_per_pathogen = 1e3 # this is really meant to be few epitopes per pathogen times many antibodies per epitope...
N_pathogens = 3

# antibody correlation matrix that defines serogroups: 
# - pathogens 1 and 2 are different strains of the same serogroup
# - pathogen 3 is a different serogroup 
antibody_correlation_matrix = diag(1, N_pathogens)
antibody_correlation_matrix[1,2] <- 0.2 -> antibody_correlation_matrix[2,1]
antibody_correlation_matrix[1,3] <- 0.1 -> antibody_correlation_matrix[3,1]
antibody_correlation_matrix[2,3] <- 0.1 -> antibody_correlation_matrix[3,2]
rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
antibody_correlation_matrix

# initialize!
set.seed(10)
pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)
immune_system = initialize_immune_system(pathogens,Duration,shape=0.72, mean_decay_time=30/21)

# tOPV-ish
  exposures = data.frame(time_exposed = c(2,3,4, 5,6,7, 8,9,10, 18,19,20, 60,61,62),
                         pathogen_exposed = c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)) # don't want to think about co-infection
  
  # run
  person = immune_system_life_history(immune_system,pathogens,exposures,Duration,gamma=0.44)
  
  ## plot some fun stuff!
  
  # serum titers over time
  gg_serum_titers(person)
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)


# bOPV-ish
  exposures = data.frame(time_exposed = c(2,3, 5,6, 8,9, 18,19, 60,61),
                         pathogen_exposed = c(1,3,1,3,1,3,1,3,1,3)) # don't want to think about co-infection
  
  # run
  person = immune_system_life_history(immune_system,pathogens,exposures,Duration,gamma=0.44)
  
  ## plot some fun stuff!
  
  # serum titers over time
  gg_serum_titers(person)
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)
  
  
# mOPV1-ish
  exposures = data.frame(time_exposed = c(2, 5, 8, 18, 60),
                         pathogen_exposed = c(1,1,1,1,1)) # don't want to think about co-infection
  
  # run
  person = immune_system_life_history(immune_system,pathogens,exposures,Duration,gamma=0.44)
  
  ## plot some fun stuff!
  
  # serum titers over time
  gg_serum_titers(person)
  
  # sampled individual antibody traces and sensitivity-weighted average by waning rate quintile
  gg_antibody_histories_by_waning_quintile(N_pathogens,Duration,pathogens,person)



  
  
  
