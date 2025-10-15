# covid-vax-like cohorts

library(tidyverse)
source('immune_system_life_history_model.R')


######### random strain antigenicity correlation model
######### pathogens 1 and 2 wuhan/WA-1 and Delta-ish
######### pathogen 3 is omicron-ish

# configure   
Duration = 26 # months
N_expected_antibodies_per_pathogen = 3e2
N_pathogens = 3

# antibody correlation matrix that defines serogroups: 
antibody_correlation_matrix = diag(1, N_pathogens)
antibody_correlation_matrix[1,2] <- 0.9 -> antibody_correlation_matrix[2,1]
antibody_correlation_matrix[1,3] <- 0.42 -> antibody_correlation_matrix[3,1]
antibody_correlation_matrix[2,3] <- 0.42 -> antibody_correlation_matrix[3,2]
rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
antibody_correlation_matrix

# original mrna-like
set.seed(10) # omicron-similarity is really sensitive to the seed

pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix,
                                alpha=1)

exposures = data.frame(time_exposed = c(11,12, 20),
                       pathogen_exposed = c(1,1,1)) 

# make plot like fig 1a of https://www.nejm.org/doi/full/10.1056/NEJMc2119912
reference_medians = expand.grid(label = factor(c('baseline','after\n2nd dose',
                                                 '7 months\nafter\n2nd dose','after boost','6 months\nafter boost'),
                                               levels = c('baseline','after\n1st dose','after\n2nd dose',
                                                          '7 months\nafter\n2nd dose','after boost','6 months\nafter boost')),
                                strain=c('WA1','Omicron')) |>
  mutate(NAb_observed = 10^c(log10(5),3.2,2.5,3.5,3,log10(5),1.5,1.25,2.9,2.2))

# run cohort
N_people = 1000
cohort=list()
# tried parallelizing using `future`, but it was slower than the simple loop
for (k in 1:N_people){
  cohort[[k]] = immune_system_life_history(pathogens,exposures,Duration,
                                            gamma=0, # IM vaccine doesn't protect from itself
                                            mu = 7.5,sigma=2.8,
                                            max_log2_NAb=24, shape=1.5,mean_decay_time=13/30)
}

sampling_strategy = data.frame(year = c(1,11,12,19,20,26)/12,
                              label = factor(c('baseline','after\n1st dose','after\n2nd dose',
                                         '7 months\nafter\n2nd dose','after boost','6 months\nafter boost'),
                                         levels = c('baseline','after\n1st dose','after\n2nd dose',
                                                    '7 months\nafter\n2nd dose','after boost','6 months\nafter boost')))

plot_dat = as.data.frame(cohort[[1]]$serum_NAb) |>
  rename(all_of(c(WA1='pathogen_1',Delta='pathogen_2',Omicron = 'pathogen_3')))|> 
  right_join(sampling_strategy, by='year') |>
  mutate(id=1)
for (k in 2:N_people){
  plot_dat = plot_dat |>
    rbind(as.data.frame(cohort[[k]]$serum_NAb) |>
            rename(all_of(c(WA1='pathogen_1',Delta='pathogen_2',Omicron = 'pathogen_3'))) |> 
            right_join(sampling_strategy, by='year') |>
            mutate(id=k))
}
plot_dat = plot_dat |> pivot_longer(-c(year,label,id),names_to = 'strain',values_to = 'NAb') |>
  mutate(strain = factor(strain,levels=c('WA1','Delta','Omicron'))) |>
  mutate(NAb_observed = pmax(10,NAb*2^rnorm(nrow(plot_dat)))) |>
  mutate(NAb_observed = if_else(NAb_observed==10,5,NAb_observed))

ggplot(plot_dat,aes(x=label,y=NAb_observed,color=strain,group=interaction(label,strain))) +
  geom_jitter(width=0.2,size=0.5) +
  stat_summary(fun = median,geom = "crossbar",width = 0.8) +
  geom_hline(aes(yintercept = 10),linetype='dashed',linewidth=0.25) +
  # stat_summary(
  #   fun.data = function(x) {
  #     m <- median(log(x))
  #     s <- sd(log(x))/sqrt(length(x))
  #     data.frame(y = exp(m), ymin = exp(m- 2*s), ymax = exp(m + 2*s))
  #   },
  #   geom = "errorbar",
  #   width = 0.8, linewidth=1
  # )  +
  geom_crossbar(data=reference_medians,aes(xmin = as.numeric(label)-0.45,xmax = as.numeric(label)+0.45 )) +
  theme_bw() +
  xlab('') +
  scale_y_continuous(trans='log10',limits=10^c(0,5),breaks=10^seq(0,5))

# tl-dr is one can kind of get the right dynamics for closely related strains, but I'm not getting 
# affinity maturation "for free". But I think the pathogen model is wrong. So will try again!




### ESCAPE PATHOGEN EVOLUTION
# original mrna-like

# antibody correlation matrix that defines serogroups: 
antibody_correlation_matrix = diag(1, N_pathogens)
antibody_correlation_matrix[1,2] <- 0.9 -> antibody_correlation_matrix[2,1]
antibody_correlation_matrix[1,3] <- 0.42 -> antibody_correlation_matrix[3,1]
antibody_correlation_matrix[2,3] <- 0.42 -> antibody_correlation_matrix[3,2]
rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')


set.seed(10) # omicron-similarity is really sensitive to the seed

pathogens = intialize_escape_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix,
                                alpha=1, w_escape = c(1,0.1,0.9))


cohort=list()
for (k in 1:N_people){
  cohort[[k]] = immune_system_life_history(pathogens,exposures,Duration,
                                           gamma=0, # IM vaccine doesn't protect from itself
                                           mu = 7.5,sigma=2.8,
                                           max_log2_NAb=24, shape=1.5,mean_decay_time=13/30)
}

plot_dat = as.data.frame(cohort[[1]]$serum_NAb) |>
  rename(all_of(c(WA1='pathogen_1',Delta='pathogen_2',Omicron = 'pathogen_3')))|> 
  right_join(sampling_strategy, by='year') |>
  mutate(id=1)
for (k in 2:N_people){
  plot_dat = plot_dat |>
    rbind(as.data.frame(cohort[[k]]$serum_NAb) |>
            rename(all_of(c(WA1='pathogen_1',Delta='pathogen_2',Omicron = 'pathogen_3'))) |> 
            right_join(sampling_strategy, by='year') |>
            mutate(id=k))
}
plot_dat = plot_dat |> pivot_longer(-c(year,label,id),names_to = 'strain',values_to = 'NAb') |>
  mutate(strain = factor(strain,levels=c('WA1','Delta','Omicron'))) |>
  mutate(NAb_observed = pmax(10,NAb*2^rnorm(nrow(plot_dat)))) |>
  mutate(NAb_observed = if_else(NAb_observed==10,5,NAb_observed))

ggplot(plot_dat,aes(x=label,y=NAb_observed,color=strain,group=interaction(label,strain))) +
  geom_jitter(width=0.2,size=0.5) +
  stat_summary(fun = median,geom = "crossbar",width = 0.8) +
  geom_hline(aes(yintercept = 10),linetype='dashed',linewidth=0.25) +
  # stat_summary(
  #   fun.data = function(x) {
  #     m <- median(log(x))
  #     s <- sd(log(x))/sqrt(length(x))
  #     data.frame(y = exp(m), ymin = exp(m- 2*s), ymax = exp(m + 2*s))
  #   },
  #   geom = "errorbar",
  #   width = 0.8, linewidth=1
  # )  +
  geom_crossbar(data=reference_medians,aes(xmin = as.numeric(label)-0.45,xmax = as.numeric(label)+0.45 )) +
  theme_bw() +
  xlab('') +
  scale_y_continuous(trans='log10',limits=10^c(0,5),breaks=10^seq(0,5))

# not making a difference but I don't know why... something wrong with the immunogenicity zeros logic maybe





