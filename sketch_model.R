#' ---
#' title: "Toy multistrain immunity model from antibodies to epitopes"
#' output:
#'   md_document
#' 
#' ---
#' 
#' 
#' # Toy multistrain immunity model from antibodies to epitopes
#' 
#' 
#' **What?** Implement a multistrain extension of the typhoid/polio/etc intrahost immunity model in the context of a constant force-of-infection cohort model.
#' 
#' **What did I learn?**
#' 
#' **What's next?** 
#'   
#' <!-- more -->
#' 
#' ## Outline
#' 
#' 
#' ## The code
#'
#'  First, let's get some boilerplate out of the way and set up our environment.
#'
#+ echo=TRUE, message=FALSE, results = 'hide'
# imports
library(tidyverse)
library(mvtnorm)
library(matrixStats)

# helper functions
fold_rise = function(log2_NAb_pre=1, weights=1,
                     max_log2_NAb = 16, mu = 16/15*7.2, sigma=16/15*2.9){
  pmax(0,rnorm(n=length(log2_NAb_pre),
               mean=pmax(0,weights * mu * (1-log2_NAb_pre/max_log2_NAb)),
               sd=pmax(0,weights*sigma * (1-log2_NAb_pre/max_log2_NAb))))
} 

rdirichlet_copula <- function(n, m, rho, total, nonzero_positions=NULL) {
  # helper to define antibody relevances for each pathogen, with approximate antibody correlations (thanks chatGPT!)
  
  # For each component k=1..m, draw an n-vector ~ N(0, rho)
  Z <- rmvnorm(n = m, sigma = rho)   # m x n
  
  # filter to non-zero number of antibodies per pathogen
  r <- t(colRanks(-Z, ties.method = "min"))    # 1 = largest
  if(is.null(nonzero_positions)){
    keep <- r <= matrix(rep(total, each = nrow(Z)),nrow=nrow(Z))
  } else {
    keep=nonzero_positions
  }
  U <- pnorm(Z)                     # m x n uniforms from Gaussian copula
  U[!keep] <- 0
  
  G <- qexp(U, rate = 1)             # m x n exponentials (Gamma(shape=1))
  X <- t(G)                          # n x m
  X <- t(X / rowMaxs(X))             # normalize to simplex

  colnames(X) = paste('pathogen_',1:n,sep = '')
  rownames(X) = paste('antibody_',1:m,sep = '')
  
  return(X)
}

# core model functions
intialize_pathogens = function(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix){
  
  # randomly draw number of actual antibodies (* antibodies per epitope) per pathogen
  n_antibodies = 1+rpois(N_expected_antibodies_per_pathogen-1,n=N_pathogens)
  
  # set the size of antibody repertoire to cover the range needed for reasonable correlations
  N_global_antibodies = 100*N_pathogens * sum(n_antibodies)
  
  # define the pathogen immunogenicity (strenth of immune system stimulation)
  immunogenicity <- rdirichlet_copula(n=N_pathogens, m=N_global_antibodies, 
                                      rho=antibody_correlation_matrix, 
                                      total = n_antibodies)
  
  # define the pathogen sensitivity (strength of neutralization)
  sensitivity <- rdirichlet_copula(n=N_pathogens, m=N_global_antibodies, 
                                   rho=antibody_correlation_matrix, 
                                   total = n_antibodies,
                                   nonzero_positions = immunogenicity>0)
  
  # get rid of unnecessary zeros
  needed_antibodies = apply(immunogenicity>0,1,any)
  immunogenicity = immunogenicity[needed_antibodies,]
  sensitivity = sensitivity[needed_antibodies,]

  # normalize sensitivity so that summing over it gives the serum titer
  sensitivity = sweep(sensitivity, 2, colSums(sensitivity), "/")
  
  return(list(immunogenicity=immunogenicity,sensitivity=sensitivity))
}

# initial immune system (antibody titer model)
initialize_immune_system = function(pathogens,
                                    Duration,
                                    shape=1, mean_decay_time=30/21){
  
  # antibody titers per antibody
  log2_NAb = matrix(rep(0,nrow(pathogens$sensitivity)*Duration),nrow=Duration)
  colnames(log2_NAb) = rownames(pathogens$immunogenicity)
  
  # waning rate per antibody
  waning_rate = rgamma(nrow(pathogens$sensitivity), shape = shape, scale = mean_decay_time/shape)
  names(waning_rate) = rownames(pathogens$immunogenicity)
  
  
  return(list(log2_NAb=log2_NAb,waning_rate=waning_rate))
}

# poisson exposure model
initialize_poisson_exposure = function(exposure_rate,Duration,N_pathogens){
  
  time_exposed = sort(1+sample(1:(Duration-1),size=floor(Duration*exposure_rate))) # 1+ so first infection is always after month 1
  pathogen_exposed = sample(1:N_pathogens,size=length(time_exposed),replace = TRUE)
  
  return(list(time_exposed=time_exposed,pathogen_exposed=pathogen_exposed))
}

# infection model, based on observed titer (could (should?) also do by individal-antibody titers but this is easier for now)
probability_infected = function(log2_NAb_pre,sensitivity=1,gamma=0.46){
  titer = 2^(log2_NAb_pre %*% sensitivity/sum(sensitivity))
  p = titer^(-gamma)
  return(p)
}

# immune system dynamics model
immune_system_life_history = function(immune_system,pathogens,exposures,Duration){
  
  exposure_counter = 0
  
  for(k in exposures$time_exposed[1]:Duration){
    
    if (k %in% exposures$time_exposed){
      exposure_counter = exposure_counter + 1
      
      p_inf = probability_infected(log2_NAb_pre = immune_system$log2_NAb[k-1,],
                                   sensitivity=pathogens$sensitivity[,exposures$pathogen_exposed[exposure_counter]])
      
      if (runif(1)<p_inf){
        idx = pathogens$immunogenicity[,exposures$pathogen_exposed[exposure_counter]] >0
        
        immune_system$log2_NAb[k,idx] = immune_system$log2_NAb[k-1,idx] + 
                                            fold_rise(log2_NAb_pre = immune_system$log2_NAb[k-1,idx], 
                                                      weights = pathogens$immunogenicity[idx,exposures$pathogen_exposed[exposure_counter]])
        immune_system$log2_NAb[k,!idx] = pmax(0,immune_system$log2_NAb[k-1,!idx] - immune_system$waning_rate[!idx]*1.44)
      } else {
        immune_system$log2_NAb[k,] = pmax(0,immune_system$log2_NAb[k-1,] - immune_system$waning_rate*1.44)
      }
      
    } else {
      immune_system$log2_NAb[k,] = pmax(0,immune_system$log2_NAb[k-1,] - immune_system$waning_rate*1.44)
    }
  }
  
  # calculate observed aggregate NAb traces
  serum_NAb=(2^(immune_system$log2_NAb)) %*% pathogens$sensitivity
  
  return(list(immune_system = immune_system,
              serum_NAb=data.frame(year = (1:Duration)/12,serum_NAb)))
}


# config!

Duration = 50*12 # months
N_expected_antibodies_per_pathogen = 1e3 # this is really few epitopes times many antibodies per epitope...
N_pathogens = 3

# antibody correlation: 
# - pathogens 1 and 2 are different strains of the same serogroup
# - pathogen 3 is a different serogroup 
antibody_correlation_matrix = m <- diag(1, N_pathogens)
antibody_correlation_matrix[1,2] <- 0.9 -> antibody_correlation_matrix[2,1]
rownames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(antibody_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
antibody_correlation_matrix

# initialize!
set.seed(10)
pathogens = intialize_pathogens(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix)
immune_system = initialize_immune_system(pathogens,Duration)
exposures = initialize_poisson_exposure(exposure_rate=1/60,Duration,N_pathogens)

# run
person = immune_system_life_history(immune_system,pathogens,exposures,Duration)

## plot some fun stuff!

# serum titers over time
serum_plot = person$serum_NAb |> pivot_longer(-year,names_to='pathogen',values_to='NAb',names_prefix='pathogen_') 
ggplot(serum_plot) +
  geom_line(aes(x=year,y=NAb,color=pathogen)) +
  theme_bw() +
  scale_y_continuous(trans='log2')


## individual antibodies
for (k in 1:N_pathogens){
  antibodies = rownames(pathogens$sensitivity)[pathogens$sensitivity[,k]>0]
  tmp_waning_rates = person$immune_system$waning_rate[antibodies]
  
  # select waning rates by quintiles
  n_per <- 40
  take <- data.frame(waning_rate = tmp_waning_rates, i = seq_along(tmp_waning_rates) ,
                 row.names = names(tmp_waning_rates)) |>
    mutate(waning_rate_quintile = ntile(waning_rate, 5)) |>
    # filter(q %% 2 == 1) |> # grab odd bins
    group_by(waning_rate_quintile) |>
    slice_sample(n = n_per, replace = FALSE) |>
    mutate(rep = 1:n_per) |>
    ungroup() |>
    mutate(antibody = sub('antibody_','',names(tmp_waning_rates)[i])) |>
    mutate(sensitivity = pathogens$sensitivity[names(tmp_waning_rates)[i],k]) |>
    mutate(waning_rate_quintile = factor(c('first','second','third','fourth','fifth')[waning_rate_quintile],
                                         levels=c('first','second','third','fourth','fifth')))
  
  plot_dat =  data.frame(year = (1:Duration)/12,
                         person$immune_system$log2_NAb[,antibodies[take$i]]) |>
    pivot_longer(-year,names_to = 'antibody',values_to = 'log2_titer',names_prefix = 'antibody_') |>
    left_join(take |> select(-i)) |>
    group_by(waning_rate_quintile,year) |>
    mutate(mean_log2_titer = log2(sum((2^log2_titer)*sensitivity)/sum(sensitivity))) |>
    mutate(pathogen=k)
  if (k==1){ 
    antibody_plot = plot_dat
  } else {
    antibody_plot = rbind(antibody_plot,plot_dat)
  }

}

ggplot(antibody_plot,aes(x=year,y=log2_titer)) +
  geom_line(aes(group=antibody),alpha=0.1) +
  geom_line(aes(y=mean_log2_titer),color='blue') +
  facet_grid('waning_rate_quintile~pathogen') + 
  theme_bw() +
  guides(color='none') + ylab('log2(titer)')
