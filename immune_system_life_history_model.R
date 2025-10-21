# imports
library(mvtnorm)
library(matrixStats)

# internal functions
fold_rise = function(log2_NAb_pre=1, weights=1,
                     max_log2_NAb = 16, mu = 16/15*7.2, sigma=16/15*2.9){
  pmax(0,rnorm(n=length(log2_NAb_pre),
               mean=pmax(0,weights * mu * (1-log2_NAb_pre/max_log2_NAb)),
               sd=pmax(0,weights*sigma * (1-log2_NAb_pre/max_log2_NAb))))
} 

# infection model, based on observed titer (could (should?) also do by individal-antibody titers but this is easier for now)
probability_infected = function(log2_NAb_pre,sensitivity=1,gamma=0.8){
  titer = 2^(log2_NAb_pre %*% sensitivity)
  p = titer^(-gamma)
  return(p)
}

rdirichlet_copula <- function(n, m, rho, total, dirichlet_alpha=1,nonzero_positions=NULL) {
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
  
  G <- qgamma(U, rate = 1,shape=dirichlet_alpha)             # m x n exponentials (Gamma(shape=1))
  X <- t(G)                          # n x m
  X <- t(X / rowMaxs(X))             # normalize to simplex

  colnames(X) = paste('pathogen_',1:n,sep = '')
  rownames(X) = paste('antibody_',1:m,sep = '')
  
  return(X)
}

sample_rates_given_immunogenicity <- function(u, rho=0.4, 
                                      dirichlet_alpha=1, shape=1, mean_time = 30/21,
                                      clip = 1e-12) {
  n <- length(u)
  
  U  <- pgamma(u,rate=1,shape=dirichlet_alpha) # uniforms via known CDF
  U  <- pmin(pmax(U, clip), 1-clip)        # avoid 0/1
  z1 <- qnorm(U)
  
  z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  w = qgamma(pnorm(z2), shape = shape, scale = mean_time/shape)
  
  return(w)
}


# core model functions
intialize_pathogens = function(N_expected_antibodies_per_pathogen,N_pathogens,antibody_correlation_matrix,
                               dirichlet_alpha=1){
  
  # randomly draw number of actual antibodies per pathogen
  n_antibodies = 1+rpois(N_expected_antibodies_per_pathogen-1,n=N_pathogens)
  
  # set the size of antibody repertoire to cover the range needed for reasonable correlations
  N_global_antibodies = 100*N_pathogens * sum(n_antibodies)
  
  # define the pathogen immunogenicity (strenth of immune system stimulation)
  immunogenicity <- rdirichlet_copula(n=N_pathogens, m=N_global_antibodies, 
                                      rho=antibody_correlation_matrix, 
                                      total = n_antibodies,
                                      dirichlet_alpha = dirichlet_alpha)
  
  # define the pathogen sensitivity (strength of neutralization)
  sensitivity <- rdirichlet_copula(n=N_pathogens, m=N_global_antibodies, 
                                   rho=antibody_correlation_matrix, 
                                   total = n_antibodies,
                                   nonzero_positions = immunogenicity>0,
                                   dirichlet_alpha = dirichlet_alpha)
  
  # get rid of unnecessary zeros
  needed_antibodies = apply(immunogenicity>0,1,any)
  immunogenicity = immunogenicity[needed_antibodies,]
  sensitivity = sensitivity[needed_antibodies,]

  # normalize sensitivity so that summing over it gives the serum titer
  sensitivity = sweep(sensitivity, 2, colSums(sensitivity), "/")
  
  return(list(immunogenicity=immunogenicity,sensitivity=sensitivity))
}

# initial immune system (antibody titer model)
initialize_immune_system = function(pathogens,Duration){
  
  # antibody titers per antibody
  log2_NAb = matrix(rep(0,nrow(pathogens$sensitivity)*Duration),nrow=Duration)
  colnames(log2_NAb) = rownames(pathogens$immunogenicity)
  
  # initialize waning rate per antibody
  waning_rate = rep(0,nrow(pathogens$immunogenicity))
  names(waning_rate) = rownames(pathogens$immunogenicity)
  
  # initialize broading rate per antibody
  broadening_rate = rep(0,nrow(pathogens$immunogenicity))
  names(broadening_rate) = rownames(pathogens$immunogenicity)
  broadening_ratio = rep(0,nrow(pathogens$immunogenicity))
  names(broadening_ratio) = rownames(pathogens$immunogenicity)
  parent = rep(NA,nrow(pathogens$immunogenicity))
  names(parent) = rownames(pathogens$immunogenicity)
  
  
  return(list(log2_NAb=log2_NAb,waning_rate=waning_rate,broadening_rate=broadening_rate,
              parent=parent))
}

# poisson exposure model
initialize_poisson_exposure = function(exposure_rate,Duration,N_pathogens){
  
  time_exposed = sort(1+sample(1:(Duration-1),size=floor(Duration*exposure_rate))) # 1+ so first infection is always after month 1
  pathogen_exposed = sample(1:N_pathogens,size=length(time_exposed),replace = TRUE)
  
  return(list(time_exposed=time_exposed,pathogen_exposed=pathogen_exposed))
}

# immune system dynamics model
immune_system_life_history = function(pathogens,exposures,Duration,
                                      gamma=0.8,
                                      max_log2_NAb = 16, mu = 16/15*7.2, sigma=16/15*2.9,
                                      dirichlet_alpha=1,
                                      shape_waning=1, mean_waning_time=21/30,  
                                      shape_broadening = 1, mean_broadening_time=2,
                                      max_mean_broadening_slots = 10, lambda_broadening = 0.1){
  
  immune_system = initialize_immune_system(pathogens,Duration)
  
  exposure_counter = 0
  infected = 0 * exposures$time_exposed
  p_inf=infected
  
  for(k in exposures$time_exposed[1]:Duration){

    # wane everything that's been initialized
    initialized_NAb_idx = immune_system$log2_NAb[k-1,]>0
    
    immune_system$log2_NAb[k,initialized_NAb_idx] = pmax(0,immune_system$log2_NAb[k-1,initialized_NAb_idx] -
                                        immune_system$waning_rate[initialized_NAb_idx]*1.44)
    
    # broadenening step
    broadened_idx = which(immune_system$parent != names(immune_system$parent))
    parent_idx = match(immune_system$parent[broadened_idx],names(immune_system$log2_NAb[k,]))

    # iterated version of log2NAb_child = log2NAb_parent * p_max * (1 - exp(-w * t))
    # log2NAb_child[k] = log2NAb_parent[k] * p_max * (1 - (1 - log2NAb_child[k-1]/(p_max*log2NAb_parent[k-1]))*exp(-w))
    w = immune_system$waning_rate[broadened_idx]
    Npk = immune_system$log2_NAb[k,parent_idx] * immune_system$broadening_ratio[broadened_idx]
    Npkm1 = immune_system$log2_NAb[k-1,parent_idx] * immune_system$broadening_ratio[broadened_idx]

    started_idx = Npkm1>0
    immune_system$log2_NAb[k,broadened_idx[started_idx]] = Npk[started_idx] *
      (1 - (1 - immune_system$log2_NAb[k-1,broadened_idx[started_idx]]/Npkm1[started_idx])*exp(-w[started_idx]))

    # check exposure
    if (k %in% exposures$time_exposed){
      exposure_counter = exposure_counter + 1
      
      p_inf[exposure_counter] = probability_infected(log2_NAb_pre = immune_system$log2_NAb[k-1,],
                                   sensitivity=pathogens$sensitivity[,exposures$pathogen_exposed[exposure_counter]],
                                   gamma=gamma)
      
      # if infected, boost
      if (runif(1)<p_inf[exposure_counter]){
        infected[exposure_counter]=1
        infection_idx = pathogens$immunogenicity[,exposures$pathogen_exposed[exposure_counter]] >0
        
        if (any(!initialized_NAb_idx & infection_idx)){
          # effects of affinity maturation are broken into two pieces:
          # first part of getting better at initializing antigen is covered by correlation of
          # waning rate with immunogenicity here
          
          immune_system$waning_rate[!initialized_NAb_idx & infection_idx] =
            sample_rates_given_immunogenicity(u = pathogens$immunogenicity[!initialized_NAb_idx & infection_idx,
                                                                                  exposures$pathogen_exposed[exposure_counter]],
                                                     dirichlet_alpha = dirichlet_alpha,
                                                     shape=shape_waning, mean_time=mean_waning_time)
          
          # parent is self
          immune_system$parent[!initialized_NAb_idx & infection_idx] = names(which(!initialized_NAb_idx & infection_idx))
            
          # second part of broadening rate here
          immune_system$broadening_rate[!initialized_NAb_idx & infection_idx] =
            sample_rates_given_immunogenicity(u = pathogens$immunogenicity[!initialized_NAb_idx & infection_idx,
                                                                                  exposures$pathogen_exposed[exposure_counter]],
                                                     dirichlet_alpha = dirichlet_alpha,
                                                     shape=shape_broadening, mean_time=mean_broadening_time)
          
        }
        
        # "fast" immune response to present antigens
        immune_system$log2_NAb[k,infection_idx] = immune_system$log2_NAb[k-1,infection_idx] + 
          fold_rise(log2_NAb_pre = immune_system$log2_NAb[k-1,infection_idx], 
                    weights = pathogens$immunogenicity[infection_idx,exposures$pathogen_exposed[exposure_counter]],
                    max_log2_NAb = max_log2_NAb, mu = mu, sigma=sigma)
        
        # second part of affinity maturation is broadening over time which needs the additional logic here
        # figure out who gets to broaden within the space of antibodies living in the simulation
        # probability of broadening correlated with immunogenicity 
        expansion_weights = (immune_system$log2_NAb[k ,infection_idx])
        expansion_weights = expansion_weights / max(expansion_weights)* max_mean_broadening_slots
        n_expansion = rpois(n=sum(infection_idx),lambda = expansion_weights)

        parents = which(infection_idx)
        parents = rep(parents, n_expansion)

        if(sum(n_expansion)>sum(immune_system$log2_NAb[k,]==0)){
          parents = sample(parents,size=sum(immune_system$log2_NAb[k,]==0))
        }

        # weight by pathogen similarity
        w = colMeans(pathogens$immunogenicity[,exposures$pathogen_exposed[exposure_counter]] *
                       pathogens$immunogenicity[,-exposures$pathogen_exposed[exposure_counter]])
        
        pathogen_idx = matrix(as.numeric(pathogens$immunogenicity[which(immune_system$log2_NAb[k,!infection_idx]==0),
                                            -exposures$pathogen_exposed[exposure_counter]]>0),ncol=ncol(pathogens$immunogenicity)-1)
        w = pathogen_idx %*% w
        
        # w=1+0*w
        
        if (any(w>0)){
          
          if(sum(w>0)<length(parents)){
            parents = sample(parents,size=sum(w>0))
          }
            
          children = sample(which(immune_system$log2_NAb[k,!infection_idx]==0), size=length(parents), prob=w)
  
          immune_system$waning_rate[children] = immune_system$waning_rate[parents]
          immune_system$broadening_ratio[children] = exp(-lambda_broadening*runif(n=length(children)))
          immune_system$broadening_rate[children] = immune_system$broadening_rate[parents]
          immune_system$parent[children] = names(parents)
  
          # I have a choice about what happens to originally cross-reactive antibodies when directly stimulated
          # this choice centers the antibody family around the newly boosted antigen and no longer the old
          boosted_broadened_idx = which(infection_idx & (immune_system$parent != names(immune_system$parent)))
          boosted_family_idx = which((immune_system$parent %in% immune_system$parent[boosted_broadened_idx]))
          boosted_family_idx = setdiff(boosted_family_idx,boosted_broadened_idx)
          matched_family_idx = match(immune_system$parent[boosted_family_idx],immune_system$parent[boosted_broadened_idx])

          # parent and child swap places
          # child becomes parent
          immune_system$parent[boosted_broadened_idx] = names(boosted_broadened_idx)
          # old parent and siblings become child of new parent
          immune_system$parent[boosted_family_idx] = names(boosted_broadened_idx)[matched_family_idx]
          
          # create broadening ratio for new children
          immune_system$broadening_ratio[boosted_family_idx[is.na(immune_system$broadening_ratio[boosted_family_idx])]] =
            exp(-lambda_broadening*runif(n=length(boosted_family_idx[is.na(immune_system$broadening_ratio[boosted_family_idx])])))
          
          # familial NAbs boost with parent boost
          immune_system$log2_NAb[k,boosted_family_idx] = pmax(immune_system$log2_NAb[k,boosted_family_idx],
            immune_system$log2_NAb[k-1,boosted_family_idx] +
            immune_system$broadening_ratio[boosted_family_idx] * 
            (immune_system$log2_NAb[k,boosted_broadened_idx[matched_family_idx]]-
               immune_system$log2_NAb[k-1,boosted_broadened_idx[matched_family_idx]]))
        }
      } 
    }
  }
  
  # calculate observed aggregate NAb traces
  serum_NAb=(2^(immune_system$log2_NAb)) %*% pathogens$sensitivity
  
  return(list(immune_system = immune_system,
              serum_NAb=data.frame(year = (1:Duration)/12,serum_NAb),
              infections = data.frame(year=exposures$time_exposed/12,pathogen=exposures$pathogen_exposed,
                                      infected=infected,p_inf=p_inf)))
}



## plot helpers
gg_serum_titers = function(person){
  serum_plot = person$serum_NAb |> pivot_longer(-year,names_to='pathogen',values_to='NAb',names_prefix='pathogen_') |>
    left_join(person$infections |> mutate(pathogen = as.character(pathogen)))
  p=ggplot(serum_plot) +
    geom_line(aes(x=year,y=NAb,color=pathogen)) +
    geom_point(data = serum_plot |> filter(!is.na(infected)),aes(x=year,y=NAb,color=pathogen)) +
    theme_bw() +
    scale_y_continuous(trans='log2')
  return(p)
}

gg_antibody_histories_by_waning_quintile = function(N_pathogens,Duration,pathogens,person){
  ## individual antibodies
  for (k in 1:N_pathogens){
    antibodies = rownames(pathogens$sensitivity)[pathogens$sensitivity[,k]>0 & person$immune_system$waning_rate>0]
    tmp_waning_rates = person$immune_system$waning_rate[antibodies]
    
    if(!is_empty(tmp_waning_rates)){
      # select waning rates by quintiles
      n_per <- pmax(1,pmin(40,floor(length(antibodies)/5)))
      take <- data.frame(waning_rate = tmp_waning_rates, i = seq_along(tmp_waning_rates) ,
                         row.names = names(tmp_waning_rates)) |>
        mutate(waning_rate_quintile = ntile(waning_rate, 5)) |>
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
        left_join(take |> select(-i), by='antibody') |>
        group_by(waning_rate_quintile,year) |>
        mutate(mean_log2_titer = log2(sum((2^log2_titer)*sensitivity)/sum(sensitivity))) |>
        mutate(pathogen=k)
      
      if (k==1){ 
        antibody_plot_dat = plot_dat
      } else {
        antibody_plot_dat = rbind(antibody_plot_dat,plot_dat)
      }
    }
  }
  
  p=ggplot(antibody_plot_dat,aes(x=year,y=log2_titer)) +
    geom_line(aes(group=antibody),alpha=0.1) +
    geom_line(aes(y=mean_log2_titer),color='blue') +
    facet_grid('waning_rate_quintile~pathogen') + 
    theme_bw() +
    guides(color='none') + ylab('log2(titer)')
  
  return(p)
}

# mab cross-reactivity over time plots
