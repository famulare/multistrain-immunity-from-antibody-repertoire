#' ---
#' title: "Toy multistrain immunity model from epitopes"
#' output:
#'   md_document
#' 
#' ---
#' 
#' 
#' # Toy multistrain immunity model from epitopes
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
rdirichlet_copula <- function(n, m, rho, total) {
  # helper to define epitope relevances for each pathogen, with approximate epitope correlations (thanks chatGPT!)
  stopifnot(all(dim(rho) == c(n, n)))
  if (min(eigen(rho, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("rho must be positive definite")
  
  # For each component k=1..m, draw an n-vector ~ N(0, rho)
  Z <- rmvnorm(n = m, sigma = rho)   # m x n
  
  # filter to non-zero number of epitopes per pathogen
  r <- t(colRanks(-Z, ties.method = "min"))    # 1 = largest
  keep <- r <= matrix(rep(total, each = nrow(Z)),nrow=nrow(Z))

  U <- pnorm(Z)                     # m x n uniforms from Gaussian copula
  U[!keep] <- 0
  
  G <- qexp(U, rate = 1)             # m x n exponentials (Gamma(shape=1))
  X <- t(G)                          # n x m
  X <- t(X / rowMaxs(X))               # normalize to simplex

  colnames(X) = paste('pathogen_',1:n,sep = '')
  rownames(X) = paste('epitope_',1:m,sep = '')
  return(X)
}

fold_rise = function(log2_NAb_pre,
                     max_log2_NAb = 15, mu = 6, sigma=2){
  pmax(0,rnorm(n=length(log2_NAb_pre),mean=mu,sd=sigma)*(1-log2_NAb_pre/max_log2_NAb))
} 


#' Next, let's define the **pathogen-epitope** model. 
# define pathogens by their epitopes

N_global_epitopes = 1e5
N_expected_epitopes_per_pathogen = 1e3
N_pathogens = 3

# epitope correlation: 
# - pathogens 1 and 2 are different strains of the same serogroup
# - pathogen 3 is a different serogroup 
epitope_correlation_matrix = m <- diag(1, N_pathogens)
epitope_correlation_matrix[1,2] <- 0.9 -> epitope_correlation_matrix[2,1]
rownames(epitope_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
colnames(epitope_correlation_matrix) = paste('pathogen_',1:N_pathogens,sep = '')
epitope_correlation_matrix


# define pathogen epitope repertoire 
n_epitopes = 1+rpois(N_expected_epitopes_per_pathogen-1,n=N_pathogens)

set.seed(1)
pathogens <- rdirichlet_copula(n=N_pathogens, m=N_global_epitopes, rho=epitope_correlation_matrix, 
                       total = n_epitopes)

colSums(pathogens)   
cor(pathogens)


#' Now, we can define the antibody repertoire and dynamics for a person
#' 
Duration = 50*12 # months

log2_NAb = matrix(rep(0,N_global_epitopes*Duration),nrow=Duration)
waning_rate = rgamma(N_global_epitopes, shape = 1, rate = 21/30)

infection_rate = 1/60
infection_times = sort(sample(2:Duration,size=floor(Duration*infection_rate)))
infecting_pathogen = sample(1:N_pathogens,size=length(infection_times),replace = TRUE)
infection_counter = 0
for(k in infection_times[1]:Duration){

  if (k %in% infection_times){
    infection_counter = infection_counter + 1
    idx = pathogens[,infecting_pathogen[infection_counter]] >0
    
    log2_NAb[k,idx] =  log2_NAb[k-1,idx] + pathogens[idx,infecting_pathogen[infection_counter]]*fold_rise(log2_NAb_pre = log2_NAb[k-1,idx])
    log2_NAb[k,!idx] = pmax(0,log2_NAb[k-1,!idx] - waning_rate[!idx]*0.69)
    
    } else {
    log2_NAb[k,] = pmax(0,log2_NAb[k-1,] - waning_rate*0.69)
  }
}



plot(log2( (2^log2_NAb) %*% (pathogens[,1])/sum(pathogens[,1])))
plot(log2( (2^log2_NAb) %*% (pathogens[,2])/sum(pathogens[,2])))
plot(log2( (2^log2_NAb) %*% (pathogens[,3])/sum(pathogens[,3])))




x1 = scale(log2( (2^log2_NAb) %*% (pathogens[,1])/sum(pathogens[,1])))
x2 = scale(log2( (2^log2_NAb) %*% (pathogens[,2])/sum(pathogens[,2])))
x3 = scale(log2( (2^log2_NAb) %*% (pathogens[,3])/sum(pathogens[,3])))

cor(x1,x2)
cor(x1,x3)


