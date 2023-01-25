#' @param y Matrix of observations
#' @param X Covariates list, each element is the design matrix for each column of y
#' @param S Number of MCMC samples to be drawn
#' @param nburn Number of MCMC samples to burn-in
#' @param initial_beta Initial value of \eqn{beta}.
#' @param initial_gamma Initital value of \eqn{gamma}.
#' @param initial_b Initital value of \eqn{b}.

#' @param prior_mean_beta Prior mean for \eqn{beta}. (Default zero vector)
#' @param prior_var_beta Prior covariance matrix for \eqn{beta}
#' @param prior_mean_gamma Prior mean for \eqn{beta}. (Default zero vector)
#' @param prior_var_gamma Prior covariance matrix for \eqn{gamma}
#' @param v_0 Prior degrees of freedom of random effects
#' @param R_0 Prior covariance matrix of random effects

#' @param intercept Logical value indicating wheter include the intercept
#' @param random_seed Random seed 
#' @param way_beta Way of sampling beta (=1: combine two betas in one vector. =2: Each vector separated) 
#' @param way_gamma Way of sampling gamma (=1: In one vector. =2: Each vector separated)  

#' @param cov_h_beta Covariance matrix of proposal distribution for beta
#' @param cov_h_gamma Covariance matrix of proposal distribution for gamma
#' @param cov_h_b Covariance matrix of proposal distribution for random effects

#' @param scale_cov_beta Scale coefficient 
#' @param scale_cov_gamma Scale coefficient 
#' @param scale_cov_b Scale coefficient 
#' @param update_jumping_rule logical TRUE Adaptive MCMC

#' @param second_kind logical if second kind movement is included - Chaniadilis (2018)
#' @param inc_burn logical: include burned samples in the return

#' @return A list:
#' \item{posterior_beta} Values of beta
#' \item{posterior_gamma} Values of gamma
#' \item{posterior_b} Values of random effects b
#' \item{posterior_D} values of covariance matrix D
#' \item{accept_rate_b} Acceptance rate of \eqn{beta}
#' \item{accept_rate_beta} Acceptance rate of \eqn{beta}
#' \item{accept_rate_gamma} Acceptance rate of \eqn{gamma}
#' \item{accept_rate_mixed} Acceptance rate of second kind movement



mcmc_cmp <- function(y, X, S = 10000, nburn = 5000, initial_beta, initial_gamma, initial_b,
                     prior_mean_beta, prior_var_beta, prior_mean_gamma, prior_var_gamma, v_0, R_0,
					 intercept = TRUE, random_seed, way_beta = 2, way_gamma = 2,
                     cov_h_beta, cov_h_gamma, cov_h_b, scale_cov_beta, 
                     scale_cov_gamma, scale_cov_b, update_jumping_rule = FALSE, 
                     second_kind = TRUE, inc_burn = FALSE,...){
  
  
  design_matrix <- X
  #Include intercept to each design_matrix
  if(intercept == TRUE)
    design_matrix <- mapply(cbind, 1, X, SIMPLIFY = F)
  
  J <- ncol(y) #Number of variables
  n <- nrow(y) #Number of observations
  k_j <- sapply(design_matrix, dim)[2,] #Dimension of each design matrix
  
  # Default priors
  if(missing(v_0))  
    v_0 <- J
  if(missing(R_0))
    R_0 <- diag(J)
  if(missing(prior_mean_beta)){
    prior_mean_beta <- list()
    for(j in 1:J) prior_mean_beta[[j]] <- rep(0, k_j[j])
  }
  if(missing(prior_var_beta)){
    prior_var_beta <- list()
    for(j in 1:J) prior_var_beta[[j]] <- diag(k_j[j])
  }
  if(missing(prior_mean_gamma)){
    prior_mean_gamma <- list()
    for(j in 1:J) prior_mean_gamma[[j]] <- rep(0, k_j[j])
  }
  if(missing(prior_var_gamma)){
    prior_var_gamma <- list()
    for(j in 1:J) prior_var_gamma[[j]] <- diag(k_j[j])
  }
  
  
  # Initial Values
  if(missing(initial_b))
    initial_b <- matrix(c(rep(0, n*J)), nrow = n, ncol = J)
	
  # GLM models for default initial values of beta and gamma - Chaniadilis (2018)
  freq_glm <- list()
  for(j in 1:J){
   freq_glm[[j]] <- stats::glm(y[,j] ~ -1 + design_matrix[[j]], family = "quasipoisson")
  }
  
  if(missing(initial_beta)){
    initial_beta <- list()
    for(j in 1:J) initial_beta[[j]] <- unname(coef(freq_glm[[j]]))
  }
  if(missing(initial_gamma)){
    initial_gamma <- list()
    for(j in 1:J) initial_gamma[[j]] <- unname(coef(freq_glm[[j]]))
  }  
  
  # Scale constants
  if(missing(scale_cov_beta))
    scale_cov_beta = 2.4/sqrt(k_j)
  if(missing(scale_cov_gamma))
    scale_cov_gamma = 2.4/sqrt(k_j)  
  if(missing(scale_cov_b))
    scale_cov_b = rep(0.1, J)
	
  # Proposal distributions
  if(missing(cov_h_b))
    cov_h_b <- scale_cov_b^2*diag(1,J)
  if(missing(cov_h_beta)){
    cov_h_beta <- list()
    for(j in 1:J) cov_h_beta[[j]] <- scale_cov_beta[j]^2*unname(summary(freq_glm[[j]])$cov.scaled)
  }
  if(missing(cov_h_gamma)){
    cov_h_gamma <- list()
    for(j in 1:J) cov_h_gamma[[j]] <- scale_cov_gamma[j]^2*unname(summary(freq_glm[[j]])$cov.scaled)
  }
  
  # Prior densities to evaluate in acceptance ratios
  prior_b <- function(x, D) mvnfast::dmvn(x, mu = rep(0,J), sigma = D)
  prior_beta <- function(x){
    prior <- c()
    if(way_beta == 1){
      prior <- mvnfast::dmvn(unlist(x), mu = unlist(prior_mean_beta), 
                             sigma = as.matrix(Matrix::bdiag(prior_var_beta)))
    } else {
      for(j in 1:J){
        prior[j] <- mvnfast::dmvn(x[[j]], mu = prior_mean_beta[[j]], 
                                  sigma = prior_var_beta[[j]])
      }
    }
    return(prior)
  }
  prior_gamma <- function(x){
    prior <- c()
    if(way_gamma == 1){
      prior <- mvnfast::dmvn(unlist(x), mu = unlist(prior_mean_gamma),
                             sigma = as.matrix(Matrix::bdiag(prior_var_gamma)))
    } else {
      for(j in 1:J){
        prior[j] <- mvnfast::dmvn(x[[j]], mu = prior_mean_gamma[[j]], 
                                  sigma = prior_var_gamma[[j]])
      }
    }
    return(prior)
  }
  
  #Auxiliary Functions
  
  #Unnormalized density COM-Poisson
  log_dcmp <- function(y, mu, nu) nu*(y*log(mu) - lfactorial(y))
  
  # Dimensions for beta_j & gamma_j
  indices <- matrix(nrow = J, ncol = 2)
  indices[1,] <- c(1, k_j[1])
  if(J > 1)
    for(j in 2:J){
      indices[j,] <- c(cumsum(k_j)[j-1]+1, cumsum(k_j)[j])   
    }
	
  #Adaptive MCMC
  update_cov <- function(posterior){
    post_l <- list()
    for(j in 1:J){
      post_l[[j]] <- cov(matrix(unlist(sapply(posterior, "[", j)), ncol = k_j[j], byrow=T))
    }
    return(post_l)
  }
  #Product of list as matrices
  prod_list = function(X, beta){
    matrix(unlist(Map("%*%", X, beta)), ncol = J)
  }
  
  #Acceptance ratios
  accept_rule = function(y, mu_new, nu_new, mu, nu, log_prior, param){
    log_qf <- log_dcmp(y, mu, nu)
    log_qf_star <- log_dcmp(y, mu_new, nu_new)
    #Auxiliary data
    y_aux <- matrix(mapply(com_sampler, c(mu_new), c(nu_new)), nrow = n, byrow = F)
    log_qf_aux <- log_dcmp(y_aux, mu, nu)
    log_qf_aux_star <- log_dcmp(y_aux, mu_new, nu_new)
	
    if(param == "r_eff"){
      loglike <- log_prior + apply(log_qf_star - log_qf, 1, sum) + apply(log_qf_aux - log_qf_aux_star, 1, sum)
      accept <- log(stats::runif(nrow(y))) < loglike
    }
    if(param == "way_1"){
      loglike <- log_prior + sum(log_qf_star - log_qf) + sum(log_qf_aux - log_qf_aux_star)
      accept <- log(stats::runif(1)) < loglike
    }
    if(param == "way_2"){
      loglike <- log_prior + apply(log_qf_star - log_qf, 2, sum) + apply(log_qf_aux - log_qf_aux_star, 2, sum)
      accept <- log(stats::runif(J)) < loglike
    }
    return(accept)
  }
  
  #First values for D, b, beta and gamma
  if (!missing(random_seed))
    set.seed(random_seed)
  D_current <- solve(drop(rWishart(n = 1, df = v_0, Sigma = R_0)))
  b_current <- t(apply(initial_b, 1, function(x) drop(mvnfast::rmvn(1, mu = x, sigma = cov_h_b))))
  beta_current <- list()
  gamma_current <- list()
  for(j in 1:J){
    beta_current[[j]] <- drop(mvnfast::rmvn(1, mu = initial_beta[[j]], sigma = cov_h_beta[[j]]))
    gamma_current[[j]] <- drop(mvnfast::rmvn(1, mu = initial_gamma[[j]], sigma = cov_h_gamma[[j]]))
  }
  mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
  nu_current <- exp(prod_list(design_matrix, gamma_current))
  
  #Progress bar to check the progress of the MCMC
  pb <- utils::txtProgressBar(min=1, max=nburn+S, style=3)
  
  #Set values for the number of acceptances
  if(way_beta == 1) accept_beta <- 0 else accept_beta <- rep(0,J)
  if(way_gamma == 1) accept_gamma <- 0 else accept_gamma <- rep(0,J)
  accept_mixed <- rep(0,J)
  accept_b <- rep(0, n)
  
  #Define matrices to save the posterior draws into
  post_b <- list()
  post_beta <- list()
  post_gamma <- list()
  post_D <- list()
  
  ##############################################################################
  ############################### MCMC - Model #################################
  
  for(s in 1:(S+nburn)){
    
    # Adaptive MCMC
    if(s == nburn & update_jumping_rule == TRUE){
      cov_h_beta <- Map("*", update_cov(post_beta), scale_cov_beta^2)
      cov_h_gamma <- Map("*", update_cov(post_gamma), scale_cov_gamma^2)
    }
    
    ############################################################################
    #                            Block 1: Sample b                                              
    ############################################################################
    #Proposal
    b_star <- t(apply(b_current, 1, function(x) drop(mvnfast::rmvn(1, mu = x, sigma = cov_h_b))))
    mu_star <- exp(prod_list(design_matrix, beta_current) + b_star)
    nu_star <- nu_current
    #Priors Evaluation
    log_prior <- log(apply(b_star, 1, prior_b, D = D_current)) - log(apply(b_current, 1, prior_b, D = D_current))
    #Which i<n accepted?
    ind_accept <- which(accept_rule(y = y, mu_new = mu_star, nu_new = nu_star, mu = mu_current, nu = nu_current, log_prior = log_prior, param = "r_eff"))
    #Accept new values
    b_current[ind_accept,] <- b_star[ind_accept,]
    mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
    accept_b[ind_accept] <- accept_b[ind_accept] + 1
    
    ############################################################################
    #               Block 2: Sample Beta & Gamma: 2nd Kind                                              
    ############################################################################
    
    if(second_kind == TRUE){
      for(j in 1:J){
        beta_star <- beta_current[[j]]
        gamma_star <- gamma_current[[j]]
        for(k in 1:k_j[j]){
          #Proposals
          beta_star[k] <- drop(mvnfast::rmvn(1, mu = beta_current[[j]][k], sigma = cov_h_beta[[j]][k,k]))
          gamma_star[k] <- drop(mvnfast::rmvn(1, mu = gamma_current[[j]][k], sigma = cov_h_gamma[[j]][k,k]))
          mu_star <- exp(design_matrix[[j]]%*%beta_star + b_current[,j])
          nu_star <- exp(design_matrix[[j]]%*%gamma_star)
          #Log Prior
          log_prior <- log(mvnfast::dmvn(gamma_star, mu = prior_mean_gamma[[j]], sigma = prior_var_gamma[[j]])) -
            log(mvnfast::dmvn(gamma_current[[j]], mu = prior_mean_gamma[[j]], sigma = prior_var_gamma[[j]])) +
            log(mvnfast::dmvn(beta_star, mu = prior_mean_beta[[j]], sigma = prior_var_beta[[j]])) - 
            log(mvnfast::dmvn(beta_current[[j]], mu = prior_mean_beta[[j]], sigma = prior_var_beta[[j]]))
          if(accept_rule(y = y[,j], mu_new = mu_star, nu_new = nu_star, mu = mu_current[,j], nu = nu_current[,j], log_prior = log_prior, param = "way_1")){
            #Accept
            beta_current[[j]][k] <- beta_star[k]
            gamma_current[[j]][k] <- gamma_star[k]
            mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
            nu_current <- exp(prod_list(design_matrix, gamma_current))
            accept_mixed[j] <- accept_mixed[j] + 1/k_j[j]
          }
        }
      }
    }
    
    ############################################################################
    #                            Block 3: Sample Beta                                              
    ############################################################################
    
    # WAY 1:
    if(way_beta == 1){
      #Proposals
      beta_star_oneblock <- mvnfast::rmvn(1, mu = unlist(beta_current), sigma = as.matrix(Matrix::bdiag(cov_h_beta)))
      beta_star <- list()
      for(j in 1:J) beta_star[[j]] <- beta_star_oneblock[indices[j,1]:indices[j,2]]
      mu_star <- exp(prod_list(design_matrix, beta_star) + b_current)
      nu_star <- nu_current
      #Log Priors
      log_prior = log(prior_beta(beta_star)) - log(prior_beta(beta_current))
      #Decision Rules
      if(accept_rule(y, mu_new = mu_star, nu_new = nu_star, mu = mu_current, nu = nu_current, log_prior = log_prior, param = "way_1")){
        beta_current <- beta_star
        mu_current <- mu_star
        accept_beta <- accept_beta + 1
      }
    } else {
      #WAY 2:
      #Proposals
      beta_star <- list()
      for(j in 1:J) beta_star[[j]] <- drop(mvnfast::rmvn(1, mu = beta_current[[j]], sigma = cov_h_beta[[j]]))
      mu_star <- exp(prod_list(design_matrix, beta_star) + b_current)
      nu_star <- nu_current
      #Log Priors
      log_prior = log(prior_beta(beta_star)) - log(prior_beta(beta_current))
      #Which variable j<J accepted?
      var_accept <- which(accept_rule(y = y, mu_new = mu_star, nu_new = nu_star, mu = mu_current, nu = nu_current, log_prior = log_prior, param = "way_2"))
      #Accept
      beta_current[var_accept] <- beta_star[var_accept]
      mu_current <- exp(prod_list(design_matrix, beta_current) + b_current)
      accept_beta[var_accept] <- accept_beta[var_accept] + 1 
    }
    
    ############################################################################
    #                         Block 4: Sample Gamma                                                  
    ############################################################################
    
    # WAY 1:
    if(way_gamma == 1){
      #Proposals
      gamma_star_oneblock <- mvnfast::rmvn(1, mu = unlist(gamma_current), sigma = as.matrix(Matrix::bdiag(cov_h_gamma))) 
      gamma_star = list()
      for(j in 1:J) gamma_star[[j]] = gamma_star_oneblock[indices[j,1]:indices[j,2]]
      mu_star <- mu_current
      nu_star <- exp(prod_list(design_matrix, gamma_star))
      #Log Priors
      log_prior = log(prior_gamma(gamma_star)) - log(prior_gamma(gamma_current))
      #Decision Rules
      if(accept_rule(y, mu_new = mu_star, nu_new = nu_star, mu = mu_current, nu = nu_current, log_prior = log_prior, param = "way_1")){
        gamma_current <- gamma_star
        nu_current <- nu_star
        accept_gamma <- accept_gamma + 1
      }
    } else {
      # WAY 2:
      #Proposals
      gamma_star = list()
      for(j in 1:J) gamma_star[[j]] <- drop(mvnfast::rmvn(1, mu = gamma_current[[j]], sigma = cov_h_gamma[[j]]))
      mu_star <- mu_current
      nu_star <- exp(prod_list(design_matrix, gamma_star))
      #Log Priors
      log_prior = log(prior_gamma(gamma_star)) - log(prior_gamma(gamma_current))
      #Which variable j<J accepted?
      var_accept <- which(accept_rule(y = y, mu_new = mu_star, nu_new = nu_star, mu = mu_current, nu = nu_current, log_prior = log_prior, param = "way_2"))
      #Accept
      gamma_current[var_accept] <- gamma_star[var_accept]
      nu_current <- exp(prod_list(design_matrix, gamma_current))
      accept_gamma[var_accept] <- accept_gamma[var_accept] + 1 
    }
    
    ############################################################################
    #                         Block 5: Sample D                                              
    ############################################################################
    diag <- apply(b_current^2, 2, sum)
    non_diag <- sum(apply(b_current, 1, prod))
    
    scale <- solve(solve(R_0) + matrix(c(diag[1], non_diag, non_diag, diag[2]), nrow = 2))
    #Sampling
    D_current <- solve(drop(rWishart(n = 1, df = n + v_0, Sigma = scale)))
    
    
    #############################  Save Simulations ############################
    
    post_b[[s]] = b_current
    post_beta[[s]] = beta_current
    post_gamma[[s]] = gamma_current
    post_D[[s]] = D_current
    
    utils::setTxtProgressBar(pb, s)
  }
  close(pb)
  
  if(inc_burn == FALSE){
    list(posterior_b = post_b[-c(1:nburn)], posterior_beta = post_beta[-c(1:nburn)], posterior_gamma = post_gamma[-c(1:nburn)], posterior_D = post_D[-c(1:nburn)],
         accept_rate_b = accept_b/(nburn+S), accept_rate_beta = accept_beta/(nburn + S), accept_rate_gamma = accept_gamma/(nburn + S), accept_mixed = accept_mixed/(nburn+S))
  } else{
    list(posterior_b = post_b, posterior_beta = post_beta, posterior_gamma = post_gamma, posterior_D = post_D,
         accept_rate_b = accept_b/(nburn+S), accept_rate_beta = accept_beta/(nburn + S), accept_rate_gamma = accept_gamma/(nburn + S), accept_mixed = accept_mixed/(nburn+S))
    
  }
}

