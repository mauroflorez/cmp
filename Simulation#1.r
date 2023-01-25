# Simulation 1: Conway-Maxwell-Poisson: Comparison performance in 3 scenarios A, B, C

#Libraries

library(mvnfast)
library(coda)

#Generate data

set.seed(12)
J <- 2; k_j <- rep(3,J)
n <- 1000
X <- list()

y1 <- rmvn(n, mu = c(0,0), sigma = matrix(c(1,0.9,0.9,1), ncol = 2))

X[[1]] = cbind(rep(1,n), y1[,1], rpois(n, lambda = 1))
X[[2]] = cbind(rep(1,n), y1[,2], rpois(n, lambda = 1))

beta = c(0.2, 0.5, -0.8)
beta_obs <- list()
gamma_obs <- list()
for(j in 1:J){
  beta_obs[[j]] <- beta
  gamma_obs[[j]] <- c(0,0.2,0)
}
gamma_obs[[2]] <- c(0,0.2,0)

Xtbeta = matrix(unlist(Map("%*%", X, beta_obs)), ncol = J)
Xtgamma = matrix(unlist(Map("%*%", X, gamma_obs)), ncol = J)

mu_obs = exp(Xtbeta)
nu_obs = exp(Xtgamma)

y = matrix(mapply(com_sampler, c(mu_obs), c(nu_obs)), nrow = n, byrow = F)
cor(y)

#Plot Data

par(mfrow = c(1,2))
for(i in 1:J){
  barplot(table(y[,i])/n)
}

#Estimate Models - 3 Chains for each scenario

ptm <- proc.time()

fita_1 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.2, 0.15),
                     scale_cov_gamma = c(0.2, 0.2),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 1, scale_cov_b = rep(0.7, J),
                     prior_var_beta = rep(list(diag(10,3)),J),
                     prior_var_gamma = rep(list(diag(10,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(0.1, J), v_0 = 2)

fita_2 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.2, 0.15),
                     scale_cov_gamma = c(0.2, 0.2),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 2, scale_cov_b = rep(0.7, J),
                     prior_var_beta = rep(list(diag(10,3)),J),
                     prior_var_gamma = rep(list(diag(10,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(0.1, J), v_0 = 2)

fita_3 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.2, 0.15),
                     scale_cov_gamma = c(0.2, 0.2),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 3, scale_cov_b = rep(0.7, J),
                     prior_var_beta = rep(list(diag(10,3)),J),
                     prior_var_gamma = rep(list(diag(10,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(0.1, J), v_0 = 2)

fitb_1 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.4, 0.35),
                     scale_cov_gamma = c(0.4, 0.3),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 1, scale_cov_b = rep(0.5, J),
                     prior_var_beta = rep(list(diag(1,3)),J),
                     prior_var_gamma = rep(list(diag(1,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(1, J), v_0 = 10)

fitb_2 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.25, 0.25),
                     scale_cov_gamma = c(0.2, 0.2),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 2, scale_cov_b = rep(0.5, J),
                     prior_var_beta = rep(list(diag(1,3)),J),
                     prior_var_gamma = rep(list(diag(1,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(1, J), v_0 = 10)


fitb_3 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.2, 0.2),
                     scale_cov_gamma = c(0.2, 0.2),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 4, scale_cov_b = rep(0.5, J),
                     prior_var_beta = rep(list(diag(1,3)),J),
                     prior_var_gamma = rep(list(diag(1,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(1, J), v_0 = 10)

fitc_1 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.4, 0.4),
                     scale_cov_gamma = c(0.3, 0.3),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 1, scale_cov_b = rep(0.4, J),
                     prior_var_beta = rep(list(diag(0.1,3)),J),
                     prior_var_gamma = rep(list(diag(0.1,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(1, J), v_0 = 50)

fitc_2 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                     scale_cov_beta = c(0.4, 0.4),
                     scale_cov_gamma = c(0.3, 0.3),
                     S = 50000, nburn = 10000,
                     update_jumping_rule = T, random_seed = 2, scale_cov_b = rep(0.4, J),
                     prior_var_beta = rep(list(diag(0.1,3)),J),
                     prior_var_gamma = rep(list(diag(0.1,3)),J),
                     initial_beta = rep(list(rnorm(3)),J),
                     initial_gamma = rep(list(rnorm(3)),J),
                     #cov_h_beta = rep(list(diag(0.1,3)),J),
                     #cov_h_gamma = rep(list(diag(0.1,3)),J),
                     R_0 = diag(1, J), v_0 = 50)


fitc_3 <-mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                    scale_cov_beta = c(0.4, 0.4),
                    scale_cov_gamma = c(0.3, 0.3),
                    S = 50000, nburn = 10000,
                    update_jumping_rule = T, random_seed = 3, scale_cov_b = rep(0.4, J),
                    prior_var_beta = rep(list(diag(0.1,3)),J),
                    prior_var_gamma = rep(list(diag(0.1,3)),J),
                    initial_beta = rep(list(rnorm(3)),J),
                    initial_gamma = rep(list(rnorm(3)),J),
                    #cov_h_beta = rep(list(diag(0.1,3)),J),
                    #cov_h_gamma = rep(list(diag(0.1,3)),J),
                    R_0 = diag(1, J), v_0 = 50)



#Check acceptance rates

fitb_1$accept_rate_beta
fitb_1$accept_rate_gamma
fitb_1$accept_rate_b
fitb_1$accept_mixed

fita_2$accept_rate_beta
fitb_2$accept_rate_gamma
fitb_2$accept_rate_b
fitb_2$accept_mixed

fitb_3$accept_rate_beta
fitb_3$accept_rate_gamma
fitb_3$accept_rate_b
fitb_3$accept_mixed


#Save values

post1_beta = fitb_1$posterior_beta
post1_gamma = fitb_1$posterior_gamma
post1_b = fita_1$posterior_b

post2_beta = fitb_2$posterior_beta
post2_gamma = fitb_2$posterior_gamma
post2_b = fitb_2$posterior_b

post3_beta = fitb_3$posterior_beta
post3_gamma = fitb_3$posterior_gamma
post3_b = fitb_3$posterior_b

post1_beta_1 = matrix(nrow = length(post1_beta), ncol = k_j[1])
post1_beta_2 = matrix(nrow = length(post1_beta), ncol = k_j[2])

post2_beta_1 = matrix(nrow = length(post2_beta), ncol = k_j[1])
post2_beta_2 = matrix(nrow = length(post2_beta), ncol = k_j[2])

post3_beta_1 = matrix(nrow = length(post3_beta), ncol = k_j[1])
post3_beta_2 = matrix(nrow = length(post3_beta), ncol = k_j[2])

post1_gamma_1 = matrix(nrow = length(post1_gamma), ncol = k_j[1])
post1_gamma_2 = matrix(nrow = length(post1_gamma), ncol = k_j[2])

post2_gamma_1 = matrix(nrow = length(post2_gamma), ncol = k_j[1])
post2_gamma_2 = matrix(nrow = length(post2_gamma), ncol = k_j[2])

post3_gamma_1 = matrix(nrow = length(post3_gamma), ncol = k_j[1])
post3_gamma_2 = matrix(nrow = length(post3_gamma), ncol = k_j[2])

post1_b1 = matrix(nrow = n, ncol = length(post1_b))
post1_b2 = matrix(nrow = n, ncol = length(post1_b))

post2_b1 = matrix(nrow = n, ncol = length(post2_b))
post2_b2 = matrix(nrow = n, ncol = length(post2_b))

post3_b1 = matrix(nrow = n, ncol = length(post3_b))
post3_b2 = matrix(nrow = n, ncol = length(post3_b))

for(i in 1:length(post1_beta)){
  post1_beta_1[i,] = post1_beta[[i]][[1]]
  post2_beta_1[i,] = post2_beta[[i]][[1]]
  post3_beta_1[i,] = post3_beta[[i]][[1]]
  post1_b1[,i] = post1_b[[i]][,1]
  
  post1_beta_2[i,] = post1_beta[[i]][[2]]
  post2_beta_2[i,] = post2_beta[[i]][[2]]
  post3_beta_2[i,] = post3_beta[[i]][[2]]
  
  post1_gamma_1[i,] = post1_gamma[[i]][[1]]
  post2_gamma_1[i,] = post2_gamma[[i]][[1]]
  post3_gamma_1[i,] = post3_gamma[[i]][[1]]
  
  post1_gamma_2[i,] = post1_gamma[[i]][[2]]
  post2_gamma_2[i,] = post2_gamma[[i]][[2]]
  post3_gamma_2[i,] = post3_gamma[[i]][[2]]
  
  
  post1_b1[,i] = post1_b[[i]][,1]
  post2_b1[,i] = post2_b[[i]][,1]
  post3_b1[,i] = post3_b[[i]][,1]
  
  post1_b2[,i] = post1_b[[i]][,2]
  post2_b2[,i] = post2_b[[i]][,2]
  post3_b2[,i] = post3_b[[i]][,2]
}

thin_MCMC <- 5
S <- 50000
idxs <- seq(0, S, by = thin_MCMC)

post_beta_1 = mcmc.list(as.mcmc(post1_beta_1[idxs,]), as.mcmc(post2_beta_1[idxs,]), as.mcmc(post3_beta_1[idxs,]))
post_beta_2 = mcmc.list(as.mcmc(post1_beta_2[idxs,]), as.mcmc(post2_beta_2[idxs,]), as.mcmc(post3_beta_2[idxs,]))

post_gamma_1 = mcmc.list(as.mcmc(post1_gamma_1[idxs,]), as.mcmc(post2_gamma_1[idxs,]), as.mcmc(post3_gamma_1[idxs,]))
post_gamma_2 = mcmc.list(as.mcmc(post1_gamma_2[idxs,]), as.mcmc(post2_gamma_2[idxs,]), as.mcmc(post3_gamma_2[idxs,]))

post_b_1 = mcmc.list(as.mcmc(t(post1_b1)[idxs,]), as.mcmc(t(post2_b1)[idxs,]), as.mcmc(t(post3_b1)[idxs,]))
post_b_2 = mcmc.list(as.mcmc(t(post1_b2)[idxs,]), as.mcmc(t(post2_b2)[idxs,]), as.mcmc(t(post3_b2)[idxs,]))

#Rhat and Effective Sample Size

Rhat = rbind(
  cbind(gelman.diag(post_beta_1)$mpsrf,
        gelman.diag(post_beta_2)$mpsrf),
  cbind(gelman.diag(post_gamma_1)$mpsrf,
        gelman.diag(post_gamma_2)$mpsrf),
  cbind(gelman.diag(post_b_1)$mpsrf,
        gelman.diag(post_b_2)$mpsrf),
  cbind(max(gelman.diag(post_b_1)$psrf),
        max(gelman.diag(post_b_2)$psrf)))

colnames(Rhat) = c("y1","y2")
rownames(Rhat) = c("beta", "gamma", "b", "max_b")

ESS = round(rbind(
  cbind(mean(effectiveSize(post_beta_1))/3,
        mean(effectiveSize(post_beta_2))/3),
  cbind(mean(effectiveSize(post_gamma_1))/3,
        mean(effectiveSize(post_gamma_2))/3),
  cbind(mean(effectiveSize(post_b_1))/3,
        mean(effectiveSize(post_b_2))/3)),3)

colnames(ESS) = c("y1","y2")
rownames(ESS) = c("beta", "gamma", "b")

# HPD's

HPD_y1 = apply(simplify2array(HPDinterval(post_beta_1)),1:2, mean)
HPD_y2 = apply(simplify2array(HPDinterval(post_beta_2)),1:2, mean)
HPDs = round(cbind(HPD_y1, HPD_y2),2)

HPD_g_y1 = apply(simplify2array(HPDinterval(post_gamma_1)),1:2, mean)
HPD_g_y2 = apply(simplify2array(HPDinterval(post_gamma_2)),1:2, mean)
HPDs_g = round(cbind(HPD_g_y1, HPD_g_y2),2)

beta1_esp = apply(simplify2array(post_beta_1),2, mean)
beta2_esp = apply(simplify2array(post_beta_2),2, mean)

gamma1_esp = apply(simplify2array(post_gamma_1),2, mean)
gamma2_esp = apply(simplify2array(post_gamma_2),2, mean)

HPD_b1 = apply(simplify2array(HPDinterval(post_b_1)),2, mean)
HPD_b2 = apply(simplify2array(HPDinterval(post_b_2)),2, mean)
HPDs_b = round(rbind(HPD_b1, HPD_b2),2)

b1_esp = apply(simplify2array(post_b_1),2, mean)
b2_esp = apply(simplify2array(post_b_2),2, mean)

#Generate Table of results

estimations = matrix(nrow = 7, ncol = 2)
row.names(estimations) <- c("beta_1", "beta_2", "beta_3", "gamma_1", "gamma_2", "gamma_3", "mean_b")
colnames(estimations) <- c("y_1", "y_2")

for(j in 1:3){
  estimations[j,1] = paste(round(beta1_esp[j],2), " (",HPDs[j,1], ", ", HPDs[j,2],")", sep = "")
  estimations[j,2] = paste(round(beta2_esp[j],2), " (",HPDs[j,3], ", ", HPDs[j,4],")", sep = "")
  estimations[j+3,1] = paste(round(gamma1_esp[j],2), " (",HPDs_g[j,1], ", ", HPDs_g[j,2],")", sep = "")
  estimations[j+3,2] = paste(round(gamma2_esp[j],2), " (",HPDs_g[j,3], ", ", HPDs_g[j,4],")", sep = "")
  estimations[7,1] = paste(round(mean(b1_esp),6), " (", round(HPD_b1[1],2),", ", round(HPD_b1[2],2), ")", sep = "")
  estimations[7,2] = paste(round(mean(b2_esp),6), " (", round(HPD_b2[1],2),", ", round(HPD_b2[2],2), ")", sep = "")
}

mu1_esp = exp(X[[1]]%*%beta1_esp + b1_esp)
nu1_esp = exp(X[[1]]%*%gamma1_esp)

mu2_esp = exp(X[[2]]%*%beta2_esp + b2_esp)
nu2_esp = exp(X[[2]]%*%gamma2_esp)

mu_esp = cbind(mu1_esp, mu2_esp)
nu_esp = cbind(nu1_esp, nu2_esp)

#MSE table

mse = round(rbind(apply((mu_obs - mu_esp)^2, 2, mean),
                  apply((nu_obs - nu_esp)^2, 2, mean)),3)

colnames(mse) = c("y1", "y2")
rownames(mse) = c("mu", "nu")

round(Rhat,2)
round(ESS,2)
mse
estimations

# Traceplots

par(mfrow = c(2,3))
for(i in 1:3){
  ts.plot(post_beta_1[[1]][,i], type = "l", ylab = paste("Beta_", i))
  lines(post_beta_1[[2]][,i], col = "red")
  lines(post_beta_1[[3]][,i], col = "blue")
}

for(i in 1:3){
  ts.plot(post_beta_2[[1]][,i], type = "l", ylab = paste("Beta_", i))
  lines(post_beta_2[[2]][,i], col = "red")
  lines(post_beta_2[[3]][,i], col = "blue")
}

par(mfrow = c(2,3))
for(i in 1:3){
  ts.plot(post_gamma_1[[1]][,i], type = "l", ylab = paste("gamma_", i))
  lines(post_gamma_1[[2]][,i], col = "red")
  lines(post_gamma_1[[3]][,i], col = "blue")
}

for(i in 1:3){
  ts.plot(post_gamma_2[[1]][,i], type = "l", ylab = paste("gamma_", i))
  lines(post_gamma_2[[2]][,i], col = "red")
  lines(post_gamma_2[[3]][,i], col = "blue")
}


