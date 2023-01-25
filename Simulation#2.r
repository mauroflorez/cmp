# Simulation 2: Comparison Fit CMP vs Poisson vs NB

#Libraries
library(rstan)
library(rstanarm)
library(cowplot)
library(gridExtra)

###########################################################
# Case 1: Over-dispersed data

#Simulate data

set.seed(1)
n <- 1000
J <- 2

int <- rep(1,n)
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
x3 <- runif(n, -1, 1)

#Chanialidis (2018) Overdispersed data
ai <- sqrt((1-x3)/2)
x4 <- runif(n, -ai, ai)

mu1 = 0.3*x3 + 2*x4
mu2 = 0.5*x3 + 2*x4

y1 <- rpois(n, lambda = exp(mu1))
y2 <- rpois(n, lambda = exp(mu2))

mu_obs <- cbind(exp(mu1), exp(mu2))

y <- data.frame(y1, y2)
cor(y)

X <- list()
k_j <- rep(4,J)

for(j in 1:J){
  X[[j]] <- cbind(int, x1, x2, x3)
}

#Plot data
par(mfrow = c(1,J))
for(i in 1:J){
  barplot(table(y[,i])/n)
}


# Fit Poisson Models

pois1 <- stan_glm(y[,1] ~ X[[1]], family = poisson,
                  iter = 1000)

pois2 <- stan_glm(y[,2] ~ X[[2]], family = poisson,
                  iter = 1000)

#Y_rep
yrep_char1 <- posterior_predict(pois1, draws = 100)
yrep_int1 <- sapply(data.frame(yrep_char1, stringsAsFactors = TRUE), as.integer)

yrep_char2 <- posterior_predict(pois2, draws = 100)
yrep_int2 <- sapply(data.frame(yrep_char2, stringsAsFactors = TRUE), as.integer)

#Rootogram

#ppc_bars(y[,2], yrep_int, prob = 0.9) +  xlim(-1,25)
p1 <- ppc_rootogram(y[,1], yrep_int1, prob = 0.9, style = "hanging")
p2 <- ppc_rootogram(y[,2], yrep_int2, prob = 0.9, style = "hanging")

# Fit Negative Binomial Models
nb1 <- stan_glm(y[,1] ~ X[[1]], family = neg_binomial_2, iter = 1000)
nb2 <- stan_glm(y[,2] ~ X[[2]], family = neg_binomial_2, iter = 1000)

# Y_rep
yrep_char3 <- posterior_predict(nb1, draws = 100)
yrep_int3 <- sapply(data.frame(yrep_char3, stringsAsFactors = TRUE), as.integer)

yrep_char4 <- posterior_predict(nb2, draws = 100)
yrep_int4 <- sapply(data.frame(yrep_char4, stringsAsFactors = TRUE), as.integer)

#Rootogram
p3 <- ppc_rootogram(y[,1], yrep_int3, style = "hanging")
p4 <- ppc_rootogram(y[,2], yrep_int4, style = "hanging")


# Our Model - CMP

thin_MCMC <- 5
S <- 50000
n = 10000
idxs <- seq(0, S, by = thin_MCMC)

fit1 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                  scale_cov_gamma = c(0.3,0.3),
                  scale_cov_beta = c(0.3,0.3), 
                  S = S, nburn = n,
                  update_jumping_rule = T, random_seed = 1, scale_cov_b = rep(0.4, J),
                  prior_var_beta = rep(list(diag(0.1,4)),2),
                  prior_var_gamma = rep(list(diag(0.1,4)),2),
                  initial_beta = rep(list(rep(0,4)),3),
                  R_0 = diag(1, 2), v_0 = 10)

#Acceptance Rates
fit1$accept_rate_beta
fit1$accept_rate_gamma
fit1$accept_rate_b
fit1$accept_mixed

#Save Values
post_beta = fit1$posterior_beta
post_gamma = fit1$posterior_gamma

post_beta_1 = matrix(nrow = length(post_beta), ncol = k_j[1])
post_beta_2 = matrix(nrow = length(post_beta), ncol = k_j[2])

post_gamma_1 = matrix(nrow = length(post_gamma), ncol = k_j[1])
post_gamma_2 = matrix(nrow = length(post_gamma), ncol = k_j[2])

for(i in 1:length(post_beta)){
  post_beta_1[i,] = post_beta[[i]][[1]]
  post_beta_2[i,] = post_beta[[i]][[2]]
  
  post_gamma_1[i,] = post_gamma[[i]][[1]]
  post_gamma_2[i,] = post_gamma[[i]][[2]]
}

#ESS
effs_beta1 = effectiveSize(post_beta_1[idxs,])
effs_beta2 = effectiveSize(post_beta_2[idxs,])

effs_gamma1 = effectiveSize(post_gamma_1[idxs,])
effs_gamma2 = effectiveSize(post_gamma_2[idxs,])

#Traceplots
par(mfrow = c(4,4))
for(i in 1:4){
  plot(post_beta_1[idxs,i], type = "l", ylab = paste("Beta_", i), main = paste("mean =", round(mean(post_beta_1[,i]),2), ", eff_s = ", round(effs_beta1[i],2)))
}

for(i in 1:4){
  plot(post_beta_2[idxs,i], type = "l", ylab = paste("Beta_", i), main = paste("mean =", round(mean(post_beta_2[,i]),2), ", eff_s = ", round(effs_beta2[i],2)))
}

for(i in 1:4){
  plot(post_gamma_1[,i], type = "l", ylab = paste("gamma_", i), main = paste("mean =", round(mean(post_gamma_1[,i]),2), ", eff_s = ", round(effs_gamma1[i],2)))
}

for(i in 1:4){
  plot(post_gamma_2[,i], type = "l", ylab = paste("gamma_", i), main = paste("mean =", round(mean(post_gamma_2[,i]),2), ", eff_s = ", round(effs_gamma2[i],2)))
}

#Estimate Random Effects
post_b <- fit1$posterior_b
post_b1 = matrix(nrow = n, ncol = length(post_b))
post_b2 = matrix(nrow = n, ncol = length(post_b))

for(i in 1:length(post_b)){
  post_b1[,i] = post_b[[i]][,1]
  post_b2[,i] = post_b[[i]][,2]
}

# Mean of random effects for each variable
par(mfrow = c(1,2))

plot.ts(apply(post_b1,2,mean), type = "l", ylab = "b_1")
plot.ts(apply(post_b2,2,mean), type = "l", ylab = "b_2")

beta1_esp = apply(post_beta_1[idxs,], 2, mean)
beta2_esp = apply(post_beta_2[idxs,], 2, mean)

gamma1_esp = apply(post_gamma_1[idxs,], 2, mean)
gamma2_esp = apply(post_gamma_2[idxs,], 2, mean) 

#Expected Values
b1_esp = apply(post_b1[,idxs],1,mean)
b2_esp = apply(post_b2[,idxs],1,mean)

mu_esp1 = exp(X[[1]]%*%beta1_esp + b1_esp)
nu_esp1 = exp(X[[1]]%*%gamma1_esp)

mu_esp2 = exp(X[[2]]%*%beta2_esp + b2_esp)
nu_esp2 = exp(X[[2]]%*%gamma2_esp)

mu_esp = cbind(mu_esp1, mu_esp2)
nu_esp = cbind(nu_esp1, nu_esp2)

#Y_rep
yrep_char5 <- matrix(nrow = S, ncol = n)
yrep_char6 <- matrix(nrow = S, ncol = n)

for(i in 1:S){
  yrep_char5[i,] = mapply(com_sampler, mu_esp1, nu_esp1)
  yrep_char6[i,] = mapply(com_sampler, mu_esp2, nu_esp2)
}

yrep_int5 <- sapply(data.frame(yrep_char5, stringsAsFactors = TRUE), as.integer)
yrep_int6 <- sapply(data.frame(yrep_char6, stringsAsFactors = TRUE), as.integer)

p5 <- ppc_rootogram(y[,1], yrep_int5, style = "hanging")
p6 <- ppc_rootogram(y[,2], yrep_int6, style = "hanging")

#Plot rootogram
plot_grid(p1 + xlim(-1,10),p2 + xlim(-1,10),p3 + xlim(-1,10),
          p4 + xlim(-1,10),p5 + xlim(-1,10),p6 + xlim(-1,10),
          nrow = 3)

#Save Values

mse_1 <- mean((p1$data$tyexp-p1$data$ty)^2)
mse_2 <- mean((p2$data$tyexp-p2$data$ty)^2)
mse_3 <- mean((p3$data$tyexp-p3$data$ty)^2)
mse_4 <- mean((p4$data$tyexp-p4$data$ty)^2)
mse_5 <- mean((p5$data$tyexp-p5$data$ty)^2)
mse_6 <- mean((p6$data$tyexp-p6$data$ty)^2)

table_1 = cbind(rbind(mse_1,mse_3,mse_5),
      rbind(mse_2,mse_4,mse_6))

colnames(table_1) = c("y_1", "y_2")
rownames(table_1) = c("Poisson", "NB", "CMP")

###########################################################
# Case 2: Equi-dispersed data

mu1 = -0.5
mu2 = 0.2

y1 <- rpois(n, lambda = exp(mu1))
y2 <- rpois(n, lambda = exp(mu2))

mu_obs <- cbind(exp(mu1), exp(mu2))

y <- data.frame(y1, y2)
cor(y)

X <- list()
k_j <- rep(4,J)

for(j in 1:J){
  X[[j]] <- cbind(int, x1, x2, x3)
}

par(mfrow = c(1,J))
for(i in 1:J){
  barplot(table(y[,i])/n)
}

# Poisson Models
pois1 <- stan_glm(y[,1] ~ X[[1]], family = poisson,
                  iter = 1000)
pois2 <- stan_glm(y[,2] ~ X[[2]], family = poisson,
                  iter = 1000)

#Y_rep
yrep_char1 <- posterior_predict(pois1, draws = 100)
yrep_int1 <- sapply(data.frame(yrep_char1, stringsAsFactors = TRUE), as.integer)

yrep_char2 <- posterior_predict(pois2, draws = 100)
yrep_int2 <- sapply(data.frame(yrep_char2, stringsAsFactors = TRUE), as.integer)

#Rootogram
p1 <- ppc_rootogram(y[,1], yrep_int1, prob = 0.9, style = "hanging")
p2 <- ppc_rootogram(y[,2], yrep_int2, prob = 0.9, style = "hanging")

# Negative Binomial Models
nb1 <- stan_glm(y[,1] ~ X[[1]], family = neg_binomial_2, iter = 1000)
nb2 <- stan_glm(y[,2] ~ X[[2]], family = neg_binomial_2, iter = 1000)

# Y_rep

yrep_char3 <- posterior_predict(nb1, draws = 100)
yrep_int3 <- sapply(data.frame(yrep_char3, stringsAsFactors = TRUE), as.integer)

yrep_char4 <- posterior_predict(nb2, draws = 100)
yrep_int4 <- sapply(data.frame(yrep_char4, stringsAsFactors = TRUE), as.integer)

#Rootogram
p3 <- ppc_rootogram(y[,1], yrep_int3, style = "hanging")
p4 <- ppc_rootogram(y[,2], yrep_int4, style = "hanging")

plot_grid(p1 + xlim(-1,10),p2 + xlim(-1,10),
          p3 + xlim(-1,10),p4 + xlim(-1,10),
          nrow = 2)

# Our Model
fit2 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                  scale_cov_gamma = c(0.5,0.5),
                  scale_cov_beta = c(0.4,0.4), 
                  S = 25000, nburn = 10000,
                  update_jumping_rule = T, random_seed = 1, scale_cov_b = rep(1/10, J),
                  prior_var_beta = rep(list(diag(0.1,4)),2),
                  prior_var_gamma = rep(list(diag(0.1,4)),2),
                  initial_beta = rep(list(rep(0,4)),3),
                  R_0 = diag(1, 2), v_0 = 10)

fit2$accept_rate_beta
fit2$accept_rate_gamma
fit2$accept_rate_b
fit2$accept_mixed

post_beta = fit2$posterior_beta
post_gamma = fit2$posterior_gamma

post_beta_1 = matrix(nrow = length(post_beta), ncol = k_j[1])
post_beta_2 = matrix(nrow = length(post_beta), ncol = k_j[2])

post_gamma_1 = matrix(nrow = length(post_gamma), ncol = k_j[1])
post_gamma_2 = matrix(nrow = length(post_gamma), ncol = k_j[2])


for(i in 1:length(post_beta)){
  post_beta_1[i,] = post_beta[[i]][[1]]
  post_beta_2[i,] = post_beta[[i]][[2]]
  
  post_gamma_1[i,] = post_gamma[[i]][[1]]
  post_gamma_2[i,] = post_gamma[[i]][[2]]
}

effs_beta1 = effectiveSize(post_beta_1)
effs_beta2 = effectiveSize(post_beta_2)

effs_gamma1 = effectiveSize(post_gamma_1)
effs_gamma2 = effectiveSize(post_gamma_2)

par(mfrow = c(4,4))
for(i in 1:4){
  plot(post_beta_1[,i], type = "l", ylab = paste("Beta_", i), main = paste("mean =", round(mean(post_beta_1[,i]),2), ", eff_s = ", round(effs_beta1[i],2)))
}

for(i in 1:4){
  plot(post_beta_2[,i], type = "l", ylab = paste("Beta_", i), main = paste("mean =", round(mean(post_beta_2[,i]),2), ", eff_s = ", round(effs_beta2[i],2)))
}

for(i in 1:4){
  plot(post_gamma_1[,i], type = "l", ylab = paste("gamma_", i), main = paste("mean =", round(mean(post_gamma_1[,i]),2), ", eff_s = ", round(effs_gamma1[i],2)))
}

for(i in 1:4){
  plot(post_gamma_2[,i], type = "l", ylab = paste("gamma_", i), main = paste("mean =", round(mean(post_gamma_2[,i]),2), ", eff_s = ", round(effs_gamma2[i],2)))
}


#Estimate Random Effects
post_b <- fit2$posterior_b
post_b1 = matrix(nrow = n, ncol = length(post_b))
post_b2 = matrix(nrow = n, ncol = length(post_b))

for(i in 1:length(post_b)){
  post_b1[,i] = post_b[[i]][,1]
  post_b2[,i] = post_b[[i]][,2]
}

# Mean of random effects for each variable
par(mfrow = c(1,3))

plot.ts(apply(post_b1,2,mean), type = "l", ylab = "b_1")
plot.ts(apply(post_b2,2,mean), type = "l", ylab = "b_2")

beta1_esp = apply(post_beta_1, 2, mean)
beta2_esp = apply(post_beta_2, 2, mean)

gamma1_esp = apply(post_gamma_1, 2, mean)
gamma2_esp = apply(post_gamma_2, 2, mean) 

#Expected Values
b1_esp = apply(post_b1,1,mean)
b2_esp = apply(post_b2,1,mean)

mu_esp1 = exp(X[[1]]%*%beta1_esp + b1_esp)
nu_esp1 = exp(X[[1]]%*%gamma1_esp)

mu_esp2 = exp(X[[2]]%*%beta2_esp + b2_esp)
nu_esp2 = exp(X[[2]]%*%gamma2_esp)

mu_esp = cbind(mu_esp1, mu_esp2)
nu_esp = cbind(nu_esp1, nu_esp2)

#Plots Observed vs Fitted
yrep_char5 <- matrix(nrow = S, ncol = n)
yrep_char6 <- matrix(nrow = S, ncol = n)

for(i in 1:S){
  yrep_char5[i,] = mapply(com_sampler, mu_esp1, nu_esp1)
  yrep_char6[i,] = mapply(com_sampler, mu_esp2, nu_esp2)
}

yrep_int5 <- sapply(data.frame(yrep_char5, stringsAsFactors = TRUE), as.integer)
yrep_int6 <- sapply(data.frame(yrep_char6, stringsAsFactors = TRUE), as.integer)

p5 <- ppc_rootogram(y[,1], yrep_int5, style = "hanging")
p6 <- ppc_rootogram(y[,2], yrep_int6, style = "hanging")

plot_grid(p1 + xlim(-1,10),p2 + xlim(-1,10),p3 + xlim(-1,10),
          p4 + xlim(-1,10),p5 + xlim(-1,10),p6 + xlim(-1,10),
          nrow = 3)

#Table values

mse_1 <- mean((p1$data$tyexp-p1$data$ty)^2)
mse_2 <- mean((p2$data$tyexp-p2$data$ty)^2)
mse_3 <- mean((p3$data$tyexp-p3$data$ty)^2)
mse_4 <- mean((p4$data$tyexp-p4$data$ty)^2)
mse_5 <- mean((p5$data$tyexp-p5$data$ty)^2)
mse_6 <- mean((p6$data$tyexp-p6$data$ty)^2)

table_2 = cbind(rbind(mse_1,mse_3,mse_5),
                rbind(mse_2,mse_4,mse_6))

colnames(table_2) = c("y_1", "y_2")
rownames(table_2) = c("Poisson", "NB", "CMP")

###########################################################
# Case 2: Under-dispersed data: CMP

mu1 = -0.5 + x1
mu2 = 0.2 - x1

y1 <- mapply(com_sampler, exp(mu1) , rep(1.4,n))
y2 <- mapply(com_sampler, exp(mu2) , rep(1.8,n))

mu_obs <- cbind(exp(mu1), exp(mu2))

y <- data.frame(y1, y2)
cor(y)

X <- list()
k_j <- rep(4,J)

for(j in 1:J){
  X[[j]] <- cbind(int, x1, x2, x3)
}

par(mfrow = c(1,J))
for(i in 1:J){
  barplot(table(y[,i])/n)
}

# Poisson Models
pois1 <- stan_glm(y[,1] ~ X[[1]], family = poisson,
                  iter = 1000)

pois2 <- stan_glm(y[,2] ~ X[[2]], family = poisson,
                  iter = 1000)

#Y_rep
yrep_char1 <- posterior_predict(pois1, draws = 100)
yrep_int1 <- sapply(data.frame(yrep_char1, stringsAsFactors = TRUE), as.integer)

yrep_char2 <- posterior_predict(pois2, draws = 100)
yrep_int2 <- sapply(data.frame(yrep_char2, stringsAsFactors = TRUE), as.integer)

#Rootogram
p1 <- ppc_rootogram(y[,1], yrep_int1, prob = 0.9, style = "hanging") + theme(legend.position = "none") 
p2 <- ppc_rootogram(y[,2], yrep_int2, prob = 0.9, style = "hanging") + theme(legend.position = "none") 
 
# Negative Binomial Models
nb1 <- stan_glm(y[,1] ~ X[[1]], family = neg_binomial_2, iter = 1000)
nb2 <- stan_glm(y[,2] ~ X[[2]], family = neg_binomial_2, iter = 1000)

# Y_rep
yrep_char3 <- posterior_predict(nb1, draws = 100)
yrep_int3 <- sapply(data.frame(yrep_char3, stringsAsFactors = TRUE), as.integer)

yrep_char4 <- posterior_predict(nb2, draws = 100)
yrep_int4 <- sapply(data.frame(yrep_char4, stringsAsFactors = TRUE), as.integer)

#Rootogram
p3 <- ppc_rootogram(y[,1], yrep_int3, style = "hanging") + theme(legend.position = "none") 
p4 <- ppc_rootogram(y[,2], yrep_int4, style = "hanging") + theme(legend.position = "none") 

plot_grid(p1 + xlim(-1,10),p2 + xlim(-1,10),
          p3 + xlim(-1,10),p4 + xlim(-1,10),
          nrow = 2)


# Our Model
fit3 <- mcmc_cmp_2(y = y, X = X, way_beta = 2, way_gamma = 2, intercept = FALSE,
                  scale_cov_gamma = c(0.4,0.4),
                  scale_cov_beta = c(0.4,0.4), 
                  S = 20000, nburn = 10000,
                  update_jumping_rule = T, random_seed = 1, scale_cov_b = rep(0.5, J),
                  prior_var_beta = rep(list(diag(0.1,4)),2),
                  prior_var_gamma = rep(list(diag(0.1,4)),2),
                  initial_beta = rep(list(rep(0,4)),3),
                  R_0 = diag(1, 2), v_0 = 10)

fit3$accept_rate_beta
fit3$accept_rate_gamma
fit3$accept_rate_b
fit3$accept_mixed

post_beta = fit3$posterior_beta
post_gamma = fit3$posterior_gamma

post_beta_1 = matrix(nrow = length(post_beta), ncol = k_j[1])
post_beta_2 = matrix(nrow = length(post_beta), ncol = k_j[2])

post_gamma_1 = matrix(nrow = length(post_gamma), ncol = k_j[1])
post_gamma_2 = matrix(nrow = length(post_gamma), ncol = k_j[2])


for(i in 1:length(post_beta)){
  post_beta_1[i,] = post_beta[[i]][[1]]
  post_beta_2[i,] = post_beta[[i]][[2]]
  
  post_gamma_1[i,] = post_gamma[[i]][[1]]
  post_gamma_2[i,] = post_gamma[[i]][[2]]
}

effs_beta1 = effectiveSize(post_beta_1)
effs_beta2 = effectiveSize(post_beta_2)

effs_gamma1 = effectiveSize(post_gamma_1)
effs_gamma2 = effectiveSize(post_gamma_2)

par(mfrow = c(4,4))
for(i in 1:4){
  plot(post_beta_1[,i], type = "l", ylab = paste("Beta_", i), main = paste("mean =", round(mean(post_beta_1[,i]),2), ", eff_s = ", round(effs_beta1[i],2)))
}

for(i in 1:4){
  plot(post_beta_2[,i], type = "l", ylab = paste("Beta_", i), main = paste("mean =", round(mean(post_beta_2[,i]),2), ", eff_s = ", round(effs_beta2[i],2)))
}

for(i in 1:4){
  plot(post_gamma_1[,i], type = "l", ylab = paste("gamma_", i), main = paste("mean =", round(mean(post_gamma_1[,i]),2), ", eff_s = ", round(effs_gamma1[i],2)))
}

for(i in 1:4){
  plot(post_gamma_2[,i], type = "l", ylab = paste("gamma_", i), main = paste("mean =", round(mean(post_gamma_2[,i]),2), ", eff_s = ", round(effs_gamma2[i],2)))
}

#Estimate Random Effects

post_b <- fit3$posterior_b
post_b1 = matrix(nrow = n, ncol = length(post_b))
post_b2 = matrix(nrow = n, ncol = length(post_b))

for(i in 1:length(post_b)){
  post_b1[,i] = post_b[[i]][,1]
  post_b2[,i] = post_b[[i]][,2]
}

# Mean of random effects for each variable
par(mfrow = c(1,3))

plot.ts(apply(post_b1,2,mean), type = "l", ylab = "b_1")
plot.ts(apply(post_b2,2,mean), type = "l", ylab = "b_2")

beta1_esp = apply(post_beta_1, 2, mean)
beta2_esp = apply(post_beta_2, 2, mean)

gamma1_esp = apply(post_gamma_1, 2, mean)
gamma2_esp = apply(post_gamma_2, 2, mean) 

#Expected Values

b1_esp = apply(post_b1,1,mean)
b2_esp = apply(post_b2,1,mean)

mu_esp1 = exp(X[[1]]%*%beta1_esp + b1_esp)
nu_esp1 = exp(X[[1]]%*%gamma1_esp)

mu_esp2 = exp(X[[2]]%*%beta2_esp + b2_esp)
nu_esp2 = exp(X[[2]]%*%gamma2_esp)

mu_esp = cbind(mu_esp1, mu_esp2)
nu_esp = cbind(nu_esp1, nu_esp2)


yrep_char5 <- matrix(nrow = S, ncol = n)
yrep_char6 <- matrix(nrow = S, ncol = n)

for(i in 1:S){
  yrep_char5[i,] = mapply(com_sampler, mu_esp1, nu_esp1)
  yrep_char6[i,] = mapply(com_sampler, mu_esp2, nu_esp2)
}

yrep_int5 <- sapply(data.frame(yrep_char5, stringsAsFactors = TRUE), as.integer)
yrep_int6 <- sapply(data.frame(yrep_char6, stringsAsFactors = TRUE), as.integer)

p5 <- ppc_rootogram(y[,1], yrep_int5, style = "hanging") + theme(legend.position = "none") 
p6 <- ppc_rootogram(y[,2], yrep_int6, style = "hanging") + theme(legend.position = "none") 

color_scheme_set("brightblue")
plot_grid(p1 + xlim(-1,10),p2 + xlim(-1,10),p3 + xlim(-1,10),
          p4 + xlim(-1,10),p5 + xlim(-1,10),p6 + xlim(-1,10),
          nrow = 3)

#General Plot

combined_plot <- list(p1, p2, p3, p4, p5, p6)

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

ggp1_legend <- ppc_rootogram(y[,1], yrep_int1) +
  theme(legend.position = "bottom")

legend <- extract_legend(ggp1_legend)

pl <- grid.arrange(combined_plot, legend, nrow = 2, heights = c(10,1))

N <- length(combined_plot)
nr <- 3
nc <- 2

combine <- rbind(tableGrob(t(c(letters[1:nc])), theme = ttheme_minimal(), rows = ""), 
                 cbind(tableGrob(LETTERS[1:nr], theme = ttheme_minimal()), 
                       arrangeGrob(grobs = combined_plot),  size = "last"), size = "last")

pl <- replicate(12, ggplot(), FALSE)

N <- length(pl)
nr <- 3
nc <- 2

combine <- rbind(tableGrob(t(c(expression("y_1"), expression("y_2"))), theme = ttheme_minimal(), rows = ""), 
                 cbind(tableGrob(c("Poisson", "NB", "CMP"), theme = ttheme_minimal()), 
                       arrangeGrob(grobs = combined_plot),  size = "last"), size = "last")

grid.arrange(combine, legend, nrow = 2, heights = c(10,1))

#Save Values
mse_1 <- mean((p1$data$tyexp-p1$data$ty)^2)
mse_2 <- mean((p2$data$tyexp-p2$data$ty)^2)
mse_3 <- mean((p3$data$tyexp-p3$data$ty)^2)
mse_4 <- mean((p4$data$tyexp-p4$data$ty)^2)
mse_5 <- mean((p5$data$tyexp-p5$data$ty)^2)
mse_6 <- mean((p6$data$tyexp-p6$data$ty)^2)

table_3 = cbind(rbind(mse_1,mse_3,mse_5),
                rbind(mse_2,mse_4,mse_6))

colnames(table_3) = c("y_1", "y_2")
rownames(table_3) = c("Poisson", "NB", "CMP")

#Over:
table_1
#Equi:
table_2
#Under
table_3

