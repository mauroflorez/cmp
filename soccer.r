
#Soccer Application
#Scores English Premier League 2016 - 2022

#Libraries
library(readxl)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggdist)
library(coda)
library(bayesplot)
library(gridExtra)
library(ggrepel)
library(latex2exp)

#Load data
data <- read_excel("C:/Users/User/OneDrive - Rice University/Conway-Maxwell-Poisson/Data/BPL_2016_2022_covid.xlsx", 
                   col_types = c("skip", "skip", "text", "text", "numeric", "numeric", "text", "text"))

data$Covid = factor(data$Covid)
data$Season = factor(data$Season)
data$HomeTeam = factor(data$HomeTeam)
data$AwayTeam = factor(data$AwayTeam)

n = nrow(data)
y = cbind(data$HG, data$AG)
data$Covid <- relevel(data$Covid, ref = "Before")

#Design matrices
X1 <- model.matrix(HG ~ 1 + Covid + HomeTeam + AwayTeam, 
                   contrasts = list(Covid = "contr.treatment",
                                    HomeTeam = "contr.sum",
                                    AwayTeam = "contr.sum"),
                   data = data)

X2 <- model.matrix(AG ~ 1 + Covid + AwayTeam + HomeTeam, 
                   contrasts = list(Covid = "contr.treatment",
                                    HomeTeam = "contr.sum",
                                    AwayTeam = "contr.sum"),
                   data = data)

k1 = ncol(X1)
k2 = ncol(X2)
k_j = c(k1, k2)
X = list(X1,X2)

#Model
S = 20000
nburn = 10000
fit <- mcmc_cmp(y = y, X = X, intercept = FALSE,
                scale_cov_b = rep(1/100,2), v_0 = 2, R_0 = diag(0.1,2),
                S = S, nburn = nburn,
                prior_var_beta = list(1/10*diag(k_j[1]), 1/10*diag(k_j[2])),
                prior_var_gamma = list(1/10*diag(k_j[1]), 1/10*diag(k_j[2])),
                scale_cov_beta = c(0.05,0.05), scale_cov_gamma = c(0.005,0.005),
                update_jumping_rule = T, random_seed = 1)

#Acceptance Ratios
fit$accept_rate_beta
fit$accept_rate_gamma
fit$accept_rate_b
fit$accept_mixed

#Thinning
thin_MCMC <- 5
S <- 20000
idxs <- seq(0, S, by = thin_MCMC)

post_beta = fit$posterior_beta
post_gamma = fit$posterior_gamma

#Save posterior values
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

#Check Effective Sizes
effs_beta1 = effectiveSize(post_beta_1[idxs,])
effs_beta2 = effectiveSize(post_beta_2[idxs,])

effs_gamma1 = effectiveSize(post_gamma_1[idxs,])
effs_gamma2 = effectiveSize(post_gamma_2[idxs,])

#Traceplots (Size 15x10)

##Home Goals:
#Page 1
par(mfrow = c(5,5))
m <- matrix(c(1:25,26,26,26,26,26), nrow = 6, ncol = 5, byrow = T)
layout(mat = m, heights = c(1.5,1.5,1.5,1.5,1.4))

for(i in 1:25){
  plot.ts(post_gamma_1[idxs,i], ylim = c(min(c(post_gamma_1[idxs,i],post_beta_1[idxs,i])), max(c(post_gamma_1[idxs,i],post_beta_1[idxs,i]))),type = "l", ylab = paste("X_", i), col = "blue")
  lines(post_beta_1[idxs,i], col = "black")
}
plot.new()
plot_colors <- c("black", "blue")
legend(x = "top",inset = 0,
       legend = c("Beta", "Gamma"), 
       col=plot_colors, lwd=5, cex=0.8, horiz = TRUE)

#Page 2
m <- matrix(c(1:25,26,27,27,27,27), nrow = 6, ncol = 5, byrow = T)
layout(mat = m, heights = c(1,1,1,1,1,1,1))

for(i in 26:51){
  plot.ts(post_gamma_1[idxs,i], ylim = c(min(c(post_gamma_1[idxs,i],post_beta_1[idxs,i])), max(c(post_gamma_1[idxs,i],post_beta_1[idxs,i]))),type = "l", ylab = paste("X_", i), col = "blue")
  lines(post_beta_1[idxs,i], col = "black")
}
plot.new()
plot_colors <- c("black", "blue")
legend(x = "center",inset = 0,
       legend = c("Beta", "Gamma"), 
       col=plot_colors, lwd=5, cex=0.8, horiz = TRUE)
	   
##Away Goals
#Page 1
m <- matrix(c(1:25,26,26,26,26,26), nrow = 6, ncol = 5, byrow = T)
layout(mat = m, heights = c(1.5,1.5,1.5,1.5,1.4))

for(i in 1:25){
  plot.ts(post_gamma_2[idxs,i], ylim = c(min(c(post_gamma_2[idxs,i],post_beta_2[idxs,i])), max(c(post_gamma_2[idxs,i],post_beta_2[idxs,i]))),type = "l", ylab = paste("X_", i), col = "blue")
  lines(post_beta_2[idxs,i], col = "black")
}
plot.new()
plot_colors <- c("black", "blue")
legend(x = "top",inset = 0,
       legend = c("Beta", "Gamma"), 
       col=plot_colors, lwd=5, cex=0.8, horiz = TRUE)
#Page 2
m <- matrix(c(1:25,26,27,27,27,27), nrow = 6, ncol = 5, byrow = T)
layout(mat = m, heights = c(1,1,1,1,1,1,1))
for(i in 26:51){
  plot.ts(post_gamma_2[idxs,i], ylim = c(min(c(post_gamma_2[idxs,i],post_beta_2[idxs,i])), max(c(post_gamma_2[idxs,i],post_beta_2[idxs,i]))),type = "l", ylab = paste("X_", i), col = "blue")
  lines(post_beta_2[idxs,i], col = "black")
}
plot.new()
plot_colors <- c("black", "blue")
legend(x = "center",inset = 0,
       legend = c("Beta", "Gamma"), 
       col=plot_colors, lwd=5, cex=0.8, horiz = TRUE)

#Estimate Random Effects

post_b <- fit$posterior_b
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

#Expected Values
beta1_esp = apply(post_beta_1[idxs,], 2, mean)
beta2_esp = apply(post_beta_2[idxs,], 2, mean)

gamma1_esp = apply(post_gamma_1[idxs,], 2, mean)
gamma2_esp = apply(post_gamma_2[idxs,], 2, mean) 

b1_esp = apply(post_b1[,10000+idxs],1,mean) # check this!!!!!!!!!!!!!!!!!!!!!
b2_esp = apply(post_b2[,10000+idxs],1,mean)

mu_esp1 = exp(X[[1]]%*%beta1_esp + b1_esp)
nu_esp1 = exp(X[[1]]%*%gamma1_esp)

mu_esp2 = exp(X[[2]]%*%beta2_esp + b2_esp)
nu_esp2 = exp(X[[2]]%*%gamma2_esp)

mu_esp = cbind(mu_esp1, mu_esp2)
nu_esp = cbind(nu_esp1, nu_esp2)

##Goodness of Fit - y^rep

yrep_char1 <- matrix(nrow = S, ncol = n)
yrep_char2 <- matrix(nrow = S, ncol = n)

for(i in 1:S){
  yrep_char1[i,] = mapply(com_sampler, mu_esp1, nu_esp1)
  yrep_char2[i,] = mapply(com_sampler, mu_esp2, nu_esp2)
}

yrep_int1 <- sapply(data.frame(yrep_char1, stringsAsFactors = TRUE), as.integer)
yrep_int2 <- sapply(data.frame(yrep_char2, stringsAsFactors = TRUE), as.integer)

color_scheme_set("brightblue")
p1 <- ppc_rootogram(y[,1], yrep_int1) + theme(legend.position = "none")  + labs(x = "Home goals scored")
p2 <- ppc_rootogram(y[,2], yrep_int2) + theme(legend.position = "none")  + labs(x = "Away goals scored")

combined_plot <- grid.arrange(p1, p2, ncol = 2)

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

ggp1_legend <- ppc_rootogram(y[,1], yrep_int1) +
  theme(legend.position = "bottom")

legend <- extract_legend(ggp1_legend)

grid.arrange(combined_plot, legend, nrow = 2, heights = c(10,1))

## Offensive and defensive strenghts
#colnames(X1)

off_home <- apply(post_beta_1[idxs,], 2, mean)[4:27]
off_home <- c(off_home, 0 - sum(off_home))

def_away <- apply(post_beta_1[idxs,], 2, mean)[28:51]
def_away <- c(def_away, 0 - sum(def_away))

def_home <- apply(post_beta_2[idxs,], 2, mean)[28:51]
def_home <- c(def_home, 0 - sum(def_home))

off_away <- apply(post_beta_2[idxs,], 2, mean)[4:27]
off_away <- c(off_away, 0 - sum(off_away))

table_teams = data.frame("Home" = levels(factor(data$HomeTeam)), "Attack Home" = exp(off_home),
                         "Defense Home" = exp(-def_home), "Attack Away" = exp(off_away), "Defense Away" = exp(-def_away))
par(mfrow = c(1,2))
plot(Defense.Home ~ Attack.Home, data = table_teams, pch = 19)
abline(v = 1)
abline(h = 1)
text(Defense.Home ~ Attack.Home, labels=Home,data=table_teams, cex=0.8, font=1,pos=1)

plot(Defense.Away ~ Attack.Away, data = table_teams, pch = 19)
abline(v = 1)
abline(h = 1)
text(Defense.Away ~ Attack.Away, labels=Home,data=table_teams, cex=0.8, font=1,pos=1)

personal_theme =   theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
  
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(3, 7, 3)]

plot_1 <- ggplot(table_teams, aes(x = Attack.Home, y = Defense.Home)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point()  + geom_label_repel(aes(label = Home),
                            box.padding = 0.25,
                            point.padding = 0.5,
                            segment.color = "grey50") +
  labs(title = "Strengths playing Home",
       x = TeX("Attack $(\\exp(\\beta^\\omega_{H_i}))$"), y = TeX("Defense $(\\exp(-\\beta^\\delta_{H_i}))$")) + 
  personal_theme

plot_2 <- ggplot(table_teams, aes(x = Attack.Away, y = Defense.Away)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point()  + geom_label_repel(aes(label = Home),
                            box.padding = 0.25,
                            point.padding = 0.5,
                            segment.color = "grey50") +
  labs(title = "Strengths playing Away",
       x = TeX("Attack $(\\exp(\\beta^\\omega_{A_i}))$"), y = TeX("Defense $(\\exp(-\\beta^\\delta_{A_i}))$")) + 
  personal_theme

grid.arrange(plot_1, plot_2, ncol = 2)

#COVID - Home Advantage

soccer_home_adv = data.frame("Home.adv" = c(exp(post_beta_1[idxs,1]) - exp(post_beta_2[idxs,1]),
                                         exp(post_beta_1[idxs,1]+post_beta_1[idxs,2]) - exp(post_beta_2[idxs,1]+post_beta_2[idxs,2]),
                                         exp(post_beta_1[idxs,1]+post_beta_1[idxs,3]) - exp(post_beta_2[idxs,1]+post_beta_2[idxs,3])),
                          "Covid" = c(rep("Before",nrow(post_beta_1[idxs,])), 
                                      rep("After",nrow(post_beta_1[idxs,])),
                                      rep("During",nrow(post_beta_1[idxs,]))))

soccer_home_adv$Covid = factor(soccer_home_adv$Covid, levels = c("After", "During", "Before"), ordered = T)

g_rides <- 
  ggplot(mlb_home_adv, aes(x = Home.adv, y = Covid, fill = Covid, height = stat(density))) +
  coord_cartesian(clip = "off") + 
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none")+
  theme_classic()

g_rides + 
  geom_density_ridges(
    quantile_lines = TRUE, quantiles = 2, 
    color = "black",
    bandwidth = 0.05,
    alpha = 0.8, size = 0.5
  ) +
  labs(x = "Home Advantage", y = "Covid season")

## Aditional plots:
# Prediction of a game Liverpool vs Manchester City

home <- "Liverpool"
away <- "Manchester City"

summary(data$HomeTeam)
#Liverpool - 13
#City - 14

mu_exp_h <- exp(sum(beta1_esp[c(1,2,16,41)]))
nu_exp_h <- exp(sum(gamma1_esp[c(1,2,16,41)]))

mu_exp_a <- exp(sum(beta2_esp[c(1,2,17,40)]))
nu_exp_a <- exp(sum(gamma2_esp[c(1,2,17,40)]))

home_goals <- c()
away_goals <- c()

for(i in 1:1000){
  home_goals[i] <- com_sampler(mu_exp_h, nu_exp_h)
  away_goals[i] <- com_sampler(mu_exp_a, nu_exp_a)
}

outcomes <- data_frame(
  club = c(home, away),
  expected_goals = c(mean(home_goals), mean(away_goals)),
  prob_win = c(length(which(home_goals > away_goals)),
               length(which(away_goals > home_goals))) / length(home_goals),
  prob_tie = rep(length(which(home_goals == away_goals)) / length(home_goals),
                 2),
  prob_loss = c(length(which(away_goals > home_goals)),
                length(which(home_goals > away_goals))) / length(home_goals)
)

heatmap <- data_frame(home_goals, away_goals) %>%
  group_by(home_goals, away_goals) %>%
  summarize(probability = n() / nrow(.)) %>%
  mutate(plot = "Goal Distribution")

plot(heatmap)

histogram <- data_frame(home_goals, away_goals) %>%
  mutate(home_mov = home_goals - away_goals) %>%
  select(home_mov) %>%
  group_by(home_mov) %>%
  summarize(probability = n() / nrow(.)) %>%
  mutate(plot = "Margin of Victory", winner = ifelse(home_mov > 0, home,
                                                     ifelse(home_mov < 0, away, "Tie"))) %>%
  mutate(winner = factor(winner, levels = c(home, "Tie", away)))

score_dist <- ggplot(heatmap, aes(x = home_goals, y = away_goals)) +
  geom_tile(aes(fill = probability)) +
  scale_fill_gradient(name = "Probability") +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  scale_y_continuous(breaks = seq(0, 100, 1)) +
  labs(x = paste0(home, " Score"), y = paste0(away, " Score")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
score_dist

mov_dist <- ggplot(histogram, aes(x = home_mov, y = probability)) +
  geom_col(aes(fill = winner)) +
  scale_fill_brewer(type = "qual", palette = 3, name = "Winner") +
  scale_x_continuous(breaks = seq(-100, 100, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05), labels = scales::percent) +
  labs(x = paste0(home, " Margin of Victory"), y = "Probability") +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )
mov_dist




