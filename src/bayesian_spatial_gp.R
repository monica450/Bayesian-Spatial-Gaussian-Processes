set.seed(360655)
library(ggplot2)
library(tidyverse)

##PUNTO 1
# Carico i dati
load("C:/Users/data_for_homework.RData")
head(data_clean)
str(data_clean)
ggplot(data_clean, aes(x = X, y = Y, color = Y)) +
  geom_point(size = 3) +
  scale_color_distiller(palette = "Spectral") +
  labs(
    title = "Distribuzione spaziale del log(PM10)",
    x = "Coordinata X",
    y = "Coordinata Y",
    color = "log(PM10)"
  ) +
  theme_minimal()


#PUNTO 2

# Risposta (già logaritmo)
data_clean$y <- data_clean$log_PM10

#SUBSET E STRUTTURE DATI

n_sub <- 150
data_sub <- data_clean %>% slice_sample(n = n_sub)

coords <- data.frame(X = data_sub$X, Y = data_sub$Y)
y <- matrix(data_sub$y, ncol = 1)

# Design matrix: (mu, beta1, beta2)
Xmat <- cbind(1, coords$X, coords$Y)
n <- nrow(Xmat)
p <- ncol(Xmat)


#DISTANZE E PRIOR PER phi

dist_mat <- as.matrix(dist(coords))
dist_vec <- as.vector(dist_mat)

distmin <- min(dist_vec[dist_vec > 0])
distmax <- max(dist_vec)

phi_min <- 3 / distmax
phi_max <- 3 / distmin


#3) PRIOR DEBOLMENTE INFORMATIVE

#beta ~ N(0, 10^4 I)
b0 <- rep(0, p)
V0 <- diag(1e4, p)
V0_inv <- solve(V0)

#varianze ~ IG(1,1)
a_tau <- 1; b_tau <- 1
a_sig <- 1; b_sig <- 1


#FUNZIONI UTILI

R_phi <- function(phi, dist_mat) exp(-phi * dist_mat)

log_pw <- function(w, sigma2, phi, dist_mat) {
  R <- R_phi(phi, dist_mat)
  L <- chol(R)
  logdetR <- 2 * sum(log(diag(L)))
  v <- backsolve(L, forwardsolve(t(L), w))
  quad <- sum(w * v)
  -0.5 * (logdetR + quad / sigma2)
}

sample_w <- function(y, Xmat, beta, sigma2, tau2, phi, dist_mat) {
  Sigma <- sigma2 * exp(-phi * dist_mat)
  cov_post <- solve(solve(Sigma) + diag(n) / tau2)
  mean_post <- cov_post %*% ((y - Xmat %*% beta) / tau2)
  t(chol(cov_post)) %*% rnorm(n) + mean_post
}

sample_beta <- function(y, Xmat, w, tau2, b0, V0_inv) {
  Vb <- solve(t(Xmat) %*% Xmat / tau2 + V0_inv)
  mb <- Vb %*% (t(Xmat) %*% (y - w) / tau2)
  t(chol(Vb)) %*% rnorm(p) + mb
}

sample_tau2 <- function(y, Xmat, beta, w, a_tau, b_tau) {
  res <- y - Xmat %*% beta - w
  1 / rgamma(1, a_tau + n/2, b_tau + sum(res^2)/2)
}

sample_sigma2 <- function(w, phi, dist_mat, a_sig, b_sig) {
  R <- R_phi(phi, dist_mat)
  L <- chol(R)
  v <- backsolve(L, forwardsolve(t(L), w))
  quad <- sum(w * v)
  1 / rgamma(1, a_sig + n/2, b_sig + quad/2)
}


#METROPOLIS ADATTIVO PER phi

phi_from_u <- function(u) {
  phi_min + (phi_max - phi_min) * exp(u) / (1 + exp(u))
}

log_jacobian <- function(u) {
  log(phi_max - phi_min) + u - 2 * log1p(exp(u))
}

log_target_u <- function(u, w, sigma2) {
  phi <- phi_from_u(u)
  log_pw(w, sigma2, phi, dist_mat) + log_jacobian(u)
}


#MCMC: GIBBS + METROPOLIS

nsim <- 20000
burn <- 5000
thin <- 10

beta <- matrix(0, p, 1)
w <- matrix(0, n, 1)
tau2 <- 0.5
sigma2 <- 1
phi <- (phi_min + phi_max) / 2
u <- log((phi - phi_min)/(phi_max - phi))

sd_prop <- 0.3
target_acc <- 0.25
nbatch <- 100
A <- 500; B <- 1000

beta_save <- matrix(NA, nsim, p)
tau2_save <- sigma2_save <- phi_save <- rep(NA, nsim)

#SALVATAGGIO W
w_save <- list()
k <- 0
acc_sum <- 0; batch_count <- 0; acc_count <- 0

for (i in 1:nsim) {
  
  w <- sample_w(y, Xmat, beta, sigma2, tau2, phi, dist_mat)
  beta <- sample_beta(y, Xmat, w, tau2, b0, V0_inv)
  tau2 <- sample_tau2(y, Xmat, beta, w, a_tau, b_tau)
  sigma2 <- sample_sigma2(w, phi, dist_mat, a_sig, b_sig)
  
  u_star <- rnorm(1, u, sd_prop)
  loga <- log_target_u(u_star, w, sigma2) - log_target_u(u, w, sigma2)
  alpha <- min(1, exp(loga))
  
  acc_sum <- acc_sum + alpha
  batch_count <- batch_count + 1
  
  if (runif(1) < alpha) {
    u <- u_star
    acc_count <- acc_count + 1
  }
  
  phi <- phi_from_u(u)
  
  if (batch_count == nbatch) {
    acc_rate <- acc_sum / nbatch
    sd_prop <- exp(log(sd_prop) + (A/(B+i)) * (acc_rate - target_acc))
    acc_sum <- 0; batch_count <- 0
  }
  
  beta_save[i,] <- beta
  tau2_save[i] <- tau2
  sigma2_save[i] <- sigma2
  phi_save[i] <- phi
  
  if (i > burn && (i - burn) %% thin == 0) {
    k <- k + 1
    w_save[[k]] <- as.vector(w)
  }
}

Wmat <- do.call(cbind, w_save)

cat("Acceptance MH(phi):", acc_count / nsim, "\n")

# 5) Burn-in + thinning

keep <- seq(burn + 1, nsim, by = thin)

post_beta <- beta_save[keep, ]
post_tau2 <- tau2_save[keep]
post_sigma2 <- sigma2_save[keep]
post_phi <- phi_save[keep]

# Traceplot
par(mfrow=c(2,2))
plot(post_phi, type="l", main="phi (post)", xlab="iter", ylab="phi")
plot(post_sigma2, type="l", main="sigma2 (post)", xlab="iter", ylab="sigma2")
plot(post_tau2, type="l", main="tau2 (post)", xlab="iter", ylab="tau2")
plot(post_beta[,1], type="l", main="mu (post)", xlab="iter", ylab="mu")
par(mfrow=c(1,1))

# ACF
acf(post_phi, main="ACF phi")


##PUNTO 3


#TRACEPLOT DELLE CATENE

par(mfrow=c(3,2))
plot(post_beta[,1], type="l", main="mu", xlab="iter", ylab="")
plot(post_beta[,2], type="l", main="beta1 (coord X)", xlab="iter", ylab="")
plot(post_beta[,3], type="l", main="beta2 (coord Y)", xlab="iter", ylab="")
plot(post_tau2, type="l", main="tau2", xlab="iter", ylab="")
plot(post_sigma2, type="l", main="sigma2", xlab="iter", ylab="")
plot(post_phi, type="l", main="phi", xlab="iter", ylab="")
par(mfrow=c(1,1))

#ACF
par(mfrow=c(2,2))
acf(post_beta[,2], main="ACF beta1")
acf(post_beta[,3], main="ACF beta2")
acf(post_sigma2, main="ACF sigma2")
acf(post_phi, main="ACF phi")
par(mfrow=c(1,1))


#EFFETTO SIGNIFICATIVO DELLE COORDINATE

beta_ci <- apply(post_beta, 2, quantile, probs=c(0.025, 0.5, 0.975))
colnames(beta_ci) <- c("mu", "beta1", "beta2")
beta_ci

beta_summary <- data.frame(
  parametro = c("mu", "beta1", "beta2"),
  mean = colMeans(post_beta),
  sd = apply(post_beta, 2, sd),
  q025 = beta_ci[1,],
  q50  = beta_ci[2,],
  q975 = beta_ci[3,]
)

print(beta_summary)


#C(h; phi): MEDIA E CI 95% DELLA COVARIANZA SPAZIALE

#Griglia di distanze
h_grid <- seq(0, distmax, length.out = 200)

# C(h) = sigma^2 * exp(-phi * h)
C_post <- sapply(h_grid, function(h)
  post_sigma2 * exp(-post_phi * h)
)

C_df <- data.frame(
  h = h_grid,
  mean = apply(C_post, 2, mean),
  q025 = apply(C_post, 2, quantile, 0.025),
  q975 = apply(C_post, 2, quantile, 0.975)
)

ggplot(C_df, aes(x=h)) +
  geom_ribbon(aes(ymin=q025, ymax=q975), alpha=0.25) +
  geom_line(aes(y=mean), linewidth=1) +
  labs(title="C(h; phi): media e intervallo di credibilità 95%",
       x="Distanza h", y="Covarianza") +
  theme_minimal()


#RAPPRESENTAZIONE DI W(s): MEDIA E CI 95%

#Wmat: matrice n x K (ottenuta dal punto 2)
#Calcolo M, Q1, Q2
M  <- rowMeans(Wmat)
Q1 <- apply(Wmat, 1, quantile, 0.025)
Q2 <- apply(Wmat, 1, quantile, 0.975)

data_w <- data.frame(
  X = coords$X,
  Y = coords$Y,
  Yobs = as.vector(y),
  M = M,
  Q1 = Q1,
  Q2 = Q2
)


#MAPPA DELLA MEDIA A POSTERIORI DI W(s)

ggplot(data_w, aes(x=X, y=Y, color=M)) +
  geom_point(size=3) +
  scale_color_distiller(palette="Spectral") +
  labs(title="W(s): media a posteriori",
       x="X", y="Y", color="M") +
  theme_minimal()


#MAPPA DELL’INCERTEZZA (AMPIEZZA CI 95%)

data_w$CI_width <- data_w$Q2 - data_w$Q1

ggplot(data_w, aes(x=X, y=Y, color=CI_width)) +
  geom_point(size=3) +
  scale_color_distiller(palette="Spectral") +
  labs(title="W(s): ampiezza intervallo di credibilità 95%",
       x="X", y="Y", color="Q2 - Q1") +
  theme_minimal()


#CONFRONTO TRA Y(s) OSSERVATO E W(s)

# Mappa di Y(s)
ggplot(data_w, aes(x=X, y=Y, color=Yobs)) +
  geom_point(size=3) +
  scale_color_distiller(palette="Spectral") +
  labs(title="Y(s) = log_PM10 osservato",
       x="X", y="Y", color="Y(s)") +
  theme_minimal()

#PUNTO 4

# Coordinate
coords_mat <- as.matrix(coords[, c("X", "Y")])
coords_grid <- as.matrix(griglia[, c("X", "Y")])

n_obs  <- nrow(coords_mat)
n_grid <- nrow(coords_grid)
n_mcmc <- length(w_save)


#Distanze osservati–griglia

D_obs_grid <- as.matrix(dist(rbind(coords_mat, coords_grid)))
D_obs_grid <- D_obs_grid[1:n_obs, (n_obs + 1):(n_obs + n_grid)]


# Storage campioni posteriori di W(s0)

W_grid <- matrix(NA, n_mcmc, n_grid)


#MCMC

for (b in 1:n_mcmc) {
  
  phi_b  <- phi_save[b]
  sig2_b <- sigma2_save[b]
  w_b <- matrix(w_save[[b]], ncol = 1)
  
  # Correlazioni
  R_obs <- exp(-phi_b * as.matrix(dist(coords_mat)))
  R_og  <- exp(-phi_b * D_obs_grid)
  
  Sigma_obs <- sig2_b * R_obs
  Sigma_og  <- sig2_b * R_og
  
  # Cholesky per stabilità numerica 
  L <- chol(Sigma_obs)
  
  # Sigma_obs^{-1} w
  a <- backsolve(L, forwardsolve(t(L), w_b))
  
  # Media condizionata
  mean_pred <- t(Sigma_og) %*% a
  
  # Salvo solo la media
  W_grid[b, ] <- as.vector(mean_pred)
}


# Riassunti posteriori

M_grid  <- apply(W_grid, 2, mean)
Q1_grid <- apply(W_grid, 2, quantile, 0.025)
Q2_grid <- apply(W_grid, 2, quantile, 0.975)

df_grid_plot <- data.frame(
  X = coords_grid[,1],
  Y = coords_grid[,2],
  M = M_grid,
  Q1 = Q1_grid,
  Q2 = Q2_grid
)


# GRAFICI (geom_tile)

ggplot(df_grid_plot, aes(x = X, y = Y, fill = M)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(title = "Media a posteriori di W(s) sulla griglia",
       fill = "E[W(s)|y]") +
  theme_minimal()

ggplot(df_grid_plot, aes(x = X, y = Y, fill = Q1)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(title = "Limite inferiore CI 95% di W(s)",
       fill = "Q1") +
  theme_minimal()

ggplot(df_grid_plot, aes(x = X, y = Y, fill = Q2)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_equal() +
  labs(title = "Limite superiore CI 95% di W(s)",
       fill = "Q2") +
  theme_minimal()

