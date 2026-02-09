library(MLBC)
set.seed(2025)
N <- 10000
n_audit <- 1000
n_sim <- 500
error_precise = error_naive = error_bca_naive = error_bca = matrix(0, n_sim, 3)
colnames(error_precise) = colnames(error_naive) = colnames(error_bca_naive) = colnames(error_bca) = c('z', 'Intercept', 'x')
for(i in 1:n_sim){
  x <- rnorm(N, 2, 1)
  z <- rbinom(N, 1, 0.7)
  z_star <- z
  z_star[z == 1] <- sample(c(0,1), size = NROW(z_star[z == 1]), prob = c(0.07, 0.93), replace = T)
  z_star[z == 0] <- sample(c(0,1), size = NROW(z_star[z == 0]), prob = c(0.93, 0.07), replace = T)
  e <- rnorm(N)
  y <- 1 + x - 2*z + e
  pop_data <- data.frame(x, z, z_star, y)
  
  s_audit <- pop_data[sample(1:nrow(pop_data), n_audit),]
  
  fpr <- sum(s_audit$z==0 & s_audit$z_star==1)/nrow(s_audit)
  errors <- s_audit$z_star - s_audit$z
  lambda <- c(mean(errors), mean(s_audit$x*errors))
  
  #precise estimator - uses data without measurement error
  ksi <- rbind(z, rep(1, N), x)
  psi_hat_precise <- solve(1/N*ksi%*%t(ksi))%*%(1/N*ksi%*%y)
  error_precise[i,] <- psi_hat_precise - c(-2, 1, 1)
  
  #naive estimator - does not correct measurement error
  ksi_hat <- rbind(z_star, rep(1, N), x)
  psi_hat_naive <- solve(1/N*ksi_hat%*%t(ksi_hat))%*%(1/N*ksi_hat%*%y)
  error_naive[i, ] <- psi_hat_naive -c(-2, 1, 1)
  
  #naive bca estimator - corrects measurement error even though assumptions are violated
  psi_bca_naive <- ols_bca(y ~ z_star + x, data = pop_data, fpr = fpr, m = n_audit)
  error_bca_naive[i, ] <- psi_bca_naive$coef[c(1, 3, 2)] - c(-2, 1, 1)
  
  #upgraded bca estimator - corrects measurement error
  psi_bca <- (diag(3)+solve(1/N*ksi_hat%*%t(ksi_hat))%*%(cbind(c(fpr, lambda), matrix(0, 3, 2))))%*%psi_hat_naive
  error_bca[i, ] <- psi_bca-c(-2, 1, 1)
}
sse_bca <- apply(error_bca, 1, function(x) sum(x^2))
sse_bca_naive <- apply(error_bca_naive, 1, function(x) sum(x^2))
sse_naive <- apply(error_naive, 1, function(x) sum(x^2))
sse_precise <- apply(error_precise, 1, function(x) sum(x^2))

results <- as.data.frame(cbind(sse_precise, sse_naive, sse_bca_naive, sse_bca))
save(results, error_bca, error_bca_naive, error_naive, error_precise, file = "Results/Upgraded bca esimator.RData")
#stargazer::stargazer(results, type='latex', title = "Results", flip=T)
