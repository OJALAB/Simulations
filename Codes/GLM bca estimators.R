set.seed(42)
N <- 100000
n_audit <- 1000
n_sim <- 50

biases_precise = biases_naive = biases_bca  = matrix(0, n_sim, 3)
colnames(biases_precise) = colnames(biases_naive) = colnames(biases_bca) = c('z', 'Intercept', 'x')

for(i in 1:n_sim){
  x <- rnorm(N, 1, 1)
  z <- rbinom(N, 1, 0.7)
  z_star <- z
  z_star[z == 1] <- sample(c(0,1), size = NROW(z_star[z == 1]), prob = c(0.07, 0.93), replace = T)
  z_star[z == 0] <- sample(c(0,1), size = NROW(z_star[z == 0]), prob = c(0.9, 0.1), replace = T)
  mu <- 1 + x - 2*z
  #y <- rbinom(N, 1, plogis(mu))
  y <- rpois(N, exp(mu))
  
  pop_data = data.frame(x, z, z_star, y)
  
  s_audit <- pop_data[sample(1:nrow(pop_data), n_audit),]
  fpr <- sum(s_audit$z==0 & s_audit$z_star==1)/nrow(s_audit)
  lambda01 = sqrt(N)*sum(s_audit$z == 0 & s_audit$z_star == 1)/sum(s_audit$z == 0)
  lambda10 = sqrt(N)*sum(s_audit$z == 1 & s_audit$z_star == 0)/sum(s_audit$z == 1)
  
  #precise <- glm(y ~ z + x, family = "binomial")
  precise <- glm(y ~ z + x, family = "poisson")
  biases_precise[i ,] <- precise$coefficients[c(2, 1, 3)] - c(-2, 1, 1)
  
  #naive <- glm(y ~ z_star + x, family = "binomial")
  naive <- glm(y~ z_star + x, family = "poisson")
  psi_hat <- naive$coefficients[c(2, 1, 3)]
  biases_naive[i, ] <- psi_hat - c(-2, 1, 1)
  
  
  ksi_hat <- cbind(z_star, rep(1, N), x)
  eta_hat <- psi_hat%*%t(ksi_hat)
  #p_hat <- plogis(t(psi_hat%*%t(ksi_hat)))
  #A <- t(sweep(ksi_hat, MARGIN = 1, STATS = p_hat*(1-p_hat), FUN = "*"))%*%ksi_hat/N
  A <- t(sweep(ksi_hat, MARGIN = 1, STATS = as.vector(exp(eta_hat)), FUN = "*"))%*%ksi_hat/N
  eta0 <- cbind(rep(1, n_audit), s_audit$x)%*%psi_hat[c(2, 3)]
  eta1 <- eta0 + rep(psi_hat[1], n_audit)
  #mu_diff <- plogis(eta0) - plogis(eta1)
  mu_diff <- exp(eta0)-exp(eta1)
  term1 <- (s_audit$z == 0)*lambda01*cbind(1, 1, s_audit$x)*as.vector(mu_diff)
  term2 <- -(s_audit$z == 1)*lambda10*cbind(0, 1, s_audit$x)*as.vector(mu_diff)
  b <- colMeans(term1 + term2)
  
  psi_bca <- psi_hat - 1/sqrt(N)*solve(A)%*%b
  biases_bca[i, ] <- psi_bca - c(-2, 1, 1)
  
}
biases <- rbind(colMeans(biases_precise), colMeans(biases_naive), colMeans(biases_bca))
sse_precise <- apply(biases_precise, 1, function(x) sum(x^2))
sse_naive <- apply(biases_naive, 1, function(x) sum(x^2))
sse_bca <- apply(biases_bca, 1, function(x) sum(x^2))

results <- as.data.frame(cbind(sse_precise, sse_naive, sse_bca))
biases
#save(results, biases, file = "Results/GLM bca estimators Poisson.RData")
#stargazer::stargazer(results, type='latex', title = "Results", flip=T)
#stargazer::stargazer(biases, type = 'latex', title = "Results", summary = FALSE)
