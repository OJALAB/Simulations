sim=function(n){
  error_bca <- numeric(50)
  error_bcm <- numeric(50)
  error_naive <- numeric(50)
  for(i in 1:50){
  gamma <- runif(10, 1, 2)
  alfa <- runif(1, -0.5, 0.5)
  p <- c(0.25, 0.2, 0.15, 0.05, 0.05, 0.1, 0.02, 0.03, 0.05, 0.1)
  misclass <- diag(10)*70+matrix(rpois(100, 3), nrow=10)
  misclass = sweep(misclass, MARGIN = 2, STATS = colSums(misclass), FUN = '/')
  q <- rbinom(n, 1, 0.7)
  theta <- sample(1:10, size=n, replace = T, prob = p)
  y <- rpois(n, exp(gamma[theta]+alfa*q))
  d <- sapply(1:n, function(x){sample(1:10, size = 1, prob = misclass[, theta[x]])})
  
  fpr <- sum(sapply(2:10, function(x){sum(misclass[x,]*p)-misclass[x, x]*p[x]}))
  omega_hat <- (diag(sapply(2:10, function(x){sum(misclass[x,]*p)}))-sweep(misclass[2:10, 2:10], MARGIN = 2, STATS = p[2:10], FUN='*'))/fpr
  omega_hat = as.matrix(Matrix::bdiag(matrix(0, 1, 1), omega_hat, matrix(0, 1, 1)))
  
  glm_naive <- glm(y ~ as.factor(d) + q, family = "poisson")
  glm_results <- coef(glm_naive)
  psi_hat <- glm_results
  d_factor <- relevel(factor(d, levels=as.character(1:10)), ref='1')
  X <- model.matrix(~d_factor +q)
  m <- solve(t(X)%*%X/n)
  
  psi_bca <- (diag(11)+fpr*m%*%omega_hat)%*%psi_hat
  psi_bca[2:10] <- psi_bca[2:10]+psi_bca[1]
  error_bca[i] <- sum((psi_bca-c(gamma, alfa))^2)
  
  psi_bcm <- solve(diag(11)-fpr*m%*%omega_hat)%*%psi_hat
  psi_bcm[2:10] <- psi_bcm[2:10]+psi_bcm[1]
  error_bcm[i] <- sum((psi_bcm-c(gamma, alfa))^2)
  
  glm_results[2:10] <- glm_results[2:10]+glm_results[1]
  error_naive[i] <- sum((glm_results-c(gamma, alfa))^2)
  }
  return(list(error_bca, error_bcm, error_naive))
}

results_1000=sim(1000)
results_10000=sim(10000)
results_100000=sim(100000)
results=as.data.frame(c(results_1000, results_10000, results_100000))
colnames(results)=c('bca1000', 'bcm1000', 'naive1000', 'bca10000', 'bcm10000', 'naive10000', 'bca100000', 'bcm100000', 'naive100000')
results=rbind(results, colMeans(results))
rownames(results)[51]='Mean Error'

save(results, file = "Results/Bias corrected estimators results.RData")
#stargazer::stargazer(wynik_bc, type='latex', title = "Results", label = "tab:summary", flip=T)
