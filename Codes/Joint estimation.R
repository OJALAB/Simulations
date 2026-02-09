library(optimx)

l_poiss = function(y, d, gamma, intercept, alpha, q, omega, p){#y - number, d - number from 1 to 10, gamma - vector of length 9, alpha- number
  ref <- alpha*q+intercept #q - 0 or 1, omega - matrix 10x10, p - vector of length 10
  result <- omega[d, 10]*p[10]*dpois(y, lambda = exp(ref))
  for(i in 1:9){
    logmi <- c(gamma %*% diag(9)[, i])+ref
    result <- result+omega[d, i]*p[i]*dpois(y, lambda = exp(logmi))
  }
  penalty <- sum(c(gamma, intercept, alpha)^2)
  return(sum(log(result))-0.5*penalty)
}

opt <- function(vector, y, d, q, omega, p){
  return(-l_poiss(y=y, d=d, gamma=vector[1:9], intercept=vector[10], alpha=vector[11], q=q, omega=omega, p=p))
}

opt2 <- function(vector, y, d, q, omega){
  gamma <- vector[1:9]
  intercept <- vector[10]
  alpha <- vector[11]
  p <- vector[12:20]
  p <- exp(c(p, 0))
  p <- p/sum(p)
  return(-l_poiss(y, d, gamma, intercept, alpha, q, omega, p))
}

sim = function(n, gamma, alpha, p, accuracy = 70, prob_q = 0.7){
  gamma <- runif(10, 1, 2)
  alpha <- runif(1, -0.5, 0.5)
  p <- c(0.25, 0.2, 0.15, 0.05, 0.05, 0.1, 0.02, 0.03, 0.05, 0.1)
  misclass <- diag(10)*accuracy + matrix(rpois(100, (100-accuracy)/10), nrow = 10)
  misclass  <-  sweep(misclass, MARGIN = 2, STATS = colSums(misclass), FUN = '/')
  q <- rbinom(n, 1, prob_q)
  theta <- sample(1:10, size = n, replace = T, prob = p)
  y <- rpois(n, exp(gamma[theta] + alpha*q))
  d <- sapply(1:n, function(x){sample(1:10, size = 1, prob = misclass[, theta[x]])})
  
  fit <- optimx(par=c(numeric(9), log(mean(y)), 0), fn=opt, y=y, d=d, q=q, omega=misclass, p=p, method = "BFGS")
  results <- as.numeric(fit[1, 1:11])
  results[1:9] <- results[1:9]+results[10]
  
  fit2 <- optimx(par=c(numeric(9), log(mean(y)), numeric(10)), fn=opt2, y=y, d=d, q=q, omega=misclass, method = "BFGS")
  results2 <- as.numeric(fit2[1, 1:11])
  results2[1:9] <- results2[1:9]+results2[10]
  
  naive <- glm(y ~ as.factor(d) + q, family = "poisson")
  glm_results <- coef(naive)
  glm_results[2:10] <- glm_results[2:10]+glm_results[1]
  
  return(list(results, results2, glm_results))
}

errors_1000 <- matrix(0, 25, 3)
#colnames(errors_1000) <- c('Joint', 'Joint with unknown probabilities', 'Naive')
for(i in 1:25){
  gamma <- runif(10, 1, 2)
  alpha <- runif(1, -0.5, 0.5)
  p <- c(0.25, 0.2, 0.15, 0.05, 0.05, 0.1, 0.02, 0.03, 0.05, 0.1)
  sim_results=sim(1000, gamma, alpha, p)
  errors_1000[i, 1] <- sum((c(gamma, alpha)-sim_results[[1]])^2)
  errors_1000[i, 2] <- sum((c(gamma, alpha)-sim_results[[2]])^2)
  errors_1000[i, 3] <- sum((c(gamma, alpha)-sim_results[[3]])^2)
}

errors_5000 <- matrix(0, 25, 3)
for(i in 1:25){
  gamma <- runif(10, 1, 3)
  alpha <- runif(1, -0.5, 0.5)
  p <- gtools::rdirichlet(1, alpha = rep(3, 10))
  accuracy <- runif(1, 60, 90)
  prob_q <- runif(1, 0.4, 0.9)
  sim_results <- sim(1000, gamma, alpha, p, accuracy, prob_q)
  errors_5000[i, 1] <- sum((c(gamma, alpha)-sim_results[[1]])^2)
  errors_5000[i, 2] <- sum((c(gamma, alpha)-sim_results[[2]])^2)
  errors_5000[i, 3] <- sum((c(gamma, alpha)-sim_results[[3]])^2)
}

results <- as.data.frame(cbind(errors_1000, errors_5000))
results <- rbind(results, colMeans(results))
rownames(results)[26] <- 'Mean Error'
colnames(results) <- c('Joint estimation 1000 observations',
                  'Joint esitmation 1000 observations with unknown probabilities',
                  'Standard GLM 1000 observations',
                  'Joint estimation 5000 observations',
                  'Joint esitmation 5000 observations with unknown probabilities',
                  'Standard GLM 5000 observations')
save(results, file = "Results/Joint estimation.RData")

#stargazer::stargazer(as.data.frame(wynik), summary.stat = c("mean", "median", "sd", "min", "max"), type='latex', title = "Results", digits=4, flip=T)
