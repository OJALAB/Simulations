pred = function(coefs, sample){
  return(coefs[1] + coefs[2]*sample[['x']] + coefs[3]*sample[['z']])
}

set.seed(42)
N <- 100000
n_sim <- 1000
n_prob <- 500
n_nonprob <- 1000
n_audit <- 1000
nonprob_ratio <- 0.7
biases <- matrix(0, n_sim, 8)
colnames(biases) <- c('gold', 'mass_imp', 'mass_impnaive', 'bca', 'bca_upgr', 'bcm', 'bcm_upgr', 'joint')

for(i in 1:n_sim){
  x <- rnorm(N, 2, 1)
  z <- rbinom(N, 1, 0.7)
  z_star <- z
  #create measurement error
  z_star[z == 1] <- sample(c(0,1), size = NROW(z_star[z == 1]), prob = c(0.1, 0.9), replace = T)
  z_star[z == 0] <- sample(c(0,1), size = NROW(z_star[z == 0]), prob = c(0.93, 0.07), replace = T)
  e <- rnorm(N)
  y <- 1 + x - 2*z + e
  y_mean <- mean(y)
  
  pop_data <- data.frame(x, z, z_star, y)
  prob_sample <- pop_data[sample(1:N, n_prob), ]
  #create a biased non-probability sample
  nonprob_sample_x_large <- pop_data[as.numeric(sample(rownames(pop_data[x>2, ]), n_nonprob*nonprob_ratio)), ]
  nonprob_sample_x_small <- pop_data[as.numeric(sample(rownames(pop_data[x<=2, ]), n_nonprob*(1-nonprob_ratio))), ]
  nonprob_sample <- rbind(nonprob_sample_x_small, nonprob_sample_x_large)
  
  #gold standard estimator - uses precise values from probability sample
  #impossible in practice but used as a benchmark
  gold_standard <- lm(y ~ x + z, data = prob_sample)
  gold_standard_coefs <- gold_standard$coefficients
  y_pred_gold_standard <- pred(gold_standard_coefs, prob_sample)
  biases[i, 1] <- mean(y_pred_gold_standard) - y_mean
  
  #mass imputation estimator - uses precise values from non-probability sample to estimate values of probability sample
  #impossible in practice but used as an benchmark
  mass_imputation <- lm(y ~ x + z,  data = nonprob_sample)
  mass_imputation_coefs <- mass_imputation$coefficients
  y_pred_mass_imputation <- pred(mass_imputation_coefs, prob_sample)
  biases[i, 2] <- mean(y_pred_mass_imputation) - y_mean
  
  #naive mass imputation estimator - uses data from non-probability sample with measurement errors but does not correct them
  mass_imputation_naive <- lm(y ~ x + z_star,  data = nonprob_sample)
  mass_imputation_naive_coefs <- mass_imputation_naive$coefficients
  y_pred_mass_imputation_naive <- pred(mass_imputation_naive_coefs, prob_sample)
  biases[i, 3] <- mean(y_pred_mass_imputation_naive) - y_mean
  
  #estimate false positive rate on an audit sample
  s_audit <- pop_data[sample(1:nrow(pop_data), n_audit),]
  fpr <- sum(s_audit$z==0 & s_audit$z_star==1)/n_audit
  errors <- s_audit$z_star - s_audit$z
  lambda <- c(mean(errors), mean(s_audit$x*errors))
  
  #bca estimator - corrects measurement error from non-probability sample and then uses mass imputation
  bca_classic <- MLBC::ols_bca(y ~ z_star + x, data = nonprob_sample, fpr = fpr, m = n_audit)
  bca_classic_coef <- bca_classic$coef[c(3, 2, 1)]
  y_pred_bca_classic <- pred(bca_classic_coef, prob_sample)
  biases[i, 4] <- mean(y_pred_bca_classic)-y_mean
  
  #upgraded bca estimator - corrects measurement error from non-probability sample and then uses mass imputation 
  ksi_hat <- rbind(z_star, rep(1, N), x)
  psi_hat <- solve(1/N*ksi_hat%*%t(ksi_hat))%*%(1/N*ksi_hat%*%y)
  bca_upgraded <- (diag(3)+solve(1/N*ksi_hat%*%t(ksi_hat))%*%(cbind(c(fpr, lambda), matrix(0, 3, 2))))%*%psi_hat
  bca_upgraded_coefs <- bca_upgraded[c(2, 3, 1)]
  y_pred_bca_upgraded <- pred(bca_upgraded_coefs, prob_sample)
  biases[i, 5] <- mean(y_pred_bca_upgraded) - y_mean
  
  #bcm estimator - corrects measurement error using asymptotical behavior of inverse matrices
  bcm <- MLBC::ols_bcm(y ~ z_star + x, data = nonprob_sample, fpr = fpr, m = n_audit)
  bcm_coef <- bcm$coef[c(3, 2, 1)]
  y_pred_bcm <- pred(bcm_coef, prob_sample)
  biases[i, 6] <- mean(y_pred_bcm) - y_mean
  
  #upgraded bcm estimator - corrects measurement error using asymptotical behavior of inverse matrices
  bcm_upgr <- solve(diag(3)-solve(1/N*ksi_hat%*%t(ksi_hat))%*%(cbind(c(fpr, lambda), matrix(0, 3, 2))))%*%psi_hat
  bcm_upgr_coefs <- bcm_upgr[c(2, 3, 1)]
  y_pred_bcm_upgr <- pred(bcm_upgr_coefs, prob_sample)
  biases[i, 7] <- mean(y_pred_bcm_upgr) - y_mean
  
  #joint estimator - uses joint estimation to correct measurement error
  joint <- MLBC::one_step(y ~ z_star + x, data = nonprob_sample)
  joint_coef <- joint$coef[c(3, 2, 1)]
  y_pred_joint <- pred(joint_coef, prob_sample)
  biases[i, 8] <- mean(y_pred_joint) - y_mean
}
squared_errors <- biases^2
results_biases <- colMeans(biases)
results_mse <- colMeans(squared_errors)
results <- as.data.frame(rbind(results_biases, results_mse))
save(results, file="Results/Combining probability and nonprobability sampling.RData")
#stargazer::stargazer(rbind(results_biases, results_mse), type='latex', digits=4, title = "Results", label = "tab:summary", flip=T)

