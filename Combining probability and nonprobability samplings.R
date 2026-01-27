set.seed(42)
N = 100000
n_sim = 1000
n_prob = 500
n_nonprob = 1000
n_audit = 1000
nonprob_ratio = 0.7
biases = matrix(0, n_sim, 5)
colnames(biases) = c('gold', 'mass_imp', 'mass_impnaive', 'bca', 'bca_upgr')
for(i in 1:n_sim){
  x = rnorm(N, 2, 1)
  z = rbinom(N, 1, 0.7)
  z_star = z
  z_star[z == 1] = sample(c(0,1), size = NROW(z_star[z == 1]), prob = c(0.10, 0.90), replace = T)
  z_star[z == 0] = sample(c(0,1), size = NROW(z_star[z == 0]), prob = c(0.93, 0.07), replace = T)
  e = rnorm(N)
  y = 1 + x - 2*z + e
  y_mean = mean(y)
  
  pop_data = data.frame(x, z, z_star, y)
  prob_sample = pop_data[sample(1:N, n_prob), ]
  nonprob_sample_x_large = pop_data[as.numeric(sample(rownames(pop_data[x>2, ]), n_nonprob*nonprob_ratio)), ]
  nonprob_sample_x_small = pop_data[as.numeric(sample(rownames(pop_data[x<=2, ]), n_nonprob*(1-nonprob_ratio))), ]
  nonprob_sample = rbind(nonprob_sample_x_small, nonprob_sample_x_large)
  
  #gold standard esitmator - uses precise values from probability sample, impossible in practice
  gold_standard = lm(y ~ x + z, data = prob_sample)
  gold_standard_coefs = gold_standard$coefficients
  y_pred_gold_standard = gold_standard_coefs[1] + gold_standard_coefs[2] * prob_sample[['x']] + gold_standard_coefs[3] * prob_sample[['z']]
  bias_gold_standard = mean(y_pred_gold_standard) - y_mean
  biases[i, 1] = bias_gold_standard
  
  #mass imputation estimator - uses precise values from nonprobability sample to estimate values of probability sample, impossible in practice
  mass_imputation = lm(y ~ x + z,  data = nonprob_sample)
  mass_imputation_coefs = mass_imputation$coefficients
  y_pred_mass_imputation = mass_imputation_coefs[1]+mass_imputation_coefs[2]*prob_sample[['x']]+mass_imputation_coefs[3]*prob_sample[['z']]
  bias_mass_imputation = mean(y_pred_mass_imputation) - y_mean
  biases[i, 2] = bias_mass_imputation
  
  #naive mass imputation estimator - uses data from nonprobability sample with measuerement errors but does not correct them
  mass_imputation_naive = lm(y ~ x + z_star,  data = nonprob_sample)
  mass_imputation_naive_coefs = mass_imputation_naive$coefficients
  y_pred_mass_imputation_naive = mass_imputation_naive_coefs[1]+mass_imputation_naive_coefs[2]*prob_sample[['x']]+mass_imputation_naive_coefs[3]*prob_sample[['z']]
  bias_mass_imputation_naive = mean(y_pred_mass_imputation_naive) - y_mean
  biases[i, 3] = bias_mass_imputation_naive
  
  #estimate false positive rate on an audit sample
  s_audit = pop_data[sample(1:nrow(pop_data), n_audit),]
  fpr=sum(s_audit$z==0 & s_audit$z_star==1)/n_audit
  errors=s_audit$z_star - s_audit$z
  lambda=c(mean(errors), mean(s_audit$x*errors))
  
  #bca estimator - correctrs measurement error from nonprobability sample and then uses mass imputation 
  bca_classic = MLBC::ols_bca(y ~ z_star + x, data = nonprob_sample, fpr = fpr, m = n_audit)
  bca_classic_coef = bca_classic$coef
  y_pred_bca_classic = bca_classic_coef[3]+bca_classic_coef[2]*prob_sample[['x']]+bca_classic_coef[1]*prob_sample[['z']]
  bias_bca_classic = mean(y_pred_bca_classic)-y_mean
  biases[i, 4] = bias_bca_classic
  
  #upgraded bca estimator - correctrs measurement error from nonprobability sample and then uses mass imputation 
  ksi_hat=rbind(z_star, rep(1, N), x)
  psi_hat=solve(1/N*ksi_hat%*%t(ksi_hat))%*%(1/N*ksi_hat%*%y)
  bca_upgraded=(diag(3)+solve(1/N*ksi_hat%*%t(ksi_hat))%*%(cbind(c(fpr, lambda), matrix(0, 3, 2))))%*%psi_hat
  y_pred_bca_upgraded = bca_upgraded[2]+bca_upgraded[3]*prob_sample[['x']]+bca_upgraded[1]*prob_sample[['z']]
  bias_bca_upgraded = mean(y_pred_bca_upgraded) - y_mean
  biases[i, 5] = bias_bca_upgraded
}
squared_errors = biases^2
results_biases = colMeans(biases)
results_mse = colMeans(squared_errors)
stargazer::stargazer(rbind(results_biases, results_mse), type='latex', title = "Results", label = "tab:summary")
