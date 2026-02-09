set.seed(42)
N = 8000
m = 1000
p = 0.7
n_sim = 1000
results = numeric(n_sim)

for(i in 1:n_sim){
  theta = rbinom(N, 1, p)
  epsilon = rnorm(N)
  y = 10 + theta + (0.3 + 0.2*theta)*epsilon
  
  kappa = 1
  fpr = kappa/sqrt(N)
  theta_star = theta
  theta_star[theta == 0] = rbinom(sum(theta == 0), 1, fpr/(1-p))
  theta_star[theta == 1] = 1 - rbinom(sum(theta == 1), 1, fpr/p)
  data = data.frame(y, theta, theta_star)
  
  bca1 = MLBC::ols_bca(y ~ theta_star, data = data, fpr = fpr, m = m)
  results[i] = bca1$coef['theta_star']
}
mean(results)