sim=function(n){
  bledy_bca=numeric(50)
  bledy_bcm=numeric(50)
  bledy_naiwne=numeric(50)
  for(i in 1:50){
  gamma=runif(10, 1, 2)
  alfa=runif(1, -0.5, 0.5)
  p=c(0.25, 0.2, 0.15, 0.05, 0.05, 0.1, 0.02, 0.03, 0.05, 0.1)
  pomylki=diag(10)*70+matrix(rpois(100, 3), nrow=10)
  pomylki = sweep(pomylki, MARGIN = 2, STATS = colSums(pomylki), FUN = '/')
  q=rbinom(n, 1, 0.7)
  thety=sample(1:10, size=n, replace = T, prob = p)
  y=rpois(n, exp(gamma[thety]+alfa*q))
  d=sapply(1:n, function(x){sample(1:10, size = 1, prob = pomylki[, thety[x]])})
  
  fpr=sum(sapply(2:10, function(x){sum(pomylki[x,]*p)-pomylki[x, x]*p[x]}))
  omega_hat=(diag(sapply(2:10, function(x){sum(pomylki[x,]*p)}))-sweep(pomylki[2:10, 2:10], MARGIN = 2, STATS = p[2:10], FUN='*'))/fpr
  omega_hat=as.matrix(Matrix::bdiag(matrix(0, 1, 1), omega_hat, matrix(0, 1, 1)))
  
  model_naiwny = glm(y ~ as.factor(d) + q, family = "poisson")
  glm_wyniki=coef(model_naiwny)
  psi_hat=glm_wyniki
  d_factor=relevel(factor(d, levels=as.character(1:10)), ref='1')
  X=model.matrix(~d_factor +q)
  m=solve(t(X)%*%X/n)
  
  psi_bca=(diag(11)+fpr*m%*%omega_hat)%*%psi_hat
  psi_bca[2:10]=psi_bca[2:10]+psi_bca[1]
  bledy_bca[i]=sum((psi_bca-c(gamma, alfa))^2)
  
  psi_bcm=solve(diag(11)-fpr*m%*%omega_hat)%*%psi_hat
  psi_bcm[2:10]=psi_bcm[2:10]+psi_bcm[1]
  bledy_bcm[i]=sum((psi_bcm-c(gamma, alfa))^2)
  
  glm_wyniki[2:10]=glm_wyniki[2:10]+glm_wyniki[1]
  bledy_naiwne[i]=sum((glm_wyniki-c(gamma, alfa))^2)
  }
  return(list(bledy_bca, bledy_bcm, bledy_naiwne))
}

wynik_1000=sim(1000)
bledy_bca_1000=wynik_1000[[1]]
bledy_bcm_1000=wynik_1000[[2]]
bledy_naiwne_1000=wynik_1000[[3]]
mean(bledy_bca_1000)
mean(bledy_bcm_1000)
mean(bledy_naiwne_1000)

wynik_10000=sim(10000)
bledy_bca_10000=wynik_10000[[1]]
bledy_bcm_10000=wynik_10000[[2]]
bledy_naiwne_10000=wynik_10000[[3]]
mean(bledy_bca_10000)
mean(bledy_bcm_10000)
mean(bledy_naiwne_10000)

wynik_100000=sim(100000)
bledy_bca_100000=wynik_100000[[1]]
bledy_bcm_100000=wynik_100000[[2]]
bledy_naiwne_100000=wynik_100000[[3]]
mean(bledy_bca_100000)
mean(bledy_bcm_100000)
mean(bledy_naiwne_100000)

wynik_bc=as.data.frame(c(wynik_1000, wynik_10000, wynik_100000))
colnames(wynik_bc)=c('bca1000', 'bcm1000', 'naive1000', 'bca10000', 'bcm10000', 'naive10000', 'bca100000', 'bcm100000', 'naive100000')
wynik_bc=rbind(wynik_bc, colMeans(wynik_bc))
rownames(wynik_bc)[51]='Mean Error'
save(wynik_bc, file = "Results/Bias corrected estimators results.RData")

#stargazer::stargazer(wynik_bc, type='latex', title = "Results", label = "tab:summary", flip=T)
