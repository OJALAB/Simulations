library(optimx)
l_poiss = function(y, d, gamma, intercept, alfa, q, omega, p){#y - liczba, d - liczba od 1 do 10, gamma - wektor dł. 9, alfa- liczba
  ref=alfa*q+intercept #q - 0 lub 1, omega - macierz 10x10, p - wektor dł. 10
  wynik=omega[d, 10]*p[10]*dpois(y, lambda = exp(ref))
  for(i in 1:9){
    logmi=c(gamma %*% diag(9)[, i])+ref
    wynik=wynik+omega[d, i]*p[i]*dpois(y, lambda = exp(logmi))
  }
  penalty=sum(c(gamma, intercept, alfa)^2)
  return(sum(log(wynik))-0.5*penalty)
}
opt = function(wektor, y, d, q, omega, p){
  return(-l_poiss(y=y, d=d, gamma=wektor[1:9], intercept=wektor[10], alfa=wektor[11], q=q, omega=omega, p=p))
}
opt2= function(wektor, y, d, q, omega){
  gamma=wektor[1:9]
  intercept=wektor[10]
  alfa=wektor[11]
  p=wektor[12:20]
  p=exp(c(p, 0))
  p=p/sum(p)
  return(-l_poiss(y, d, gamma, intercept, alfa, q, omega, p))
}

bledy_1000=numeric(25)
bledy2_1000=numeric(25)
bledy_uproszczone_1000=numeric(25)
for(i in 1:25){
  gamma=runif(10, 1, 2)
  alfa=runif(1, -0.5, 0.5)
  p=c(0.25, 0.2, 0.15, 0.05, 0.05, 0.1, 0.02, 0.03, 0.05, 0.1)
  pomylki=diag(10)*70+matrix(rpois(100, 3), nrow=10)
  pomylki = sweep(pomylki, MARGIN = 2, STATS = colSums(pomylki), FUN = '/')
  q=rbinom(1000, 1, 0.7)
  thety=sample(1:10, size=1000, replace = T, prob = p)
  y=rpois(1000, exp(gamma[thety]+alfa*q))
  d=sapply(1:1000, function(x){sample(1:10, size = 1, prob = pomylki[, thety[x]])})

  fit=optimx(par=c(numeric(9), log(mean(y)), 0), fn=opt, y=y, d=d, q=q, omega=pomylki, p=p, method = "BFGS")
  wyniki=as.numeric(fit[1, 1:11])
  wyniki[1:9]=wyniki[1:9]+wyniki[10]
  bledy_1000[i]=sum((c(gamma, alfa)-wyniki)^2)
  
  fit2=optimx(par=c(numeric(9), log(mean(y)), numeric(10)), fn=opt2, y=y, d=d, q=q, omega=pomylki, method = "BFGS")
  wyniki2=as.numeric(fit2[1, 1:11])
  wyniki2[1:9]=wyniki2[1:9]+wyniki2[10]
  bledy2_1000[i]=sum((c(gamma, alfa)-wyniki2)^2)

  model_naiwny = glm(y ~ as.factor(d) + q, family = "poisson")
  glm_wyniki=coef(model_naiwny)
  glm_wyniki[2:10]=glm_wyniki[2:10]+glm_wyniki[1]
  bledy_uproszczone_1000[i]=sum((c(gamma, alfa)-glm_wyniki)^2)
}
mean(bledy_1000)
mean(bledy2_1000)
mean(bledy_uproszczone_1000)

bledy_5000=numeric(25)
bledy_uproszczone_5000=numeric(25)
bledy2_5000=numeric(25)
for(i in 1:25){
  gamma=runif(10, 1, 3)
  alfa=runif(1, -0.5, 0.5)
  p=gtools::rdirichlet(1, alpha = rep(3, 10))
  dokladnosc=runif(1, 60, 90)
  pomylki=diag(10)*dokladnosc+matrix(rpois(100, (100-dokladnosc)/10), nrow=10)
  pomylki = sweep(pomylki, MARGIN = 2, STATS = colSums(pomylki), FUN = '/')
  prawdopodobienstwo_q=runif(1, 0.4, 0.9)
  q=rbinom(5000, 1, prawdopodobienstwo_q)
  thety=sample(1:10, size=5000, replace = T, prob = p)
  y=rpois(5000, exp(gamma[thety]+alfa*q))
  d=sapply(1:5000, function(x){sample(1:10, size = 1, prob = pomylki[, thety[x]])})
  
  fit=optimx(par=c(numeric(9), log(mean(y)), 0), fn=opt, y=y, d=d, q=q, omega=pomylki, p=p, method = "BFGS")
  wyniki=as.numeric(fit[1, 1:11])
  wyniki[1:9]=wyniki[1:9]+wyniki[10]
  bledy_5000[i]=sum((c(gamma, alfa)-wyniki)^2)
  
  fit2=optimx(par=c(numeric(9), log(mean(y)), numeric(10)), fn=opt2, y=y, d=d, q=q, omega=pomylki, method = "BFGS")
  wyniki2=as.numeric(fit2[1, 1:11])
  wyniki2[1:9]=wyniki2[1:9]+wyniki2[10]
  bledy2_5000[i]=sum((c(gamma, alfa)-wyniki2)^2)

  model_naiwny = glm(y ~ as.factor(d) + q, family = "poisson")
  glm_wyniki=coef(model_naiwny)
  glm_wyniki[2:10]=glm_wyniki[2:10]+glm_wyniki[1]
  bledy_uproszczone_5000[i]=sum((c(gamma, alfa)-glm_wyniki)^2)
}
mean(bledy_5000)
mean(bledy2_5000)
mean(bledy_uproszczone_5000)

wynik=data.frame(c1=bledy_1000, c2=bledy2_1000, c3=bledy_uproszczone_1000, c4=bledy_5000, c5=bledy2_5000, c6=bledy_uproszczone_5000)
wynik=rbind(wynik, colMeans(wynik))
rownames(wynik)[26]='Mean Error'
colnames(wynik)=c('Joint estimation 1000 observations',
                  'Joint esitmation 1000 observations with unknown probabilities',
                  'Standard GLM 1000 observations',
                  'Joint estimation 5000 observations',
                  'Joint esitmation 5000 observations with unknown probabilities',
                  'Standard GLM 5000 observations')
save(wynik, file = "Wyniki symulacji.RData")

stargazer::stargazer(as.data.frame(wynik), summary.stat = c("mean", "median", "sd", "min", "max"), type='latex', title = "Results", digits=4, flip=T)
