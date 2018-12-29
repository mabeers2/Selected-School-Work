#hw4

####################################################################
######## PROBLEM 3              ####################################
####################################################################
setwd("./Documents/school/AMS 274/hw3/hw4")
mouse = read.csv("mouse.csv")




# two likelihoods to be maximized using optim 

#likelihood 1 : Bin(y_1 | m, pi1)
# y1 = dead
# y2 = malformation
# y3 = normal
ll1 = function(v, y, m, x){
  #alpha = v[1], beta = v[2]
  pi1 = exp(v[1] + v[2]*x)/(1 + exp(v[1] + v[2]*x))
  #print(pi1)
  return (sum(dbinom(x = y, size = m, prob = pi1, log = TRUE)))
}

ll1.optim = optim(c(0,0), function(x)ll1(x, mouse$dead, mouse$total, mouse$concentration), 
                  control = list(fnscale = -1), hessian = TRUE)

ll1.optim
MLEs.ll1 = ll1.optim$par
SEs.ll1 =  sqrt(diag(solve(-1*ll1.optim$hessian)))
MLEs.ll1
SEs.ll1


#likelihood 2 : Bin(y_2 | m-y1, pi2/(1-pi1))
# y1 = dead
# y2 = malformation
# y3 = normal
ll2 = function(v, y, m, x, ab1){
  #alpha = v[1], beta = v[2]
  pi1 = exp(ab1[1] + ab1[2]*x)/(1 + exp(ab1[1] + ab1[2]*x))
  pi2 = exp(v[1] + v[2]*x)/((1 + exp(v[1] + v[2]*x))*(1 + exp(ab1[1] + ab1[2]*x)))
  pi = pi2/(1-pi1)
  #print(pi1)
  return (sum(dbinom(x = y, size = m, prob = pi, log = TRUE)))
}

ll2.2 = function(v,y,m,x){
  #alpha = v[1], beta = v[2]
  pi1 = exp(v[1] + v[2]*x)/(1 + exp(v[1] + v[2]*x))
  #print(pi1)
  return (sum(dbinom(x = y, size = m, prob = pi1, log = TRUE)))
  
}

ll2.2.optim = optim(c(0,0), function(x)ll2.2(x, mouse$malformation, mouse$total- mouse$dead, 
                                           mouse$concentration), 
                    control = list(fnscale = -1), hessian = TRUE)
ll2.2.optim


ll2.optim = optim(c(0,0), function(x)ll2(x, mouse$malformation, mouse$total- mouse$dead, 
                                         mouse$concentration, ll1.optim$par), 
                  control = list(fnscale = -1), hessian = TRUE)

ll2.optim

MLEs.ll2 = ll2.2.optim$par
SEs.ll2 =  sqrt(diag(solve(-1*ll2.2.optim$hessian)))
MLEs.ll2
SEs.ll2

#plot estimated response curves. 
x.grid = seq(0,500, length.out = 501)
pi1.hat = exp(MLEs.ll1[1] + MLEs.ll1[2]*x.grid)/(1 + exp(MLEs.ll1[1] + MLEs.ll1[2]*x.grid))
pi2.hat = exp(MLEs.ll2[1] + MLEs.ll2[2]*x.grid)/((1 + exp(MLEs.ll1[1] + MLEs.ll1[2]*x.grid))*(1 + exp(MLEs.ll2[1] + MLEs.ll2[2]*x.grid)))
pi3.hat = 1/((1 + exp(MLEs.ll1[1] + MLEs.ll1[2]*x.grid))*(1 + exp(MLEs.ll2[1] + MLEs.ll2[2]*x.grid)))

pdf("plot13.pdf", width = 8, height = 4)
par(mfrow =c(1,1))
plot(x.grid, pi1.hat, type= "l", ylim = c(0,1), 
     xlab = "Concentration (mg/kg per day)", ylab = "pi(x)", main = "Estimated Response Curves (MLEs)") #dead 
lines(pi2.hat, col = 'red') #malformed
lines(pi3.hat, col = "blue") #normal
points(mouse$concentration, mouse$dead/mouse$total)
points(mouse$concentration, mouse$malformation/mouse$total, col = "red")
points(mouse$concentration, mouse$normal/mouse$total, col = "blue")
legend(x = "topright", legend = c("Dead", "Malformation", "Normal"), 
       col = c("Black", "red", "blue"), lty = 1, pch=1)
dev.off()

#discuss that we don't include the error associated with pi1 in these estimates, but despite this, they
# seem to work pretty well following the empirical proportions fairly closely. 




############################################################################
######## part c               ##############################################
############################################################################

#binomial MH log likelihood1

library(MASS)
n.iterations= 10000
n.burn = 5000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)

coefs = matrix(0,ncol = 2, nrow = n.iterations)

for (i in 1:(n.iterations-1)){
  
  #generate candidates
  candidate_vec = mvrnorm(n = 1, mu = c(coefs[i, ]), Sigma = -1*solve(ll1.optim$hessian))
  
  #calculate log posterior at candidate point and at current point. 
  #ll1(x, mouse$dead, mouse$total, mouse$concentration)
  ll_current = ll1(coefs[i, ], mouse$dead, mouse$total, mouse$concentration)
  ll_proposal = ll1(candidate_vec, mouse$dead, mouse$total, mouse$concentration)
  
  priors_current = sum(dnorm(coefs[i, ], 0, 100, log = TRUE))
  priors_proposal = sum(dnorm(candidate_vec, 0, 100, log = TRUE))
  
  log_post_current = ll_current + priors_current
  log_post_proposal = ll_proposal + priors_proposal
  
  q = min(1, exp(log_post_proposal - log_post_current))
  
  if (runif(1) < q){
    coefs[i+1, ] = candidate_vec
  } else{
    coefs[i + 1, ] = coefs[i, ]
  }
  print (i)
}

par(mfrow = c(1,1))
coefs.ll1 = coefs[keeps, ]
ts.plot(coefs.ll1[, 1])
ts.plot(coefs.ll1[, 2])
colMeans(coefs.ll1)
ll1.optim$par

mcmc_plots = function(v, t, MLE, t2 = t){
  par(mfrow = c(1,2), mar=c(5.1,4.1,4.1,2.1))
  m = mean(v)
  ts.plot(v, ylab = t, main = t)
  hist(v, main = t2, xlab = t2, sub = paste("Mean = ", round(m, 4)))
  abline(v = m, col = "blue")
  abline(v = MLE, col = "green")
  abline(v = quantile(v, c(.025, .975)), col = "red", lty = 2)
}

pdf("plot14.pdf", width = 8, height = 4)
mcmc_plots(coefs.ll1[, 1], "alpha1", ll1.optim$par[1])
dev.off()
pdf("plot15.pdf", width = 8, height = 4)
mcmc_plots(coefs.ll1[, 2], "beta1", ll1.optim$par[2])
dev.off()





#binomial MH log likelihood2

ll2 = function(v, y, m, x, ab1){
  #alpha = v[1], beta = v[2]
  pi1 = exp(ab1[1] + ab1[2]*x)/(1 + exp(ab1[1] + ab1[2]*x))
  pi2 = exp(v[1] + v[2]*x)/((1 + exp(v[1] + v[2]*x))*(1 + exp(ab1[1] + ab1[2]*x)))
  pi = pi2/(1-pi1)
  #print(pi1)
  return (sum(dbinom(x = y, size = m, prob = pi, log = TRUE)))
}

ll2.c = function(v, y, m, x, ab1){
  #alpha = v[1], beta = v[2]
  n = nrow(ab1)
  lls = numeric(n)
  for (i in 1:n){
    lls[i] = ll2(v, y, m ,x, ab1[i, ])
  }
  #hist(lls, sub = paste("mean = ", mean(lls)))

  return (mean(lls))
}
#ll2(x, mouse$malformation, mouse$total- mouse$dead, mouse$concentration, ll1.optim$par)
ll2.c(MLEs.ll2, mouse$malformation, mouse$total - mouse$dead, mouse$concentration, coefs.ll1)
ll2(MLEs.ll2, mouse$malformation, mouse$total - mouse$dead, mouse$concentration, MLEs.ll2)

n.iterations= 10000
n.burn = 5000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)

coefs = matrix(0,ncol = 2, nrow = n.iterations)

for (i in 1:(n.iterations-1)){
  
  #generate candidates
  candidate_vec = mvrnorm(n = 1, mu = c(coefs[i, ]), Sigma = -1*solve(ll2.optim$hessian))
  
  #calculate log posterior at candidate point and at current point. 
  #ll1(x, mouse$dead, mouse$total, mouse$concentration)
  ll_current = ll2.2(coefs[i, ], mouse$malformation, mouse$total - mouse$dead, mouse$concentration)
  ll_proposal = ll2.2(candidate_vec, mouse$malformation, mouse$total - mouse$dead, mouse$concentration)
  
  priors_current = sum(dnorm(coefs[i, ], 0, 100, log = TRUE))
  priors_proposal = sum(dnorm(candidate_vec, 0, 100, log = TRUE))
  
  log_post_current = ll_current + priors_current
  log_post_proposal = ll_proposal + priors_proposal
  
  q = min(1, exp(log_post_proposal - log_post_current))
  
  if (runif(1) < q){
    coefs[i+1, ] = candidate_vec
  } else{
    coefs[i + 1, ] = coefs[i, ]
  }
  print (i)
}

par(mfrow = c(1,1))
coefs.ll2 = coefs[keeps, ]
ts.plot(coefs.ll2[, 1])
ts.plot(coefs.ll2[, 2])
colMeans(coefs.ll2)
ll2.optim$par


pdf("plot16.pdf", width = 8, height = 4)
mcmc_plots(coefs.ll2[, 1], "alpha2", ll2.optim$par[1])
dev.off()
pdf("plot17.pdf", width = 8, height = 4)
mcmc_plots(coefs.ll2[, 2], "beta2", ll2.optim$par[2])
dev.off()



# ok so you have decent posterior samples, now do the response probabilities. 
x.grid = seq(0,500, length.out = 501)
n = length(x.grid)
M = matrix(0,nrow = n, ncol = 9)


for (i in 1:n){
  xit.beta.ll1 = coefs.ll1[, 1] + coefs.ll1[,2]*x.grid[i]
  xit.beta.ll2 = coefs.ll2[, 1] + coefs.ll2[,2]*x.grid[i]
  pi1.hat.xi.samples = exp(xit.beta.ll1)/(1 + exp(xit.beta.ll1))
  pi2.hat.xi.samples = exp(xit.beta.ll2)/((1 + exp(xit.beta.ll1))*(1 + exp(xit.beta.ll2)))
  pi3.hat.xi.samples = 1/((1 + exp(xit.beta.ll1))*(1 + exp(xit.beta.ll2)))
  M[i, 1] = quantile(pi1.hat.xi.samples, .025)
  M[i, 2] = quantile(pi2.hat.xi.samples, .025)
  M[i, 3] = quantile(pi3.hat.xi.samples, .025)
  M[i, 4] = mean(pi1.hat.xi.samples)
  M[i, 5] = mean(pi2.hat.xi.samples)
  M[i, 6] = mean(pi3.hat.xi.samples)
  M[i, 7] = quantile(pi1.hat.xi.samples, .975)
  M[i, 8] = quantile(pi2.hat.xi.samples, .975)
  M[i, 9] = quantile(pi3.hat.xi.samples, .975)
}


pdf("plot18.pdf", width = 8, height = 8)
par(mfrow = c(3,1))
plot(x.grid, M[, 4], type = "l", ylab = "pi1(x)", xlab = "Concentration (mg/kg per day)", 
     main = "Estimated Response: Dead", ylim = c(0,1))
lines(x.grid, M[,1], col = "red", lty = 2)
lines(x.grid, M[,7], col = "red", lty = 2)
points(mouse$concentration, mouse$dead/mouse$total, col = "blue")

plot(x.grid, M[, 5], type = "l", ylab = "pi2(x)", xlab = "Concentration (mg/kg per day)", 
     main = "Estimated Response: Malformation", ylim = c(0,1))
lines(x.grid, M[,2], col = "red", lty = 2)
lines(x.grid, M[,8], col = "red", lty = 2)
points(mouse$concentration, mouse$malformation/mouse$total, col = "blue")

plot(x.grid, M[, 6], type = "l", ylab = "pi3(x)", xlab = "Concentration (mg/kg per day)", 
     main = "Estimated Response: Normal", ylim = c(0,1))
lines(x.grid, M[,3], col = "red", lty = 2)
lines(x.grid, M[,9], col = "red", lty = 2)
points(mouse$concentration, mouse$normal/mouse$total, col = "blue")
dev.off()

