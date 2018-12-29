library(MASS)
setwd("./documents/school/ams 274/hw5")
fabric = read.csv("fabric.csv")


#MH step to get betas 
# necessary arguments are current beta
beta.ll = function(beta, X, lambda, MU){
  V = -1*lambda*(beta[1] + beta[2] * X) - (lambda * MU)/exp(beta[1] + beta[2]*X)
  return (sum(V))
}


get.beta = function(beta.current, X, lambda, MU){
  #generate candidates
  #S = diag(2)*c(.001,.001) # maybe fit maximum likelihood and figure somethin good out. 
  # but for now forget it. 
  Jinv = matrix(c(2.620501e-02, -3.561562e-05, -3.561562e-05,  5.576763e-08), ncol = 2, nrow = 2)
  S = 3*Jinv
  candidate_vec = mvrnorm(n = 1, mu = beta.current, Sigma = S)
  
  #calculate log posterior at candidate point and at current point. 
  ll_current = beta.ll(beta.current, X, lambda, MU)
  ll_proposal = beta.ll(candidate_vec, X, lambda, MU)
  
  priors_current = 0#sum(dnorm(coefs[i, ], 0, 100, log = TRUE))
  priors_proposal = 0#sum(dnorm(candidate_vec, 0, 100, log = TRUE))

  log_post_current = ll_current + priors_current
  log_post_proposal = ll_proposal + priors_proposal
  #print(log_post_proposal - log_post_current)
  q = min(1, exp(log_post_proposal - log_post_current))
  
  if (runif(1) < q){
     return(candidate_vec)
  } else{
    return (beta.current)
  }
}


#MH step to get lambda 
lambda.ll = function(lambda, X, beta, MU){
  V = lambda*log(lambda) - lambda*(beta[1] + beta[2]*X) - lgamma(lambda) + 
    (lambda - 1)*log(MU) - (lambda*MU)/exp(beta[1] + beta[2]*X)
  lambda.prior = -2 * log(1+lambda)
  return (sum(V) + lambda.prior)
}

get.lambda = function(lambda.current, X, beta, MU){
  #generate candidates
  # but for now forget it. 
  k = 1
  candidate_vec = exp(rnorm(1, log(lambda.current), 1))
  
  #calculate log posterior at candidate point and at current point. 
  ll_current = lambda.ll(lambda.current, X, beta, MU)
  ll_proposal = lambda.ll(candidate_vec, X, beta, MU)
  
  priors_current = 0#sum(dnorm(coefs[i, ], 0, 100, log = TRUE))
  priors_proposal = 0#sum(dnorm(candidate_vec, 0, 100, log = TRUE))
  
  log_post_current = ll_current + priors_current +log(lambda.current)
  log_post_proposal = ll_proposal + priors_proposal +log(candidate_vec)
  
  #density.ratio = dgamma(lambda.current, k*candidate_vec, k, log = TRUE) - dgamma(candidate_vec, k*lambda.current, k, log = TRUE)
  
  #print (exp(log_post_proposal - log_post_current))
  q = min(1, exp(log_post_proposal - log_post_current))
  
  if (runif(1) < q){
    return(candidate_vec)
  } else{
    return (lambda.current)
  }
}


#step to get muis
get.MU = function(Y, X, lambda, beta){
  alpha = Y + lambda
  beta = 1 + lambda/exp(beta[1] + beta[2]*X)
  return (rgamma(32, alpha, beta))
}



################################################################
# GIBBS SAMPLER !!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################

n.iterations= 50000
n.burn = 10000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)

betas = matrix(0,ncol = 2, nrow = n.iterations)
#betas[1, ] = c(.9711, .0019)

lambdas = rep(2, n.iterations)

MUS = matrix(0,ncol = 32, nrow = n.iterations)
#MUS[1, ] = fabric$faults


for (i in 1:(n.iterations-1)){
  #Mu update
  MUS[i+1,] = get.MU(fabric$faults, 
                     fabric$length, 
                     lambdas[i], 
                     betas[i, ])
  
  #beta update
  betas[i+1, ] = get.beta(betas[i, ], fabric$length, lambdas[i], MUS[i+1, ])
  
  #lambda update
  #lambda.current, X, beta, MU
  lambdas[i+1] = get.lambda(lambdas[i], fabric$length, betas[i+1, ], MUS[i+1,])
  #print (i)
}


betas.f = betas[keeps, ]
lambdas.f = lambdas[keeps]
MUS.f = MUS[keeps, ]

ts.plot(betas.f[,1], main = paste("mean = ", mean(betas.f[,1])))
ts.plot(betas.f[,2],  main = paste("mean = ", mean(betas.f[,2])))
ts.plot(lambdas.f, main = paste('mean = ', mean(lambdas.f)))

################################################################
# GIBBS SAMPLER FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################

gibbs = function(n.iterations, n.burn, n.thin, fabric){
  keeps = seq(n.burn, n.iterations, by = n.thin)
  
  #sample storage 
  betas = matrix(0,ncol = 2, nrow = n.iterations)
  lambdas = rep(2, n.iterations)
  MUS = matrix(0,ncol = 32, nrow = n.iterations)
  
  for (i in 1:(n.iterations-1)){
    #Mu update
    MUS[i+1,] = get.MU(fabric$faults, 
                       fabric$length, 
                       lambdas[i], 
                       betas[i, ])
    
    #beta update
    betas[i+1, ] = get.beta(betas[i, ], fabric$length, lambdas[i], MUS[i+1, ])
    
    #lambda update
    #lambda.current, X, beta, MU
    lambdas[i+1] = get.lambda(lambdas[i], fabric$length, betas[i+1, ], MUS[i+1,])
    #print (i)
  }
  betas.f = betas[keeps, ]
  lambdas.f = lambdas[keeps]
  MUS.f = MUS[keeps, ]
  
  return (list(betas = betas.f, lambdas = lambdas.f, MUS = MUS.f))
}

post.samples = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric)

O = optim(c(0,0), function(x)beta.ll(x, fabric$length, mean(lambdas.f), colMeans(MUS.f)), 
      method = "BFGS", control = list(fnscale = -1), hessian = T)

Jinv = solve(-1 * O$hessian)

pdf("plot7.pdf", width = 8, height = 4)
mcmc_plots(betas.f[, 1], "beta1", O$par[1])
dev.off()
pdf("plot8.pdf", width = 8, height = 4)
mcmc_plots(betas.f[, 2], "beta2", O$par[2])
dev.off()
pdf("plot9.pdf", width = 8, height = 4)
mcmc_plots(lambdas.f, "lambda", NA)
dev.off()



####################################################################################
#ok so now what you want to do is generate point and interval estimates for the response mean. 
lengths.grid = seq(100, 1000, length.out = 1000)
fill = matrix(0, nrow = length(lengths.grid), ncol = 3)
for (i in 1:length(lengths.grid)){
  responses = exp(betas.f[,1] + betas.f[,2]*lengths.grid[i])
  #responses = rgamma(length(gamma.i), lambdas.f, lambdas.f/gamma.i)
  fill[i, 1] = quantile(responses, 0.025)
  fill[i, 2] = mean(responses)
  fill[i, 3] = quantile(responses, 0.975) 
}

pdf("plot5.pdf", width = 8, height = 4)
par(mfrow =  c(1,1))
plot(lengths.grid, fill[, 2], type = "l", main = "Response Mean: Hierarchical", 
     xlab = "Length of Roll", ylab = "Number of Faults", ylim = c(0,30))
lines(lengths.grid, fill[,1], col = "red", lty = 2)
lines(lengths.grid, fill[,3], col = "red", lty = 2)
points(fabric$length, fabric$faults, col = "blue")
dev.off()
#######################################################################################



#######################################################################################
# Now do posterior predictive residuals. 

pdf("plot6.pdf", width = 8, height = 4)
n = nrow(fabric)
ppr = matrix(0, nrow = n, ncol = 2)
for (i in 1:n){
  #mus = exp(betas.f[,1] + betas.f[,2]*fabric[i, "length"])
  mus2 = rpois(nrow(MUS.f), MUS.f[, i])#sapply(mus, function(x)rpois(1,x))
  ris = fabric$faults[i] - mus2
  ppr[i, 1] = quantile(ris, 0.025)
  ppr[i, 2] = quantile(ris, 0.975) 
}



#residuals plot, ok this can't be right. Can it? 
plot(1, type="n", xlab="Length of Roll", ylab="Posterior Predictive Residuals", 
     xlim=c(100, 1000), ylim=c(-20, 20), main = "Posterior Predictive Residuals: Hierarchical")
abline(h=0)
for (j in 1:n){
  if (ppr[j,1] <0 & ppr[j,2] >0){
    c = "red"
  } else{
    c = "blue"
  }
  print(ppr[j,1])
  print(ppr[j,2])
  segments(x0 = fabric$length[j], y0 = ppr[j, 1], y1 = ppr[j, 2], col = c)
}
dev.off()
###########################################################################################



###########################################################################################
pdf("plot10.pdf", width = 8, height = 4)

prior.powers = c(1.01, 2, 10)
lambdas.grid = seq(0,10, length.out = 10000)

prior.provided = function(lambda, p){
  return ((1 + lambda)^(-1*p))
}
plot(lambdas.grid, prior.provided(lambdas.grid, prior.powers[1]), main = "Prior Choices",
     xlab = "lambda", ylab = "Density", type = "l", ylim = c(0,1))
lines(lambdas.grid, prior.provided(lambdas.grid, prior.powers[2]), col = "red")
lines(lambdas.grid, prior.provided(lambdas.grid, prior.powers[3]), col = "blue")
lines(lambdas.grid, dgamma(lambdas.grid, .01, .01), col = "green")
legend(x = "topright", legend = c("p = 1.01", "p = 2", "p = 10", "Gamma(.01, .01)"), 
       col = c("Black", "red", "blue", "green"), lty = 1)
dev.off()













