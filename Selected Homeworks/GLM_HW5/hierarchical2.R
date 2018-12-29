#hierarchical 2, actually decent code. 

library(MASS)
setwd("./documents/school/ams 274/hw5")
fabric = read.csv("fabric.csv")

#MH step to get betas 
# necessary arguments are current beta
beta.ll = function(beta, X, lambda, MU){
  V = -1*lambda*(beta[1] + beta[2] * X) - (lambda * MU)/exp(beta[1] + beta[2]*X)
  return (sum(V))
}


mcmc_plots = function(v, t, MLE, t2 = t){
  par(mfrow = c(1,2), mar=c(5.1,4.1,4.1,2.1))
  m = mean(v)
  ts.plot(v, ylab = t, main = t)
  hist(v, main = t2, xlab = t2, sub = paste("Mean = ", round(m, 4)))
  abline(v = m, col = "blue")
  abline(v = MLE, col = "green")
  abline(v = quantile(v, c(.025, .975)), col = "red", lty = 2)
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
  #lambda.prior = -2 * log(1+lambda)
  return (sum(V))
}

lambda.prior.p = function(lambda, p){
  return (-p*log(1 + lambda))
}

lambda.prior.gamma = function(lambda, alpha = .01, beta = .01){
  return (dgamma(lambda, alpha, beta, log = TRUE))
}

get.lambda = function(lambda.current, X, beta, MU, f.lambda.prior){
  #generate candidates
  # but for now forget it. 
  k = 1
  candidate_vec = exp(rnorm(1, log(lambda.current), 1))
  
  #calculate log posterior at candidate point and at current point. 
  ll_current = lambda.ll(lambda.current, X, beta, MU)
  ll_proposal = lambda.ll(candidate_vec, X, beta, MU)
  
  priors_current = f.lambda.prior(lambda.current)#sum(dnorm(coefs[i, ], 0, 100, log = TRUE))
  priors_proposal = f.lambda.prior(candidate_vec)#sum(dnorm(candidate_vec, 0, 100, log = TRUE))
  
  log_post_current = ll_current + priors_current +log(lambda.current)
  log_post_proposal = ll_proposal + priors_proposal +log(candidate_vec)
  
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


gibbs = function(n.iterations, n.burn, n.thin, fabric, lprior){
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
    lambdas[i+1] = get.lambda(lambdas[i], fabric$length, betas[i+1, ], MUS[i+1,], lprior)
    #print (i)
  }
  betas.f = betas[keeps, ]
  lambdas.f = lambdas[keeps]
  MUS.f = MUS[keeps, ]
  
  return (list(betas = betas.f, lambdas = lambdas.f, MUS = MUS.f))
}




post.samples.1.01 = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                     lprior = function(x)lambda.prior.p(x, 1.01))

post.samples.2 = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                          lprior = function(x)lambda.prior.p(x, 2))

post.samples.10 = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                       lprior = function(x)lambda.prior.p(x, 10))

post.samples.20 = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                        lprior = function(x)lambda.prior.p(x, 20))

post.samples.9 = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                        lprior = function(x)lambda.prior.p(x, 9))



post.samples.gamma = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                        lprior = function(x)lambda.prior.gamma(x))

pdf("plot11.pdf", width = 8, height = 4)
mcmc_plots(post.samples.1.01$lambdas, "lambda, p = 1.01", NA)
dev.off()
pdf("plot12.pdf", width = 8, height = 4)
mcmc_plots(post.samples.2$lambdas, "lambda, p = 2", NA)
dev.off()
pdf("plot13.pdf", width = 8, height = 4)
mcmc_plots(post.samples.10$lambdas, "lambda, p = 10", NA)
dev.off()
pdf("plot14.pdf", width = 8, height = 4)
mcmc_plots(post.samples.gamma$lambdas, "lambda, Gamma prior", NA)
dev.off()


# Ok so now the question is how to compare the effect of prior. Maybe Gelfand and Ghosh? 
# see if there's a difference there. What is the argument k in gelfand and ghosh? 

get.posterior.predictive = function(posterior_samples, x0){
  gamma_i = exp(posterior_samples[["betas"]][,1] + posterior_samples[["betas"]][,2] * x0)
  n = length(gamma_i)
  mu0 = rgamma(n, posterior_samples[['lambdas']], posterior_samples[['lambdas']]/gamma_i)
  y0 = rpois(n, mu0)
  return (y0)
}
test = get.posterior.predictive(post.samples.1.01, 500)

x0.grid = seq(min(fabric$length)-20, max(fabric$length)+20, length.out = 400)

get.post.pred.summary = function(samples, x.grid){
  summary.df = data.frame(x = x.grid)
  n = length(x.grid)
  q025 = numeric(n)
  means = numeric(n)
  q975 = numeric(n)
  for (i in 1:n){
    post.pred = get.posterior.predictive(samples, x.grid[i])
    q025[i] = quantile(post.pred, 0.025)
    q975[i] = quantile(post.pred, 0.975)
    means[i] = mean(post.pred)
  }
  summary.df$q025 = q025
  summary.df$q975 = q975
  summary.df$means = means
  return (summary.df)
}

summary.1.01  = get.post.pred.summary(post.samples.1.01, x0.grid)
summary.2  = get.post.pred.summary(post.samples.2, x0.grid)
summary.10 = get.post.pred.summary(post.samples.10, x0.grid)
summary.gamma = get.post.pred.summary(post.samples.gamma, x0.grid)



pdf("plot15.pdf", width = 8, height = 4)
par(mfrow = c(1,1))
plot(summary.1.01$x, summary.1.01$means, type = "l", ylim = c(0, 45), xlab = "Length", 
     ylab = "faults", main = "Posterior Predictive Intervals Under different Priors")
lines(summary.1.01$x, summary.1.01$q025, lty = 2)
lines(summary.1.01$x, summary.1.01$q975, lty = 2)

lines(summary.1.01$x, summary.2$means, col = "red")
lines(summary.1.01$x, summary.2$q025, lty = 2, col = "red")
lines(summary.1.01$x, summary.2$q975, lty = 2, col = "red")

lines(summary.1.01$x, summary.10$means, col = "blue")
lines(summary.1.01$x, summary.10$q025, lty = 2, col = "blue")
lines(summary.1.01$x, summary.10$q975, lty = 2, col = "blue")

lines(summary.1.01$x, summary.gamma$means, col = "green")
lines(summary.1.01$x, summary.gamma$q025, lty = 2, col = "green")
lines(summary.1.01$x, summary.gamma$q975, lty = 2, col = "green")

legend(x = "topleft", legend = c("p = 1.01", "p = 2", "p = 10", "Gamma(.01, .01)"), 
       col = c("Black", "red", "blue", "green"), lty = 1)
dev.off()


#gelfand and ghosh, posterior predictive replicates for the hierarchical model. 
# in the other file, compute g&g for the standard glm. 
#probably be concerned about numerical stability of the variance computation. I really don't want to 
# deal with that shit. 

post.pred.samples.1.01 = apply(post.samples.1.01$MUS, 2, function(x)rpois(nrow(post.samples.1.01$MUS), x))
post.pred.samples.2 = apply(post.samples.2$MUS, 2, function(x)rpois(nrow(post.samples.1.01$MUS), x))
post.pred.samples.10 = apply(post.samples.10$MUS, 2, function(x)rpois(nrow(post.samples.1.01$MUS), x))
post.pred.samples.gamma = apply(post.samples.gamma$MUS, 2, function(x)rpois(nrow(post.samples.1.01$MUS), x))
post.pred.samples.20 = apply(post.samples.20$MUS, 2, function(x)rpois(nrow(post.samples.1.01$MUS), x))


computeGG = function(post.pred.samples, fabric){
  mu.ls = colMeans(post.pred.samples)
  G = sum((fabric$faults - mu.ls)^2)
  P = sum(apply(post.pred.samples, 2, var))
  return (c(G, P, G+P))
}

GG.results = matrix(0, nrow = 5, ncol = 3)
GG.results[1, ] = computeGG(post.pred.samples.1.01, fabric)
GG.results[2, ] = computeGG(post.pred.samples.2, fabric)
GG.results[3, ] = computeGG(post.pred.samples.10, fabric)
GG.results[4, ] = computeGG(post.pred.samples.gamma, fabric)
GG.results[5, ] = computeGG(post.pred.samples.20, fabric)
GG.results



GGPLOTS = function(fabric, post_pred_replicates){
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
  
}


prior.powers = 5:15
GG.results = matrix(0, nrow = length(prior.powers), ncol = 4)
GG.results[,1] = prior.powers
for (i in 1:length(prior.powers)){
  ps = gibbs(n.iterations = 50000, n.burn = 10000, n.thin = 1, fabric=fabric,
                         lprior = function(x)lambda.prior.p(x, prior.powers[i]))
  post.pred.samples = apply(ps$MUS, 2, function(x)rpois(nrow(ps$MUS), x))
  GG.results[i, ] = c(prior.powers[i], computeGG(post.pred.samples, fabric))
  print(i)
}


