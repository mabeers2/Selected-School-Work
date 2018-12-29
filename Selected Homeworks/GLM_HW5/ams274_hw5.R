# ams 274 hw5 
library(MASS)
setwd("./documents/school/ams 274/hw5")
fabric = read.csv("fabric.csv")


ll.poisson = function(v, y, x){
  l = exp(v[1] + v[2]*x)
  return (sum(dpois(x = y, lambda = l, log = TRUE)))
}


optim.poisson = optim(c(0,0), function(v) ll.poisson(v, fabric$faults, fabric$length), 
                      control = list(fnscale = -1), method = "BFGS", hessian = TRUE)
  
optim.poisson
#blah= glm(faults~length, data=fabric, family = poisson(link = "log"))

n.iterations= 50000
n.burn = 10000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)

coefs = matrix(0,ncol = 2, nrow = n.iterations)

for (i in 1:(n.iterations-1)){
  
  #generate candidates
  candidate_vec = mvrnorm(n = 1, mu = coefs[i, ], Sigma = -1*solve(optim.poisson$hessian))
  
  #calculate log posterior at candidate point and at current point. 
  #ll1(x, mouse$dead, mouse$total, mouse$concentration)
  ll_current = ll.poisson(coefs[i, ], fabric$faults, fabric$length)
  ll_proposal = ll.poisson(candidate_vec, fabric$faults, fabric$length)
  
  priors_current = 0#sum(dnorm(coefs[i, ], 0, 100, log = TRUE))
  priors_proposal = 0#sum(dnorm(candidate_vec, 0, 100, log = TRUE))
  
  log_post_current = ll_current + priors_current
  log_post_proposal = ll_proposal + priors_proposal
  
  q = min(1, exp(log_post_proposal - log_post_current))
  
  if (runif(1) < q){
    coefs[i+1, ] = candidate_vec
  } else{
    coefs[i + 1, ] = coefs[i, ]
  }
  #print (i)
}


par(mfrow = c(1,1))
coefs.ll1 = coefs[keeps, ]
ts.plot(coefs.ll1[, 1])
ts.plot(coefs.ll1[, 2])
colMeans(coefs.ll1)
optim.poisson$par


mcmc_plots = function(v, t, MLE, t2 = t){
  par(mfrow = c(1,2), mar=c(5.1,4.1,4.1,2.1))
  m = mean(v)
  ts.plot(v, ylab = t, main = t)
  hist(v, main = t2, xlab = t2, sub = paste("Mean = ", round(m, 4)))
  abline(v = m, col = "blue")
  abline(v = MLE, col = "green")
  abline(v = quantile(v, c(.025, .975)), col = "red", lty = 2)
}

pdf("plot1.pdf", width = 8, height = 4)
mcmc_plots(coefs.ll1[, 1], "beta1", optim.poisson$par[1])
dev.off()
pdf("plot2.pdf", width = 8, height = 4)
mcmc_plots(coefs.ll1[, 2], "beta2", optim.poisson$par[2])
dev.off()


# point and interval estimates for response mean as a function of the covariate. 
lengths.grid = seq(100, 1000, length.out = 1000)
fill = matrix(0, nrow = length(lengths.grid), ncol = 3)
for (i in 1:length(lengths.grid)){
  responses = exp(coefs.ll1[,1] + coefs.ll1[,2]*lengths.grid[i])
  fill[i, 1] = quantile(responses, 0.025)
  fill[i, 2] = mean(responses)
  fill[i, 3] = quantile(responses, 0.975) 
}

pdf("plot3.pdf", width = 8, height = 4)
par(mfrow =  c(1,1))
plot(lengths.grid, fill[, 2], type = "l", main = "Response Mean", 
     xlab = "Length of Roll", ylab = "Number of Faults", ylim = c(0,30))
lines(lengths.grid, fill[,1], col = "red", lty = 2)
lines(lengths.grid, fill[,3], col = "red", lty = 2)
points(fabric$length, fabric$faults, col = "blue")
dev.off()


#posterior predictive residuals, check to make sure this is right. 
pdf("plot4.pdf", width = 8, height = 4)
n = nrow(fabric)
ppr = matrix(0, nrow = n, ncol = 2)
for (i in 1:n){
  mus = exp(coefs.ll1[,1] + coefs.ll1[,2]*fabric[i, "length"])
  mus2 = sapply(mus, function(x)rpois(1,x))
  ris = fabric$faults[i] - mus2
  ppr[i, 1] = quantile(ris, 0.025)
  ppr[i, 2] = quantile(ris, 0.975) 
}

#residuals plot, ok this can't be right. Can it? 
plot(1, type="n", xlab="Length of Roll", ylab="Posterior Predictive Residuals", 
     xlim=c(100, 1000), ylim=c(-15, 20), main = "Posterior Predictive Residuals")
abline(h=0)
for (j in 1:n){
  if (ppr[j,1] <0 & ppr[j,2] >0){
    c = "red"
  } else{
    c = "blue"
  }
  segments(x0 = fabric$length[j], y0 = ppr[j, 1], y1 = ppr[j, 2], col = c)
}
dev.off()




########################################################################
########################################################################

n = nrow(fabric)
ppr = matrix(0, nrow = n, ncol = 2)
G = 0
P = 0
for (i in 1:n){
  mus = exp(coefs.ll1[,1] + coefs.ll1[,2]*fabric[i, "length"])
  mus2 = sapply(mus, function(x)rpois(1,x))
  G =  G + (fabric$faults[i] - mean(mus2))^2
  P = P + var(mus2)
}
G 
P
G+P











