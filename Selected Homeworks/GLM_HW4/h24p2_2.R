ll5 = function(coef, Lengths, food, gender){
  n = length(Lengths)
  S = 0
  #say O corresponds to positions one and two in coef vector
  for (i in 1:n){
    v1 = coef[1] + coef[2]*Lengths[i] + coef[3]*gender[i]
    v2 = coef[4] + coef[5]*Lengths[i] + coef[6]*gender[i]
    if (food[i] == "O"){
      #print ((coef[1] + coef[2]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S + v1 - log(1 + exp(v1) + exp(v2))
    } else if (food[i] == "I"){
      #print((coef[3] + coef[4]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S + v2 - log(1 + exp(v1) + exp(v2))
    } else {
      #print(- log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S - log(1 + exp(v1) + exp(v2))
    }
  }
  return (S)
}

ll5.optim = optim(rep(1, 6), function(x) ll5(x, gator$Length, gator$Food, gator$Gender), 
                  control = list(fnscale = -1), hessian=T, method = "BFGS")
ll5.optim$par



n.iterations = 20000
n.burn = 5000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)
coefs = matrix(3, nrow = n.iterations, ncol = 6)
d = 1

system.time(for (i in 1:(n.iterations-1)){
  
  #generate candidates
  candidate_vec = mvrnorm(n = 1, mu = c(coefs[i, ]), Sigma = -1*solve(ll5.optim$hessian))
  
  #calculate log posterior at candidate point and at current point. 
  ll_current = ll5(coefs[i, ], gator$Length, gator$Food, gator$Gender)
  ll_proposal = ll5(candidate_vec, gator$Length, gator$Food, gator$Gender)
  
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
)


par(mfrow = c(1,1))
coefs = coefs[keeps, ]
ts.plot(coefs[, 1])
ts.plot(coefs[, 2])
ts.plot(coefs[, 3])
ts.plot(coefs[, 4])
ts.plot(coefs[, 5])
colMeans(coefs)
maxlik$par


pdf("plot7.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 1], "Intercept: Other", ll5.optim$par[1])
dev.off()

pdf("plot8.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 2], "Coefficient: Other", ll5.optim$par[2])
dev.off()

pdf("plot211.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 3], "Coefficient: Other, Gender", ll5.optim$par[3])
dev.off()

pdf("plot9.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 4], "Intercept: Invertebrates", ll5.optim$par[4])
dev.off()

pdf("plot10.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 5], "Coefficient: Invertebrates", ll5.optim$par[5])
dev.off()

pdf("plot11.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 6], "Coefficient: Invertebrate, Gender", ll5.optim$par[6])
dev.off()

#coefficient invertebrate -3.279318
#intercept invertebrate 6.387922 
#-0.1318982 -1.151515 -1.32701 0.2075341




length.grid = seq(min(gator$Length), max(gator$Length), length.out = 200)
x.male = cbind(1, length.grid, 0)
x.female = cbind(1, length.grid, 1)
#what's the difference. The difference is in the numerator and denominator. We have to pass in 
# coefs and somehow specify which columns to use. 

plot_2b_helper = function(x, mins, means, maxs, t){
  plot(x, means, type = "l", ylim = c(0,1), 
       ylab = "Pi(x)", main = t, xlab = "Length (m)")
  lines(x, mins, col= "red", lty=2)
  lines(x, maxs, col = "red", lty = 2)
}

plot_2b = function(mcmc_samples, X, mf){
  n.grid.points = nrow(X)
  summary_df = matrix(0,nrow = n.grid.points, ncol = 9)
  for (i in 1:n.grid.points){
    xi = X[i, ]
    e1 = exp(mcmc_samples[, c(1,2,3)] %*% xi)
    e2 = exp(mcmc_samples[, c(4,5,6)] %*% xi)
    den = 1 + e1 + e2
    p1 = e1/den
    p2 = e2/den
    p3 = 1/den
    summary_df[i, ] = c(quantile(p1, .025), 
                        quantile(p2, .025),
                        quantile(p3, .025),
                        mean(p1), 
                        mean(p2), 
                        mean(p3), 
                        quantile(p1, .975), 
                        quantile(p2, .975),
                        quantile(p3, .975))
  }
  
  plot_2b_helper(X[,2], summary_df[, 1], summary_df[, 4], summary_df[, 7], paste("Other: ", mf))
  plot_2b_helper(X[,2], summary_df[, 2], summary_df[, 5], summary_df[, 8], paste("Invertebrates: ", mf))
  plot_2b_helper(X[,2], summary_df[, 3], summary_df[, 6], summary_df[, 9], paste("Fish: ", mf))
}

pdf("plot12.pdf", width = 8, height = 8)
layout(matrix(1:6, 3, 2, byrow = FALSE))
plot_2b(coefs, x.male, "Male")
plot_2b(coefs, x.female, "Female")
dev.off()

layout(matrix(1:6, 3, 2, byrow = FALSE))
plot_2b(coefs, x.male, "Male")
plot_2b(coefs, x.female, "Female")


