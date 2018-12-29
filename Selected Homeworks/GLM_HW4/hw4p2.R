# ams 274 hw 4

####################################################
###        P 2                     @###############
####################################################

#read in data. 
setwd("./documents/school/AMS 274/hw3/hw4")
gator = read.csv('alligator.csv')
gator2a = gator[order(gator$Length), ]
#length in meters
#food choice: I invertebrates, F fish, O other. 
# Gender 0 : male, 1: female

#quick EDA
par(mfrow = c(1,1))
plot(density(gator$Length[gator$Food == "F"], bw = .1), col = "black", ylim = c(0,1.7))
lines(density(gator$Length[gator$Food == "I"], bw = .1), col = "red")
lines(density(gator$Length[gator$Food == "O"], bw = .1), col = "green")

hist(gator$Length[gator$Food == "F"], col=rgb(151,255,255,100,maxColorValue = 255))
hist(gator$Length[gator$Food == "I"], add = TRUE, col=rgb(204,0,102,150, maxColorValue = 255))
hist(gator$Length[gator$Food == "O"], add = TRUE, col=rgb(255,128,0,100, maxColorValue = 255))


pdf("plot1.pdf", width = 8, height = 4)
hist(gator$Length[gator$Food == "F"], border = "darkturquoise", 
     main = "Length by Food Preference", xlab = "Length (m)")
hist(gator$Length[gator$Food == "I"], add = TRUE, border = "red")
hist(gator$Length[gator$Food == "O"], add = TRUE, border = "blue")
legend(x = "topright", legend = c("Fish", "Invertebrates", "Other"), 
       col = c("darkturquoise", "red", "blue"), lty = 1)
dev.off()

# here's the approach you take: get posterior samples using mcmcpack, then do metropolis hastings 
# yourself. 
library(MCMCpack)
post_samples = MCMCmnl(Food ~ Length, baseline = "F", mcmc.method = "IndMH", 
                       mcmc = 100000, thin = 10, tune = 1.0, burnin = 5000, 
                       b0 = 0, B0 = 0, data = gator2a)

post_samples = as.data.frame(post_samples)


get_pi_2a = function(mcmc_results, length){
  v = c(1, length)
  eI = exp(cbind(mcmc_results[, "(Intercept).I"], mcmc_results[, "Length.I"]) %*% v)
  eO = exp(cbind(mcmc_results[, "(Intercept).O"], mcmc_results[, "Length.O"]) %*% v)
  #cat("eI = ", eI, "\n")
  #cat("eO = ", eO, "\n")
  piI = eI/(1 + eI + eO)
  piO = eO/(1 + eI + eO)
  piF = 1/(1 + eI + eO)
  df = data.frame(cbind(piI = piI, piO = piO, piF = piF))
  names(df) = c("piI", "piO", "piF")
  return (df)
}

get_summary = function(v, qs){
  r = c(mean(v), quantile(v, qs))
  return (r)
}

qs = c(.025, .975)
n = nrow(gator)

piI = data.frame(matrix(0, nrow = n, ncol = 3))
names(piI) = c("mean", "q025", "q975")

piO = data.frame(matrix(0, nrow = n, ncol = 3))
names(piO) = c("mean", "q025", "q975")

piF = data.frame(matrix(0, nrow = n, ncol = 3))
names(piF) = c("mean", "q025", "q975")

for (i in 1:n){
  pis = get_pi_2a(post_samples, gator2a$Length[i])
  piI[i, ] = get_summary(pis$piI, qs = qs)
  piO[i, ] = get_summary(pis$piO, qs = qs)
  piF[i, ] = get_summary(pis$piF, qs = qs)
}

pdf("plot2.pdf", width = 8, height = 8)
par(mfrow = c(3,1), mar = c(4,4,4,4))
eps = .05

#plot piI
plot(gator2a$Length, piI$mean, type = "l", xlab = "Length (m)", ylab = "pi(x)", main = "Invertebrates",
     ylim = c(min(piI$q025) - eps, max(piI$q975) + eps))
lines(gator2a$Length, piI$q025, col = "red", lty = 2)
lines(gator2a$Length, piI$q975, col = "red", lty = 2)
#plot piO
plot(gator2a$Length, piO$mean, type = "l", xlab = "Length (m)", ylab = "pi(x)", main = "Other",
     ylim = c(min(piO$q025) - eps, max(piO$q975) + eps))
lines(gator2a$Length, piO$q025, col = "red", lty = 2)
lines(gator2a$Length, piO$q975, col = "red", lty = 2)
#plot piF
plot(gator2a$Length, piF$mean, type = "l", xlab = "Length (m)", ylab = "pi(x)", main = "Fish", 
     ylim = c(min(piF$q025) - eps, max(piF$q975) + eps))
lines(gator2a$Length, piF$q025, col = "red", lty = 2)
lines(gator2a$Length, piF$q975, col = "red", lty = 2)
dev.off()


################################################################################################
#alright so you have something functional with mcmc pack. Lets see if you can get something functional
# with a custom MH algorithm.  Good luck 
################################################################################################
library(MASS)

ll2 = function(coef, Lengths, food){
  n = length(Lengths)
  S = 0
  #say O corresponds to positions one and two in coef vector
  for (i in 1:n){
    if (food[i] == "O"){
      #print ((coef[1] + coef[2]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S + (coef[1] + coef[2]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i]))
    } else if (food[i] == "I"){
      #print((coef[3] + coef[4]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S + (coef[3] + coef[4]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i]))
    } else {
      #print(- log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i]))
    }
  }
  return (S)
}

get_pi = function(coef, l){
  X = cbind(1, l)
  v1 = coef[1:2]
  v2 = coef[3:4]
  e1 = exp(X %*% v1)
  e2 = exp(X %*% v2)
  den = 1 + e1 + e2
  num = c(1, e1, e2)
  return (num/den)
}

ll3 = function(coef, Lengths, food){
  M = model.matrix(~0 + food)
  n = length(Lengths)
  S = numeric(n)
  for (i in 1:n){
    #print(get_pi(coef, Lengths[i]))
    S[i] = dmultinom(M[i, ], prob = get_pi(coef, Lengths[i]), log = T)
  }
  #print (S)
  return (sum(S))
}

ll2(c(1,1,1,1), gator$Length, gator$Food)
ll3(c(1,1,1,1), gator$Length, gator$Food)

colMeans(post_samples)
maxlik = optim(c(1,1,1,1), function(x) ll2(x, gator$Length, gator$Food), control = list(fnscale = -1), hessian=T)
# so your log likelihood function appears to be working properly. That's good. 

n.iterations = 10000
n.burn = 5000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)
coefs = matrix(3, nrow = n.iterations, ncol = 4)
d = 1

system.time(for (i in 1:(n.iterations-1)){
  
  #generate candidates
  candidate_vec = mvrnorm(n = 1, mu = c(coefs[i, ]), Sigma = -1*solve(maxlik$hessian))
  
  #calculate log posterior at candidate point and at current point. 
  ll_current = ll2(coefs[i, ], gator$Length, gator$Food)
  ll_proposal = ll2(candidate_vec, gator$Length, gator$Food)

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
  #print (i)
}
)


par(mfrow = c(1,1))
coefs = coefs[keeps, ]
ts.plot(coefs[, 1])
ts.plot(coefs[, 2])
ts.plot(coefs[, 3])
ts.plot(coefs[, 4])
colMeans(coefs)
colMeans(post_samples)


mcmc_plots = function(v, t, MLE, t2 = t){
  par(mfrow = c(1,2), mar=c(5.1,4.1,4.1,2.1))
  m = mean(v)
  ts.plot(v, ylab = t, main = t)
  hist(v, main = t2, xlab = t2, sub = paste("Mean = ", round(m, 4)))
  abline(v = m, col = "blue")
  abline(v = MLE, col = "green")
  abline(v = quantile(v, c(.025, .975)), col = "red", lty = 2)
}


pdf("plot3.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 1], "Intercept: Other", maxlik$par[1])
dev.off()
pdf("plot4.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 2], "Coefficient: Other", maxlik$par[2])
dev.off()
pdf("plot5.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 3], "Intercept: Invertebrates", maxlik$par[3])
dev.off()
pdf("plot6.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 4], "Coefficient: Invertebrates", maxlik$par[4])
dev.off()

############################################################################################

get_pi_2a = function(mcmc_results, length){
  v = c(1, length)
  eI = exp(cbind(mcmc_results[, 3], mcmc_results[, 4]) %*% v)
  eO = exp(cbind(mcmc_results[, 1], mcmc_results[, 2]) %*% v)
  #cat("eI = ", eI, "\n")
  #cat("eO = ", eO, "\n")
  piI = eI/(1 + eI + eO)
  piO = eO/(1 + eI + eO)
  piF = 1/(1 + eI + eO)
  df = data.frame(cbind(piI = piI, piO = piO, piF = piF))
  names(df) = c("piI", "piO", "piF")
  return (df)
}

get_summary = function(v, qs){
  r = c(mean(v), quantile(v, qs))
  return (r)
}

qs = c(.025, .975)
n = nrow(gator)

piI = data.frame(matrix(0, nrow = n, ncol = 3))
names(piI) = c("mean", "q025", "q975")

piO = data.frame(matrix(0, nrow = n, ncol = 3))
names(piO) = c("mean", "q025", "q975")

piF = data.frame(matrix(0, nrow = n, ncol = 3))
names(piF) = c("mean", "q025", "q975")

for (i in 1:n){
  pis = get_pi_2a(coefs, gator2a$Length[i])
  piI[i, ] = get_summary(pis$piI, qs = qs)
  piO[i, ] = get_summary(pis$piO, qs = qs)
  piF[i, ] = get_summary(pis$piF, qs = qs)
}

pdf("plot2.pdf", width = 8, height = 8)
par(mfrow = c(3,1), mar = c(4,4,4,4))
eps = .05

#plot piI
plot(gator2a$Length, piI$mean, type = "l", xlab = "Length (m)", ylab = "pi(x)", main = "Invertebrates",
     ylim = c(0, 1))
lines(gator2a$Length, piI$q025, col = "red", lty = 2)
lines(gator2a$Length, piI$q975, col = "red", lty = 2)
#plot piO
plot(gator2a$Length, piO$mean, type = "l", xlab = "Length (m)", ylab = "pi(x)", main = "Other",
     ylim = c(0, 1))
lines(gator2a$Length, piO$q025, col = "red", lty = 2)
lines(gator2a$Length, piO$q975, col = "red", lty = 2)
#plot piF
plot(gator2a$Length, piF$mean, type = "l", xlab = "Length (m)", ylab = "pi(x)", main = "Fish", 
     ylim = c(0, 1))
lines(gator2a$Length, piF$q025, col = "red", lty = 2)
lines(gator2a$Length, piF$q975, col = "red", lty = 2)
dev.off()


#######################################################################################
###### part b                 ######################################################
####################################################################


ll4 = function(coef, Lengths, food, gender){
  n = length(Lengths)
  S = 0
  #say O corresponds to positions one and two in coef vector
  for (i in 1:n){
    if (food[i] == "O"){
      #print ((coef[1] + coef[2]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S + (coef[1] + coef[2]*Lengths[i] + coef[5]*gender[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]+ coef[5]*gender[i]) + exp(coef[3] + coef[4]*Lengths[i]+ coef[5]*gender[i]))
    } else if (food[i] == "I"){
      #print((coef[3] + coef[4]*Lengths[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S + (coef[3] + coef[4]*Lengths[i]+ coef[5]*gender[i]) - log(1 + exp(coef[1] + coef[2]*Lengths[i]+ coef[5]*gender[i]) + exp(coef[3] + coef[4]*Lengths[i]+ coef[5]*gender[i]))
    } else {
      #print(- log(1 + exp(coef[1] + coef[2]*Lengths[i]) + exp(coef[3] + coef[4]*Lengths[i])))
      S = S - log(1 + exp(coef[1] + coef[2]*Lengths[i]+ coef[5]*gender[i]) + exp(coef[3] + coef[4]*Lengths[i]+ coef[5]*gender[i]))
    }
  }
  return (S)
}

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

ll5.optim = optim(c(1,1,1,1,1,1), function(x) ll5(x, gator$Length, gator$Food, gator$Gender), 
               control = list(fnscale = -1), hessian=T)
ll5.optim$par



maxlik = optim(rep(0,5), function(x) ll4(x, gator$Length, gator$Food, gator$Gender), 
               control = list(fnscale = -1), hessian=T)
maxlik$par



n.iterations = 10000
n.burn = 5000
n.thin = 1
keeps = seq(n.burn, n.iterations, by = n.thin)
coefs = matrix(3, nrow = n.iterations, ncol = 5)
d = 1

system.time(for (i in 1:(n.iterations-1)){
  
  #generate candidates
  candidate_vec = mvrnorm(n = 1, mu = c(coefs[i, ]), Sigma = -1*solve(maxlik$hessian))
  
  #calculate log posterior at candidate point and at current point. 
  ll_current = ll4(coefs[i, ], gator$Length, gator$Food, gator$Gender)
  ll_proposal = ll4(candidate_vec, gator$Length, gator$Food, gator$Gender)
  
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
mcmc_plots(coefs[, 1], "Intercept: Other", maxlik$par[1])
dev.off()
pdf("plot8.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 2], "Coefficient: Other", maxlik$par[2])
dev.off()
pdf("plot9.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 3], "Intercept: Invertebrates", maxlik$par[3])
dev.off()
pdf("plot10.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 4], "Coefficient: Invertebrates", maxlik$par[4])
dev.off()
pdf("plot11.pdf", width = 8, height = 4)
mcmc_plots(coefs[, 5], "Coefficient: Gender", maxlik$par[5])
dev.off()






# now get response probabilities from the appropriate mcmc samples. 
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
    e1 = exp(mcmc_samples[, c(1,2,5)] %*% xi)
    e2 = exp(mcmc_samples[, c(3,4,5)] %*% xi)
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











