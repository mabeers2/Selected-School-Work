############################################################################
########            PROBLEM 3                       ########################
############################################################################

get_data3 = function(n.samples){
  set.seed(127)
  x = runif(n.samples, min = -100, max = 100)
  df = data.frame(intercept = rep(1, n.samples), x = x, x2 = x^2, x3 = x^3)
  fx = .2 + 2*x + x^2 + 3*x^3
  t = fx +rnorm(n.samples, mean = 0, sd = 1)
  return (cbind(t=t, df))
}

get_data5 = function(n.samples){
  set.seed(127)
  x = runif(n.samples, min = -100, max = 100)
  df = data.frame(intercept = rep(1, n.samples), x = x, x2 = x^2, x3 = x^3, x4 = x^4, x5 = x^5)
  fx = .2 + 2*x + x^2 + 3*x^3
  t = fx +rnorm(n.samples, mean = 0, sd = 1)
  return (cbind(t=t, df))
}


n.samples = 10000
dat = get_data3(n.samples)
fit = lm(t~0+. , data =dat)
coef(fit)
sum((fit$residuals)^2)/n.samples


n.samples = 10000
dat = get_data5(n.samples)
fit = lm(t~0+. , data =dat)
coef(fit)
sum((fit$residuals)^2)/n.samples


# PART C
n.samples = 100
dat = get_data5(n.samples)
Phi = as.matrix(dat[, 2:7])
w.hat = solve(diag(6) + t(Phi)%*% Phi) %*% t(Phi) %*% as.matrix(dat[, 1])


n.samples = 100
dat = get_data3(n.samples)
Phi = as.matrix(dat[, 2:5])
w.hat.100 = solve(diag(4) + t(Phi)%*% Phi) %*% t(Phi) %*% as.matrix(dat[, 1])

n.samples = 1000
dat = get_data3(n.samples)
Phi = as.matrix(dat[, 2:5])
w.hat.1000 = solve(diag(4) + t(Phi)%*% Phi) %*% t(Phi) %*% as.matrix(dat[, 1])

n.samples = 10000
dat = get_data3(n.samples)
Phi = as.matrix(dat[, 2:5])
w.hat.10000 = solve(diag(4) + t(Phi)%*% Phi) %*% t(Phi) %*% as.matrix(dat[, 1])


#3data plot
pdf("third_degree_bayesian_fits.pdf")
x.space = seq(-100, 100, length.out = 5000)
X = cbind(1, x.space, x.space^2, x.space^3)
xw100 = X %*% w.hat.100
xw1000 = X %*% w.hat.1000
xw10000 = X %*% w.hat.10000
true = X %*% as.matrix(c(.2, 2, 1, 3))
plot(x.space, xw100, type = "l", main = "Third Degree Bayesian Linear Regression Fits")
lines(x.space, xw1000, col = "red")
lines(x.space, xw10000, col  = "blue")
lines(x.space, true, col = "green", lty = 2)
legend(x = "topleft", legend = c("100", "1000", "10000", "True"), 
       col = c("black", "red", "blue", "green"), lty = c(1,1, 1, 2))
dev.off()
