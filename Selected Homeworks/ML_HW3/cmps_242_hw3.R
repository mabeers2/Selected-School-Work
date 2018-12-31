# cmps 242 hw3

####################################################
###         PROBLEM 2               ################
####################################################

#PART A
setwd("./documents/school/cmps242")
#define your data and labels
m = matrix(0, ncol = 15, nrow = 300)
for (i in 1:300){
  m[i, ] = sample(x = c(-1, 1), size = 15, replace = TRUE, prob = c(.5, .5))
}

t = m[, 1]

#implement the perceptron algorithm. 
w = rep(0,15)
epochs = 0
prediction_errors = 0
another_iteration_required = T

while (another_iteration_required){
  another_iteration_required = F
  for (i in 1:300){
    eta.i = 1
    if (sum(w * m[i, ])*t[i] <= 0) {
      #mistake made
      w = w + eta.i*t[i]*m[i, ]
      another_iteration_required =  T
      prediction_errors = prediction_errors + 1
    }
  }
  epochs = epochs + 1
  cat("total epochs = ", epochs, "\n")
  
}
w
prediction_errors
#################################################################################
# PART B
t = sign(rowSums(m))
table(t)
#implement the perceptron algorithm. 
w = rep(0,15)
epochs = 0
prediction_errors = 0
another_iteration_required = T

while (another_iteration_required){
  another_iteration_required = F
  for (i in 1:300){
    eta.i = 1
    if (sum(w * m[i, ])*t[i] <= 0) {
      #mistake made
      w = w + eta.i*t[i]*m[i, ]
      another_iteration_required =  T
      prediction_errors = prediction_errors + 1
    }
  }
  epochs = epochs + 1
  cat("total epochs = ", epochs, "\n")
  
}
w
prediction_errors


#############################################################
# PART C

#run 3 epochs
m.c = matrix(0, ncol = 15, nrow = 300)
for (i in 1:300){
  m.c[i, ] = sample(x = c(-1, 1), size = 15, replace = TRUE, prob = c(.5, .5))
}

#generate labels in this exceptionally complicated way. 
t = numeric(300)
for (i in 1:300){
  r = sample(-4:4, size = 1)
  q = r + sum(m.c[i, 1:11])
  if (q <= 0){
    t[i] = -1
  }else{
    t[i] = 1
  }
}


w = rep(0, 15)
w.matrix = matrix(0, ncol = 15, nrow = 900)
eta.i = 1
for (ep in 1:3){
  for (i in 1:300){
    if (sum(w * m.c[i, ])*t[i] <= 0) {
      #mistake made
      w = w + eta.i*t[i]*m.c[i, ]
    }
    w.matrix[300*(ep - 1) + i, ] = w
  }
  #shuffle m.c
  m.c = m.c[sample(1:300, size = 300, replace = FALSE), ]
}


#generate test set + random noise on labels. 
m.test = matrix(0, ncol = 15, nrow = 500)
for (i in 1:500){
  m.test[i, ] = sample(x = c(-1, 1), size = 15, replace = TRUE, prob = c(.5, .5))
}

t.test = numeric(500)
for (i in 1:500){
  r = sample(-4:4, size = 1)
  q = r + sum(m.test[i, 1:11])
  if (q <= 0){
    t.test[i] = -1
  }else{
    t.test[i] = 1
  }
}


#5 hypotheses

#1
last_hypothesis_accuracy = function(w.matrix, x, t){
   w = w.matrix[900, ]
   xw = x %*% w #500x15 X 15x1
   predictions = sapply(xw, function(x){ifelse(x <= 0, -1, 1)})
   return(mean(predictions == t))
}
last_hypothesis_accuracy(w.matrix, m.test, t.test)

#2
average_hypothesis_accuracy = function(w.matrix, x, t){
   w = colMeans(w.matrix)
   xw = x %*% w #500x15 X 15x1
   predictions = sapply(xw, function(x){ifelse(x <= 0, -1, 1)})
   return(mean(predictions == t)) 
}
average_hypothesis_accuracy(w.matrix, m.test, t.test)

#3
voted_hypothesis_accuracy = function(w.matrix, x, t){
   pred_matrix = matrix(0, nrow = nrow(w.matrix), ncol = nrow(x))
   for (i in 1:nrow(w.matrix)){
     w = w.matrix[i, ]
     xw = x %*% w #500x15 X 15x1
     predictions = sapply(xw, function(x){ifelse(x <= 0, -1, 1)})
     pred_matrix[i, ] = predictions
   }
   p = colSums(pred_matrix)
   predictions = sapply(p, function(x){ifelse(x <= 0, -1, 1)})
   return (mean(predictions == t))
}
voted_hypothesis_accuracy(w.matrix, m.test, t.test)
#4
last_epoch_avg_accuracy = function(w.matrix, x, t){
  w = colMeans(w.matrix[601:900, ])
  xw = x %*% w #500x15 X 15x1
  predictions = sapply(xw, function(x){ifelse(x <= 0, -1, 1)})
  return(mean(predictions == t)) 
}
last_epoch_avg_accuracy(w.matrix, m.test, t.test)
#5
last_epoch_vote_accuracy = function(w.matrix, x, t){
  pred_matrix = matrix(0, nrow = 300, ncol = nrow(x))
  for (i in 601:900){
    w = w.matrix[i, ]
    xw = x %*% w #500x15 X 15x1
    predictions = sapply(xw, function(x){ifelse(x <= 0, -1, 1)})
    pred_matrix[i-600, ] = predictions
  }
  p = colSums(pred_matrix)
  predictions = sapply(p, function(x){ifelse(x <= 0, -1, 1)})
  return (mean(predictions == t))
}
last_epoch_vote_accuracy(w.matrix, m.test, t.test)


#repeat experiment 10 times 
get_m = function(n){
  m.test = matrix(0, ncol = 15, nrow = n)
  for (i in 1:n){
    m.test[i, ] = sample(x = c(-1, 1), size = 15, replace = TRUE, prob = c(.5, .5))
  }
  return(m.test)
}

get_t = function(m.test){
  n = nrow(m.test)
  t.test = numeric(n)
  for (i in 1:n){
    r = sample(-4:4, size = 1)
    q = r + sum(m.test[i, 1:11])
    if (q <= 0){
      t.test[i] = -1
    }else{
      t.test[i] = 1
    }
  }
  return (t.test)
}

get_w.matrix = function(m.c, t){
  w = rep(0, 15)
  w.matrix = matrix(0, ncol = 15, nrow = 900)
  eta.i = 1
  for (ep in 1:3){
    for (i in 1:300){
      if (sum(w * m.c[i, ])*t[i] <= 0) {
        #mistake made
        w = w + eta.i*t[i]*m.c[i, ]
      }
      w.matrix[300*(ep - 1) + i, ] = w
    }
    #shuffle m.c
    m.c = m.c[sample(1:300, size = 300, replace = FALSE), ]
  }
  return(w.matrix)
}

n.iter = 10
accuracy_matrix = matrix(0,ncol = 5, nrow = n.iter)
for (i in 1:n.iter){
  m.train = get_m(300)
  m.test = get_m(500)
  t.train = get_t(m.train)
  t.test = get_m(m.test)
  w.matrix = get_w.matrix(m.train, t.train)
  accuracy_matrix[i, 1] = last_hypothesis_accuracy(w.matrix, m.test, t.test)
  accuracy_matrix[i, 2] = average_hypothesis_accuracy(w.matrix, m.test, t.test)
  accuracy_matrix[i, 3] = voted_hypothesis_accuracy(w.matrix, m.test, t.test)
  accuracy_matrix[i, 4] = last_epoch_avg_accuracy(w.matrix, m.test, t.test)
  accuracy_matrix[i, 5] = last_epoch_vote_accuracy(w.matrix, m.test, t.test)
  cat("iteration ", i, " complete\n")
}

colMeans(accuracy_matrix)




#switch over from python was that horrible or what. 
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


#5 data plot
# so here we can't directly invert the matrix so we have to use some kind of gradient descent approach.
# this is completed in python. 


ham_dir = paste0("./enron", 1:5, "/ham")
spam_dir = paste0("./enron", 1:5, "/spam")
ham_files = c()
spam_files = c()

for (dir in ham_dir){
  ham_files = c(ham_files, paste0(dir, "/", list.files(dir)))
}

for (dir in spam_dir){
  spam_files = c(spam_files, paste0(dir, "/", list.files(dir)))
}

gsub("\xa9", "", "ajsdklfjdskl\xa9fjkewlhf")






