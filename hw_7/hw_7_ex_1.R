###########################################################
# Course:       Modern econometric and statistical learning
#               methods forquantitative asset management
#
# Instructor:   Prof. Dr. Marc Paolella,
#               Urban Ulrych,
#               Simon Hediger,
#               University of Zurich
#
# Author:       Rino Beeli
#
# Date:         May 25th, 2020
# 
# Topic:        Homework 7 - Exercise 1.
###########################################################

library(rpart)        # regression trees
library(rpart.plot)
library(ranger)       # random forests
library(MASS)         # multivariate Normal distribution
library(glmnet)       # LASSO regression

rm(list=ls())
set.seed(42)



# -------------
# parameters
# -------------
S <- 100  # number of simulations
N <- 500  # T training samples from AR(2) process
p <- 100  # number of predictor columns in D



# -------------
# functions
# -------------

sim.AR2 <- function(N) {
  # Simulates N observations of a stationary AR(2) process given by:
  #   X_t = 0.3 X_{t-1} - 0.4 X_{t-2} + e_t
  # with
  #   e_t ~ iid N(0,1)
  #   X_0 = X_1 = 0
  e <- rnorm(N, mean=0, sd=1)
  X <- rep(0, N)
  for (t in 3:N) {
    X[t] <- 0.3*X[t-1] - 0.4*X[t-2] + e[t]
  }
  return(X)
}

make.D <- function(X, p) {
  # Creates the data matrix D as specified
  # in exercise 1. First column represents
  # the response variable X_{p+1}, ..., X_{T+1}.
  D <- matrix(0, nrow=length(X)-p, ncol=p+1)
  for (t in 1:nrow(D)) {
    D[t,] <- X[(p+t):t]
  }
  return(D)
}

# verify AR(2) simulation function by plotting PACF
par(mfrow=c(1,1), mar=c(5, 5, 4, 3))
pacf(sim.AR2(5000), lag.max=10, main='PACF simulated AR(2) process with N=5000')



# -------------
# simulation
# -------------

# prediction errors vectors
err.rnd <- rep(0, S)
err.ar <- rep(0, S)
err.rtree <- rep(0, S)
err.rtree.pruned <- rep(0, S)
err.rforest <- rep(0, S)
err.rforest.oob <- rep(0, S)
err.rforest.cv <- rep(0, S)
err.ols <- rep(0, S)
err.lasso <- rep(0, S)

# e) keep track if first two variables have been selected
#    by the variable importance measure of the random forest
rforest.var.imp.sel <- matrix(0, nrow=S, ncol=2)

for (s in 1:S) {
  cat(sprintf('Simulation %i of %i \n', s, S))
  
  # simulate AR(2) process
  X <- sim.AR2(N+1)
  X.train <- X[-(N+1)]
  X.test <- X[N+1]
  
  # data matrix D
  D <- as.data.frame(make.D(X, p))
  colnames(D) <- c('Y', paste0('X', 1:p))
  
  # split into training (n=500) and test data (n=1)
  D.train <- D[-nrow(D), ]
  D.test <- D[nrow(D), ]
  
  # squared test error function
  test.serr <- function(pred) (D.test[1,'Y'] - pred)^2
  
  
  
  # ---------------- AR(2) benchmark ----------------
  ar.fit <- ar(X.train, aic=F, method='mle', order.max=2)
  ar.pred <- predict(ar.fit, newdata=X.train[(N-1):N], n.ahead=1)$pred[1]
  err.ar[s] <- test.serr(ar.pred)
  
  
  # ---------- Gaussian random prediction -----------
  err.rnd[s] <- test.serr(rnorm(1, mean=mean(X), sd=sd(X)))
  
  
  # ---------------- regression tree ----------------
  rtree.fit <- rpart(Y ~ ., data=D.train, method='anova', control=rpart.control(cp=0))
  
  # b) prune the tree based on minimum cost complexity criterion
  rtree.fit.pruned <- prune(rtree.fit, cp=rtree.fit$cptable[which.min(rtree.fit$cptable[,'xerror']),'CP'])
  
  if (s == 1) {
    # plot regression trees
    prp(rtree.fit, type=2, extra=1)
    prp(rtree.fit.pruned, type=2, extra=1)
  }
  
  # uses last row of D for test prediction
  rtree.pred <- predict(rtree.fit, newdata=D.test)
  rtree.pruned.pred <- predict(rtree.fit.pruned, newdata=D.test)
  err.rtree[s] <- test.serr(rtree.pred)
  err.rtree.pruned[s] <- test.serr(rtree.pruned.pred)
  
  
  # ----------------- random forest -----------------
  rforest.fit <- ranger(Y ~ ., data=D.train, mtry=floor(sqrt(p)), num.trees=100, importance='impurity')
  rforest.pred <- predict(rforest.fit, data=D.test)$predictions[1]
  err.rforest[s] <- test.serr(rforest.pred)
  
  # c) out-of-bag error
  err.rforest.oob[s] <- rforest.fit$prediction.error
  
  # d) BONUS: compare the OOB-error with the leave-one-out CV error
  tmp.rforest.cv <- numeric()
  for(i in 300:400) {
    i.rffit <- ranger(Y ~ ., data=D[(i-300+1):i, ], mtry=floor(sqrt(p)), num.trees=100)
    i.pred <- predict(i.rffit, data=D[i+1, ])
    tmp.rforest.cv[i] <- (D[i+1, 'Y'] - i.pred$predictions[1])^2
  }
  err.rforest.cv[s] <- mean(tmp.rforest.cv, na.rm=T)
  
  # e) selected variables by variable importance measure
  rforest.var.imp.sel[s,] <- names(sort(rforest.fit$variable.importance, decreasing=T))[1:2] %in% names(D)[2:3]
  
  
  # --------------- OLS regression ----------------
  ols.fit <- lm(Y ~ -1 + ., data=D.train)
  ols.pred <- predict(ols.fit, newdata=D.test[2:ncol(D)])
  err.ols[s] <- test.serr(ols.pred)
  if (s == 1) {
    print(summary(ols.fit))
  }
  
  
  # -------------- LASSO regression ---------------
  # (not part of exercise)
  lasso.fit <- cv.glmnet(as.matrix(D.train[, 2:ncol(D)]), as.matrix(D.train[,1]), alpha=1, intercept=F)
  lasso.pred <- predict(lasso.fit, newx=as.matrix(D.test[2:ncol(D)]), s='lambda.1se')
  err.lasso[s] <- test.serr(lasso.pred)
}



# print mean squared errors
cat(sprintf('
   MSE Gaussian random:         %.3f
   MSE AR(2) fitted:            %.3f
   MSE regression tree:         %.3f
   MSE regression tree pruned:  %.3f
   MSE random forest:           %.3f (OOB=%.3f, %.0f%% of times greater than MSE)
   MSE OLS:                     %.3f
   MSE LASSO:                   %.3f
', mean(err.rnd), mean(err.ar), mean(err.rtree),
   mean(err.rtree.pruned), mean(err.rforest),
   mean(err.rforest.oob), 100*mean(err.rforest.oob > err.rforest),
   mean(err.ols), mean(err.lasso)))

# barplot of MSE by method
err.df <- data.frame(Method=c(
  'Gaussian random', 'AR(2) fitted', 'Regression tree',
  'Regression tree pruned', 'Random forest', 'OLS', 'LASSO'),
  MSE=c(mean(err.rnd), mean(err.ar), mean(err.rtree),
       mean(err.rtree.pruned), mean(err.rforest),
       mean(err.ols), mean(err.lasso)))
err.df <- err.df[order(err.df$MSE), ] # sort by MSE
par(mar=c(5, 10, 5, 2))
b <- barplot(err.df$MSE, names.arg=err.df$Method, main='MSE test prediction error',
        horiz=T, xlab='MSE',xlim=c(0, max(err.df$MSE)*1.4), las=2)
text(x=err.df$MSE + 0.1*max(err.df$MSE), y=b, label=format(err.df$MSE, digits=2), cex=0.75)


# e) show how many times variable of lag 1 and lag 2 have been selected
#    by the variable importance measure of the random forest
cat(sprintf('
   Selected lag 1 and 2:  %.0f%%
   Selected lag 1:        %.0f%%
   Selected lag 2:        %.0f%%
', 100*mean(rowSums(rforest.var.imp.sel) == 2),
   100*mean(rforest.var.imp.sel[,1]),
   100*mean(rforest.var.imp.sel[,2])))


cat('Random forest error summary: \n')
summary(err.rforest)


# d) BONUS: compare the OOB-error with the leave-one-out CV error
cat('Random forest OOB summary: \n')
summary(err.rforest.oob)

cat('Random forest CV summary: \n')
summary(err.rforest.cv)




