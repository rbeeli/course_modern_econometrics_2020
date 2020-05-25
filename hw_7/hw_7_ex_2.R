############################################################
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
# Topic:        Homework 7 - Exercise 2.
###########################################################

library(rpart)        # regression trees
library(rpart.plot)
library(ranger)       # random forests
library(MASS)         # simulate Multivariate Normal distribution
library(glmnet)       # LASSO regression
library(foreach)      # parallel processing
library(doParallel)   # parallel processing

rm(list=ls())
set.seed(42)




# -------------
# functions
# -------------
sim.mv.AR2 <- function(N, d, rho) {
  # Simulates N observations of a multivariate stationary AR(2)
  # process given by:
  #   X_{t,j} = 0.3 X_{t-1,j} - 0.4 X_{t-2,j} + e_{t,j}
  # with
  #   j=1,...,d
  #   e_{t,j} ~ N(0, Sigma)
  #   X_{0,j} = X_{1,j} = 0
  #   Sigma_{i,i} = 1
  #   Sigma_{i,j} = rho for j != i
  #   -1/(d-1) < rho < 1
  
  # ensure positive definite
  stopifnot(-1/(d-1) < rho)
  stopifnot(rho < 1)
  
  # create equi-correlated matrix
  Sigma <- matrix(rho, d, d)
  diag(Sigma) <- 1
  
  # create matrix
  X <- matrix(0, ncol=d, nrow=N)
  for (t in 3:N) {
    e_t <- mvrnorm(1, mu=rep(0,d), Sigma=Sigma)
    X[t,] <- 0.3*X[t-1,] - 0.4*X[t-2,] + e_t
  }
  
  return(X)
}

dgp <- function(X) {
  return(-sin(2*X[1]) + (X[2]^2 - (25/12)) + X[3] + (exp(-X[4]) - (2/5)*sinh(5/2)))
}

make.D <- function(X, d) {
  # Creates the data matrix D as specified
  # in exercise 2. First column represents
  # the response variable Y_t.
  u <- rnorm(nrow(X), mean=0, sd=1)
  Y <- matrix(0, nrow=nrow(X), ncol=1)
  for (t in 1:nrow(X)) {
    Y[t] <- dgp(X[t,]) + u[t]
  }
  D <- cbind(Y, X)
  return(D)
}


run.sim <- function(S, N, d, rho, mtry) {
  # prediction errors vectors
  err.rnd <- rep(0, S)
  err.dgp <- rep(0, S)
  err.ar <- rep(0, S)
  err.rtree <- rep(0, S)
  err.rtree.pruned <- rep(0, S)
  err.rforest <- rep(0, S)
  err.rforest.oob <- rep(0, S)
  err.ols <- rep(0, S)
  err.lasso <- rep(0, S)
  
  # e) keep track if first four variables have been selected
  #    by the variable importance measure of the random forest
  rforest.var.imp.sel <- matrix(0, nrow=S, ncol=4)
  
  for (s in 1:S) {
    cat(sprintf('Simulation %i of %i \n', s, S))
    
    # simulate multivariate AR(2) process
    X <- sim.mv.AR2(N+1, d, rho)
    
    # data matrix D
    D <- as.data.frame(make.D(X, d))
    colnames(D) <- c('Y', paste0('X', 1:d))
    
    # split into training (n=500) and test data (n=1)
    D.train <- D[-nrow(D), ]
    D.test <- D[nrow(D), ]
    
    # squared test error function
    test.serr <- function(pred) (D.test[1,'Y'] - pred)^2
    
    
    
    # ---------------- AR(2) benchmark ----------------
    ar.fit <- ar(D.train[,1], aic=F, method='mle', order.max=2)
    ar.pred <- predict(ar.fit, newdata=tail(D.train[,1], n=2), n.ahead=1)$pred[1]
    err.ar[s] <- test.serr(ar.pred)
    
    
    # ---------- Gaussian random prediction -----------
    err.rnd[s] <- test.serr(rnorm(1, mean=mean(X), sd=sd(X)))
    
    
    # --------- True data generating process ----------
    dgp.pred <- dgp(unlist(D.test[,2:ncol(D.test)]))
    err.dgp[s] <- test.serr(dgp.pred)
    
    
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
    rforest.fit <- ranger(Y ~ ., data=D.train, mtry=mtry, num.trees=100, importance='impurity')
    rforest.pred <- predict(rforest.fit, data=D.test)$predictions[1]
    err.rforest[s] <- test.serr(rforest.pred)
    
    # c) Out-of-bag error
    err.rforest.oob[s] <- rforest.fit$prediction.error
    
    # e) selected variables by variable importance measure
    rforest.var.imp.sel[s,] <- names(sort(rforest.fit$variable.importance, decreasing=T))[1:4] %in% names(D)[2:5]
    
    
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
  out <- sprintf('
     rho=%.1f
     MSE Gaussian random:         %.3f
     MSE AR(2) fitted:            %.3f
     MSE DGP:                     %.3f
     MSE regression tree:         %.3f
     MSE regression tree pruned:  %.3f
     MSE random forest:           %.3f (OOB=%.3f, %.0f%% of times greater than MSE)
     MSE OLS:                     %.3f
     MSE LASSO:                   %.3f
  ', rho, mean(err.rnd), mean(err.ar), mean(err.dgp),
     mean(err.rtree), mean(err.rtree.pruned),
       mean(err.rforest), mean(err.rforest.oob),
     100*mean(err.rforest.oob > err.rforest),
     mean(err.ols), mean(err.lasso))
  cat(out)
  
  # save to txt file
  hdl <- file(sprintf('output/2_MSE_rho=%.1f_mtry=%i_d=%i.txt', rho, mtry, d))
  writeLines(c(out), hdl)
  close(hdl)
  
  
  
  
  
  
  # barplot of MSE by method
  pdf(file=sprintf('output/2_MSE_rho=%.1f_mtry=%i_d=%i.pdf', rho, mtry, d), width=8, height=8)
  
  err.df <- data.frame(Method=c(
    'Gaussian random', 'DGP', 'AR(2) fitted', 'Regression tree',
    'Regression tree pruned', 'Random forest', 'OLS', 'LASSO'),
    MSE=c(mean(err.rnd), mean(err.dgp), mean(err.ar), mean(err.rtree),
          mean(err.rtree.pruned), mean(err.rforest),
          mean(err.ols), mean(err.lasso)))
  err.df <- err.df[order(err.df$MSE), ] # sort by MSE
  par(mfrow=c(1,1), mar=c(5, 10, 5, 2))
  b <- barplot(err.df$MSE, names.arg=err.df$Method,
               main=sprintf('MSE test prediction error (rho=%.1f)', rho),
               horiz=T, xlab='MSE',
               xlim=c(0, max(err.df$MSE)*1.4), las=2)
  text(x=err.df$MSE + 0.1*max(err.df$MSE), y=b, label=format(err.df$MSE, digits=2), cex=0.75)
  
  dev.off()
  
  
  # e) show how many times variable of lag 1-4 have been selected
  #    by the variable importance measure of the random forest
  out <- sprintf('
     rho=%.1f
     Selected lag 1,2,3,4:  %.0f%%
     Selected lag=1:        %.0f%%
     Selected lag=2:        %.0f%%
     Selected lag=3:        %.0f%%
     Selected lag=4:        %.0f%%
     Overall:               %.0f%%
  ', rho,
     100*mean(rowSums(rforest.var.imp.sel) == 4),
     100*mean(rforest.var.imp.sel[,1]),
     100*mean(rforest.var.imp.sel[,2]),
     100*mean(rforest.var.imp.sel[,3]),
     100*mean(rforest.var.imp.sel[,4]),
     100*mean(rforest.var.imp.sel))
  cat(out)
  
  # save to txt file
  hdl <- file(sprintf('output/2_var_imp_rho=%.1f_mtry=%i_d=%i.txt', rho, mtry, d))
  writeLines(c(out), hdl)
  close(hdl)
}

# ------------------------------------------------------------------------------------------



# verify multivariate AR(2) simulation function by plotting PACFs
pdf(file='output/2_PACF.pdf', width=8, height=8)
par(mfrow=c(2,2), mar=c(5, 5, 4, 1))
ts <- sim.mv.AR2(5000, 4, 0.45)
for (i in 1:4) {
  pacf(ts[,i], lag.max=10, main='PACF simulated AR(2) with N=5000')
}
dev.off()

# print correlation matrix
cor(ts[,1:4])



# ------------- 
# simulations
# -------------

# parameters
S <- 100  # number of simulations
N <- 200  # T training samples from AR(2) process
d <- 150  # number of predictor columns in D


# parameters matrix
params <- matrix(0, nrow=10, ncol=2)

# increasing rho
params[1,] = c(0.0, floor(sqrt(d)))
params[2,] = c(0.1, floor(sqrt(d)))
params[3,] = c(0.5, floor(sqrt(d)))
params[4,] = c(0.9, floor(sqrt(d)))

# increasing mtry
params[5,] = c(0.1, floor(sqrt(d)))
params[6,] = c(0.1, 50)
params[7,] = c(0.1, 120)

params[8,] = c(0.9, floor(sqrt(d)))
params[9,] = c(0.9, 50)
params[10,] = c(0.9, 120)



# create parallel computation processes
cl <- makeCluster(detectCores())
registerDoParallel(cl)
              
# run simulations using different rho and mtry
# values, see params matrix
pkgs <- c('rpart', 'rpart.plot', 'ranger', 'MASS', 'glmnet')
foreach(i=1:nrow(params), .packages=pkgs) %dopar% {
  run.sim(S, N, d, params[i,1], params[i,2])
}

# dispose parallel computation processes
stopCluster(cl)

    









