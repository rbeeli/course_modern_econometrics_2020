############################################################
# Course:       Modern econometric and statistical learning
#               methods forquantitative asset management
#
# Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
#               University of Zurich
#
# Author:       Rino Beeli
#
# Date:         March 31st, 2020
# 
# Topic:        Homework 2 - Exercise 5.
###########################################################

library(fields)
library(glmnet)
library(numbers)
library(matrixStats)

# clear environment
rm(list=ls())

# set RNG seed for reproducibility
set.seed(5)

# load dataset
load('JanFeb.rda')
X = X
Y = as.matrix(Y)

# print dimensions
dim(X)
dim(Y)

# show data summary
summary(as.vector(X))
summary(as.vector(Y))

# plot first observation
plot(X[1,], type='l', xlab='x', ylab='y')
title('First row of matrix X')

# spatial data - what's the stride length?
plot_act = acf(X[1,], lag.max=250, plot=F)
plot(plot_act, main='Autocorrelation of 1st row of matrix X')

# we see an autocorrelation peak at around ~145.
# let's see which values around 145 divide the number of columns of X without rest
divisors(ncol(X))

# 144 divides ncol(X) = 10368 without rest (=72)

# plot 2D image function
plot_2d = function(mat, image_plot_format) {
  if (is.vector(mat))
    mat = as.matrix(mat)
  if (!image_plot_format) {
    mat = as.matrix(mat)
    mat = t(apply(mat, 2, rev))
  }
  image.plot(mat)
}

# plot first row of X as 2D image
plot_2d(matrix(X[1,], nrow=144), image_plot_format=T)

# plot correlation of each area to target variable
plot_2d(matrix(cor(X, Y), nrow=144), image_plot_format=T)
plot_2d(matrix(cor(X, rowMeans(X)), nrow=144), image_plot_format=T)

# correlation between Y and rowMeans(X)
cor(Y, rowMeans(X))



# perform k-fold cross-validation for Ridge (alpha=0), LASSO (alpha=1)
# and Elastic Net regressions (alpha in 12.5% steps)
alphas = seq(0, 1, by=0.125)
cv_models = list()
cv_coefs = list()
cv_min_mse = matrix(NA, length(alphas))
cv_lambda = matrix(NA, length(alphas))

for (alpha in alphas) {
  i = length(cv_models)+1
  
  # fit model using cross validation
  cv_fit = cv.glmnet(X, Y, alpha=alpha)
  cv_models[[i]] = cv_fit

  # coefficients at lambda with minimum MSE
  cv_coefs[[i]] = coef(cv_fit, s=cv_fit$lambda.min)

  # minimum MSE
  cv_min_mse[i] = cv_fit$cvm[which(cv_fit$lambda == cv_fit$lambda.min)]
  
  # optimal lambda
  cv_lambda[i] = cv_fit$lambda.min
}




# plot MSE
par(mfrow=c(3,3), mar=c(4,4.5,5,1))
for (i in 1:length(alphas)) {
  cv_model = cv_models[[i]]
  plot(cv_model)
  name = 'Elastic Net'
  if (alphas[i] == 0) {
    name = 'Ridge'
  }
  if (alphas[i] == 1) {
    name = 'LASSO'
  }
  title(bquote(bold(.(name) ~ 'regression' ~ alpha==.(alphas[i]))), line=3.3)
}


# plot MSE depending on alpha
par(mfrow=c(1,1))
plot(alphas, cv_min_mse, type='l',
     xlab=expression(bold(alpha)), ylab='MSE',
     main=expression(bold('Elastic Net MSE at different' ~ alpha)))


# plot coefficients
par(mfrow=c(3,3), mar=c(2,3,3,6))
for (i in 1:length(alphas)) {
  coefs = cv_coefs[[i]]
  image.plot(matrix(coefs[-1], nrow=144), main=bquote(bold('Elastic Net regression' ~ alpha == .(alphas[i]))))
}


# optimal lambdas
cbind(alphas, cv_lambda)

# minimum MSEs
mean((Y - rowMeans(X))^2)
cbind(alphas, cv_min_mse, cv_min_mse == min(cv_min_mse))







