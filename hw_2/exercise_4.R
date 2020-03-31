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
# Topic:        Homework 2 - Exercise 4.
###########################################################

library(glmnet)

# clear environment
rm(list=ls())

# set RNG seed for reproducibility.
# the choice of seed is influential, because some
# samples with n=100 randomly show significant correlations,
# which influences the cross-validation result for Ridge and LASSO.
set.seed(5)



# create dataset of standard normal and independent variables
n = 100
data = data.frame(matrix(rnorm(n * 11), n, 11))
colnames(data)[11] = 'Y'



# ------------------------------------------------------------------------------
# a)

# split data into training and test datasets 50:50
data_train = data[1:(n/2), ]
data_test  = data[(n/2+1):n, ]


# fit models onto training data (first one is the empty model)
models = list()
models[[1]] = lm(Y ~ 1, data=data_train)

# add and fit one predictor after the other, up to full model
for (x in 1:10)
  models[[x+1]] = lm(paste("Y ~ 1 + ", paste(paste0('X', 1:x), collapse="+")), data=data_train)


# extract MSE, R-squared, Adj. R-squared for all fitted models
train_mse = sapply(models, function(m) mean(m$residuals^2))
train_r2 = sapply(models, function(m) summary(m)$r.squared)
train_adj_r2 = sapply(models, function(m) summary(m)$adj.r.squared)


# plot MSE, R-squared, Adj. R-squared
par(mfrow=c(1,2))
plot(train_mse, col='blue', type='l', xlab='Model complexity', ylab='MSE',
     xaxt='n', ylim=c(-0.1*(max(train_mse) - min(train_mse)) + min(train_mse), max(train_mse)))
axis(side=1, at=1:11, labels=0:10)
title('OLS MSE on training data')
legend('bottomleft', inset=c(0.01, 0.01), legend=c('MSE'), cex=1,
       col=c('blue'), lty=1, box.lty=0, horiz=T, bg='transparent')

plot(train_r2, col='blue', type='l', xlab='Model complexity',
     ylab=as.expression(R^2 ~ ''), lty=1, xaxt='n',
     ylim=c(-0.2*(max(train_r2, train_adj_r2) - min(train_r2, train_adj_r2)) + min(train_r2, train_adj_r2), max(train_r2, train_adj_r2)))
axis(side=1, at=1:11, labels=0:10)
title(expression(bold('OLS' ~ R^2 ~ 'and Adjusted' ~ R^2 ~ 'on training data')))
lines(train_adj_r2, col='blue', lty=5)
legend('bottomleft', inset=c(0.01, 0.01), cex=1, col=c('blue', 'blue'),
       lty=c(1, 5), box.lty=0, bg='transparent',
       legend=c(as.expression(R^2 ~ ''), as.expression('Adjusted' ~ R^2)), horiz=T)



# ------------------------------------------------------------------------------
# b)

# compute MSE on test dataset
test_mse = sapply(models, function(m) mean((data_test$Y - predict.lm(m, data_test))^2))


# plot training MSE and test MSE
par(mfrow=c(1,1))
plot(train_mse, col='blue', type='l', xlab='Model complexity', ylab='MSE', xaxt='n',
     ylim=c(-0.1*(max(train_mse, test_mse) - min(train_mse, test_mse)) + min(train_mse, test_mse), max(train_mse, test_mse)))
title('OLS MSE training vs. test data')
lines(test_mse, col='red', lty=5)
axis(side=1, at=1:11, labels=0:10)
legend('bottomleft', inset=c(0.01, 0.01), legend=c('MSE training', 'MSE test'),
       cex=1, col=c('blue', 'red'), lty=c(1, 5), box.lty=0, horiz=T, bg='transparent')





# ------------------------------------------------------------------------------
# c)

# create vector of lambdas
lambdas = seq(0, 4, length=500)
x = as.matrix(data[, -ncol(data)])
y = as.matrix(data[, ncol(data)])

# perform k-fold cross-validation for Ridge and LASSO regression
k = 10
foldid = sample(1:k, size=length(y), replace=T) # use same k-fold CV for both approaches
cv_ridge = cv.glmnet(x, y, foldid=foldid, lambda=lambdas, alpha=0, standardize=F)
cv_lasso = cv.glmnet(x, y, foldid=foldid, lambda=lambdas, alpha=1, standardize=F)

# find minimum MSE and corresponding lambda
cv_ridge_min_mse = cv_ridge$cvm[which(cv_ridge$cvm == min(cv_ridge$cvm))]
cv_ridge_min_mse_lambda = cv_ridge$lambda.min

  sprintf('Ridge min. MSE: %4f', cv_ridge_min_mse)
  sprintf('Ridge min. MSE lambda: %.4f', cv_ridge_min_mse_lambda)

cv_lasso_min_idx = which(cv_lasso$cvm == min(cv_lasso$cvm))
cv_lasso_min_mse = min(cv_lasso$cvm[cv_lasso_min_idx])
cv_lasso_min_lambda = cv_lasso$lambda[cv_lasso_min_idx]

print('LASSO min. MSE')
cv_lasso_min_mse

print('LASSO min. MSE lambdas:')
cv_lasso_min_lambda

print('Ridge coefficients at min. MSE:')
coef(cv_ridge, s=cv_ridge$lambda.min)

print('LASSO coefficients at min. MSE:')
coef(cv_lasso, s=cv_lasso$lambda.min)

mean(y)



# plot MSE against lambda
par(mfrow=c(1,1))
plot(log(cv_ridge$lambda), cv_ridge$cvm, type='l', col='red',
     xlab=expression(log(lambda)), ylab='MSE', ylim=c(min(cv_ridge$cvm, cv_lasso$cvm), max(cv_ridge$cvm, cv_lasso$cvm)))
title('MSE of Ridge and LASSO regression')
points(log(cv_lasso$lambda), cv_lasso$cvm, type='l', col='blue')
legend('topright', inset=c(0.03, 0.03),
       legend=c(expression('Ridge' ~ (alpha==0.0)), expression('LASSO' ~ (alpha==1.0))),
       col=c('red', 'blue'), lty=1, box.lty=0, bg='transparent')



# plot Ridge coefficients against log of lambda
par(mfrow=c(1,2))
plot(cv_ridge$glmnet.fit, 'lambda')
title(expression(bold('Ridge coefficients against penalty' ~ log(lambda))), line=2.8)
vnat = coef(cv_ridge$glmnet.fit)
vnat = vnat[-1, ncol(vnat)]
axis(4, at=vnat, line=-.5, label=labels(vnat), las=1, tick=F, cex.axis=0.7)

# plot LASSO coefficients against log of lambda
plot(cv_lasso$glmnet.fit, 'lambda')
title(expression(bold('LASSO coefficients against penalty' ~ log(lambda))), line=2.8)
vnat = coef(cv_lasso$glmnet.fit)
vnat = vnat[-1, ncol(vnat)]
axis(4, at=vnat, line=-.5, label=labels(vnat), las=1, tick=F, cex.axis=0.7)







