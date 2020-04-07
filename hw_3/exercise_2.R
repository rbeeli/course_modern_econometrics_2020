############################################################
# Course:       Modern econometric and statistical learning
#               methods forquantitative asset management
#
# Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
#               University of Zurich
#
# Author:       Rino Beeli
#
# Date:         April 7th, 2020
# 
# Topic:        Homework 3 - Exercise 2.
###########################################################

library(car)

rm(list=ls())     # clear environment
set.seed(42)      # set RNG seed for reproducibility
Sys.setlocale(locale="en_US.UTF-8") # english locale



# S&P 500 (log-returns) and annualized risk-free returns
market = read.csv('data/m_sp500ret_3mtcm.txt', sep='\t', skip=1)
market = market[c('Date', 'sp500', 'X3mTCM')]
colnames(market) = c('Date', 'Mkt', 'Rf')

# parse date and extract month and year
market$Date = as.Date(gsub("^", "01-", market$Date), format="%d-%b-%y")

# convert annual to monthly returns (log-returns)
# convert percentage to decimal returns
market$Rf = market$Rf / 12 / 100
market$Mkt.ex = market$Mkt - market$Rf
  
# ten stock returns matrix (log-returns)
stocks = read.csv('data/m_logret_10stocks.txt', sep='\t')
stocks = na.omit(stocks) # remove empty rows
stocks$Date = as.Date(strftime(as.Date(stocks$Date, '%m/%d/%Y'), '%Y-%m-01'))

# sanity check for matching dates
stopifnot(all.equal(market$Date, stocks$Date))


N = nrow(market)
K = 10




# --------------------------------------------------
# b) CAPM on 10 stocks using linear regression
# --------------------------------------------------

alphas = matrix(0, K, 1)
betas = matrix(0, K, 1)
SRs = matrix(0, K, 1) # annualized Sharpe ratios
TRs = matrix(0, K, 1) # annualized Treynor ratios

for (i in 1:K) {
  rets.ex = stocks[,i+1] - market$Rf
  mkt.ex = market$Mkt.ex
  
  mu_i = mean(rets.ex)
  std_m = sd(mkt.ex)
  std_i = sd(rets.ex)
  cov_im = cov(rets.ex, mkt.ex)
  corr_im = cor(rets.ex, mkt.ex)
  
  fit = lm(rets.ex ~ 1 + mkt.ex)
  fit_sum = summary(fit)
  
  # statistics of interest
  alpha = alphas[i] = fit_sum$coef[1,1]
  beta = betas[i] = fit_sum$coef[2,1]
  SR = mu_i / std_i
  TR = (mu_i*std_m^2)/cov_im
  SR.ann = SRs[i] = SR * sqrt(12) # annualized
  TR.ann = TRs[i] = TR * 12 # annualized
  
  # 5% confidence intervals for coefficients
  ci = confint(fit, level=0.95)
  ci_alpha = ci[1,]
  ci_beta = ci[2,]
  
    # CIs taken from:
    #   Performance Hypothesis Testing with the Sharpe and Treynor Measures
    #      by J. D. Jobson and Bob M. Korkie (1981)
    #   https://www.jstor.org/stable/2327554?seq=1
  ci_SR = SR + c(-1,1)*qnorm(0.975)*sqrt(1/N*(1 + 1/2*SR^2))
  ci_TR = TR + c(-1,1)*qnorm(0.975)*sqrt(1/N*(std_m^4/cov_im^2)*(std_i^2 + mu_i^2*(1-1/corr_im^2)))
  
  # annualize
  ci_SR = ci_SR * sqrt(12)
  ci_TR = ci_TR * 12
  
  if (i == 1) {  
    cat("Stock  Alpha            Beta             SR                TR\n")
  }
  
  cat(sprintf("%02.0f     %.3f            %.3f            %.3f             %.3f
         [%.3f, %.3f]   [%.3f, %.3f]   [%.3f, %.3f]   [%.3f, %.3f]\n",
          i, alpha, beta, SR.ann, TR.ann,
          ci_alpha[1], ci_alpha[2],
          ci_beta[1], ci_beta[2],
          ci_SR[1], ci_SR[2],
          ci_TR[1], ci_TR[2]))
}



# --------------------------------------------------
# c) Confidence intervals of alpha, beta, SR and TR
# based on bootstrap method
# --------------------------------------------------

T = 10000

b_alphas = matrix(0, T, K)
b_betas = matrix(0, T, K)
b_SRs = matrix(0, T, K)
b_TRs = matrix(0, T, K)

for (i in 1:K) {
  rets.ex = as.matrix(stocks[,i+1] - market$Rf)
  mkt.ex = as.matrix(market$Mkt.ex)
  
  cat(paste('Bootstrapping stock', i, '/', K,'... \n'))
  
  for (t in 1:T) {  
    # sample with replacement
    idxs = sample(length(rets.ex), length(rets.ex), replace=T)
    b_rets.ex = rets.ex[idxs]
    b_mkt.ex = mkt.ex[idxs]
    
    fit = lm(b_rets.ex ~ 1 + b_mkt.ex)
    fit_sum = summary(fit)
    
    # statistics of interest
    alpha = b_alphas[t, i] = fit_sum$coef[1,1]
    beta = b_betas[t, i] = fit_sum$coef[2,1]
    SR = b_SRs[t, i] = mean(b_rets.ex) / sd(b_rets.ex) * sqrt(12) # annualize
    TR = b_TRs[t, i] = mean(b_rets.ex) / beta * 12 # annualize
  }
}



cat(sprintf('Bootstrap results with T=%i replications\n', T))

for (i in 1:K) {
  # bootstrapped 5% confidence intervals for coefficients
  ci_alpha = quantile(b_alphas[,i], probs=c(0.025, 0.975))
  ci_beta = quantile(b_betas[,i], probs=c(0.025, 0.975))
  ci_SR = quantile(b_SRs[,i], probs=c(0.025, 0.975))
  ci_TR = quantile(b_TRs[,i], probs=c(0.025, 0.975))
  
  if (i == 1) {  
    cat("Stock  Alpha             Beta             SR                TR\n")
  }
  
  cat(sprintf("%02.0f     [%.3f, %.3f]   [%.3f, %.3f]   [%.3f, %.3f]   [%.3f, %.3f]\n",
              i,
              ci_alpha[1], ci_alpha[2],
              ci_beta[1], ci_beta[2],
              ci_SR[1], ci_SR[2],
              ci_TR[1], ci_TR[2]))
  
    # histograms
    par(mfrow=c(1,4))
    
    hist(b_alphas[,i], main=sprintf('Alpha of stock %i', i), xlab=expression(alpha))
    abline(v=alphas[i], col='red')
    
    hist(b_betas[,i], main=sprintf('Beta of stock %i', i), xlab=expression(beta))
    abline(v=betas[i], col='red')
    
    hist(b_SRs[,i], main=sprintf('Sharpe Ratio of stock %i', i), xlab='Sharpe Ratio')
    abline(v=SRs[i], col='red')
    
    hist(b_TRs[,i], main=sprintf('Treynor Ratio of stock %i', i), xlab='Treynor Ratio')
    abline(v=TRs[i], col='red')
}




# --------------------------------------------------
# f) Segmented/piecewise regression in order to test
#    time-dependence of CAPM beta coefficient
# --------------------------------------------------

t_0 = which(market$Date == '2001-02-01')

# indicator variables
ind_1 = as.integer(1:N < t_0)
ind_2 = as.integer(1:N >= t_0)

for (i in 1:K) {
  rets.ex = stocks[,i+1] - market$Rf
  mkt.ex = market$Mkt.ex
  
  segment1 = ind_1*mkt.ex
  segment2 = ind_2*mkt.ex
  
  # no intercept!
  fit = lm(rets.ex ~ -1 + segment1 + segment2)
  
  fit_sum = summary(fit)
  print(fit_sum)
  
  # test hypothesis beta_1 = beta_2
  test_eq = linearHypothesis(fit, 'segment1 = segment2')
  cat(sprintf('%02.0f   H_0: beta_1 = beta_2, p-value=%.4f\n', i, test_eq[2, 'Pr(>F)']))
}




# --------------------------------------------------
# f) Find best segmentation point t_0 for
#    time dependent CAPM beta
# --------------------------------------------------

MSEs = matrix(0, N, 1)

for (t_0 in 1:N) {
  # indicator variables
  ind_1 = as.integer(1:N < t_0)
  ind_2 = as.integer(1:N >= t_0)
  
  for (i in 1:K) {
    rets.ex = stocks[,i+1] - market$Rf
    mkt.ex = market$Mkt.ex
    
    segment1 = ind_1*mkt.ex
    segment2 = ind_2*mkt.ex
    
    # no intercept!
    fit = lm(rets.ex ~ -1 + segment1 + segment2)
    fit_sum = summary(fit)
    MSEs[t_0] = MSEs[t_0] + mean(fit_sum$residuals^2)
  }
}

# find minimum MSE
MSEs.min = which.min(MSEs)
MSEs.min.date = strftime(stocks$Date[MSEs.min], '%d.%m.%Y')

par(mfrow=c(1,1))
plot(MSEs, type='l', main='Total MSE in dependenc of t_0',
     xlab=expression(t_0), ylab='Average MSE')
abline(v=MSEs.min, col='red')

cat(sprintf('Minimum MSE at t_0=%i, %s\n', MSEs.min, stocks$Date[MSEs.min]))



# plot stock and market returns and t_0 which minimizes average MSE
cols = c('black', 'black', 'blue', 'green', 'cyan', 'red',
         'chartreuse', 'cadetblue3', 'burlywood3', 'darkorchid4', 'deeppink')
cumrets = 1 + cumsum(cbind(market$Mkt, stocks[,-1]))
colnames(cumrets)[1] = 'Market'
matplot(1:nrow(cumrets), cumrets, type='l',
        lty=c(1, rep(6, 10)), lwd=c(2, rep(1, 10)), xlim=c(0, 170),
        col=cols, main='Stock returns', xlab='Date index', ylab='Cumulative returns')
abline(v=MSEs.min, col='red', lwd=3)
text(nrow(cumrets) + 15, cumrets[nrow(cumrets),][seq(2, 10, 2)],
     labels=colnames(cumrets)[seq(2, 10, 2)], cex=0.7,
     col=cols[seq(2, 10, 2)])
text(nrow(cumrets) + 5, cumrets[nrow(cumrets),][seq(1, 10, 2)],
     labels=colnames(cumrets)[seq(1, 10, 2)], cex=0.7,
     col=cols[seq(1, 10, 2)])
text(MSEs.min+18, 0.5, labels=bquote(t[0]^"*"==.(MSEs.min.date)), cex=1.1)



  




  