  ############################################################
  # Course:       Modern econometric and statistical learning
  #               methods forquantitative asset management
  #
  # Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
  #               University of Zurich
  #
  # Author:       Rino Beeli
  #
  # Date:         April 16th, 2020
  # 
  # Topic:        Homework 4 - Exercise 1.
  ###########################################################
  
  library(purrr)
  library(glmnet)
  
  rm(list=ls())     # clear environment
  set.seed(42)      # set RNG seed for reproducibility
  Sys.setlocale(locale="en_US.UTF-8") # english locale
  
  
  # read NASDAQ data file
  data =  read.csv('data/d_nasdaq_82stocks.txt', sep='\t')
  
  # remove empty column
  data = data[,-which(colnames(data) == 'X')]
  
  # sanity check for missing values
  stopifnot(which(is.na(data)) == 0)
  
  dates = data$Date
  nasdaq = data$Nasdaq
  stocks = data[,-c(1,2)]
  
  
  
  # function for plotting NASDAQ vs. tracking portoflio
  plot.nasdaq.vs.pf = function(filename, pf.rets, stocks, plot.title) {
    # save PDF of plot to filesystem
    pdf(paste0('../report/', filename), width=6, height=6*9/16, pointsize=8)
    
    plot.data = cbind(1 + cumsum(nasdaq), 1 + cumsum(pf.rets))
    tracking_error = sqrt(252) * sd(pf.rets - nasdaq)
    
    par(mfrow=c(1,1))
    matplot(plot.data, type='l', col=c('red', 'blue'), lty=1,
            main=plot.title,
            xlab='Date', xaxt='n', ylab='Cumulative returns')
    axis_idx = seq(1, nrow(data), by=floor(nrow(data) / 12))
    axis(1, at=axis_idx, labels=dates[axis_idx])
    legend('topleft', legend=c('NASDAQ', sprintf('Tracking portfolio (N=%i)', length(stocks))),
           col=c('red', 'blue'), lty=1, cex=0.8)
    text(length(dates)*0.87, 0.9, labels=sprintf('Tracking Error: %.3f%%', tracking_error*100))
    
    # write PDF to filesystem
    dev.off()
  }
  
  
  # -------------------------------------------------------------
  # a) Tracking benchmark index using only significant
  #    constutients returns contributors (backward selection)
  # -------------------------------------------------------------
  
  if (T) {
    #   ------------------ Backward selection  ------------------
    
    # full model
    fit.full = lm(nasdaq ~ -1 + ., data=stocks)
    
    #alpha = 0.1
    #alpha = 0.05
    alpha = 0.01
    
    # start with full model
    bw.stocks = stocks
    bw.fit = fit.full
    bw.fit.sum = summary(bw.fit)
    
    print(bw.fit.sum)
    
    for (pass in 1:ncol(stocks)) {
      # check if all regressors have significant p-values
      insig.num = sum(coef(bw.fit.sum)[,'Pr(>|t|)'] > alpha)
      if (insig.num == 0) {
        # all regressors are significant, stop loop
        break
      }
      
      cat(sprintf('Testing model with %i of %i variables (%i insignificant regressors present) \n', 
                  ncol(bw.stocks) - 1, ncol(stocks), insig.num))
      
      # fit all one-term deletion models
      fw.partial = list()
      fw.fstats = rep(0, ncol(bw.stocks))
      for (i in 1:ncol(bw.stocks)) {
        # fit model with removed stock
        fit2 = lm(nasdaq ~ -1 + ., data=bw.stocks[,-i])
        fw.partial[[i]] = fit2
        
        # calc partial F-statistic
        RSS0 = sum(fit2$residuals^2)
        RSS1 = sum(bw.fit$residuals^2)
        n = nrow(stocks)
        k = length(coef(bw.fit))
        fw.fstats[i] = (RSS0 - RSS1)/(RSS1/(n - k))
      }
      
      # remove regressor/stock with least significant OLS coefficient
      insig.idx = which.min(fw.fstats)
      insig.stock = colnames(bw.stocks)[insig.idx]
      cat(sprintf('Removed stock %s (coef. p-value %.4f, partial F-stat %.4f)\n',
                  insig.stock, coef(bw.fit.sum)[insig.idx, 'Pr(>|t|)'], fw.fstats[insig.idx]))
      bw.stocks = bw.stocks[,-insig.idx]
      
      # choose the new model as current best model
      bw.fit = lm(nasdaq ~ -1 + ., data=bw.stocks)
      bw.fit.sum = summary(bw.fit)
    }
    
    bw.fit
    bw.fit.stocks = names(coef(bw.fit))
    bw.fit.idxs = which(colnames(stocks) %in% bw.fit.stocks)
    
    cat('------------------------------------------------\n')
    cat(sprintf('Backward-selection resulted in %i out of %i stocks for tracking NASDAQ at alpha=%.2f\n',
                length(bw.fit.stocks), ncol(stocks), alpha))
    cat(sprintf('Selected stocks:\n%s\n', paste(bw.fit.stocks, collapse=' ')))
    cat('------------------------------------------------\n')
    
    
    # plot NASDAQ vs. tracking portfolio
    nasdaq.tracked = as.matrix(stocks[, bw.fit.idxs]) %*% coef(bw.fit)
    plot.nasdaq.vs.pf(sprintf('1_a_BW_alpha=%.2f.pdf', alpha), nasdaq.tracked, bw.fit.stocks,
                      sprintf('Tracking portfolio using backward-selection (alpha=%.2f)', alpha))
  
  }
  
  # -------------------------------------------------------------
  # b) Tracking benchmark index using predefined
  #    number of stocks N (forward selection)
  # -------------------------------------------------------------
  
  Ns = c(3, 5, 10)
  
  if (T) {
    #   ------------------ Forward selection  ------------------
    
    for (N in Ns) {
      
      fw.stocks.in = data.frame(matrix(nrow=nrow(data), ncol=0)) # empty
      fw.stocks.out = stocks
      
      # start with empty model
      fw.fit = lm(nasdaq ~ -1)
      
      for (pass in 1:N) {
        cat(sprintf('Fitting %i out of %i variables...\n', pass, N))
        
        # fit all one-term addition models
        fit.partial = list()
        fit.fstats = rep(0, ncol(fw.stocks.out))
        for (i in 1:ncol(fw.stocks.out)) {
          # fit model with additional stock
          fit.data = cbind(fw.stocks.in, fw.stocks.out[,i])
          colnames(fit.data)[ncol(fit.data)] = colnames(fw.stocks.out)[i]
          fit2 = lm(nasdaq ~ ., data=fit.data)
          fit.partial[[i]] = fit2
          
          # calc partial F-statistic
          RSS0 = sum(fw.fit$residuals^2)
          RSS1 = sum(fit2$residuals^2)
          n = nrow(stocks)
          k = length(coef(fit2))
          fit.fstats[i] = (RSS0 - RSS1)/(RSS1/(n - k))
        }
        
        # add regressor/stock with most significant OLS coefficient
        sig.idx = which.max(fit.fstats)
        sig.stock = colnames(fw.stocks.out)[sig.idx]
        cat(sprintf('Added stock %s\n', sig.stock))
        fw.stocks.in = cbind(fw.stocks.in, fw.stocks.out[,sig.idx])
        colnames(fw.stocks.in)[ncol(fw.stocks.in)] = colnames(fw.stocks.out)[sig.idx]
        fw.stocks.out = fw.stocks.out[,-sig.idx]
        
        # choose the new model as current best model
        fw.fit = lm(nasdaq ~ -1 + ., data=fw.stocks.in)
      }
      
      fw.fit
      fw.fit.stocks = names(coef(fw.fit))
      
      cat(sprintf('Forward-selection selected %i out of %i stocks for tracking NASDAQ\n',
                  length(fw.fit.stocks), ncol(stocks)))
      cat(sprintf('Selected stocks:  %s\n', paste(sort(fw.fit.stocks), collapse=' ')))
      
      
      # plot NASDAQ vs. tracking portfolio
      nasdaq.tracked = as.matrix(stocks[, fw.fit.stocks]) %*% coef(fw.fit)
      plot.nasdaq.vs.pf(sprintf('1_b_FW_N=%i.pdf', N), nasdaq.tracked, fw.fit.stocks,
                        sprintf('Tracking portfolio using forward-selection with N=%i', N))
      
    }
  }
  
  
  if (T) {
    #   ------------------ LASSO  ------------------
    
    for (N in Ns) {
      
      lasso.all = glmnet(as.matrix(stocks), y=nasdaq, alpha=1, intercept=F)
      lasso.all.coefs = coef(lasso.all)
      
      # extract model with N regressors
      coefs.idx = which(colSums(abs(lasso.all.coefs) > 0.000001) == N)[1]
      lambda = lasso.all$lambda[coefs.idx]
      
      # fit best LASSO model again
      lasso.fit = glmnet(as.matrix(stocks), y=nasdaq, alpha=1, lambda=lambda, intercept=F)
      coefs.lasso = coef(lasso.fit)[,1]
      
      coefs.lasso = coefs.lasso[which(abs(coefs.lasso) > 0.000001)]
      lasso.stocks = names(coefs.lasso)
      
      # fit OLS using LASSO selected stocks (no regularization)
      lasso.ols = lm(nasdaq ~ -1 + ., data=stocks[, lasso.stocks])
      coefs.ols = coef(lasso.ols)
      
      cat(sprintf('LASSO stocks at N=%i:  %s\n', N, paste(sort(lasso.stocks), collapse=' ')))
      
      
      # plot NASDAQ vs. tracking portfolio using LASSO
      nasdaq.tracked = as.matrix(stocks[, lasso.stocks]) %*% coefs.lasso
      plot.nasdaq.vs.pf(sprintf('1_b_LASSO_N=%i.pdf', N), nasdaq.tracked, lasso.stocks,
                        sprintf('Tracking portfolio using LASSO with N=%i', length(lasso.stocks)))
      
      # plot NASDAQ vs. tracking portfolio (OLS using LASSO stocks)
      nasdaq.tracked = as.matrix(stocks[, lasso.stocks]) %*% coefs.ols
      plot.nasdaq.vs.pf(sprintf('1_b_LASSO_OLS_refit_N=%i.pdf', N), nasdaq.tracked, lasso.stocks,
                        sprintf('Tracking portfolio using LASSO stocks with N=%i (OLS refit)', length(lasso.stocks)))
     
      # plot coefficients of LASSO and OLS refit as barplot
      pdf(sprintf('../report/1_b_coefs_N=%i.pdf', N), width=6, height=6*9/16, pointsize=8)
      
      barplot(rbind(coefs.lasso, coefs.ols), beside=T, col=c('blue', 'red'),
              main='LASSO vs. OLS refit coefficients weights')
      legend('topleft', fill=c('blue', 'red'), legend=c('LASSO', 'OLS'))
       
      dev.off()
    }
  }
  
