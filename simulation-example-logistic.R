# Load Packages
library(MASS)
library(fdapace)
library(tidyverse)
library(fcr)
# source("../uncond_imp.R")
# source("../cond_imp.R")
# source("../param_est_logistic.R")

## Simulation Parameters
sims <- 1
M <- 100 # grid size
w <- 1

## Data set-up
var_delt0 <- 0.5
var_delt1 <- 0.5
grid <- seq(from=0,to=1,length.out = M)
nfpc <- 2 # number of fpcs used for mu1
p <- 0.5
mu0 <- rep(0,M)
C0_f<-function(t,s,sig2=1,rho=0.5){ # Matern covariance function with nu = 5/2
  d <- abs(outer(t,s,"-"))
  tmp2 <- sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
C0 <- C0_f(grid,grid)
lam <- eigen(C0,symmetric = T)$values/M
# lam <- eigen(C0,symmetric = T)$values
phi <- eigen(C0,symmetric = T)$vectors*sqrt(M)
# phi <- eigen(C0,symmetric = T)$vectors
if(nfpc==1){
  mu1 <- phi[,1]*w
}else{
  mu1 <- rowSums(phi[,1:nfpc])*w
}
C1 <- C0

## Slope Function
if(nfpc==1){
  beta <- phi[,1]*(1/lam[1])*w
}else{
  beta <- phi[,1:nfpc]%*%(1/lam[1:nfpc])*w
  # beta <- phi%*%(c(t(mu1)%*%phi/M)/lam)
}

v <- c("N","m","J")
for(mm in 1:length(v)){

  set.seed(123)
  ## Simulation Parameters
  N <- 400; m <-2; J <- 2; k <- 10
  N_all <- NULL; m_all <- NULL
  J_all <- NULL; k_all <- NULL
  if(v[mm]=="N"){
    N_all <- c(100,200,400,800)
    rm(N)
    lv <- length(N_all)
  }else if(v[mm]=="m"){
    m_all <- c(2,5,10,20)
    rm(m)
    lv <- length(m_all)
  } else if(v[mm]=="J"){
    J_all <- 1:6
    rm(J)
    lv <- length(J_all)
  }else if(v[mm]=="k"){
    k_all <- c(5,10,20,50,100)
    rm(k)
    lv <- length(k_all)
  }
  #-------------------------------
  ## Comparison metrics
  comp_names <- c("MeC","MuC","MeU","MuU","Pace","True")
  ln <- length(comp_names)
  name_list <- list(NULL,comp_names)
  mse_b <- matrix(NA,sims,ln,dimnames = name_list)
  mse_x <- matrix(NA,sims,ln,dimnames = name_list)
  betas <- array(NA,dim = c(M,sims,ln),dimnames = list(NULL,NULL,comp_names))
  mse_x_med <- matrix(NA,lv,ln,dimnames = name_list)
  mse_x_mean <- matrix(NA,lv,ln,dimnames = name_list)
  mse_b_med <- matrix(NA,lv,ln,dimnames = name_list)
  mse_b_mean <- matrix(NA,lv,ln,dimnames = name_list)
  mse_b_tmean <- matrix(NA,lv,ln,dimnames = name_list)
  bias <- matrix(NA,lv,ln,dimnames = name_list)
  variance <- matrix(NA,lv,ln,dimnames = name_list)
  #-----------------------------------------
  for(jj in 1:lv){
    if(!is.null(N_all)){N <- N_all[jj]}
    if(!is.null(m_all)){m <- m_all[jj]}
    if(!is.null(J_all)){J <- J_all[jj]}
    if(!is.null(k_all)){k <- k_all[jj]}

    for(ii in 1:sims){
      ## Simulate Data
      y <- sort(rbinom(N,1,p))
      N_0 <- sum(y==0)
      N_1 <- sum(y==1)
      Xs_0 <- mvrnorm(N_0,mu0,C0)
      Xs_1 <- mvrnorm(N_1,mu1,C1)
      Xn_0 <- Xs_0 + rnorm(N_0*M,sd = sqrt(var_delt0))
      Xn_1 <- Xs_1 + rnorm(N_1*M,sd = sqrt(var_delt1))
      X_s <- rbind(Xs_0,Xs_1)
      X_comp <- rbind(Xn_0,Xn_1)

      Xi0 <- (Xs_0 - mu0)%*%phi/M
      Xi1 <- (Xs_1 - mu1)%*%phi/M
      Xi <- rbind(Xi0,Xi1)
      #---------------------------------------------------------------------------------------------
      X_mat<-matrix(nrow=N,ncol=m)
      T_mat<-matrix(nrow=N,ncol=m)
      ind_obs<-matrix(nrow=N,ncol=m)

      for(i in 1:N){
        ind_obs[i,]<-sort(sample(1:M,m,replace=FALSE))
        X_mat[i,]<-X_comp[i,ind_obs[i,]]
        T_mat[i,]<-grid[ind_obs[i,]]
      }

      spt<-c(1,N_0+1)
      ind_obs[spt,1] = 1; ind_obs[spt,m] = M
      X_mat[spt,]<-X_comp[spt,ind_obs[spt[1],]]
      T_mat[spt,]<-rbind(grid[ind_obs[spt[1],]],grid[ind_obs[spt[1],]])

      ## Create data frame for observed data
      obsdf <- data.frame("X" = c(t(X_mat)),"argvals" = c(t(T_mat)),
                          "y" = rep(y,each = m),"subj" = rep(1:N,each = m))
      #---------------------------------------------------------------------------------------------
      ## Estimate distribution parameters
      phat <- mean(y)
      np <- min(floor((nrow(obsdf) - 1)/N),J)
      cpar <- param_est(obsdf,grid,cond.y = T,p = phat,nPhi = np) # parameters for distribution of scores conditional on x_ij and Y_i
      upar <- param_est(obsdf,grid,cond.y = F) # parameters for distribution of scores conditional on x_ij but not Y_i
      #------------------------------------------------------------------------------------------------------
      # PACE (Mean imputation, unconditional on outcome):
      X_obs <- split(t(X_mat), c(col(t(X_mat))))
      T_obs <- split(t(T_mat), c(col(t(T_mat))))
      pace <- FPCA(Ly = X_obs,Lt = T_obs,
                   optns = list(dataType = "Sparse",methodSelectK=J,nRegGrid = M,error=TRUE))
      X_pace <- fitted.values(pace)
      mse_x[ii,"Pace"] <- mean(rowMeans((X_s-X_pace)^2))

      # Beta estimate
      fit_pace <- if(J==1){
        glm(y~c(pace$xiEst[,1]),family = "binomial")
      }else{
        glm(y~pace$xiEst[,1:J],family = "binomial")
      }
      b.hat.pace <- coef(fit_pace)[-1]
      beta.hat.pace <- if(J==1){pace$phi[,1]*b.hat.pace}else{pace$phi[,1:J]%*%b.hat.pace}
      betas[,ii,"Pace"] <- beta.hat.pace
      mse_b[ii,"Pace"] <- sum((beta.hat.pace-beta)^2)/M

      #---------------------------------------------------------------------------------------------
      ## Mean Imputation, Unconditional on outcome
      Xi_MeU_all <- uncond_imp(obsdf,workGrid = grid,k = k,impute_type = "Mean",Cx = upar$params$Cx,
                               phi = upar$params$phi,lam = upar$params$lam,mux = upar$params$mux,
                               var_delt = upar$params$var_delt,seed = ii)

      # Beta estimate
      Xi_MeU <- Xi_MeU_all[,1:J]
      fit_meu <- if(J==1){glm(y~c(Xi_MeU),family = "binomial")}else{glm(y~Xi_MeU,family = "binomial")}
      b.hat.meu <- coef(fit_meu)[-1]
      beta.hat.meu <- if(J==1){upar$params$phi[,1]*b.hat.meu}else{upar$params$phi[,1:J]%*%b.hat.meu}
      betas[,ii,"MeU"] <- beta.hat.meu
      mse_b[ii,"MeU"] <- sum((beta.hat.meu-beta)^2)/M

      # X estimate
      X_MeU <- if(J==1){
        sweep(matrix(rep(upar$params$phi[,1],N),N,M,byrow = T)*c(Xi_MeU),2,upar$params$mux,FUN = "+")
      }else{
        t(upar$params$mux + upar$params$phi[,1:J]%*%t(Xi_MeU))
      }
      mse_x[ii,"MeU"] <- mean(rowMeans((X_s-X_MeU)^2))

      #---------------------------------------------------------------------------------------------------------
      ## Mean Imputation, Conditional on outcome
      Xi_MeC_all <- cond_imp(obsdf,workGrid = grid,k = k,seed = ii,impute_type = "Mean",
                             mu0 = cpar$params$mu0,mu1 = cpar$params$mu1,var_delt = cpar$params$var_delt,
                             Cx = cpar$params$Cx,phi = cpar$params$phi,lam = cpar$params$lam)
      Xi_MeC <- Xi_MeC_all[,1:J]
      X_MeC <- Xi_MeC%*%t(cpar$params$phi[,1:J])
      for(i in 1:N){
        mu_y <- if(y[i]==0){cpar$params$mu0}else{cpar$params$mu1}
        X_MeC[i,] <- X_MeC[i,] + mu_y
      }
      mse_x[ii,"MeC"] <- mean(rowMeans((X_s - X_MeC)^2))

      # Add back the means
      Xi_MeC_hat <- Xi_MeC
      if(J==1){
        mean0j <- mean(cpar$params$mu0*cpar$params$phi[,1])
        mean1j <- mean(cpar$params$mu1*cpar$params$phi[,1])
        Xi_MeC_hat[which(y==0)] <- Xi_MeC_hat[which(y==0)] + mean0j
        Xi_MeC_hat[which(y==1)] <- Xi_MeC_hat[which(y==1)] + mean1j
      }else{
        mean0j <- colMeans(cpar$params$mu0*cpar$params$phi[,1:J])
        mean1j <- colMeans(cpar$params$mu1*cpar$params$phi[,1:J])
        Xi_MeC_hat[which(y==0),] <- Xi_MeC_hat[which(y==0),] + matrix(mean0j,sum(y==0),J,byrow = T)
        Xi_MeC_hat[which(y==1),] <- Xi_MeC_hat[which(y==1),] + matrix(mean1j,sum(y==1),J,byrow = T)
      }

      # Estimate Beta
      fit_mec <- glm(y~Xi_MeC_hat,family = "binomial")
      b.hat.mec <- coef(fit_mec)[-1]
      beta.hat.mec <- if(J==1){cpar$params$phi[,1]*b.hat.mec}else{cpar$params$phi[,1:J]%*%b.hat.mec}
      betas[,ii,"MeC"] <- beta.hat.mec
      mse_b[ii,"MeC"] <- sum((beta.hat.mec-beta)^2)/M

      #---------------------------------------------------------------------------------------------
      ## Multitple Imputation, Conditional on outcome
      Xi_MuC_all <- cond_imp(obsdf,workGrid = grid,k = k,seed = ii,impute_type = "Multiple",
                             mu0 = cpar$params$mu0,mu1 = cpar$params$mu1,var_delt = cpar$params$var_delt,
                             Cx = cpar$params$Cx,phi = cpar$params$phi,lam = cpar$params$lam)
      Xi_MuC_imp <- Xi_MuC_all[,1:J,]

      # Estimate X's from imputed Xi
      X_MuC_all <- array(NA,c(N,M,k))
      if(J==1){
        for(i in 1:k){
          X_MuC_all[,,i] <- Xi_MuC_imp[,i]%*%t(cpar$params$phi[,1])
          for(j in 1:N){
            mu_y <- if(y[j]==0){cpar$params$mu0}else{cpar$params$mu1}
            X_MuC_all[j,,i] <- X_MuC_all[j,,i] + mu_y
          }
        }
      }else{
        for(i in 1:k){
          X_MuC_all[,,i] <- Xi_MuC_imp[,,i]%*%t(cpar$params$phi[,1:J])
          for(j in 1:N){
            mu_y <- if(y[j]==0){cpar$params$mu0}else{cpar$params$mu1}
            X_MuC_all[j,,i] <- X_MuC_all[j,,i] + mu_y
          }
        }
      }
      X_MuC <- apply(X_MuC_all,c(1,2),mean)
      mse_x[ii,"MuC"] <- mean(rowMeans((X_s-X_MuC)^2))

      # Add back the means
      Xi_MuC_hat <- Xi_MuC_imp
      if(J==1){
        mean0j <- mean(cpar$params$mu0*cpar$params$phi[,1:J])
        mean1j <- mean(cpar$params$mu1*cpar$params$phi[,1:J])
        Xi_MuC_hat[which(y==0),] <- Xi_MuC_hat[which(y==0),] + mean0j
        Xi_MuC_hat[which(y==1),] <- Xi_MuC_hat[which(y==1),] + mean1j
      }else{
        mean0j <- colMeans(cpar$params$mu0*cpar$params$phi[,1:J])
        mean1j <- colMeans(cpar$params$mu1*cpar$params$phi[,1:J])
        Xi_MuC_hat[which(y==0),,] <- aperm(apply(Xi_MuC_hat[which(y==0),,],c(1,3),function(x) x + mean0j),c(2,1,3))
        Xi_MuC_hat[which(y==1),,] <- aperm(apply(Xi_MuC_hat[which(y==1),,],c(1,3),function(x) x + mean1j),c(2,1,3))
      }

      # Estimate Beta
      b.hat.muc <- matrix(NA,J,k)
      beta.hat.mat.muc <- matrix(NA,M,k)
      for(i in 1:k){
        fit_muc <- if(J==1){glm(y~c(Xi_MuC_hat[,i]),family = "binomial")
        }else{glm(y~Xi_MuC_hat[,,i],family = "binomial")}
        b.hat.muc[,i] <- coef(fit_muc)[-1]
        beta.hat.mat.muc[,i] <- if(J==1){cpar$params$phi[,1]*b.hat.muc[,i]}else{cpar$params$phi[,1:J]%*%b.hat.muc[,i]}
      }
      beta.hat.muc <- rowMeans(beta.hat.mat.muc)
      betas[,ii,"MuC"] <- beta.hat.muc
      mse_b[ii,"MuC"] <- sum((beta.hat.muc-beta)^2)/M

      #----------------------------------------------------------------------------------------
      ## Multiple Imputation, Unconditional on outcome
      Xi_MuU_all <- uncond_imp(obsdf,workGrid = grid,k = k,impute_type = "Multiple",Cx = upar$params$Cx,
                               phi = upar$params$phi,lam = upar$params$lam,mux = upar$params$mux,
                               var_delt = upar$params$var_delt,seed = ii)

      # Beta estimate
      Xi_MuU <- Xi_MuU_all[,1:J,]
      b.hat.muu <- matrix(NA,nrow = J,ncol = k)
      beta.hat.mat.muu <- matrix(NA,nrow = M,ncol = k)
      for(i in 1:k){
        fit_muu <- if(J==1){glm(y~c(Xi_MuU[,i]),family = "binomial")
        }else{glm(y~Xi_MuU[,,i],family = "binomial")}
        b.hat.muu[,i] <- coef(fit_muu)[-1]
        beta.hat.mat.muu[,i] <- if(J==1){upar$params$phi[,1]*b.hat.muu[,i]}else{upar$params$phi[,1:J]%*%b.hat.muu[,i]}
      }
      beta.hat.muu <- rowMeans(beta.hat.mat.muu)
      betas[,ii,"MuU"] <- beta.hat.muu
      mse_b[ii,"MuU"] <- sum((beta.hat.muu-beta)^2)/M

      # X estimate
      X_MuU_all <- array(NA,c(N,M,k))
      for(i in 1:k){
        X_MuU_all[,,i] <- if(J==1){
          sweep(matrix(rep(upar$params$phi[,1],N),N,M,byrow = T)*c(Xi_MuU[,i]),2,upar$params$mux,FUN = "+")
        }else{
          t(upar$params$mux + upar$params$phi[,1:J]%*%t(Xi_MuU[,,i]))
        }
      }
      X_MuU <- apply(X_MuU_all,c(1,2),FUN = mean)
      mse_x[ii,"MuU"] <- mean(rowMeans((X_s-X_MuU)^2))

      #-----------------------------------------------------------------------------------------
      ## Beta True (approximated by J FPCs)
      fit_true <- if(J==1){glm(y~c(Xi[,1]),family = "binomial")}else{glm(y~Xi[,1:J],family = "binomial")}
      b.hat.true <- coef(fit_true)[-1]
      beta.hat.true <- if(J==1){phi[,1]*b.hat.true}else{phi[,1:J]%*%b.hat.true}
      betas[,ii,"True"] <- beta.hat.true
      mse_b[ii,"True"] <- sum((beta.hat.true - beta)^2)/M
    }
    mse_x_med[jj,] <- apply(mse_x,2,median)
    mse_x_mean[jj,] <- apply(mse_x,2,mean)
    mse_b_med[jj,] <- apply(mse_b,2,median)
    mse_b_mean[jj,] <- apply(mse_b,2,mean)
    mse_b_tmean[jj,] <- apply(mse_b,2,mean,trim = 0.05)
    bias[jj,] <- apply(betas,3,function(x){
      sum((rowMeans(x) - beta)^2)/M
    })
    variance[jj,] <- apply(betas,3,function(x){
      mean(colSums((x - rowMeans(x))^2)/M)
    })
    cat(jj)
  }
  if(!is.null(N_all)){
    Nout <- list(N_all = N_all,m = m,J = J,k = k,w = w,bias = bias,variance = variance,
                 mse_b_mean = mse_b_mean,mse_b_med = mse_b_med,
                 mse_x_mean = mse_x_mean,mse_x_med = mse_x_med,beta = beta,mu1 = mu1,mu0 = mu0,
                 C0 = C0,grid = grid,sims = sims)
  }else if(!is.null(m_all)){
    mout <- list(m_all = m_all,N = N,J = J,k = k,w = w,bias = bias,variance = variance,
                 mse_b_mean = mse_b_mean,mse_b_med = mse_b_med,
                 mse_x_mean = mse_x_mean,mse_x_med = mse_x_med,beta = beta,mu1 = mu1,mu0 = mu0,
                 C0 = C0,grid = grid,sims = sims)
  }else if(!is.null(J_all)){
    Jout <- list(J_all = J_all,N = N,m = m,k = k,w = w,bias = bias,variance = variance,
                 mse_b_mean = mse_b_mean,mse_b_med = mse_b_med,
                 mse_x_mean = mse_x_mean,mse_x_med = mse_x_med,beta = beta,mu1 = mu1,mu0 = mu0,
                 C0 = C0,grid = grid,sims = sims)
  }else if(!is.null(k_all)){
    kout <- list(k_all = k_all,N = N,m = m,J = J,w = w,bias = bias,variance = variance,
                 mse_b_mean = mse_b_mean,mse_b_med = mse_b_med,
                 mse_x_mean = mse_x_mean,mse_x_med = mse_x_med,beta = beta,mu1 = mu1,mu0 = mu0,
                 C0 = C0,grid = grid,sims = sims)
  }
}
logistic_est <- list(N = Nout,m = mout,J = Jout)
save(logistic_est,file = "logistic_est.RData")
