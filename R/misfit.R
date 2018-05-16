#' MISFIT (Multiple Imputation for Sparsely-sampled Functions at Irregular Times)
#'
#'Performs MISFIT for either linear (\code{family="Gaussian"}) or logistic
#'(\code{family="Binomial"}) regression.
#'
#'@param dat
#'@param grid
#'@param K
#'@param J
#'@param family
#'@param seed
#'@param ret_allxi
#'@param user_params
#'@param fcr.args
#'@param k
#'@param nPhi
#'@details
#'@return
#'@author
#'@references
#'@example
#'@export

misfit <- function(dat,grid,K=10,J,family="Gaussian",seed=NULL,ret_allxi = F,user_params = NULL,
                   fcr.args = list(use_bam = T,niter = 1),k = 15, nPhi = NULL,
                   face.args=list(knots = 12, pve = 0.95)){
  M <- length(grid)
  N <- length(unique(dat$subj))

  if(family=="Gaussian"){
    # Estimate imputation parameters
    if(is.null(user_params)){
      ipars <- param_est_linear(obsdf,grid,T,fcr.args = fcr.args,k = k,nPhi = nPhi,face.args = face.args)
      muy <- ipars$params$muy;  var_y <- ipars$params$var_y
      Cxy <- ipars$params$Cxy;  Cx <- ipars$params$Cx
      phi <- ipars$params$phi;  lam <- ipars$params$lam
      mux <- ipars$params$mux;  var_delt <- ipars$params$var_delt
    }else{
      ipars = list(params = user_params)
      muy <- user_params$muy;  var_y <- user_params$var_y
      Cxy <- user_params$Cxy;  Cx <- user_params$Cx
      phi <- user_params$phi;  lam <- user_params$lam
      mux <- user_params$mux;  var_delt <- user_params$var_delt
    }

    ## Multitple Imputation, Conditional on outcome
    xi_all <- cond_imp_lm(dat,workGrid = grid,k = K,impute_type = "Multiple",
                          muy = muy,var_y = var_y,Cxy = Cxy,var_delt = var_delt,
                          Cx = Cx,mux = mux,phi = phi,lam = lam,seed = seed)
    xihat <- xi_all[,1:J,]

    # Estimate X's from imputed Xi
    Xall <- array(0,c(N,M,K))
    Xall <- sweep(Xall,MARGIN = c(1,2),FUN = "+",STATS = matrix(mux,N,M,byrow = T))
    if(J==1){
      for(i in 1:K){  ################# probably a better way to do this loop; could possibly do without the loop?
        Xall <- sweep(Xall,1:2,STATS = xihat[,i]%*%t(phi[,1])/K,FUN = "+")
      }
    }else{
      for(i in 1:K){
        Xall[,,i] <- Xall[,,i] + xihat[,,i]%*%t(phi[,1:J])
      }
    }
    Xhat <- apply(Xall,c(1,2),mean)

    # Estimate Beta
    bhat <- matrix(NA,nrow = J,ncol = K)
    beta.hat.mat <- matrix(NA,nrow = M,ncol = K)
    beta.var <- array(NA,dim = c(M,M,K))
    alpha <- numeric(K)
    for(i in 1:K){
      fit <- if(J==1){lm(y~xihat[,i])}else{lm(y~xihat[,,i])}
      veps <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
      bhat[,i] <- coef(fit)[-1]
      if(J==1){
        beta.hat.mat[,i] <- phi[,1]*bhat[,i]
        beta.var[,,i] <- phi[,1]%*%t(phi[,1])*veps/c(t(xihat[,i])%*%xihat[,i])
        alpha[i] <- coef(fit)[1] - mean(phi[,1]*mux)*bhat[i]
      }else{
        beta.hat.mat[,i] <- phi[,1:J]%*%bhat[,i]
        beta.var[,,i] <- phi[,1:J]%*%solve(t(xihat[,,i])%*%xihat[,,i])%*%t(phi[,1:J])*veps
        alpha[i] <- coef(fit)[1] - sum((phi[,1:J]*mux)%*%bhat[,i])/M
      }
    }
    beta.hat <- rowMeans(beta.hat.mat)
    alpha.hat <- mean(alpha)
    W <- apply(beta.var,c(1,2),mean)
    B <- (beta.hat.mat - beta.hat)%*%t(beta.hat.mat - beta.hat)/(K-1)
    Cbeta <- W + ((K+1)/K)*B
    ebeta <- eigen(Cbeta)$values
    beta_phi <- eigen(Cbeta)$vectors

    # p-value
    Tb <- sum(beta.hat^2)
    # sum(ebeta[ebeta>0]*qchisq(0.95,1,lower.tail = T))
    pvnorm <- imhof(Tb,ebeta[ebeta>0])[[1]]
    pvnorm <- ifelse(pvnorm <0 ,0,pvnorm)

    if(!ret_allxi){
      xi_all <- apply(xi_all,c(1,2),mean)
    }

  }else if(family=="Binomial"){
    sum_y <- dat %>% group_by(subj) %>% summarise(y = first(y)) %>% summarise(vy = var(y),my = mean(y))
    muy <- sum_y[['my']]
    vary <- sum_y[['vy']]

    # Estimate imputation parameters
    if(is.null(user_params)){
      ipars <- param_est_logistic(obsdf,grid,cond.y = T,p = muy,fcr.args = fcr.args,
                                  k = k,nPhi = nPhi,face.args = face.args)
      mu0 <- ipars$params$mu0;  mu1 <- ipars$params$mu1
      var_delt <- ipars$params$var_delt;  Cx <- ipars$params$Cx
      phi <- ipars$params$phi;  lam <- ipars$params$lam
    }else{
      ipars = list(params = user_params)
      mu0 <- user_params$mu0;  mu1 <- user_params$mu1
      var_delt <- user_params$var_delt;  Cx <- user_params$Cx
      phi <- user_params$phi;  lam <- user_params$lam
    }

    ## Multitple Imputation, Conditional on outcome
    xi_all <- cond_imp_logistic(dat,workGrid = grid,K = K,seed = seed,impute_type = "Multiple",
                                mu0 = mu0,mu1 = mu1,var_delt = var_delt,Cx = Cx,phi = phi,lam = lam)
    xihat <- xi_all[,1:J,]

    # Estimate X's from imputed Xi
    Xall <- array(NA,c(N,M,K))
    if(J==1){
      for(i in 1:K){
        Xall[,,i] <- xihat[,i]%*%t(phi[,1])
        for(j in 1:N){
          mu_y <- if(y[j]==0){mu0}else{mu1}
          Xall[j,,i] <- Xall[j,,i] + mu_y
        }
      }
    }else{
      for(i in 1:K){
        Xall[,,i] <- xihat[,,i]%*%t(phi[,1:J])
        for(j in 1:N){
          mu_y <- if(y[j]==0){mu0}else{mu1}
          Xall[j,,i] <- Xall[j,,i] + mu_y
        }
      }
    }
    Xhat <- apply(Xall,c(1,2),mean)

    # Add back the means to xihat
    Xitilde <- xihat
    if(J==1){
      mean0j <- mean(mu0*phi[,1])
      mean1j <- mean(mu1*phi[,1])
      Xitilde[which(y==0),] <- Xitilde[which(y==0),] + mean0j
      Xitilde[which(y==1),] <- Xitilde[which(y==1),] + mean1j
    }else{
      mean0j <- colMeans(mu0*phi[,1:J])
      mean1j <- colMeans(mu1*phi[,1:J])
      Xitilde[which(y==0),,] <- aperm(apply(Xitilde[which(y==0),,],c(1,3),function(x) x + mean0j),c(2,1,3))
      Xitilde[which(y==1),,] <- aperm(apply(Xitilde[which(y==1),,],c(1,3),function(x) x + mean1j),c(2,1,3))
    }

    # Estimate Beta
    bhat <- matrix(NA,J,K)
    beta.hat.mat <- matrix(NA,M,K)
    beta.var <- array(NA,dim = c(M,M,K))
    alpha <- numeric(K)
    for(i in 1:K){
      if(J==1){
        fit <- glm(y~c(Xitilde[,i]),family = "binomial")
        beta.var[,,i] <- phi[,1]%*%solve(t(Xitilde[,i])%*%Xitilde[,i])%*%t(phi[,1])*vary
      }else{
        fit <- glm(y~Xitilde[,,i],family = "binomial")
        beta.var[,,i] <- phi[,1:J]%*%solve(t(Xitilde[,,i])%*%Xitilde[,,i])%*%t(phi[,1:J])*vary
      }
      bhat[,i] <- coef(fit)[-1]
      alpha[i] <- log(muy) - log(1 - muy) - sum(((t(phi[,1:J])%*%(mu1 - mu0)/M)^2)/(lam[1:J]^2))/2
      beta.hat.mat[,i] <- if(J==1){phi[,1]*bhat[,i]}else{phi[,1:J]%*%bhat[,i]}
    }
    beta.hat <- rowMeans(beta.hat.mat)
    alpha.hat <- mean(alpha)

    # confidence interval
    W <- apply(beta.var,c(1,2),mean)
    B <- (beta.hat.mat - beta.hat)%*%t(beta.hat.mat - beta.hat)/(K-1)
    Cbeta <- W + ((K+1)/K)*B
    nu <- (K - 1)*(1+(diag(W))/((1+1/K)*diag(B)))^2 # adjusted d.f. for MI
    tbeta <- qt(p = c(.975),df = nu)
    # beta_muc_lwr <- beta.hat - tbeta*sqrt(diag(Cbeta))
    # beta_muc_upr <- beta.hat + tbeta*sqrt(diag(Cbeta))
    # plot(grid,beta.hat,type = 'l',xlab = "Beta",main = "Multiple Conditional")
    # lines(grid,beta_muc_lwr,lty = 2)
    # lines(grid,beta_muc_upr,lty = 2)

    # P-value
    # Imhof approximation
    ebeta <- eigen(Cbeta)$values
    Tmuc <- sum(beta.hat^2)
    pvnorm <- imhof(Tmuc,ebeta[ebeta > 0])[[1]]
    pvnorm <- ifelse(pvnorm <0 ,0,pvnorm)

    if(!ret_allxi){
      if(J==1){
        xi_all <- apply(xi_all,1,mean)
      }else
        xi_all <- apply(xi_all,c(1,2),mean)
    }
  }

  out <- list(params = ipars[['params']], xiest = xi_all, Xest = Xhat, pvnorm = pvnorm,
              beta.hat = beta.hat, alpha.hat = alpha.hat, Cbeta = Cbeta, W = W, B = B)
  return(out)
}
