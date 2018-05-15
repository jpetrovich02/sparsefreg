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

misfit <- function(dat,grid,K=10,J,family="Gaussian",seed,ret_allxi = F,user_params = NULL,
                   fcr.args = list(use_bam = T,niter = 1),k=15,nPhi = NULL,
                   face.args=list(knots = 12, pve = 0.95)){
  M <- length(grid)
  N <- length(unique(dat$subj))

  if(family=="Gaussian"){
    # Estimate imputation parameters
    if(is.null(user_params)){
      ipars <- param_est_linear(obsdf,grid,T,fcr.args = fcr.args,
                                k = k,nPhi = nPhi,face.args = face.args)
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

    # p-value
    Tb <- sum(beta.hat^2)
    # sum(ebeta[ebeta>0]*qchisq(0.95,1,lower.tail = T))
    pv <- imhof(Tb,ebeta[ebeta>0])[[1]]
    pv <- ifelse(pv<0,0,pv)

    if(!ret_allxi){
      xi_all <- apply(xi_all,c(1,2),mean)
    }

  }else if(family=="Binomial"){
    sum_y <- dat %>% group_by(subj) %>% summarise(y = first(y)) %>% summarise(vy = var(y),my = mean(y))
    var_y <- sum_y[['vy']]

    # Estimate imputation parameters
    if(is.null(user_params)){
      ipars <- param_est_logistic(obsdf,grid,cond.y = T,p = sum_y[['my']],
                                  fcr.args = fcr.args,k = k,nPhi = nPhi)
      mu0 <- ipars$params$mu0;  mu1 <- ipars$params$mu1
      var_delt <- ipars$params$var_delt;  Cx <- ipars$params$Cb
      phi <- ipars$params$phi;  lam <- ipars$params$lam
    }else{
      ipars = list(params = user_params)
      mu0 <- user_params$mu0;  mu1 <- user_params$mu1
      var_delt <- user_params$var_delt;  Cx <- user_params$Cb
      phi <- user_params$phi;  lam <- user_params$lam
    }

    ## Multitple Imputation, Conditional on outcome
    xi_all <- cond_imp(dat,workGrid = grid,k = K,seed = seed,impute_type = "Multiple",
                       mu0 = mu0,mu1 = mu1,var_delt = var_delt,Cx = Cx,phi = phi,lam = lam)
    if(J==1){
      xihat <- xi_all[,1:J]
    }else{
      xihat <- xi_all[,1:J,]
    }

    # Estimate X's from imputed Xi
    Xall <- array(NA,c(N,M,K))
    for(i in 1:K){
      Xall[,,i] <- xihat[,,i]%*%t(phi[,1:J])
      for(j in 1:N){
        mu_y <- if(y[i]==0){mu0}else{mu1}
        Xall[j,,i] <- Xall[j,,i] + mu_y
      }
    }
    Xhat <- apply(Xall,c(1,2),mean)

    # Add back the means to xihat
    Xitilde <- xihat
    mean0j <- colMeans(mu0*phi[,1:J])
    mean1j <- colMeans(mu1*phi[,1:J])
    if(J==1){stop("J==1")}else{
      Xitilde[which(y==0),,] <- aperm(apply(Xitilde[which(y==0),,],c(1,3),function(x) x + mean0j),c(2,1,3))
      Xitilde[which(y==1),,] <- aperm(apply(Xitilde[which(y==1),,],c(1,3),function(x) x + mean1j),c(2,1,3))
    }
    # if(J==1){
    # Xitilde <- Xitilde + matrix(meanj,dim(Xitilde)[1],dim(Xitilde)[2])
    # }else{
    # Xitilde <- aperm(apply(Xitilde,c(1,3),function(x) x + meanj),c(2,1,3))
    # }

    # Estimate Beta
    bhat <- matrix(NA,J,K)
    beta.hat.mat <- matrix(NA,M,K)
    beta.var <- array(NA,dim = c(M,M,K))
    alpha <- numeric(K)
    for(i in 1:K){
      fit <- if(J==1){glm(y~c(Xitilde[,i]),family = "binomial")
      }else{glm(y~Xitilde[,,i],family = "binomial")}
      bhat[,i] <- coef(fit)[-1]
      beta.var[,,i] <- phi[,1:J]%*%solve(t(Xitilde[,,i])%*%Xitilde[,,i])%*%t(phi[,1:J])*var_y
      alpha[i] <- coef(fit)[1] - sum((phi[,1:J]*mux)%*%bhat[,i])/M
      beta.hat.mat[,i] <- if(J==1){phi[,1]*bhat[,i]}else{phi[,1:J]%*%bhat[,i]}
    }
    beta.hat <- rowMeans(beta.hat.mat)
    alpha.hat <- mean(alpha)

    # confidence interval
    W <- apply(beta.var,c(1,2),mean)
    B <- (beta.hat.mat - beta.hat.muc)%*%t(beta.hat.mat - beta.hat.muc)/(K-1)
    Cbeta <- W + ((K+1)/K)*B
    nu <- (K - 1)*(1+(diag(W))/((1+1/K)*diag(B)))^2 # adjusted d.f. for MI
    tbeta <- qt(p = c(.975),df = nu)
    # beta_muc_lwr <- beta.hat.muc - tbeta*sqrt(diag(Cbeta))
    # beta_muc_upr <- beta.hat.muc + tbeta*sqrt(diag(Cbeta))
    # plot(grid,beta.hat.muc,type = 'l',xlab = "Beta",main = "Multiple Conditional")
    # lines(grid,beta_muc_lwr,lty = 2)
    # lines(grid,beta_muc_upr,lty = 2)

    # P-value
    # Imhof approximation
    ebeta <- eigen(Cbeta)$values
    Tmuc <- sum(beta.hat^2)
    pv <- imhof(Tmuc,ebeta[ebeta > 0])[[1]]

    # Satterthwaite approximation
    sc <- sum(ebeta)
    sc2 <- sum(ebeta^2)
    pv_sat <- pgamma(Tmuc,shape = sc^2/(2*sc2),rate = sc/(2*sc2),lower.tail = F)

    if(!ret_allxi){
      xi_all <- apply(xi_all,c(1,2),mean)
    }
  }

  out <- list(params = ipars[['params']], xiest = xi_all, Xest = Xhat, pv = pv,
              beta.hat = beta.hat, alpha.hat = alpha.hat, Cbeta = Cbeta, W = W, B = B)
  return(out)
}
