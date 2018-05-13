#' MISFIT (Multiple Imputation for Sparsely-sampled Functions at Irregular Times)
#'
#'

misfit <- function(dat,grid,K,J,family="Gaussian"){
  if(family=="Gaussian"){
    ipars <- param_est_linear(obsdf,J,grid,T)
    muy <- ipars$params$muy;  var_y <- ipars$params$var_y
    Cxy <- ipars$params$Cxy;  Cx <- ipars$params$Cx
    phi <- ipars$params$phi;  lam <- ipars$params$lam
    mux <- ipars$params$mux;  var_delt <- ipars$params$var_delt

    ## Multitple Imputation, Conditional on outcome
    xi_all <- cond_imp_lm(obsdf,workGrid = grid,k = K,impute_type = "Multiple",
                              muy = muy,var_y = var_y,Cxy = Cxy,var_delt = var_delt,
                              Cx = Cx,mux = mux,phi = phi,lam = lam)
    xihat <- xi_all[,1:J,]

    # Estimate Beta
    bhat <- matrix(NA,nrow = J,ncol = K)
    beta.hat.mat <- matrix(NA,nrow = M,ncol = K)
    beta.var <- array(NA,dim = c(M,M,K))
    for(i in 1:K){
      fit <- if(J==1){lm(y~xihat[,i])}else{lm(y~xihat[,,i])}
      veps_muc <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
      bhat[,i] <- coef(fit)[-1]
      if(J==1){
        beta.hat.mat[,i] <- phi[,1]*bhat[,i]
        beta.var[,,i] <- phi[,1]%*%t(phi[,1])*veps_muc/c(t(xihat[,,i])%*%xihat[,,i])
      }else{
        beta.hat.mat[,i] <- phi[,1:J]%*%bhat[,i]
        beta.var[,,i] <- phi[,1:J]%*%solve(t(xihat[,,i])%*%xihat[,,i])%*%t(phi[,1:J])*veps_muc
      }
    }
    beta.hat <- rowMeans(beta.hat.mat)
    v.w <- apply(beta.var,c(1,2),mean)
    v.b <- (beta.hat.mat - beta.hat)%*%t(beta.hat.mat - beta.hat)/(K-1)
    v.t <- v.w + ((K+1)/K)*v.b
    v.evals <- eigen(v.t)$values

    # p-value
    Tb <- sum(beta.hat^2)
    # sum(v.evals[v.evals>0]*qchisq(0.95,1,lower.tail = T))
    pv <- imhof(Tb,v.evals[v.evals>0])[[1]]
    pv <- ifelse(pv<0,0,pv)
  }else if(family=="Binomial"){

  }
}
