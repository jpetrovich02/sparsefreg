#' MISFIT (Multiple Imputation for Sparsely-sampled Functions at Irregular Times)
#'
#'Performs MISFIT for either linear (\code{family="Gaussian"}) or logistic
#'(\code{family="Binomial"}) regression.
#'
#'@param dat A data frame with \eqn{n} rows (where \eqn{N} is the number of subjects,
#'  each with \eqn{m_i} observations, so that \eqn{\sum_{i=1}^N m_i = n})
#'  expected to have either 3 or 4. If \code{cond.y} is TRUE, should include
#'  4 columns, with variables 'X','y','subj', and 'argvals'. If \code{cond.y}
#'  is FALSE, only 3 columns are needed (no 'y' variable is used).
#'@param grid A length \eqn{M} vector of the unique desired grid points on which to evaluate the function.
#'@param K An integer specifying the number of desired imputations, if \code{impute_type} is "Multiple".
#'@param J An integer specifying the number of FPCs to include in the regression model.
#'@param family A string indicating the family of the response variable. Currently only "Gaussian"
#'  (linear regression) and "Binomial" (logistic regression) are supported.
#'@param seed An integer used to specify the seed. Optional, but useful for making results reproducible in the
#'  Multiple Imputation step.
#'@param ret_allxi logical-valued.. Indicates whether or not to return all \code{K}
#'  imputed sets of scores. If FALSE (default), returns the average scores across the \code{K} imputations.
#'@param user_params An optional list of user-defined imputation parameters. Currently, the user must provide
#'  either all necessary imputation parameters, or none. See 'Details'.
#'@param fcr.args A list of arguments which can be passed to \code{fcr} (for the estimation of imputaion
#'  parameters). Default is to use \code{use_bam} = T and \code{niter} = 1. The list must not contain the
#'  formula, which is constructed within \code{misfit}. See \code{fcr} for more details.
#'@param k Dimension of the smooth terms used in \code{fcr}. Default is 15.
#'@param nPhi
#'@param face.args A list of arguments to be passed to the underlying function \code{face.sparse}.
#'  Currently defaults to setting \code{knots} = 12 and \code{pve} = 0.95. See \code{face.sparse} for
#'  more details.
#'@details
#'  When using the \code{user_params} argument, the user must supply a list containing the
#'  following elements.
#'
#'  \bold{Linear Regression}:
#'  \itemize{
#'   \item 'Cx': An \eqn{M\times M} matrix representing the covariance function of \eqn{X(t)},
#'     evaluated on \code{grid}. Should not be missing any values.
#'   \item 'mux': A length \eqn{M} numeric vector representing the mean function of \eqn{X(t)},
#'     evaluated on \code{grid}. Should not be missing any values.
#'   \item 'var_delt': A single numeric value representing the variance of \eqn{\delta}, the
#'     measurement error associated with \eqn{X(t)}.
#'   \item 'muy': A single numeric value representing the mean of \eqn{Y}.
#'   \item 'lam': A numeric vector of length \emph{at least} \code{J}, representing the eigenvalues
#'     of \eqn{C_X(t,s)}, the covariance function of \eqn{X(t)}.
#'   \item 'phi': A matrix with \eqn{M} rows and \emph{at least} \code{J} columns, representing the
#'     eigenfunctions of \eqn{C_X(t,s)} (one per column) evaluated on \code{grid}. Should not be missing
#'     any values.
#'   \item 'Cxy': A numeric vector of length \eqn{M}, representing the cross-covariance \eqn{C_{XY}(t)}
#'     evaluated on \code{grid}. Should not be missing any values.
#'   \item 'var_y': A single numeric value representing the varinace of \eqn{Y}.
#'  }
#'  \bold{Logistic Regression}:
#'  \itemize{
#'   \item 'Cx': An \eqn{M\times M} matrix representing the covariance function of \eqn{X(t)},
#'     evaluated on \code{grid}. Should not be missing any values.
#'   \item 'mu0': A length \eqn{M} numeric vector representing the mean function of \eqn{X(t)|Y = 1},
#'     evaluated on \code{grid}. Should not be missing any values.
#'   \item 'mu1': A length \eqn{M} numeric vector representing the mean function of \eqn{X(t)|Y = 0},
#'     evaluated on \code{grid}. Should not be missing any values.
#'   \item 'var_delt': A single numeric value representing the variance of \eqn{\delta}, the
#'     measurement error associated with \eqn{X(t)}.
#'   \item 'lam': A numeric vector of length \emph{at least} \code{J}, representing the eigenvalues
#'     of \eqn{C_X(t,s)}, the covariance function of \eqn{X(t)}.
#'   \item 'phi': A matrix with \eqn{M} rows and \emph{at least} \code{J} columns, representing the
#'     eigenfunctions of \eqn{C_X(t,s)} (one per column) evaluated on \code{grid}. Should not be missing
#'     any values.
#'  }
#'@return
#'@references
#'@examples
#'
#'\dontrun{
#'
#'###################################################################
#'#------- Example Using MISFIT for a Linear SoF Model -------------#
#'###################################################################
#'
#'set.seed(123)
#'
#'## Data generation
#'M <- 100 # grid size
#'N <- 400 # sample size
#'m <- 2 # observations per subject
#'J <- 5 # number of FPCs to use
#'K <- 10 # number of imputations
#'var_eps <- 1 # variance of model error
#'var_delt <- 0.5 # variance of measurement error
#'grid <- seq(from=0,to=1,length.out = M)
#'mux <- rep(0,M)
#'Cx_f<-function(t,s,sig2=1,rho=0.5){ # Matern covariance function with nu = 5/2
#'  d <- abs(outer(t,s,"-"))
#'  tmp2 <- sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
#'Cx <- Cx_f(grid,grid)
#'lam <- eigen(Cx,symmetric = T)$values/M
#'phi <- eigen(Cx,symmetric = T)$vectors*sqrt(M)
#'
#'beta <- 10*(sin(2*pi*grid)+1)
#'alpha <- 0
#'
#'X_s <- mvrnorm(N,mux,Cx)
#'X_comp <- X_s + rnorm(N*M,sd = sqrt(var_delt))
#'Xi <- (X_s-mux)%*%phi/M
#'eps <- rnorm(N,0,sd = sqrt(var_eps))
#'y <- c(alpha + X_s%*%beta/M + eps)
#'
#'X_mat<-matrix(nrow=N,ncol=m)
#'T_mat<-matrix(nrow=N,ncol=m)
#'ind_obs<-matrix(nrow=N,ncol=m)
#'
#'for(i in 1:N){
#'  ind_obs[i,]<-sort(sample(1:M,m,replace=FALSE))
#'  X_mat[i,]<-X_comp[i,ind_obs[i,]]
#'  T_mat[i,]<-grid[ind_obs[i,]]
#'}
#'
#'spt<-1
#'ind_obs[spt,1] = 1; ind_obs[spt,m] = M
#'X_mat[spt,]<-X_comp[spt,ind_obs[spt,]]
#'T_mat[spt,]<-grid[ind_obs[spt,]]
#'
#'## Create data frame for observed data
#'obsdf <- data.frame("X" = c(t(X_mat)),"argvals" = c(t(T_mat)),
#'                    "y" = rep(y,each = m),"subj" = rep(1:N,each = m))
#'
#'misfit_out <- misfit(obsdf,grid = grid,K = K,J = J,family = "Gaussian",user_params = NULL,k = 12)
#'
#'}
#'
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
