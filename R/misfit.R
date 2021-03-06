#' MISFIT (Multiple Imputation for Sparsely-sampled Functions at Irregular Times)
#'
#'Performs MISFIT for either linear (\code{family="gaussian"}) or logistic
#'(\code{family="binomial"}) regression.
#'
#'@param dat A data frame with \eqn{n} rows (where \eqn{N} is the number of subjects,
#'  each with \eqn{m_i} observations, so that \eqn{\sum_{i=1}^N m_i = n})
#'  expected to have either 3 or 4. If \code{cond.y} is TRUE, should include
#'  4 columns, with variables 'X','y','subj', and 'argvals'. If \code{cond.y}
#'  is FALSE, only 3 columns are needed (no 'y' variable is used).
#'@param grid A length \eqn{M} vector of the unique desired grid points on which to evaluate the function.
#'@param nimps An integer specifying the number of desired imputations, if \code{impute_type} is "Multiple".
#'@param J An integer specifying the number of FPCs to include in the regression model. By default (NULL),
#'  J will be chosen as the minimum number of FPCs required to explain a given percentage of variance.
#'@param pve The desired percentage of variance to be explained by the FPCs.
#'  Only used if \code{J} is not supplied. Defaults to 0.95.
#'@param family A string indicating the family of the response variable. Currently only "gaussian"
#'  (linear regression) and "binomial" (logistic regression) are supported.
#'@param seed An integer used to specify the seed. Optional, but useful for making results reproducible in the
#'  Multiple Imputation step.
#'@param impute_type A string indicating whether to use mean or multiple imputation.
#'  Only accepts "Mean" or "Multiple". Defaults to "Multiple".
#'@param cond.y A boolean indicating whehter to condition on the response variable when imputing.
#'  Defaults to TRUE.
#'@param user_params An optional list of user-defined imputation parameters. Currently, the user must provide
#'  either all necessary imputation parameters, or none. See 'Details'.
#'@param use_fcr A boolean indicating whether to use \code{\link[fcr]{fcr}} or \code{\link[fdapace]{FPCA}} when estimating the necessary imputation
#'  parameters. TRUE indicates fcr, FALSE indicates pace. Default is TRUE. See 'Details' for more discussion.
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
#'
#'  By default, use_fcr is TRUE, meaning that fcr is used to estimate imputation parameters.
#'  Using FPCA (i.e. use_fcr = FALSE) is roughly 10 times faster, at least for small to moderate data sets.
#'  For a single use of the function, this difference is not meaningful as both methods complete in under a minute.
#'  But when performing simulations, this speed difference is significant.
#'  More testing is needed to determine which method more accuartely estimates the imputation parameters. See
#'  'References' below for details on the methods used in fcr and FPCA.
#'
#'
#'@return
#'@references
#'\cite{Leroux, A., Xiao, L., Crainiceanu, C., & Checkley, W. (2018). Dynamic prediction in functional concurrent regression with an application to child growth. Statistics in medicine, 37(8), 1376-1388.}
#'
#'\cite{Yao, Fang, Hans-Georg Mueller, and Jane-Ling Wang. "Functional data analysis for sparse longitudinal data." Journal of the American Statistical Association 100, no. 470 (2005): 577-590.}
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
#'nimps <- 10 # number of imputations
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
#'misfit_out <- misfit(obsdf,grid = grid,nimps = nimps,J = J)
#'
#'}
#'
#'@export

misfit <- function(dat,grid,nimps=10,J = NULL,pve = 0.95,
                   family="gaussian",link = NULL,
                   impute_type = "Multiple",cond.y = T,seed=NULL,
                   user_params = NULL,use_fcr = TRUE,k = -1, #nPhi = NULL,
                   fcr.args = list(use_bam = T,niter = 1),
                   face.args=list(knots = 12, lower = -3, pve = 0.95)){

  ## Check arguments
  if(!is.null(user_params)){
    # issue error if any necessary params are missing

    # issue error if cond.y is TRUE and user_params is non-NULL valued but missing muy and var_y
  }
  if(!(impute_type %in% c("Multiple","Mean"))){
    # issue error if impute_type is something other than "Multiple" or "Mean"
  }
  if(!(family %in% c("gaussian","binomial"))){
    # issue error if family is something other than "gaussian" or "binomial"
    stop(paste(family,"family not recognized. Must be one of: gaussian or binomial"))
  }

  ## Set up outputs
  Xitilde <- NULL
  run.time <- list(est = NULL,imp = NULL)

  ## Data values
  M <- length(grid)
  N <- length(unique(dat$subj))
  y <- (dat %>% group_by(subj) %>% summarise(y = first(y)))$y

  ## Begin imputation
  if(family=="gaussian"){

    ## Estimate imputation parameters
    if(is.null(user_params)){
      par.est <- param_est_linear(dat,y,grid,M,cond.y,use_fcr,k = k,#nPhi = nPhi,
                                  fcr.args = fcr.args,face.args = face.args)
      ipars <- par.est[["params"]]
      run.time[["est"]] <- par.est[["runtime"]]
    }else{
      ipars <- user_params
    }
    if(is.null(J)){
      J <- which.max(cumsum(ipars$lam)/sum(ipars$lam) > pve)[1] #determine number of PCs
    }

    ## Impute Scores
    imp.start <- proc.time()
    if(cond.y){
      scores_all <- cond_imp_lm(dat,workGrid = grid,nimps = nimps,seed = seed,impute_type = impute_type,
                            muy = ipars[["muy"]],var_y = ipars[["var_y"]],Cxy = ipars[["Cxy"]],
                            var_delt = ipars[["var_delt"]],Cx = ipars[["Cx"]],mux = ipars[["mux"]],
                            phi = ipars[["phi"]],lam = ipars[["lam"]])
    }else if(!cond.y){
      scores_all <- uncond_imp(dat,workGrid = grid,nimps = nimps,seed = seed,impute_type = impute_type,
                           var_delt = ipars[["var_delt"]],Cx = ipars[["Cx"]],
                           mux = ipars[["mux"]],phi = ipars[["phi"]],lam = ipars[["lam"]])
    }
    run.time[["imp"]] <- proc.time() - imp.start

    ## Obtain regression estimates using imputed scores
    if(impute_type=="Multiple"){
      scores_imp <- scores_all[,1:J,]
      Xiest <- scores_imp

      ## Estimate X's from imputed scores
      Xall <- array(NA,c(N,M,nimps))
      for(i in 1:nimps){
        # Xall[,,i] <- t(ipars[["mux"]] + ipars[["phi"]][,1:J]%*%t(Xitilde[,,i]))
        if(J==1){
          Xall[,,i] <- sweep(Xiest[,i]%*%t(ipars[["phi"]][,1]),
                             MARGIN = 2,STATS = ipars[["mux"]],FUN = "+")
        }else{
          Xall[,,i] <- t(ipars[["mux"]] + ipars[["phi"]][,1:J]%*%t(Xiest[,,i]))
        }
      }
      Xhat <- apply(Xall,c(1,2),mean)

      ## Estimate Beta
      b.hat.mat <- matrix(NA,nrow = J,ncol = nimps)
      beta.hat.mat <- matrix(NA,nrow = M,ncol = nimps)
      beta.var <- array(NA,dim = c(M,M,nimps))
      # alpha <- numeric(nimps)
      veps <- numeric(nimps)
      for(i in 1:nimps){
        if(J==1){
          fit <- lm(y~scores_imp[,i])
          veps[i] <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
          b.hat.mat[,i] <- coef(fit)[-1]
          beta.hat.mat[,i] <- ipars[["phi"]][,1]*b.hat.mat[,i]
          # alpha[i] <- coef(fit)[1] - mean(ipars[["phi"]][,1]*ipars[["mux"]])*b.hat.mat[,i]
          beta.var[,,i] <- ipars[["phi"]][,1:J]%*%as.matrix(vcov(fit)[-1,-1])%*%t(ipars[["phi"]][,1:J])
        }else{
          fit <- lm(y~scores_imp[,,i])
          veps[i] <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
          b.hat.mat[,i] <- coef(fit)[-1]
          beta.hat.mat[,i] <- ipars[["phi"]][,1:J]%*%b.hat.mat[,i]
          # alpha[i] <- coef(fit)[1] - sum((ipars[["phi"]][,1:J]*ipars[["mux"]])%*%b.hat.mat[,i])/M
          beta.var[,,i] <- ipars[["phi"]][,1:J]%*%vcov(fit)[-1,-1]%*%t(ipars[["phi"]][,1:J])
        }
      }
      ipars[["var_eps"]] <- mean(veps)
      beta.hat <- rowMeans(beta.hat.mat)
      # alpha.hat <- mean(alpha)

      var.w <- apply(beta.var,c(1,2),mean)
      var.b <- (beta.hat.mat - beta.hat)%*%t(beta.hat.mat - beta.hat)/(nimps-1)
      var.t <- var.w + ((nimps+1)/nimps)*var.b
      Cbeta <- list(var.w = var.w,var.b = var.b,var.t = var.t)

    }else if(impute_type=="Mean"){
      Xiest <- scores_all[,1:J]

      ## Estimate X's from imputed scores
      Xhat <- if(J==1){
        sweep(matrix(rep(ipars[["phi"]][,1],N),N,M,byrow = T)*c(Xiest),2,ipars[["mux"]],FUN = "+")
      }else{
        t(c(ipars[["mux"]]) + ipars[["phi"]][,1:J]%*%t(Xiest))
      }

      ## Beta estimate
      fit <- lm(y~Xiest)
      veps <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
      ipars[["var_eps"]] <- veps
      b.hat <- coef(fit)[-1]
      if(J==1){
        beta.hat <- ipars[["phi"]][,1]*b.hat
        # alpha.hat <- coef(fit)[1] - mean(ipars[["phi"]][,1]*ipars[["mux"]])*b.hat
        beta.var <- ipars[["phi"]][,1:J]%*%as.matrix(vcov(fit)[-1,-1])%*%t(ipars[["phi"]][,1:J])
      }else{
        beta.hat <- ipars[["phi"]][,1:J]%*%b.hat
        # alpha.hat <- coef(fit)[1] - sum((ipars[["phi"]][,1:J]*ipars[["mux"]])%*%b.hat)/M
        beta.var <- ipars[["phi"]][,1:J]%*%vcov(fit)[-1,-1]%*%t(ipars[["phi"]][,1:J])
      }
      Cbeta <- beta.var

    }


  }else if(family=="binomial"){
    if(is.null(link)){
      glm.fam <- binomial(link = "logit")
    }else{
      glm.fam <- binomial(link = link)
    }
    muy <-  mean(y)

    ## Estimate imputation parameters
    if(is.null(user_params)){
      par.est <- param_est_logistic(obsdf,grid,cond.y = cond.y,p = muy,fcr.args = fcr.args,
                                  k = k,face.args = face.args)#,nPhi = nPhi)
      ipars <- par.est[["params"]]
      run.time[["est"]] <- par.est[["runtime"]]
    }else{
      ipars = user_params
    }
    if(is.null(J)){
      J <- which.max(cumsum(ipars$lam)/sum(ipars$lam) > pve)[1] #determine number of PCs
    }

    ## Impute Scores
    imp.start <- proc.time()
    if(cond.y){
      scores_all <- cond_imp_logistic(dat,workGrid = grid,nimps = nimps,seed = seed,impute_type = impute_type,
                                      mu0 = ipars[["mu0"]],mu1 = ipars[["mu1"]],var_delt = ipars[["var_delt"]],
                                      Cx = ipars[["Cx"]],phi = ipars[["phi"]],lam = ipars[["lam"]])
    }else if(!cond.y){
      scores_all <- uncond_imp(dat,workGrid = grid,nimps = nimps,seed = seed,impute_type = impute_type,
                               var_delt = ipars[["var_delt"]],Cx = ipars[["Cx"]],
                               mux = ipars[["mux"]],phi = ipars[["phi"]],lam = ipars[["lam"]])
    }
    run.time[["imp"]] <- proc.time() - imp.start

    ## Obtain regression estimates using imputed scores
    if(impute_type=="Multiple"){
      scores_imp <- scores_all[,1:J,]
      Xiest <- scores_imp

      if(cond.y){
        ## Add back the means to scores_imp
        Xitilde <- scores_imp
        if(J==1){
          mean0j <- mean(ipars[["mu0"]]*ipars[["phi"]][,1])
          mean1j <- mean(ipars[["mu1"]]*ipars[["phi"]][,1])
          Xitilde[which(y==0),] <- Xitilde[which(y==0),] + mean0j
          Xitilde[which(y==1),] <- Xitilde[which(y==1),] + mean1j
        }else{
          mean0j <- colMeans(ipars[["mu0"]]*ipars[["phi"]][,1:J])
          mean1j <- colMeans(ipars[["mu1"]]*ipars[["phi"]][,1:J])
          Xitilde[which(y==0),,] <- aperm(apply(Xitilde[which(y==0),,],c(1,3),function(x) x + mean0j),c(2,1,3))
          Xitilde[which(y==1),,] <- aperm(apply(Xitilde[which(y==1),,],c(1,3),function(x) x + mean1j),c(2,1,3))
        }

        ## Estimate X's from imputed scores
        Xall <- array(NA,c(N,M,nimps))
        if(J==1){
          for(i in 1:nimps){
            Xall[,,i] <- Xiest[,i]%*%t(ipars[["phi"]][,1])
            for(j in 1:N){
              if(y[j]==0){
                mu_y <- ipars[["mu0"]]
              }else if(y[j]==1){
                mu_y <- ipars[["mu1"]]
              }
              Xall[j,,i] <- Xall[j,,i] + mu_y
            }
          }
        }else{
          for(i in 1:nimps){
            Xall[,,i] <- Xiest[,,i]%*%t(ipars[["phi"]][,1:J])
            for(j in 1:N){
              if(y[j]==0){
                mu_y <- ipars[["mu0"]]
              }else if(y[j]==1){
                mu_y <- ipars[["mu1"]]
              }
              Xall[j,,i] <- Xall[j,,i] + mu_y
            }
          }
        }
        # for(i in 1:nimps){
        #   # Xall[,,i] <- t(ipars[["mux"]] + ipars[["phi"]][,1:J]%*%t(Xitilde[,,i]))
        #   if(J==1){
        #     Xall[,,i] <- Xitilde[,i]%*%t(ipars[["phi"]][,1])
        #   }else{
        #     Xall[,,i] <- Xitilde[,,i]%*%t(ipars[["phi"]][,1:J])
        #   }
        # }
        Xhat <- apply(Xall,c(1,2),mean)

        ## Estimate Beta
        bhat <- matrix(NA,J,nimps)
        beta.hat.mat <- matrix(NA,M,nimps)
        beta.var <- array(NA,dim = c(M,M,nimps))
        # alpha <- numeric(nimps)
        for(i in 1:nimps){
          if(J==1){
            fit <- glm(y~c(Xitilde[,i]),family = glm.fam)
            bhat[,i] <- coef(fit)[-1]
            beta.hat.mat[,i] <- ipars[["phi"]][,1]*bhat[,i]
            beta.var[,,i] <- ipars[["phi"]][,1]%*%as.matrix(vcov(fit)[-1,-1])%*%t(ipars[["phi"]][,1])
            # alpha[i] <-
          }else{
            fit <- glm(y~Xitilde[,,i],family = glm.fam)
            bhat[,i] <- coef(fit)[-1]
            beta.hat.mat[,i] <- ipars[["phi"]][,1:J]%*%bhat[,i]
            beta.var[,,i] <- ipars[["phi"]][,1:J]%*%vcov(fit)[-1,-1]%*%t(ipars[["phi"]][,1:J])
            # alpha[i] <-
          }
        }

      }else if(!cond.y){
        ## Estimate X's from imputed scores
        Xall <- array(NA,c(N,M,nimps))
        for(i in 1:nimps){
          if(J==1){
            Xall[,,i] <- t(ipars[["mux"]] +ipars[["phi"]][,1]%*%t(scores_imp[,i]))
            # Xall[,,i] <- sweep(scores_imp[,i]%*%t(ipars[["phi"]][,1]),
            #                    MARGIN = 2,STATS = ipars[["mux"]],FUN = "+")
          }else{
            Xall[,,i] <- t(ipars[["mux"]] + ipars[["phi"]][,1:J]%*%t(scores_imp[,,i]))
          }
        }
        Xhat <- apply(Xall,c(1,2),FUN = mean)

        ## Estimate Beta
        bhat <- matrix(NA,nrow = J,ncol = nimps)
        beta.hat.mat <- matrix(NA,nrow = M,ncol = nimps)
        beta.var <- array(NA,dim = c(M,M,nimps))
        # alpha <- numeric(nimps)
        for(i in 1:nimps){
          if(J==1){
            fit <- glm(y ~ scores_imp[,i],family = glm.fam)
            bhat[,i] <- coef(fit)[-1]
            beta.hat.mat[,i] <- ipars[["phi"]][,1]*bhat[,i]
            beta.var[,,i] <- ipars[["phi"]][,1]%*%as.matrix(vcov(fit)[-1,-1])%*%t(ipars[["phi"]][,1])
            # alpha[i] <-
          }else{
            fit <- glm(y ~ scores_imp[,,i],family = glm.fam)
            bhat[,i] <- coef(fit)[-1]
            beta.hat.mat[,i] <- ipars[["phi"]][,1:J]%*%bhat[,i]
            beta.var[,,i] <- ipars[["phi"]][,1:J]%*%vcov(fit)[-1,-1]%*%t(ipars[["phi"]][,1:J])
            # alpha[i] <-
          }
        }
      }

      beta.hat <- rowMeans(beta.hat.mat)
      # alpha.hat <- mean(alpha)

      var.w <- apply(beta.var,c(1,2),mean)
      var.b <- (beta.hat.mat - beta.hat)%*%t(beta.hat.mat - beta.hat)/(nimps-1)
      var.t <- var.w + ((nimps+1)/nimps)*var.b
      Cbeta <- list(var.w = var.w,var.b = var.b,var.t = var.t)

    }else if(impute_type=="Mean"){

      Xiest <- scores_all[,1:J]

      if(cond.y){
        ## Add back means
        Xitilde <- Xiest
        if(J==1){
          mean0j <- mean(ipars[["mu0"]]*ipars[["phi"]][,1])
          mean1j <- mean(ipars[["mu1"]]*ipars[["phi"]][,1])
          Xitilde[which(y==0)] <- Xitilde[which(y==0)] + mean0j
          Xitilde[which(y==1)] <- Xitilde[which(y==1)] + mean1j
        }else{
          mean0j <- colMeans(ipars[["mu0"]]*ipars[["phi"]][,1:J])
          mean1j <- colMeans(ipars[["mu1"]]*ipars[["phi"]][,1:J])
          Xitilde[which(y==0),] <- Xitilde[which(y==0),] + matrix(mean0j,sum(y==0),J,byrow = T)
          Xitilde[which(y==1),] <- Xitilde[which(y==1),] + matrix(mean1j,sum(y==1),J,byrow = T)
        }

        ## Estimate X's from imputed scores
        # Xhat <- Xiest%*%t(ipars[["phi"]][,1:J])
        Xhat <- Xiest%*%t(ipars[["phi"]][,1:J])
        for(i in 1:N){
          if(y[i]==0){
            mu_y <- ipars[["mu0"]]
          }else if(y[i]==1){
            mu_y <- ipars[["mu1"]]
          }
          Xhat[i,] <- Xhat[i,] + mu_y
        }

        ## Estimate Beta
        fit <- glm(y ~ Xitilde,family = glm.fam)
        bhat <- coef(fit)[-1]
        if(J==1){
          beta.var <- ipars[["phi"]][,1]%*%as.matrix(vcov(fit)[-1,-1])%*%t(ipars[["phi"]][,1])
          beta.hat <- ipars[["phi"]][,1]*bhat
          # alpha.hat <-
        }else{
          beta.var <- ipars[["phi"]][,1:J]%*%vcov(fit)[-1,-1]%*%t(ipars[["phi"]][,1:J])
          beta.hat <- ipars[["phi"]][,1:J]%*%bhat
          # alpha.hat <-
        }

        Cbeta <- beta.var

      }else if(!cond.y){
        ## Estimate X's from imputed scores
        Xhat <- t(ipars[["mux"]] + ipars[["phi"]][,1:J]%*%t(Xiest))

        ## Beta estimate
        fit <- glm(y~Xiest,family = glm.fam)
        bhat <- coef(fit)[-1]
        if(J==1){
          beta.hat <- ipars[["phi"]][,1]*bhat
          beta.var <- ipars[["phi"]][,1]%*%as.matrix(vcov(fit)[-1,-1])%*%t(ipars[["phi"]][,1])
          # alpha.hat <-
        }else{
          beta.hat <- ipars[["phi"]][,1:J]%*%bhat
          beta.var <- ipars[["phi"]][,1:J]%*%vcov(fit)[-1,-1]%*%t(ipars[["phi"]][,1:J])
          # alpha.hat <-
        }
        Cbeta <- beta.var
      }

    }

  }

  out <- list(params = ipars, Xiest = Xiest, Xhat = Xhat,
              beta.hat = beta.hat, #alpha.hat = alpha.hat,
              Cbeta = Cbeta,J = J,pve = pve, Xitilde = Xitilde,
              run.time = run.time)
  return(out)
}
