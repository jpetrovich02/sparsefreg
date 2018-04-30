uncond_imp <- function(dat,workGrid,k=10,seed=NULL,impute_type="Multiple",
                       var_delt=NULL,Cx=NULL,mux=NULL,phi=NULL,lam=NULL,
                       tol=1e-05){
  if(!is.null(seed)){set.seed(seed)}
  J <- ncol(phi)
  N <- length(unique(dat[,"subj"]))
  subjs <- sort(unique(dat$subj))
  a <- phi%*%diag(lam)

  if(impute_type=="Multiple"){
    out <- array(NA,c(N,J,k))
    rownames(out) <- subjs
  }
  if(impute_type=="Mean"){
    out <- matrix(NA,N,J)
    rownames(out) <- subjs
  }
  for(i in 1:N){
    rows <- which(dat[,"subj"]==subjs[i])
    t_obs <- dat[rows,"argvals"]
    ind <- match(t_obs,workGrid)

    # Observed values of X, pooled estimate of the covariance
    x_obs <- dat[rows,"X"]

    Bi <- Cx[ind,ind] + diag(var_delt,nrow = length(x_obs))
    di <- x_obs - mux[ind]
    # ai <- phi[ind,1:J]%*%diag(lam[1:J],nrow = J)
    ai <- a[ind,]

    # mu_star <- t(ai)%*%solve(Bi)%*%di
    if(length(x_obs)==1){
      mu_star <- ai/c(Bi)*di
    }else{
      mu_star <- t(ai)%*%solve(Bi)%*%di
    }
    if(impute_type=="Mean"){
      out[i,] <- mu_star
    }

    if(impute_type=="Multiple"){
      if(length(x_obs)==1){
        sig_star <- diag(lam[1:J],nrow = J) - ai%*%t(ai)/c(Bi)
      }else{
        sig_star <- diag(lam[1:J],nrow = J) - t(ai)%*%solve(Bi)%*%ai
      }

      eig <- eigen(sig_star,symmetric = TRUE)
      eval <- eig$values
      U <- if(J==1){c(eig$vectors)}else{eig$vectors}
      if(!all(eval >= tol * abs(eval[1]))){
        eval <- Re(eval)
        eval[eval<0] <- 0
      }
      Lam <- if(J==1){eval}else{diag(eval)}
      Z <- matrix(rnorm(J*k),J,k)
      out[i,,] <- if(J==1){c(mu_star)+U*sqrt(Lam)*Z}else{c(mu_star) + U%*%sqrt(Lam)%*%Z}
      # out[i,,] <- t(mvrnorm(k,mu_star,sig_star,tol = tol))
    }
  }
  return(out)
}
