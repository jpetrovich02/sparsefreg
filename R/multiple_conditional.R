multiple_conditional <- function(df,params = NULL,workGrid,fullGrid,J,nimps,seed = NULL){
  ## Multitple Imputation, Conditional on outcome

  # param_names <- c("muy","var_y","Cxy","var_delt","Cx","mux","phi","lam")
  # N <- length(unique(df$subj))
  # M <- length(unique(fullGrid))
  #
  # if(is.null(params)){
  #   cpar <- param_est(obsdf,J,grid,T) # parameters for distribution of scores conditional on x_ij and Y_i
  #   params <- cpar$params
  # }else if(!all(param_names %in% names(params))){
  #   stop("At least one needed imputation parameter is missing")
  # }

  # Should also check the values/data type of the imputation parameters themselves


  # Impute Scores
  scores_all <- cond_imp_lm(df,workGrid = workGrid,nimps = nimps,seed = seed,impute_type = "Multiple",
                            muy = params[["muy"]],var_y = params[["var_y"]],Cxy = params[["Cxy"]],
                            var_delt = params[["var_delt"]],Cx = params[["Cx"]],mux = params[["mux"]],
                            phi = params[["phi"]],lam = params[["lam"]])

  scores_imp <- scores_all[,1:J,]
  Xiest <- apply(scores_imp,MARGIN = c(1,2),mean)

  # Estimate X's from imputed scores
  X <- if(J==1){
    sweep(matrix(rep(params[["phi"]][,1],N),N,M,byrow = T)*c(Xiest),2,params[["mux"]],FUN = "+")
  }else{
    t(c(params[["mux"]]) + params[["phi"]][,1:J]%*%t(Xiest))
  }

  # Estimate Beta
  b.hat.mat <- matrix(NA,nrow = J,ncol = nimps)
  beta.hat.mat <- matrix(NA,nrow = M,ncol = nimps)
  beta.var <- array(NA,dim = c(M,M,nimps))
  veps <- numeric(nimps)
  for(i in 1:nimps){
    if(J==1){
      fit <- lm(y~scores_imp[,i])
      veps[i] <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
      b.hat.mat[,i] <- coef(fit)[-1]
      beta.hat.mat[,i] <- params[["phi"]][,1]*b.hat.mat[,i]
      beta.var[,,i] <-
        params[["phi"]][,1:J]%*%solve(t(scores_imp[,i])%*%scores_imp[,i])%*%t(params[["phi"]][,1:J])*veps[i]
    }else{
      fit <- lm(y~scores_imp[,,i])
      veps[i] <- sum((fit$residuals)^2)/(N-length(fit$coefficients))
      b.hat.mat[,i] <- coef(fit)[-1]
      beta.hat.mat[,i] <- params[["phi"]][,1:J]%*%b.hat.mat[,i]
      beta.var[,,i] <-
        params[["phi"]][,1:J]%*%solve(t(scores_imp[,,i])%*%scores_imp[,,i])%*%t(params[["phi"]][,1:J])*veps[i]
    }
  }
  params[["var_eps"]] <- mean(veps)
  beta.hat <- rowMeans(beta.hat.mat)

  v.w <- apply(beta.var,c(1,2),mean)
  v.b <- (beta.hat.mat - beta.hat)%*%t(beta.hat.mat - beta.hat)/(nimps-1)
  v.t <- v.w + ((nimps+1)/nimps)*v.b

  ret <- list(params = params, Xiest = Xiest, X = X,
              v.w = v.w, v.b = v.b, v.t = v.t, beta = beta.hat)
  ret
}
