#'@return
#'@author
#'@references
#'@example
#'@export

param_est_pace <- function(dat,M,cond.y,FPCA.args = NULL){
  if(is.null(FPCA.args)){
    # M <- length(workGrid)
    FPCA.args <- list(dataType = "Sparse",error = TRUE,nRegGrid = M)
  }

  X_obs <- split(dat[["X"]], dat[["subj"]])
  T_obs <- split(dat[["argvals"]], dat[["subj"]])
  pace <- FPCA(Ly = X_obs,Lt = T_obs, optns = FPCA.args)

  if(cond.y){
    dat[["XY"]] <- dat[["X"]]*dat[["y"]]
    XY_obs <- split(dat[["XY"]], dat[["subj"]])

    xy_pace <- FPCA(Ly = XY_obs,Lt = T_obs, optns = FPCA.args)
    Cxy_hat <- xy_pace$mu - pace$mu*muy_hat

    list(mux = pace$mu,Cx = pace$fittedCov,phi = pace$phi,lam = pace$lambda,
         Cxy = Cxy_hat,var_delt = pace$sigma2)
  }else{
    list(mux = pace$mu,Cx = pace$fittedCov,phi = pace$phi,
         lam = pace$lambda,var_delt = pace$sigma2)
  }
}
