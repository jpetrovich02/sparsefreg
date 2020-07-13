#'@return
#'@author
#'@references
#'@example
#'@export

param_est_fcr <- function(dat,workGrid,muy,var_y,fcr.args,k,face.args){#,nPhi){
  rhs <- paste("~ ", "s(argvals, k =", k,", bs = \"ps\") + s(argvals, by = y, k =", k,", bs = \"ps\")")
  model <- update.formula(rhs, "X ~ .")
  fit <- do.call("fcr",c(list(formula = model, data = dat, subj = "subj",
                              argvals = "argvals", face.args = face.args,
                              argvals.new = workGrid),#,nPhi = nPhi),
                         fcr.args))
  Cb <- fit$face.object$Chat.new
  ci <- match(workGrid,fit$face.object$argvals.new)
  Cb <- Cb[ci,ci]
  var_delt <- fit$face.object$var.error.hat[1]

  predcoefs <- predict(fit,newdata = data.frame("argvals" = workGrid,"y" = 1,"subj"=1),type='iterms',se.fit=TRUE)
  predcoefs <- predcoefs$insample_predictions
  predcoefs$fit[,1] <- predcoefs$fit[,1] + attributes(predcoefs)$constant ## add back in \hat{\beta_0}
  f1 <- predcoefs$fit[,1]
  f2 <- predcoefs$fit[,2]
  mux <- f1 + muy*f2
  Cx <- var_y*f2%*%t(f2) + Cb
  Cxeig <- eigen(Cx)
  Cxy <- var_y*f2
  lam <- Cxeig$values/length(workGrid)
  phi <- Cxeig$vectors*sqrt(length(workGrid))

  # Output parameters
  list(mux = mux,Cx = Cx,lam = lam,phi = phi,var_delt = var_delt,Cxy = Cxy)
}
