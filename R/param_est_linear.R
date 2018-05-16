#'Imputation Parameter Estimation for Linear Regression
#'
#'Given the data, this function estimates the imputation parameters.
#'This function is a wrapper for \code{fcr} if \code{cond.y} is set to TRUE,
#'and \code{face.sparse} if \code{cond.y} is set to FALSE.
#'
#'@param dat A data frame with \eqn{n} rows (where \eqn{N} is the number of subjects,
#'       each with \eqn{m_i} observations, so that \eqn{\sum_{i=1}^N m_i = n})
#'       expected to have either 3 or 4. If \code{cond.y} is TRUE, should include
#'       4 columns, with variables 'X','y','subj', and 'argvals'. If \code{cond.y}
#'       is FALSE, only 3 columns are needed (no 'y' variable is used).
#'@param workGrid A length \eqn{M} vector of the unique desired grid points on which to evaluate the function.
#'@param cond.y A logical indicator to determine whether imputation will be done conditional on \eqn{Y}.
#'@param fcr.args A list of arguments to be passed to the underlying function \code{fcr}.
#'@param k Dimension of the smooth terms used in \code{fcr}. Only needs to be specified if \code{cond.y} is TRUE.
#'@param nPhi An integer value, indicating the number of random effects to include in the model. This is
#'passed to \code{fcr} and is only used when \code{cond.y} is TRUE. See \code{fcr} for more details.
#'@param face.args A list of arguments to be passed to the underlying function \code{face.sparse}.
#'@details
#'@return
#'@author
#'@references
#'@example
#'@export


param_est_linear <- function(dat,workGrid,cond.y=TRUE,fcr.args = list(use_bam = T,niter = 1),k = 15,nPhi = NULL,
                             face.args=list(knots = 12, pve = 0.95)){
  N <- length(unique(dat[,"subj"]))
  # cat("k is ",deparse(substitute(k)), " deparsed")
  # fit <- NULL
  if(cond.y){
    # nPhi <- min(c(floor((nrow(dat) - 2*k)/N),J))
    ks <- deparse(substitute(k))
    if(!is.null(nPhi)){fcr.args['nPhi'] <- nPhi}
    rhs <- paste("~ ", "s(argvals, k =", ks,", bs = \"ps\") + s(argvals, by = y, k =", ks,", bs = \"ps\")")
    # rhs <- paste("s(argvals, k =", ks,", bs = \"ps\") + s(argvals, by = y, k =", ks,", bs = \"ps\")")
    # model <- reformulate(response = "X",termlabels = rhs)
    model <- update.formula(rhs, "X ~ .")
    # model <- update(~ s(argavls, k = k, bs = "ps") + s(argvals, by = y, k = k, bs = "ps"),X ~ .)
    fit <- do.call("fcr",c(list(formula = model, data = substitute(dat), subj = "subj",
                              argvals = "argvals", face.args = face.args,
                              argvals.new = workGrid),
                   fcr.args))
    muy <- (dat %>% group_by(subj) %>% summarise(y = first(y)) %>% summarise(mean(y)))[[1]]
    var_y <- (dat %>% group_by(subj) %>% summarise(y = first(y)) %>% summarise(var(y)))[[1]]
    Cb <- fit$face.object$Chat.new
    ci <- match(workGrid,fit$face.object$argvals.new)
    Cb <- Cb[ci,ci]
    var_delt <- fit$face.object$var.error.hat[1]
    # pd <- plot.gam.invisible(fit)
    # f1 <- pd[[1]]$fit + fit$fit$coefficients[1]
    # f2 <- pd[[2]]$fit
    # alternatively:
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
    params <- list(mux = mux,Cx = Cx,lam = lam,phi = phi,var_delt = var_delt,
                   muy = muy,var_y = var_y,Cxy = Cxy)
    pve = fit$face.object$pve
  }else{
    ## Parameters for distribution of Xi_i|{x_ij}:
    facedf <- dat[,c("X","argvals","subj")]
    colnames(facedf) <- c("y","argvals","subj")
    start_time <- proc.time()
    fit <- do.call("face.sparse",c(list(data = facedf,argvals.new = workGrid),face.args))
    run_time <- proc.time() - start_time
    mux <- fit$mu.new
    Cx <- fit$Chat.new
    Cxeig <- eigen(Cx)
    lam <- Cxeig$values/length(workGrid)
    phi <- Cxeig$vectors*sqrt(length(workGrid))
    var_delt <- fit$var.error.new[1]
    fit$runtime <- run_time
    params <- list(mux = mux,Cx = Cx,lam = lam,phi = phi,var_delt = var_delt)
    pve = fit$pve
  }
  return(list(params = params,runtime = fit$runtime,pve = pve))
}
