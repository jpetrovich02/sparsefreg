#' Conditional Imputation for Logistic Scalar-on-Function Regression
#'
#'Given the imputation parameters, this function imputes the scores using
#'their distribution conditional on the observed values of the curves and the
#'response variable. In the case of logistic regression, this essentially means
#'that imputation is performed separately within each group of the response
#'variable. Muliple or mean imputation can be chosen.
#'
#'@param dat An \eqn{n \times 4} data frame (where \eqn{N} is the number of subjects,
#'       each with \eqn{m_i} observations, so that \eqn{\sum_{i=1}^N m_i = n})
#'       expected to have variables 'X','y','subj', and 'argvals'.
#'@param workGrid A vector of the unique desired grid points on which to evaluate the function.
#'       The length of this vector will be called \eqn{M}.
#'@param nimps An integer specifying the number of desired imputations, if \code{impute_type} is "Multiple".
#'@param seed A numeric value used to set the seed for reproducibility (useful for multiple imputation).
#'@param impute_type A string used to choose between mean and multiple imputation. Should be one of "Mean"
#'       or "Multiple".
#'@param mu0 A numeric vector of length \code{M} specifying the mean function for the 0 group, evaluated at \code{workGrid}.
#'@param mu1 A numeric vector of length \code{M} specifying the mean function for the 1 group, evaluated at \code{workGrid}.
#'@param var_delt A number representing \eqn{\sigma^2_\delta}, the variance of the noise.
#'@param Cx An \eqn{M \times M} matrix for the covariance function of \eqn{X}, evaluated at \code{workGrid}.
#'@param phi An \eqn{M \times J} matrix whose columns are the \eqn{J} eigenfunctions of \eqn{C_X}, each evaluate at \code{workGrid}.
#'@param lam A length-\eqn{J} numeric vector containing the eigenvalues of \eqn{C_X}.
#'@param tol A (small) numerical value that sets the tolerance for trimming the eigenvalues of the conditional covariance.
#'@details The variables of \code{dat} should be specified as follows: 'X' specifies the observed values of the curves
#'(no missing values here); 'y' should be a vector of 1's and 0's such that the values of 'y' are the same for each subject;
#''subj' should contain unique numeric identifiers for each subject, and 'argvals' should indicate the time point at which
#'each observation was made (note that these values should be a subset of \code{workGrid}).
#'@return Either a \eqn{N\times J} matrix of imputed scores if \code{impute_type} is set to "Mean", or a 3-dimensional array
#'of dimension \eqn{N\times J\times nimps} if \code{impute_type} is set to "Multiple".
#'@author Jusitn Petrovich, \email{jpetrovich02@@gmail.com}
#'@references
#'@example
#'@export

cond_imp_logistic <- function(dat,workGrid,nimps=10,seed=NULL,impute_type="Multiple",
                          mu0=NULL,mu1=NULL,var_delt=NULL,Cx=NULL,
                          phi=NULL,lam=NULL,tol=1e-05){
  if(!is.null(seed)){set.seed(seed)}
  J <- ncol(phi)
  N <- length(unique(dat[,"subj"]))
  subjs <- sort(unique(dat$subj))
  y <- dat[,'y']
  a <- phi%*%diag(lam)
  if(impute_type=="Multiple"){
    out <- array(NA,c(N,J,nimps))
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
    mu_y <- if(y[rows[1]]==0){mu0}else{mu1}

    Bi <- Cx[ind,ind] + diag(var_delt,nrow = length(x_obs))
    di <- x_obs - mu_y[ind]
    # ai <- phi[ind,1:J]%*%diag(lam[1:J],nrow = J)
    ai <- a[ind,]

    # mu_star <- t(ai)%*%solve(Bi)%*%di
    if(length(x_obs)==1){
      mu_star <- ai%*%solve(Bi)%*%di
    }else{
      mu_star <- t(ai)%*%solve(Bi)%*%di
    }
    if(impute_type=="Mean"){
      out[i,] <- mu_star
    }

    if(impute_type=="Multiple"){
      # sig_star <- if(J==1){lam[1] - t(ai)%*%solve(Bi)%*%ai}else{diag(lam[1:J]) - t(ai)%*%solve(Bi)%*%ai}
      if(length(x_obs)==1){
        sig_star <- diag(lam[1:J],nrow = J) - ai%*%solve(Bi)%*%t(ai)
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
      Z <- matrix(rnorm(J*nimps),J,nimps)
      out[i,,] <- if(J==1){c(mu_star)+U*sqrt(Lam)*Z}else{c(mu_star) + U%*%sqrt(Lam)%*%Z}
      # out[i,,] <- t(mvrnorm(nimps,mu_star,sig_star,tol = tol))
    }
  }
  return(out)
}
