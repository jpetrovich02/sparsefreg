# Load Packages
library(MASS)
# library(fcr)
library(fdapace)

#---------------------------------------------------------------------------------------------
set.seed(123)

## Data generation
M <- 100 # grid size
w <- 5 # weight on coefficient function
sims <- 1
var_eps <- 1
var_delt <- 0.5
grid <- seq(from=0,to=1,length.out = M)
mux <- rep(0,M)
Cx_f<-function(t,s,sig2=1,rho=0.5){ # Matern covariance function with nu = 5/2
  d <- abs(outer(t,s,"-"))
  tmp2 <- sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
Cx <- Cx_f(grid,grid)
# lam <- eigen(Cx,symmetric = T)$values/M
phi <- eigen(Cx,symmetric = T)$vectors*sqrt(M)

# beta <- w*sin(2*pi*grid)
beta <- w*(sin(2*pi*grid)+1)
# beta = -dnorm(grid, mean=.2, sd=.03)+3*dnorm(grid, mean=.5, sd=.04)+dnorm(grid, mean=.75, sd=.05)

Cxy <- Cx%*%beta/M
muy <- c(t(mux)%*%beta/M)
var_y <- c(t(beta)%*%Cx%*%beta/(M^2)) + var_eps

## Simulation Parameters
N <- 200; J <- 2; k <- 10
m_all <- c(2,5,10,20,50)
lm <- length(m_all)

## True values of all imputation parameters
# user_params <- list(Cx = Cx, mux = mux, var_delt = var_delt,
#                     muy = muy,lam = lam, phi = phi, Cxy = Cxy,
#                     var_y = var_y)

#-------------------------------
# ## Comparison metrics
comp_names <- c("MeC","MuC","MeU","MuU","Pace","True")
ln <- length(comp_names)
name_list <- list(NULL,comp_names)
mse_b <- matrix(NA,sims,ln,dimnames = name_list)
mse_x <- matrix(NA,sims,ln,dimnames = name_list)
betas <- array(NA,dim = c(M,sims,ln),dimnames = list(NULL,NULL,comp_names))
mse_x_med <- matrix(NA,lm,ln,dimnames = name_list)
mse_x_mean <- matrix(NA,lm,ln,dimnames = name_list)
mse_b_med <- matrix(NA,lm,ln,dimnames = name_list)
mse_b_mean <- matrix(NA,lm,ln,dimnames = name_list)
mse_b_tmean <- matrix(NA,lm,ln,dimnames = name_list)
bias <- matrix(NA,lm,ln,dimnames = name_list)
variance <- matrix(NA,lm,ln,dimnames = name_list)

#-----------------------------------------
for(jj in 1:lm){

  for(ii in 1:sims){
    ## Simulate Data
    X_s <- mvrnorm(N,mux,Cx)
    X_comp <- X_s + rnorm(N*M,sd = sqrt(var_delt))
    Xi <- (X_s-mux)%*%phi/M
    eps <- rnorm(N,0,sd = sqrt(var_eps))
    y <- c(X_s%*%beta/M + eps)

    #---------------------------------------------------------------------------------------------
    # Sample mi from a discrete uniform distribution s.t. E(x) = (a + b)/2 = m
    m <- m_all[jj]
    a <- 1
    b <- m*2 - a
    mi <- sample(a:b,size = N,replace = T)

    # Observed values
    Xl <- vector(mode = "list",length = N)
    Tl <- vector(mode = "list",length = N)
    ind_obs <- vector(mode = "list",length = N)

    for(i in 1:N){
      ind_obs[[i]] <- sort(sample(1:M,mi[i],replace=FALSE))
      Xl[[i]] <- X_comp[i,ind_obs[[i]]]
      Tl[[i]] <- grid[ind_obs[[i]]]
    }

    # Guarantee that the boundaries of the grid have been sampled
    spt <- which(mi > 1)[1]
    ind_obs[[spt]][1] <- 1
    ind_obs[[spt]][mi[spt]] = M
    Xl[[spt]] <- X_comp[spt,ind_obs[[spt]]]
    Tl[[spt]] <- grid[ind_obs[[spt]]]

    ## Create data frame for observed data
    obsdf <- data.frame("X" = unlist(Xl),
                        "argvals" = unlist(Tl),
                        "y" = rep(y,times = lengths(Xl)),
                        "subj" = rep(1:N,times = lengths(Xl)))

    #---------------------------------------------------------------------------------------------
    ## PACE (using its own estimation of distribution parameters)
    pace <- FPCA(Ly = Xl,Lt = Tl,
                 optns = list(dataType = "Sparse",methodSelectK=J,nRegGrid = M,error=TRUE))
    X.pace <- fitted.values(pace)
    mse_x[ii,"Pace"] <- mean(rowMeans((X_s-X.pace)^2))

    # Beta estimate
    fit_pace <- lm(y~pace$xiEst)
    b.hat.pace <- coef(fit_pace)[-1]
    beta.hat.pace <- if(J==1){pace$phi[,1]*b.hat.pace}else{pace$phi[,1:J]%*%b.hat.pace}
    betas[,ii,"Pace"] <- beta.hat.pace
    mse_b[ii,"Pace"] <- sum((beta.hat.pace-beta)^2)/M
    #---------------------------------------------------------------------------------------------
    ## Mean Imputation, Unconditional on outcome (PACE)
    meu <- misfit(obsdf,grid = grid,nimps = 10,J = J,user_params = NULL,seed = ii,
                  family = "Gaussian",
                  impute_type = "Mean",
                  cond.y = FALSE)

    mse_x[ii,"MeU"] <- mean(rowMeans((X_s - meu$Xhat)^2))
    betas[,ii,"MeU"] <- meu$beta.hat
    mse_b[ii,"MeU"] <- sum((meu$beta.hat - beta)^2)/M
    #---------------------------------------------------------------------------------------------
    ## Multitple Imputation, Conditional on outcome
    muc <- misfit(obsdf,grid = grid,nimps = 10,J = J,user_params = NULL,seed = ii,
                  family = "Gaussian",
                  impute_type = "Multiple",
                  cond.y = TRUE)

    mse_x[ii,"MuC"] <- mean(rowMeans((X_s - muc$Xhat)^2))
    betas[,ii,"MuC"] <- muc$beta.hat
    mse_b[ii,"MuC"] <- sum((muc$beta.hat - beta)^2)/M
    #----------------------------------------------------------------------------------------
    ## Mean Imputation, Conditional on outcome
    mec <- misfit(obsdf,grid = grid,nimps = 10,J = J,user_params = muc$params,seed = ii,
                  family = "Gaussian",
                  impute_type = "Mean",
                  cond.y = TRUE)

    mse_x[ii,"MeC"] <- mean(rowMeans((X_s - mec$Xhat)^2))
    betas[,ii,"MeC"] <- mec$beta.hat
    mse_b[ii,"MeC"] <- sum((mec$beta.hat - beta)^2)/M
    #----------------------------------------------------------------------------------------
    ## Multiple Imputation, Unconditional on outcome
    muu <- misfit(obsdf,grid = grid,nimps = 10,J = J,user_params = meu$params,seed = ii,
                  family = "Gaussian",
                  impute_type = "Multiple",
                  cond.y = FALSE)
    mse_x[ii,"MuU"] <- mean(rowMeans((X_s - muu$Xhat)^2))
    betas[,ii,"MuU"] <- muu$beta.hat
    mse_b[ii,"MuU"] <- sum((muu$beta.hat - beta)^2)/M

    #-----------------------------------------------------------------------------------------
    ## Beta True (approximated by J FPCs)
    fit_true <- lm(y~Xi[,1:J])
    b.hat.true <- coef(fit_true)[-1]
    beta.hat.true <- if(J==1){phi[,1]*b.hat.true}else{phi[,1:J]%*%b.hat.true}
    betas[,ii,"True"] <- beta.hat.true
    mse_b[ii,"True"] <- sum((beta.hat.true - beta)^2)/M
  }
  mse_x_med[jj,] <- apply(mse_x,2,median)
  mse_x_mean[jj,] <- apply(mse_x,2,mean)
  mse_b_med[jj,] <- apply(mse_b,2,median)
  mse_b_mean[jj,] <- apply(mse_b,2,mean)
  mse_b_tmean[jj,] <- apply(mse_b,2,mean,trim = 0.05)
  bias[jj,] <- apply(betas,3,function(x){
    sum((rowMeans(x) - beta)^2)/M
  })
  variance[jj,] <- apply(betas,3,function(x){
    mean(colSums((x - rowMeans(x))^2)/M)
  })

  cat(jj)
}

random.m <- list(m_all = m_all,N = N,J = J,k = k,w = w,bias = bias,variance = variance,
                 mse_b_mean = mse_b_mean,mse_b_med = mse_b_med,mse_x_mean = mse_x_mean,
                 mse_x_med = mse_x_med,beta = beta,grid = grid)
