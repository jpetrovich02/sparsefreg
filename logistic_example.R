#####################################################################
#------- Example Using MISFIT for a Logistic SoF Model -------------#
#####################################################################

set.seed(123)

## Load packages
library(MASS)
library(fcr)
library(dplyr)
library(CompQuadForm)
library(sparsefreg)

## Data generation/imputation parameters
M <- 100 # grid size
N <- 200
m <-10
J <- 2
nimps <- 10
w <- 1
var_delt <- 0.5
grid <- seq(from=0,to=1,length.out = M)
nfpc <- 2 # number of fpcs used for mu1
p <- 0.5
mu0 <- rep(0,M)
Cx_f<-function(t,s,sig2=1,rho=0.5){ # Matern covariance function with nu = 5/2
  d <- abs(outer(t,s,"-"))
  tmp2 <- sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
Cx <- Cx_f(grid,grid)
lam <- eigen(Cx,symmetric = T)$values/M
# lam <- eigen(C0,symmetric = T)$values
phi <- eigen(Cx,symmetric = T)$vectors*sqrt(M)
# phi <- eigen(C0,symmetric = T)$vectors
if(nfpc==1){
  mu1 <- phi[,1]*w
}else{
  mu1 <- rowSums(phi[,1:nfpc])*w
}
mux <- (mu1 - mu0)*p

## Slope Function
if(nfpc==1){
  beta <- phi[,1]*(1/lam[1])*w
}else{
  beta <- phi[,1:nfpc]%*%(1/lam[1:nfpc])*w
  # beta <- phi%*%(c(t(mu1)%*%phi/M)/lam)
}

## Intercept
alpha <- log(p) - log(1 - p) - sum(((t(phi)%*%(mu1 - mu0)/M)^2)/(lam^2))/2

## Simulate Data
y <- sort(rbinom(N,1,p))
N_0 <- sum(y==0)
N_1 <- sum(y==1)
Xs_0 <- mvrnorm(N_0,mu0,Cx)
Xs_1 <- mvrnorm(N_1,mu1,Cx)
Xn_0 <- Xs_0 + rnorm(N_0*M,sd = sqrt(var_delt))
Xn_1 <- Xs_1 + rnorm(N_1*M,sd = sqrt(var_delt))
X_s <- rbind(Xs_0,Xs_1)
X_comp <- rbind(Xn_0,Xn_1)

Xi0 <- (Xs_0 - mu0)%*%phi/M
Xi1 <- (Xs_1 - mu1)%*%phi/M
Xi <- rbind(Xi0,Xi1)

X_mat<-matrix(nrow=N,ncol=m)
T_mat<-matrix(nrow=N,ncol=m)
ind_obs<-matrix(nrow=N,ncol=m)

for(i in 1:N){
  ind_obs[i,]<-sort(sample(1:M,m,replace=FALSE))
  X_mat[i,]<-X_comp[i,ind_obs[i,]]
  T_mat[i,]<-grid[ind_obs[i,]]
}

spt<-c(1,N_0+1)
ind_obs[spt,1] = 1; ind_obs[spt,m] = M
X_mat[spt,]<-X_comp[spt,ind_obs[spt[1],]]
T_mat[spt,]<-rbind(grid[ind_obs[spt[1],]],grid[ind_obs[spt[1],]])

## Create data frame for observed data
obsdf <- data.frame("X" = c(t(X_mat)),"argvals" = c(t(T_mat)),
                    "y" = rep(y,each = m),"subj" = rep(1:N,each = m))

# # To use the true imputation parameters, uncomment below and input these in the user_params argument
# unconditional_params <- list(Cx = Cx, mux = mux,var_delt = var_delt,
#                              lam = lam, phi = phi)
#
# conditional_params <- list(Cx = Cx, mu0 = mu0, mu1 = mu1, mux = mux,
#                            var_delt = var_delt, lam = lam, phi = phi)

####################################################
## Mean Imputation, Unconditional on the response ##
####################################################
meu <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = NULL,
              family = "binomial",
              cond.y = F,
              impute_type = "Mean")

########################################################
## Multiple Imputation, Unconditional on the response ##
########################################################
muu <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = meu$params,
              family = "binomial",
              cond.y = F,
              impute_type = "Multiple")

##################################################
## Mean Imputation, Conditional on the response ##
##################################################
mec <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = NULL,
              family = "binomial",
              cond.y = T,
              impute_type = "Mean")

######################################################
## Multiple Imputation, Conditional on the response ##
######################################################
muc <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = mec$params,
              family = "binomial",
              cond.y = T,
              impute_type = "Multiple")

#####################
## Alpha Estimates ##
#####################
alpha
meu$alpha.hat
muu$alpha.hat
mec$alpha.hat
muc$alpha.hat

####################
## Beta Estimates ##
####################
plot(grid,beta,type = 'l')
lines(grid,meu$beta.hat,col = 'red')
lines(grid,muu$beta.hat,col = 'orange')
lines(grid,mec$beta.hat,col = 'green')
lines(grid,muc$beta.hat,col = 'blue')

################################
## Pointwise Confidence Bands ##
################################
par(mfrow = c(2,2))
ylim = c(-8,8)
plot(grid,meu$beta.hat,type = 'l',ylim = ylim,main = "Mean Unconditional")
lines(grid,meu$beta.hat + 1.96*sqrt(diag(meu$Cbeta)),lty = 2)
lines(grid,meu$beta.hat - 1.96*sqrt(diag(meu$Cbeta)),lty = 2)
plot(grid,muu$beta.hat,type = 'l',ylim = ylim,main = "Multiple Unconditional")
lines(grid,muu$beta.hat + 1.96*sqrt(diag(muu$Cbeta$var.t)),lty = 2)
lines(grid,muu$beta.hat - 1.96*sqrt(diag(muu$Cbeta$var.t)),lty = 2)
plot(grid,mec$beta.hat,type = 'l',ylim = ylim,main = "Mean Conditional")
lines(grid,mec$beta.hat + 1.96*sqrt(diag(mec$Cbeta)),lty = 2)
lines(grid,mec$beta.hat - 1.96*sqrt(diag(mec$Cbeta)),lty = 2)
plot(grid,muc$beta.hat,type = 'l',ylim = ylim,main = "Multiple Conditional")
lines(grid,muc$beta.hat + 1.96*sqrt(diag(muc$Cbeta$var.t)),lty = 2)
lines(grid,muc$beta.hat - 1.96*sqrt(diag(muc$Cbeta$var.t)),lty = 2)

############################
## Estimated X's (curves) ##
############################
par(mfrow = c(1,1))
ylim = c(-5,4)
matplot(grid,t(X_s[which(y==0),]),
        type = 'l',col = 'black',ylim = ylim,main = "True")
matplot(grid,t(X_s[which(y==1),]),type = 'l',col = 'red',add = T)

matplot(grid,t(meu$Xhat[which(y==0),]),
        type = 'l',col = 'black',ylim = ylim,main = "Mean Unconditional")
matplot(grid,t(meu$Xhat[which(y==1),]),type = 'l',col = 'red',add = T)

matplot(grid,t(muu$Xhat[which(y==0),]),
        type = 'l',col = 'black',ylim = ylim,main = "Multiple Unconditional")
matplot(grid,t(muu$Xhat[which(y==1),]),type = 'l',col = 'red',add = T)

matplot(grid,t(mec$Xhat[which(y==0),]),
        type = 'l',col = 'black',ylim = ylim,main = "Mean Conditional")
matplot(grid,t(mec$Xhat[which(y==1),]),type = 'l',col = 'red',add = T)

matplot(grid,t(muc$Xhat[which(y==0),]),
        type = 'l',col = 'black',ylim = ylim,main = "Multiple Conditional")
matplot(grid,t(muc$Xhat[which(y==1),]),type = 'l',col = 'red',add = T)



## Individual Trajectories
id <- sample(1:N,size = 1)
plot(grid,X_s[id,],lwd = 3,type = 'l',ylim = c(-4,4))
lines(grid,X_comp[id,],lty = 2)
points(T_mat[id,],X_mat[id,],lwd = 2)
lines(grid,muc$Xhat[id,],col = 'blue')
lines(grid,mec$Xhat[id,],col = "green3")
lines(grid,meu$Xhat[id,],col = 'red')
lines(grid,muu$Xhat[id,],col = 'orange3')
