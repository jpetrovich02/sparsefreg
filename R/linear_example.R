###################################################################
#------- Example Using MISFIT for a Linear SoF Model -------------#
###################################################################


set.seed(123)

library(MASS)
library(fcr)
# library(tidyverse)
library(dplyr)
library(CompQuadForm)
library(sparsefreg)

## Data generation
M <- 100 # grid size
N <- 200
m <- 10
J <- 4
K <- 10
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
lam <- eigen(Cx,symmetric = T)$values/M
phi <- eigen(Cx,symmetric = T)$vectors*sqrt(M)

# beta <- w*sin(2*pi*grid)
beta <- w*(sin(2*pi*grid)+1)
# beta = -dnorm(grid, mean=.2, sd=.03)+3*dnorm(grid, mean=.5, sd=.04)+dnorm(grid, mean=.75, sd=.05)
alpha <- 0

X_s <- mvrnorm(N,mux,Cx)
X_comp <- X_s + rnorm(N*M,sd = sqrt(var_delt))
Xi <- (X_s-mux)%*%phi/M
eps <- rnorm(N,0,sd = sqrt(var_eps))
y <- c(alpha + X_comp%*%beta/M + eps)


X_mat<-matrix(nrow=N,ncol=m)
T_mat<-matrix(nrow=N,ncol=m)
ind_obs<-matrix(nrow=N,ncol=m)

for(i in 1:N){
  ind_obs[i,]<-sort(sample(1:M,m,replace=FALSE))
  X_mat[i,]<-X_comp[i,ind_obs[i,]]
  T_mat[i,]<-grid[ind_obs[i,]]
}

spt<-1
ind_obs[spt,1] = 1; ind_obs[spt,m] = M
X_mat[spt,]<-X_comp[spt,ind_obs[spt,]]
T_mat[spt,]<-grid[ind_obs[spt,]]

## Create data frame for observed data
obsdf <- data.frame("X" = c(t(X_mat)),"argvals" = c(t(T_mat)),
                    "y" = rep(y,each = m),"subj" = rep(1:N,each = m))

user_params <- list(Cx = Cx, mux = mux, var_delt = var_delt,
                    muy = mean(mux*beta),lam = lam, phi = phi, Cxy = Cx%*%beta/M,
                    var_y = t(beta)%*%Cx%*%beta/(M^2) - (mean(mux*beta))^2 + var_eps)
check <- misfit(obsdf,grid = grid,K = K,J = J,family = "Gaussian",user_params = NULL,k = 20)
check$pvnorm
check$alpha.hat



plot(grid,mux,type = 'l')
lines(grid,check$params$mux,lty = 2)

plot(grid,beta,type = 'l')
lines(grid,check$beta.hat,lty = 2)

check$params$var_delt
var_delt

alpha
check$alpha.hat

mean(mux*beta)
check$params$muy

t(beta)%*%Cx%*%beta/(M^2) - (mean(mux*beta))^2 + var_eps
check$params$var_y

plot(grid,Cx%*%beta/M,type = 'l')
lines(grid,check$params$Cxy,lty = 2)

library(plot3D)
persp3D(grid,grid,Cx)
persp3D(grid,grid,check$params$Cx)

par(mfrow = c(1,2))
matplot(t(check$Xest),type = 'l')
matplot(t(X_s),type = 'l')
par(mfrow = c(1,1))
