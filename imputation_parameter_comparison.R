##########################################################
#------- Comparison of Imputation Estiamtes -------------#
##########################################################

set.seed(123)

## Load packages
library(MASS)
library(fcr)
library(fdapace)
library(dplyr)
library(CompQuadForm)
library(sparsefreg)
library(plot3D)

## Data generation/imputation parameters
M <- 100 # grid size
N <- 200
m <- 10
J <- 5
nimps <- 10
w <- 10
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

## Slope function
beta <- w*sin(2*pi*grid)
# beta <- w*(sin(2*pi*grid)+1)
# beta = -dnorm(grid, mean=.2, sd=.03)+3*dnorm(grid, mean=.5, sd=.04)+dnorm(grid, mean=.75, sd=.05)

## Intercept
alpha <- 0

Cxy <- Cx%*%beta/M
muy <- c(t(mux)%*%beta/M)
var_y <- c(t(beta)%*%Cx%*%beta/(M^2)) + var_eps

## Simulate data
X_s <- mvrnorm(N,mux,Cx)
X_comp <- X_s + rnorm(N*M,sd = sqrt(var_delt))
Xi <- (X_s-mux)%*%phi/M
eps <- rnorm(N,0,sd = sqrt(var_eps))
y <- c(alpha + X_s%*%beta/M + eps)

# Sample mi from a discrete uniform distribution s.t. E(x) = (a + b)/2 = m
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

##############################################
## Use fcr to estimate impuation parmaeters ##
##############################################
fcr_est <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = NULL,
                  use_fcr = T,
                  family = "Gaussian",
                  cond.y = T,
                  impute_type = "Multiple")

###############################################
## Use FPCA to estimate impuation parmaeters ##
###############################################
fpca_est <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = NULL,
                   use_fcr = F,
                   family = "Gaussian",
                   cond.y = T,
                   impute_type = "Multiple")

############################################
## Compare imputation parameter estimates ##
############################################

## Cx
par(mfrow = c(1,3))
persp3D(grid,grid,Cx,main = "True")
persp3D(grid,grid,fcr_est$params$Cx,main = "fcr")
persp3D(grid,grid,fpca_est$params$Cx,main = "fpca")

## lam
cbind(lam,fcr_est$params$lam,fpca_est$params$lam)

## phi
ylim = c(-3,3)
plot(grid,phi[,1],type = 'l',ylim = ylim,main = "True")
lines(grid,phi[,2],type = 'l',col = 'blue')
lines(grid,phi[,3],type = 'l',col = 'red')
lines(grid,phi[,4],type = 'l',col = 'green')
lines(grid,phi[,5],type = 'l',col = 'orange')

plot(grid,fcr_est$params$phi[,1],type = 'l',ylim = ylim,main = "True")
lines(grid,fcr_est$params$phi[,2],type = 'l',col = 'blue')
lines(grid,fcr_est$params$phi[,3],type = 'l',col = 'red')
lines(grid,fcr_est$params$phi[,4],type = 'l',col = 'green')
lines(grid,fcr_est$params$phi[,5],type = 'l',col = 'orange')

plot(grid,fpca_est$params$phi[,1],type = 'l',ylim = ylim,main = "True")
lines(grid,fpca_est$params$phi[,2],type = 'l',col = 'blue')
lines(grid,fpca_est$phi[,3],type = 'l',col = 'red')
lines(grid,fpca_est$params$phi[,4],type = 'l',col = 'green')
lines(grid,fpca_est$params$phi[,5],type = 'l',col = 'orange')

## mux
ylim = c(-1,1)
plot(grid,mux,type = 'l',ylim = ylim,main = "True")
plot(grid,fcr_est$params$mux,type = 'l',ylim = ylim,main = "fcr")
plot(grid,fpca_est$params$mux,type = 'l',ylim = ylim,main = "fpca")

## var_delt
var_delt
fcr_est$params$var_delt
fpca_est$var_delt

## var_eps
var_eps
fcr_est$params$var_eps
fpca_est$params$var_eps

## Cxy
plot(grid,Cxy,type = 'l',main = 'True')
plot(grid,fcr_est$params$Cxy,type = 'l',main = "fcr")
plot(grid,fpca_est$params$Cxy,type = 'l',main = "fpca")
