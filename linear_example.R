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
N <- 400
m <- 20
J <- 5
K <- 10
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

beta <- w*sin(2*pi*grid)
# beta <- w*(sin(2*pi*grid)+1)
# beta = -dnorm(grid, mean=.2, sd=.03)+3*dnorm(grid, mean=.5, sd=.04)+dnorm(grid, mean=.75, sd=.05)
alpha <- 0

X_s <- mvrnorm(N,mux,Cx)
X_comp <- X_s + rnorm(N*M,sd = sqrt(var_delt))
Xi <- (X_s-mux)%*%phi/M
eps <- rnorm(N,0,sd = sqrt(var_eps))
y <- c(alpha + X_s%*%beta/M + eps)

Cxy <- Cx%*%beta/M
muy <- c(t(mux)%*%beta/M)
var_y <- c(t(beta)%*%Cx%*%beta/(M^2)) + var_eps


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

user_params <- list(Cx = Cx, mux = mux, var_delt = var_delt,
                    muy = muy,lam = lam, phi = phi, Cxy = Cxy,
                    var_y = var_y)

# nPhi <- min(c(floor((nrow(obsdf) - 2*as.numeric(ks))/N),J))
check <- misfit(obsdf,grid = grid,nimps = 10,J = J,family = "Gaussian",user_params = NULL,k = -1)
check$pvnorm
check$alpha.hat
sum((check$beta.hat-beta)^2)/M
mean(rowMeans((X_s-check$Xhat)^2))


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
matplot(t(check$Xhat),type = 'l')
matplot(t(X_s),type = 'l')
par(mfrow = c(1,1))


{ids <- c(1,10,100,300)
  # ids <- sample(1:N,size = 4,replace = F)
par(mfrow = c(2,2))
for(i in 1:4){
  sid <- ids[i]
  ylim <- range(c(X_s[sid,],Xl[[sid]],check$Xhat[sid,]))
  plot(grid,X_s[sid,],type = 'l',ylim = ylim,main = paste("Subject Number ",sid))
  points(Tl[[sid]],Xl[[sid]])
  lines(grid,check$Xhat[sid,],lty = 2)
}
}
