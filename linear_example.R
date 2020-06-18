###################################################################
#------- Example Using MISFIT for a Linear SoF Model -------------#
###################################################################

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

# user_params <- list(Cx = Cx, mux = mux, var_delt = var_delt,
#                     muy = muy,lam = lam, phi = phi, Cxy = Cxy,
#                     var_y = var_y)

####################################################
## Mean Imputation, Unconditional on the response ##
####################################################
meu <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = NULL,k = -1,
              family = "Gaussian",
              cond.y = F,
              impute_type = "Mean")

########################################################
## Multiple Imputation, Unconditional on the response ##
########################################################
muu <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = meu$params,
              family = "Gaussian",
              cond.y = F,
              impute_type = "Multiple")

##################################################
## Mean Imputation, Conditional on the response ##
##################################################
mec <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = NULL,k = -1,
              family = "Gaussian",
              cond.y = T,
              impute_type = "Mean")

######################################################
## Multiple Imputation, Conditional on the response ##
######################################################
muc <- misfit(obsdf,grid,J = J,nimps = nimps,user_params = mec$params,
              family = "Gaussian",
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
ylim = c(-4,3)
matplot(grid,t(X_s),type = 'l',ylim = ylim,
        main = "True")

matplot(grid,t(meu$Xhat),type = 'l',ylim = ylim,
        main = "Mean Unconditional")

matplot(grid,t(muu$Xhat),type = 'l',ylim = ylim,
        main = "Multiple Unconditional")

matplot(grid,t(mec$Xhat),type = 'l',ylim = ylim,
        main = "Mean Conditional")

matplot(grid,t(muc$Xhat),type = 'l',ylim = ylim,
        main = "Multiple Conditional")
