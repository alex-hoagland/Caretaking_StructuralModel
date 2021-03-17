########## mainEstimation.R ####
# Creator: Alex Hoagland
# Created: 1/26/2021
# Last modified: 2/3/2021
#
# PURPOSE:
#   Main file for estimation
#
# NOTES: 
#    - Adapted from Matlab code in Einav et al., "Selection on Moral Hazard"
#
# MAJOR EDITS:  
#    - TBA: 
################################################################################

##### LAST RUN INFO:
# Families: 100
# Convergence: 1e-2
# Method: Nelder-Mead 
# Run Time: 19 hours
# Memory: 129.830G
# Number of iterations: 923
#############################


##### 0. Packages and functions ####
# If data has been pre-organized, load now 
# print("Loading pre-formatted data for the cluster and updating functions.")
load("SpeedyData_2021-03-13.Rdata") # set initial parameters to last iteration's output
load("FinalData_20210311_SCC.RData") # main data set
newstart <- as.vector(thetastar$estimate)
newstart[40] <- 0.9 # Reset pi0 to be .9
newstart[44] <- 1 # This one hardly changes outcomes, keep small? Need to look into that

library(data.table)
library(tidyr)
library(dplyr)
library(gmp) # for handling large numbers -- note that this masks matrix multiplication, will need to work with that? 
library(Rmpfr) # used in connection with gmp
library(nloptr)
library(doMC)
library(doParallel)
library(parallel)
library(foreach)
library(matrixStats) # for colMax()
library(matrixcalc) # for is.positive.definite
library(maxLik) # new Maximum Likelihood optimization
library(Rcpp)

### functions from other files
rm(list = lsf.str()) # remove all functions
Sys.setenv("PKG_LIBS" = "-lmpfr -lgmp") # Need to use MPFR for large exponentiation
sourceCpp("findChoices_cpp.cpp")
source("spendingDensities.R")
source("calculateLikelihood.R")

### functions for Gauss-Hermite ####
### 0. Gauss-Hermite nodes and weights for integration (single and multivariate)
# UNIVARIATE ref: http://www.aims.ac.za/resources/courses/na/gauss_quad.pdf
# n: order of the hermite polynomila
# x : nodes
# w : weights
gqzero <- function(n) {
   d <- sqrt((1:n)/2)
   
   # Assign d to offset diagnoals of matrix
   T <- matrix(0, n+1, n+1)
   T[row(T)-1==col(T)] <- d
   T[row(T)+1==col(T)] <- d
   V <- eigen(T)$vectors
   x <- eigen(T)$values
   w <- V[1,]^2
   w <- w*sqrt(pi)
   
   ind <- order(x)
   w <- t(w[ind])
   
   returnlist <- list("n"=n,"x"=x,"w"=w)
   return(returnlist)
} 

# MULTIVARIATE ref: https://www.r-bloggers.com/2015/09/notes-on-multivariate-gaussian-quadrature-with-r-code/ 
# n: number of points in each dimension (before pruning)
# mu: mean vector
# sigma: covariance matrix
# prune - NULL = no pruning; [0-1] is the fraction to prune
hermite <- function (points, z) {
   p1 <- 1/pi^0.4
   p2 <- 0
   for (j in 1:points) {
      p3 <- p2
      p2 <- p1
      p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
   }
   pp <- sqrt(2 * points) * p2
   c(p1, pp)
} # startup function
gauss.hermite <- function (points, iterlim = 50) { # startup function
   x <- w <- rep(0, points)
   m <- (points + 1)/2
   for (i in 1:m) {
      z <- if (i == 1) 
         sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
      else if (i == 2) 
         z - sqrt(points)/z
      else if (i == 3 || i == 4) 
         1.9 * z - 0.9 * x[i - 2]
      else 2 * z - x[i - 2]
      for (j in 1:iterlim) {
         z1 <- z
         p <- hermite(points, z)
         z <- z1 - p[1]/p[2]
         if (abs(z - z1) <= 1e-15) 
            break
      }
      if (j == iterlim) 
         warning("iteration limit exceeded")
      x[points + 1 - i] <- -(x[i] <- z)
      w[i] <- w[points + 1 - i] <- 2/p[2]^2
   }
   r <- cbind(x * sqrt(2), w/sum(w))
   colnames(r) <- c("Points", "Weights")
   r
} # startup function
mgauss.hermite <- function(n, mu, prod, prune=NULL) {
   # note: to improve speed, I've pulled some of this out into calculateLikelihood.R
   if(!all(dim(sigma) == length(mu)))
      stop("mu and sigma have nonconformable dimensions")
   
   # dm  <- length(mu)
   # gh  <- gauss.hermite(n)
   # #idx grows exponentially in n and dm
   # idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
   # pts <- matrix(gh[idx,1],nrow(idx),dm)
   # wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
   
   ## prune
   if(!is.null(prune)) {
      qwt <- quantile(wts, probs=prune)
      pts <- pts[wts > qwt,]
      wts <- wts[wts > qwt]
   }
   
   ## rotate, scale, translate points
   # eig <- eigen(sigma) 
   # rot <- eig$vectors %*% diag(sqrt(eig$values))
   # pts <- t(rot %*% t(pts) + mu)
   out <- t(prod + mu)
   #return(list(points=pts, weights=wts))
   return(out)
}
########


### set up parallelization on the cluster
nCores <- as.numeric(Sys.getenv("NSLOTS"))
if (is.na(nCores))nCores = 1
paste("Running code using",nCores,"cores") 
#registerDoMC(nCores)
set.seed(03262020) # This is repeated for each function and done on the cluster later
################################################################################


##### 1. Load data (if I need to start over from raw data files) ####
#file <- "/project/caretaking/Outputs/StructuralModel/SMData_20210223.csv" # Location of data 
#file <- "C:\\Users\\alexh\\Dropbox\\Caretaker_Spending\\2_Data\\7.StructuralModel\\SMData_20210223.csv"

# # Covariates: note variable list is included in loadData.m
# # lx{1} keeps columns of x's for p (transition probability)
# # lx{2} keeps columns of x's for psi (risk aversion)
# # lx{3} keeps columns of x's for lambda 
# # lx{4} keeps columns of x's for kappa
lx1 <- c(88:90,94) # % intercept, age, age2, female, rs, any pe in family, # of pe in family, parent has pe, child has pe, sibling has pe, spouse has pe
lx2 <- c(120,91,94) # % intercept, average age, family size, average risk score, any pe in family, # of pe in family, any adult has pe, any child has pe
lx3 <- c(88:90,92) # % intercept, age, age2, female, risk score, has pe
lx4 <- c(88:90,92) # % intercept, age, age2, female, risk score, has pe
# lx <- list(lx1, lx2, lx3, lx4)
# # NOTE: CHANGE THE INDICES HERE ONCE FULL SAMPLE IS INCLUDED
# 
# print('Loading and formatting data')
# source("loadData.R")
# print("Data formatted. Saving workspace so you don't have to do this every time")
# save.image("FinalData_20210223_SCC.RData")
######################################################################################


##### 2. Set theta0 (TODO: Check 20-50 different values here; start priors with results from previous estimation) ####
# recall that theta1 = (beta_p, beta_psi, beta_lambda, beta_kappa,
                       # sigma_p, sigma_psi, sigma_mu, sigma_kappa, sigma_p_psi, sigma_p_mu, sigma_psi_mu,
                       # gamma1, gamma2) (these two are paramters for sigma_lambda)
# recall that theta2 = (p_it, lambda_i, mu_lambda_i, sigma_lambda_i, psi_ft, kappa_lambda_i, omega)
# beta is a 12x1 vector in the most broad use (NOTE: NEED TO UPDATE THIS ONCE FULL SAMPLE IS INCORPORATED)

# start with initial values of beta (mean shifters)
beta_p <- rep(.1,length(lx1)+1) # intercept, age, age2, female, rs, any pe in family, # of pe in family, parent has pe, child has pe, sibling has pe, spouse has pe
beta_psi <- rep(.1,length(lx2)+1) # intercept, average age, family size, average risk score, any pe in family, # of pe in family, any adult has pe, any child has pe
beta_lambda <- rep(.2,length(lx3)+1) # intercept, age, age2, female, family size, risk score, has pe
beta_kappa <- rep(-.1,length(lx4)+1) # intercept, age, age2, female, family size, risk score, has pe
# 
beta <- list(beta_p, beta_psi, beta_lambda, beta_kappa)

# Sigma (the variance-covariance matrix for p_i0, psi_f_pre, mu_lambda_i) has dimension 3x3
sig2_p <- 1
sig_p_lam <- .2
sig_p_psi <- .4
sig2_lam <- 6
sig_lam_psi <- .2
sig2_psi <- 2

# Gamma paramters for sigma_lambda
gamma1 <- 1
gamma2 <-  1

# Variance of kappa
sig_kappa <- 1 

sigma_E <- 10 # prior for size of idiosyncratic shock 

# State-dependent utility parameters NOTE: for now, assume that alpha is the same across families -- eventually allow heterogeneity
alf1 <- .5 
alf2 <- .5
 
omega <- 202.3502 # average log(omega) = 5.31 in Einav et al.
sc <- 1.729 # average switching cost is $1,729 in Handel (2013)

# Evolution of psi (dp for delta_psi)
dp0 <- 1 # coefficient on psi_{t-1}
dp1 <- .5 # coefficient on any major health event dummy
dp2 <- .5 # any major health event x total diagnostic cost
dp3 <- .5 # any major health event x OOP spending
dp4 <- .5 # any major health event x hospitalization
var_dp <- 1 # variance of psi shock in each period

# Evolution of p (pi)
pi0 <- 1 # coefficient on p_{t-1}
pi1 <- .5 # coefficient on major health event dummy
pi2 <- .5 # coefficient on own acute health shock
pi3 <- .5 # coefficient on any family acute health shock
var_pi <- 1 # variance of p shock in each period

theta0 <- c(beta_p, beta_psi, beta_lambda, beta_kappa,sig2_p,sig_p_lam,sig_p_psi,sig2_lam,sig_lam_psi,sig2_psi,
            gamma1,gamma2,sig_kappa,sigma_E,alf1,alf2,omega,sc,dp0,dp1,dp2,dp3,dp4,var_dp,pi0,pi1,pi2,pi3,var_pi)
rm(lx1,lx2,lx3,lx4,sig2_p,sig_p_lam,sig_p_psi,sig2_lam,sig_lam_psi,sig2_psi,
   gamma1,gamma2,sig_kappa,sigma_E,alf1,alf2,omega,sc,dp0,dp1,dp2,dp3,dp4,var_dp,pi0,pi1,pi2,pi3,var_pi)
#####################################################################################


##### 4. Maximizing likelihood function ####
# Gauss-Hermite nodes and weights for UNIVARIATE normal (used in solving utility-maximization problem)
# I am only using 9 points in GQ across lambda, rather than original 24 (given double integration)
out <- gqzero(9)
order <- length(out$x)
xint <- t(out$x)*sqrt(2)
wint <- out$w/sqrt(pi)
nint <- length(xint)
rm(out)
###############


### Run likelihood once, to time ####
# Date (for filenames)
currentDate <- Sys.Date()
# file1 <- paste("FinalData_InitialLikelihood_",currentDate,".Rdata",sep="")
file2 <- paste("SpeedyData500_",currentDate,".Rdata",sep="")

samplefams <- sample.int(nfam,5)
ptm <- proc.time() # Start the clock
test <- calculateLikelihood(theta0)
proc.time()-ptm
print(paste0("Calculated likelihood for theta0 to test that it works: log-likelihood is ", sum(test)))
# save.image(file=file1)
###############


# Optimize likelihood function ####
# print("Now choosing theta to optimize likelihood function")
# TODO: Decide whether to use full sample or random subset

# First step: simplex estimator until some level of convergence:
# Note: skipping this, appears BFGS is faster even from the beginning
length_allb <- length(beta_p) + length(beta_lambda) + length(beta_psi) + length(beta_kappa)
# thetastar <- maxNM(calculateLikelihood,start=newstart,finalHessian=FALSE, 
#                      constraints=list(ineqA=rbind(c(rep(0, length_allb),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
#                                                   c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)),
#                                       ineqB=rep(0,10)), # require that (i) variances, gamma dist. parameters, switching costs are nonnegative
#                      control=list(reltol=1e-2,printLevel=1,iterlim=1000))

# Second step: use BFGS to estimate until full convergence: 
thetastar <- maxBFGS(calculateLikelihood,start=newstart,finalHessian=FALSE,
                     constraints=list(ineqA=rbind(c(rep(0, length_allb),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                                                c(rep(0, length_allb),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)),
                                    ineqB=rep(0,10)), # require that (i) variances, gamma dist. parameters, switching costs are nonnegative
                   control=list(printLevel=1,iterlim=100))
# #######################################################################################


# ##### 5. Any saving you want to do ####
# ### For Rmpi, if you ever get it to work: Concluding the parallelization
# # stopCluster(cluster)
# # mpi.exit()
# 
save.image(file2) # saves output to use in next loop
########################################################################################
