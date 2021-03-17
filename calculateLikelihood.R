########## calculateLikelihood.R
# Creator: Alex Hoagland
# Created: 1/26/2021
# Last modified: 2/3/2021
#
# PURPOSE:
#   Constructs likelihood function given a guess for theta and my data
#
# NOTES: 
#    - Adapted from Matlab code in Einav et al., "Selection on Moral Hazard" (testGibbs.R)
#    - Also borrows heavily from Marone and Sabety (2021) and Handel (2013)
#    - See https://www.r-bloggers.com/2015/09/notes-on-multivariate-gaussian-quadrature-with-r-code/ for notes on Gaussian quadrature
#
# MAJOR EDITS:  
#    - TBA
################################################################################


##### 0. Needed functions (move to aloca.R?)
findFamilyLL <- function(j,newfamcw,loopfam,myfam_sigL,myfam_typepoints,myfam_kappa,myfam_weights,
                         icprime,ispend,imch,iconvol,alpha_fn1,alpha_fn2,sigE,omega_fn,sc,choice_noswitch,
                         pi0,pi1,pi_vec,p_shock,dp0,dp_vec,psi_shock){
  # Print statement to track progress
  if (loopfam %% 10000 == 0) {
    print(paste0("Working on family ",loopfam))
  }
  
  # start with building objects with all 27 simulation points at once (below, will iterate across simulation points)
  #myfam_sigL <- ind_sigL[which(famcw==loopfam)]
  myfam_typepoints <- array(myfam_typepoints, dim=c(27,3,length(myfam_sigL)))
  # myfam_kappa <- ind_mean_kappa[which(famcw==loopfam)]
  
  # construct family-level weights (product of weights across individuals for each simulation, then reweight)
  # myfam_weights <- apply(matrix(allind_typeweights[,which(famcw==loopfam)]),MARGIN=1,function(x) prod(x))
  # myfam_weights <- allind_typeweights # myfam_weights/sum(myfam_weights)
  
  # Since risk aversion depends only on family-specific variables, just keep 
  #   first individual's draw as family-level psi. 
  myfam_psi <- myfam_typepoints[,3,1]
  
  # organize outputs
  myfam_spendDens <- matrix(NA,27,nyr) # 1 density for each simulation-year
  myfam_choiceProbs <- array(NA,dim=c(27,nyr,6)) # 6 choice probabilities for each simulation-year (if there are 6 plans)
  
  # Loop over all years
  for (t in 1:nyr) {
    # organize data
    
    # Family spending
    # myfam_spend <- ispend[,t]
    
    if (length(ispend[,t])==0) {
      # Family has no spending in year t, skip this one
      next
    }
    
    # Family plan choices
    myfam_d <- as.double(allind_chosen_ded %>% filter(newfam==loopfam) %>% select(paste("y",2005+t,sep="")) %>% na.omit() %>% distinct())
    myfam_c <- as.double(allind_chosen_copay %>% filter(newfam==loopfam) %>% select(paste("y",2005+t,sep="")) %>% na.omit() %>% distinct())
    myfam_m <- as.double(allind_chosen_moop %>% filter(newfam==loopfam) %>% select(paste("y",2005+t,sep="")) %>% na.omit() %>% distinct())
    
    if(is.na(myfam_d)) {
      # print(paste0("Family ",loopfam, " has spending information but no plan choices in year ",t,". Look into this!"))
      next
    }
    myfam_alld <- allp_deduct[,t,loopfam]
    myfam_allc <- allp_copay[,t,loopfam]
    myfam_allm <- allp_maxoop[,t,loopfam]
    myfam_allp <- allp_prem[,t,loopfam]
    
    if (sum(is.na(myfam_alld))==length(myfam_alld)) {
      # household year has no available plan information, skip
      next
    }
    
    # Chronic care costs -- select appropriate points/weights
    itype <- as.double(allfam_types %>% ungroup() %>% filter(newfam==loopfam & year == 2005+t) %>% select(type))
    if(itype == 0) {
      imch_p <- dxc_points
      imch_w <- dxc_weights
    } else {
      # Identify all HCCs
      ihccs <- unlist(allind_types %>% filter(newfam==loopfam & year == 2005+t) %>% select(hcc1:hcc10),use.names=F) 
      ihccs <- ihccs[which(!is.na(ihccs))]
      
      # If only one HCC, draw from that distribution 
      if (length(ihccs) == 1) {
        imch_p <- get(paste0("mch",ihccs[1],"_points",sep=""))
        imch_w <- get(paste0("mch",ihccs[1],"_weights",sep=""))
      } else {
        # If multiple, average quadrature points and weights
        # Combine points and weights into matrices
        allpoints <- matrix(nrow=length(dxc_weights))
        allweights <- allpoints
        for (h in 1:length(ihccs)) {
          allpoints <- cbind(allpoints,get(paste0("mch",ihccs[h],"_points",sep="")))
          allweights <- cbind(allweights,get(paste0("mch",ihccs[h],"_weights",sep="")))
        }
        imch_p <- rowMeans(allpoints,na.rm=T)
        imch_w <- rowMeans(allweights,na.rm=T)
      }
    }
    
    for (si in 1:27) {
      # print(paste0("Family ",j," and si ",si))
      # first, calculate each conditional density (one per family-year-simulation draw)
      # note: there will be warnings about NaNs, but the function spendingDensities does take care of them
      myfam_spendDens[si,t] <- spendingDensities(t(myfam_typepoints[si,,]),myfam_sigL,myfam_kappa,
                                                 icprime[,t],ispend[,t],imch[,t],iconvol[,t],
                                                 omega_fn,alpha_fn1,alpha_fn2)
      # Function args are (myfam_typepoints,myfam_sigL,myfam_kappa,icprime,ispend,imch,convol)
      
      # second, calculate choice probabilities at plan level
      myfam_choiceProbs[si,t,] <- findChoices_cpp(t(myfam_typepoints[si,,])[,1],t(myfam_typepoints[si,,])[,2],myfam_sigL,myfam_kappa,
                                                  alpha_fn1,alpha_fn2,imch_p,imch_w,xint,wint,
                                                  myfam_alld,myfam_allm,myfam_allc,myfam_allp,sigE,myfam_psi[si],
                                                  omega_fn,ispend[,t],sc,as.double(choice_noswitch[t]))
      # findChoices_cpp(t(myfam_typepoints[si,,])[,1],t(myfam_typepoints[si,,])[,2],myfam_sigL,myfam_kappa,
      #                 alpha_fn1,alpha_fn2,imch_p,imch_w,xint,wint,
      #                 myfam_alld,myfam_allm,myfam_allc,myfam_allp,sigE,myfam_psi[si],
      #                 omega_fn,ispend[,t],sc,as.double(choice_noswitch[t]))
      # Function args are (myfam_p, myfam_muL,,myfam_sigL,myfam_kappa,alpha_fn1,alpha_fn2,imch_p,imch_w,xint,wint,
      # allp_deduct, allp_maxoop, allp_copay, allp_prem, sigma_E, psi,omega_fn,ispend,sc,choice_noswitch)
      
      # Update psi for the next period
      myfam_psi[si] <- exp(log(myfam_psi[si])*dp0 + 
                             as.matrix(covariates.psi[which(covariates.psi$newfam==loopfam & covariates.psi$year==2005+t),3:6])%*%dp_vec + 
                             psi_shock[j,t])
      
      # Update p for next period (move to 1 if has chronic; update otherwise)
      myfam_typepoints[si,1,which(myfam_typepoints[si,1,]==1)] <- .99 
        # during updating of beliefs, change 1 to .99 (changed back later; for log odds ratio)
      logodds <- log(myfam_typepoints[si,1,]/(1-myfam_typepoints[si,1,]+1e-25))
      newlogodds <- logodds*pi0+as.double(pi1*covariates.psi[which(covariates.psi$newfam==loopfam & covariates.psi$year==2005+t),3]+
                                            as.matrix(covariates.p[which(covariates.p$newfam==loopfam & covariates.p$year==2005+t),4:5])%*%pi_vec+
                                            p_shock[which(newfamcw==loopfam),t])
      myfam_typepoints[si,1,] <- as.double(mpfr(exp(newlogodds)/(1+exp(newlogodds)),precBits=512))
      
      # Require probabilities to become 1 if the individual was affected by the event
      test <- allind_types %>% ungroup() %>% filter(newfam==loopfam & year == 2005+t) %>% select(type)
      myfam_typepoints[si,1,which(test$type==1)] <- 1
      rm(test)
      
      # Require new probabilities to be between 0 and 1
      myfam_typepoints[si,1,][which(myfam_typepoints[si,1,]>1)] <- 1
      myfam_typepoints[si,1,][which(myfam_typepoints[si,1,]<0)] <- 0
      
    }
  }
  
  ### Aggregate into likelihood function for that family
  if(length(myfam_choiceProbs) == sum(is.na(myfam_choiceProbs))) {
    print(paste0("(findFamilyLL) No choice probabilities returned for family ",loopfam))
    return(1e-323)
  } else if (length(myfam_spendDens) == sum(is.na(myfam_spendDens))) {
    print(paste0("(findFamilyLL) No spending densities returned for family ",loopfam))
    return(1e-323)
  } else {
    # for each year, need to take product of spendDens and actual chosen plan's choiceProbs (take product over years as well)
    choiceprod <- rep(1,27)
    for (t in 1:nyr) {
      choice <- as.double(allfams_choiceid[loopfam,t])
      if (!is.na(choice) & length(myfam_choiceProbs[,t,choice]) != sum(is.na(myfam_choiceProbs[,t,choice]))) { 
        # If the family made a valid choice in that year, update choiceprod vector
        choiceprod <- choiceprod * mpfr(myfam_spendDens[,t]*myfam_choiceProbs[,t,choice],precBits = 512)
      }
    }
    
    hh_likely <- sum(myfam_weights*choiceprod, na.rm=T)
    # sum over all simulation points
    
    if (hh_likely == 0) {
      print(paste0("(findFamilyLL): Family ", loopfam, " has a 0 likelihood function"))
      return(1e-323)
    }
    
    return(hh_likely)
    #out <- list(hh_likely,choiceprod,myfam_spendDens,myfam_choiceProbs) 
    #return(out) # for debugging only
  }
}
################################################################################


##### 1. Main function
calculateLikelihood <- function(theta) {
  set.seed(03262020) 
  
  ### 0. Assign values from guess of theta ####
  beta_fn1 <- theta[1:length(beta_p)]
  beta_fn2 <- theta[length(beta_p)+1:length(beta_psi)]
  beta_fn3 <- theta[length(beta_p)+length(beta_psi)+1:length(beta_lambda)]
  beta_fn4 <- theta[length(beta_p)+length(beta_psi)+length(beta_lambda)+1:length(beta_kappa)]
  allb <- length(beta_p)+length(beta_psi)+length(beta_lambda)+length(beta_kappa)
  
  sig2_p_fn <- theta[allb+1]
  sig_p_lam_fn <- theta[allb+2]
  sig_p_psi_fn <- theta[allb+3]
  sig2_lam_fn <- theta[allb+4]
  sig_lam_psi_fn <- theta[allb+5]
  sig2_psi_fn <- theta[allb+6]
  sigma_GQ <- rbind(c(sig2_p_fn,sig_p_lam_fn,sig_p_psi_fn),
                    c(sig_p_lam_fn,sig2_lam_fn,sig_lam_psi_fn),
                    c(sig_p_psi_fn,sig_lam_psi_fn,sig2_psi_fn))
  
  ### First, check that covariance matrix works
  if (sig_p_lam_fn^2-sig2_p_fn*sig2_lam_fn >= 0 | 
      sig_p_psi_fn^2-sig2_p_fn*sig2_psi_fn >= 0 | 
      sig_lam_psi_fn^2-sig2_psi_fn*sig2_lam_fn >= 0) {
    print("Invalid covariance matrix")
    return(-Inf)
  }
  
  gamma1_fn <- theta[allb+7]
  gamma2_fn <- theta[allb+8]
  sigK <- theta[allb+9]
  sigE <- theta[allb+10]
  alpha_fn1 <- theta[allb+11]
  alpha_fn2 <- theta[allb+12]
  omega_fn <- theta[allb+13]
  sc <- theta[allb+14]
  
  dp0 <- theta[allb+15]
  dp_vec <- theta[allb+16:19]
  var_dp <- theta[allb+20]
  
  pi0 <- theta[allb+21]
  pi1 <- theta[allb+22]
  pi_vec <- theta[allb+23:24]
  var_pi <- theta[allb+25]
  
  ### 1. Based on guess for theta (prior values), sample across distribution of types (p_it,lambda_it,logpsi_it) ####
  
  # note: this builds a i x 1 vector, each row of which is the "mean" vector for that individual. 
  mu_GQ <- cbind(as.matrix(covariates.x1)%*%beta_fn1,
                 as.matrix(covariates.x3)%*%beta_fn3, 
                 as.matrix(covariates.x2)%*%beta_fn2)
  # mean vector for each individual's p (to be logit-ed), muL, and log(psi)
  
  # Limit size of mu according to samplefams
  key <- famcw %>% mutate(test = newfam %in% samplefams)
  newfamcw <- key %>% filter(test == 1) %>% select(newfam)
  mu_GQ <- mu_GQ[which(key$test == 1),]
  
  # Assign GQ weights from sigma (same for all individuals--only need to run once)
  # Note: to assign weights, variance/covariance matrix must be positive definite: 
  if (is.positive.semi.definite(sigma_GQ)) {
    # note: dm = 3, the dimensions of mu; n = 3, the number of support points in each dimension
    gh  <- gauss.hermite(3)
    #idx grows exponentially in n and dm
    idx <- as.matrix(expand.grid(rep(list(1:3),3)))
    pts <- matrix(gh[idx,1],nrow(idx),3)
    eig <- eigen(sigma_GQ) 
    rot <- eig$vectors %*% diag(sqrt(eig$values))
    prodmat <- rot %*% t(pts) # Each individual's support points are just this vector shifted by their mean Xbeta
    
    # allind_typepoints <- array(apply(mu_GQ,MARGIN=1,function(x) t(prodmat+x)), dim=c(27,3,nind))
    allind_typepoints <- array(apply(mu_GQ,MARGIN=1,function(x) t(prodmat+x)), dim=c(27,3,nrow(mu_GQ)))
    allind_typeweights <- apply(matrix(gh[idx,2],nrow(idx),3), 1, prod) 
    rm(gh,idx,pts,eig,rot,prodmat)
    # Note: all individual-level weights are the same, since sigma is the same
  } else {
    print("Sigma_GQ is not positive definite")
    return(-Inf)
  }
  ###############################################################################
  
  
  ### 1a. Require that each individual's type points make sense ####
  # first, convert normal draws for p into a probability
  allind_typepoints[,1,] <- exp(allind_typepoints[,1,])/(1+exp(allind_typepoints[,1,]))
  
  # also require probabilities to be 1 if person has a chronic illness in their first year (type == 1)
  newpe <- pe[which(key$test == 1),]
  allind_typepoints[,1,which(newpe$test==1)] <- 1
  
  # second, exponentiate psi
  allind_typepoints[,3,] <- exp(allind_typepoints[,3,])
  
  # Draw individual sigma_lambdas from the gamma distribution (note: these are constant over time)
  ind_sigL <- 1/sqrt(matrix(rgamma(nrow(mu_GQ),shape=gamma1_fn,scale=gamma2_fn),nrow(mu_GQ),1))
  ind_sigL[which(ind_sigL>100)] <- 100 # Truncate the distribution
  
  # Draw individual kappas from a normal distribution with mean (as.matrix(covariates.x4)%*%beta[[4]]) and standard deviation sigK
  ind_mean_kappa <- rnorm(nrow(mu_GQ),mean = as.matrix(covariates.x4)%*%beta_fn4,sd=sigK)
  ind_mean_kappa[which(ind_mean_kappa>0)] <- 0 # require kappa to be non-positive
  
  # Draw individual-year shocks for p
  p_shock <- matrix(rnorm(nrow(mu_GQ)*nyr,mean=0,sd=sqrt(var_pi)),nrow=nrow(mu_GQ),ncol=nyr)
  
  # Draw family-year shocks for psi
  psi_shock <- matrix(rnorm(length(samplefams)*nyr,mean=0,sd=sqrt(var_dp)),nrow=length(samplefams),ncol=nyr)
  ##################################################################################
  
  
  ### 2. Construct family-level objects and weights ####
  # Loop through all families
  #print("Calculating likelihood function for each household")
  
  cl <- parallel::makeCluster(nCores,type="FORK")#,outfile="Debug_HHLikelihood.txt")
  parallel::clusterEvalQ(cl,set.seed(3262020))

  doParallel::registerDoParallel(cl)
  allfams_ll <- foreach(j = seq_along(samplefams),.combine=rbind) %dopar% {
    findFamilyLL(j,newfamcw,samplefams[j],ind_sigL[which(newfamcw==samplefams[j])],allind_typepoints[,,which(newfamcw==samplefams[j])],
                 ind_mean_kappa[which(newfamcw==samplefams[j])],allind_typeweights,
                 allind_cprime[which(famcw==samplefams[j]),],stripped_indspend[which(famcw==samplefams[j]),],stripped_indmch[which(famcw==samplefams[j]),],convol[which(famcw==samplefams[j]),],
                 alpha_fn1,alpha_fn2,sigE,omega_fn,sc,allfam_choicenoswitch[samplefams[j],],pi0,pi1,pi_vec,p_shock,dp0,dp_vec,psi_shock) 
    # function args are (loopfam,myfam_sigL,myfam_typepoints,myfam_kappa,myfam_weights,
    #   icprime,ispend,imch,iconvol,alpha_fn1,alpha_fn2,sigE,omega_fn)
  }
  parallel::stopCluster(cl)
  
  # Once done looping through families, construct log-likelihood function!
  allfams_ll <- allfams_ll[!is.na(allfams_ll)]
  if (grepl("mpfr", class(allfams_ll ), fixed = TRUE)) {
    allfams_ll <- new("mpfr", unlist(allfams_ll))
  }
  allfams_ll[allfams_ll <= 0] <- 1e-323 # To make sure there are no NaNs introduced by logs of negatives (from the convolutions? TODO CHECK)
  
  ll <- as.double(sum(log(allfams_ll),na.rm=T)) 
  #print(paste0("The log-likelihood for theta0 is ",ll))
  return(as.numeric(log(allfams_ll)))
  #return(ll)
}
##################################################################################