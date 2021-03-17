########## spendingDensities.R
# Creator: Alex Hoagland
# Created: 2/2/2021
# Last modified:
#
# PURPOSE:
#   For a given family-year-simulation number, calculates the implied conditional pdf of household spending
#
# NOTES: 
#    - Adapted from Marone and Sabety, 2021
#    - Nf refers to the # of individuals in each family-year
#    - Includes reference to omega just in case that changes over individuals later. 
#
# Input:
#   - myfam_typepoints_fn : indvidual-level types (p, muL, psi), Nf x 3
#   - myfam_sigL : indvidual-level sigma_lambda, Nf x 1
#   - myfam_kappa : indvidual-level kappa, Nf x 1
#   DATA INPUTS
#     - base : simple data set for keeping track of enrolids and newfamids
#     - allind_spend: individuals' total spending
#     - d,c,maxoop: family's plan characteristics (each 1x1)
# Output:
#   - single parameter for conditional probability of observing m given theta, data, and family type points
#
# MAJOR EDITS:  
#    - TBA: 
################################################################################


spendingDensities <- function(myfam_typepoints_fn,myfam_sigL,myfam_kappa,icprime_fn,ispend_fn,imch_fn,iconvol_fn,omega_fn,alpha_fn1,alpha_fn2) {
    # This starts at the individual level and is currently aggregated to household level by assuming spending choices are conditionally independent
  
    # a. Calculate implied lambda, based on allind_cprime (built in loadData.R)
  # now create the implied_lambda, which is -999 if spending is 0 and the solution to equation (8) solved for lambda otherwise
  implam <- ispend_fn*(1+myfam_typepoints_fn[,1]*(alpha_fn1-1))+
    myfam_typepoints_fn[,1]*alpha_fn2*imch_fn-
    omega_fn*(1+myfam_typepoints_fn[,1]*(alpha_fn1-1)-icprime_fn)
  implam[ispend_fn==0] <- -999
  
    # b. if medical spending is 0, then f_m(m) = NormCDF((log(-kappa)-mu_lambda)/sigma_lambda)
    # note: none of these parameters change over time
    nospend_vec <- pnorm((log(-1*myfam_kappa)-myfam_typepoints_fn[,2])/myfam_sigL)+iconvol_fn
    
    # c. if medical spending is >0, then f_m(m) = NormPDF((log(implied_lambda-kappa)-mu_lambda)/sigma_lambda)
    allind_implied_spenddist <- dnorm((log(implam-myfam_kappa)-myfam_typepoints_fn[,2])/myfam_sigL) +iconvol_fn
    allind_implied_spenddist[is.na(allind_implied_spenddist)] <- nospend_vec[is.na(allind_implied_spenddist)] 
    
    # d. Aggregate to family level (by assuming independence)
    # don't count household members not present in a given year. 
    allind_implied_spenddist[is.na(ispend_fn)] <- NA
    return(prod(allind_implied_spenddist,na.rm=T))
}