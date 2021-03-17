#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double solveSpend_cpp(NumericVector p, NumericVector lambda, double omega, double alpha_fn1, double alpha_fn2, 
                      double deduct, double copay, double maxoop, double mch) {
  
  // Keep array of possible spending solutions (nf rows, 4 columns)
  // Case 1: no spending
  // Case 2: M < d, so that c'(m) = 1
  // Case 3: M >= d, so that c'(m) = c
  // Case 4: M >= x, so that c'(m) = 0
  int NROW = lambda.size(); 
  NumericMatrix allspend (NROW , 4); 
  int i,j; 
  
  //int myprint=1/(1+p(0)*(alpha_fn1-1))*(lambda(0)+omega*p(0)*(alpha_fn1-1)-alpha_fn2*p(0)*mch);
  //Rcout << myprint << "\n";
  for (j = 0; j < NROW; ++j){
    allspend(j,0) = 0;
    allspend(j,1) = 1/(1+p(j)*(alpha_fn1-1))*(lambda(j)+omega*p(j)*(alpha_fn1-1)-alpha_fn2*p(j)*mch);
    allspend(j,2) = 1/(1+p(j)*(alpha_fn1-1))*(lambda(j)+omega*(1-copay+p(j)*(alpha_fn1-1))-alpha_fn2*p(j)*mch);
    allspend(j,3) = 1/(1+p(j)*(alpha_fn1-1))*(lambda(j)+omega*(1+p(j)*(alpha_fn1-1))-alpha_fn2*p(j)*mch);
  }
  //Rf_PrintValue(allspend);
  
  // Convert each possible spending solution to OOP spending, given plan characteristics and mch's cost
  NumericMatrix M (NROW , 4);
  NumericMatrix oopspend (NROW, 4);
  
  for (i =0; i < 4; ++i  ) {
    for (j = 0; j < NROW; ++j){
      M(j,i) = allspend(j, i) + mch;
      // Individual non-chronic spending + family chronic spending
    }
  }
  
  // Assign cost of chosen spending
  for (i =0; i < 4; ++i  )
    for (j = 0; j < NROW; ++j)
      
      if ( (M(j,i) < deduct) && (mch < deduct)) {
        oopspend(j, i) = allspend(j, i);
      }else if( (mch < deduct) && ( M(j, i) >= deduct ) && (deduct + copay*(M(j, i)-deduct) < maxoop)) {
        oopspend(j, i) = deduct - mch + copay* ( allspend(j, i) -deduct );
      } else if( mch >= deduct && deduct + copay*(mch-deduct)<maxoop && deduct + copay*(M(j,i)-deduct)<maxoop ){
        oopspend(j, i) = copay* allspend(j, i);
      } else if(mch >= deduct && deduct+copay*(mch-deduct)<maxoop && deduct+copay*(M( j, i) - deduct) >= maxoop ){
        oopspend(j, i) = maxoop - deduct-copay *(mch - deduct);
      } else if (deduct+copay*(mch-deduct) >= maxoop ){
        oopspend(j, i) = 0;
      } else if(mch < deduct & deduct + copay*(M (j, i)-deduct) >= maxoop){
        oopspend(j, i) = maxoop - mch;
        // This is the case where transitory shocks are so huge that you skip the middle case (unlikely in eqbm)
      }
      
      // Construct expected utility for each possible solution
      NumericMatrix util (NROW , 4);
      
      //Rf_PrintValue(oopspend);
      for (i =0; i < 4; ++i  ) {
        for (j = 0; j < NROW; ++j){
          if (lambda(j)<0) { // if lambda is negative, make sure utility is 0. 
            util(j,i) = 0; 
          } else{
            util(j,i) = p(j)*(allspend(j,i)-lambda(j)-1/(2*omega)*(allspend(j,i)-lambda(j))*(allspend(j,i)-lambda(j))) + 
              (1-p(j))*(alpha_fn1*allspend(j,i)+alpha_fn2*mch-lambda(j)-
              1/(2*omega)*(alpha_fn1*allspend(j,i)+alpha_fn2*mch-lambda(j))*(alpha_fn1*allspend(j,i)+alpha_fn2*mch-lambda(j))) -
              - oopspend(j,i);
          }
        }
      }
      //Rf_PrintValue(util);
      
      // Return utility maximizing amount (in 1000s of dollars)
      NumericVector utildestar(NROW);
      for (j = 0; j < NROW; ++j){
        utildestar(j) = util(j,0);
        for (i = 1; i < 4; ++i) {
          if (util(j,i) > utildestar(j)) {
            utildestar(j) = util(j,i);
          }
        }
      }
      
      double ss_out;
      for (j = 0; j < NROW; ++j) {
        ss_out += utildestar(j)/1000;
        // Report utility in 1000s of dollars
      }
      //Rcout << ss_out << "\n";
      
      return(ss_out); 
}

// [[Rcpp::export]]
NumericVector findChoices_cpp(NumericVector myfam_p, NumericVector myfam_muL, NumericVector myfam_sigL, NumericVector myfam_kappa,
                              double alpha_fn1, double alpha_fn2, NumericVector imch_p, NumericVector imch_w, NumericVector xint, NumericVector wint,
                              NumericVector allp_deduct, NumericVector allp_maxoop, NumericVector allp_copay, NumericVector allp_prem, 
                              double sigma_E, double psi, double omega_fn, NumericVector ispend, double sc, double choice_noswitch) {
  
  int i,j,k,ii; // loop through all available plans
  int myNROW = xint.size();
  int myNCOL = imch_p.size();
  int isize = ispend.size();
  double mch,c_mch;
  NumericVector lambda(isize);
  NumericVector exVal(6);
  NumericVector exPayoff(6);
  NumericMatrix mywt = (myNROW,myNCOL);
  NumericMatrix u = (myNROW,myNCOL);
  NumericMatrix u_ce = (myNROW,myNCOL); // certainty equivalent
  
  // MPFR objects
  mpfr_t lamwork;
  mpfr_init2(lamwork, 512);
  
  mpfr_t myce,working,weight,myexp;
  mpfr_init2(myce, 512); // sum of certainty equivalent (in mpfr)
  mpfr_init2(working,512); // working MPFR variable
  double myce_out; 
  
  mpfr_t cp1,cp2,cp3,cp4,cp5,cp6,sum;
  mpfr_init2(cp1,512);
  mpfr_init2(cp2,512);
  mpfr_init2(cp3,512);
  mpfr_init2(cp4,512);
  mpfr_init2(cp5,512);
  mpfr_init2(cp6,512);
  mpfr_init2(sum,512);
  mpfr_set_d(sum,0,GMP_RNDN);
  NumericVector choiceProbs(6);
  
  // Count number of plans
  int np; 
  np = 0;
  for (int i=0; i<6; ++i) {
    if (!isnan(allp_deduct(i))) {
      np = np + 1;
    }
  }
  if (np == 0) { 
    NumericVector NanVec(6);
    for (int j=0; j<6; ++j) NanVec(j) = NAN;
    return(NanVec);
  }
  
  // Initialize exVal and exPayoff to NaNs, not 0s
  for (int i=0; i<6; ++i) exVal(i) = NAN;
  for (int i=0; i<6; ++i) exPayoff(i) = NAN;
  
  // Build weight matrix
  for (j=0; j < myNROW; ++j) {
      for (i=0; i < myNCOL; ++i) {
        mywt(j,i) = wint(j)*imch_w(i); 
      }
    }
  
  for(k = 0; k < np; ++k) {
      if (allp_deduct(k) < 0 || isnan(allp_deduct(k))) {
         continue; // Skip to next period in iteration if no deductible information
      }
          
    exVal(k) = 0; // Move from NA to 0 if not missing information
    exPayoff(k) = 0;
                
    // Numeric (double) integration; loop over quadrature points for (i) mch at family level and (j) lambda
     for(j=0; j<myNROW; ++j) { 
        // j fixes the lambda amount (acute health shocks)
        
        for (i=0; i<myNCOL; ++i) {
          // i fixes the MCH amount and calculates related OOP (spending on chronic care)
          
          // Skip if weight is smaller than 1e-6. 
          if (mywt(j,i) < 1.0e-6) {
            continue; 
          }
          else { 
            mch = imch_p(i);
            if (mch < allp_deduct(k)) {
              c_mch = mch;
            } else if (mch >= allp_deduct(k) && allp_deduct(k)+allp_copay(k)*(mch-allp_deduct(k))<allp_maxoop(k)) {
              c_mch = allp_deduct(k)+allp_copay(k)*(mch-allp_deduct(k));
            } else {
              c_mch = allp_maxoop(k);
            }
            
            for(ii=0; ii<isize; ++ii) {
              mpfr_set_d(lamwork,xint(j)*myfam_sigL(ii) + myfam_muL(ii) + myfam_kappa(ii),MPFR_RNDN);
              mpfr_exp(lamwork,lamwork,GMP_RNDN); // Calculate lambda for each family member at once
              lambda(ii) = mpfr_get_d(lamwork,GMP_RNDN);
              
              //Rcout << lambda(ii) << "\n"; 
              if (isnan(ispend(ii))) {
                lambda(ii) = -9999; // Revert lambda to a large negative number if a household member is missing (so that they are assigned 0 spending)
              }
            }
            
            //Rf_PrintValue(lambda);
            
            //u(j,i) = solveSpend_cpp(myfam_p,lambda,omega_fn,alpha_fn1,alpha_fn2,allp_deduct(k),allp_copay(k),allp_maxoop(k),mch);
            //Rcout << u(j,i) << "\n";
            exPayoff(k) = exPayoff(k) + u(j,i)*mywt(j,i);
            //Rcout << exPayoff(k) << "\n";
          }
       }
     }
                
    // Certainty equivalent expected value for plan k (after calculating the full exPayoff)
    mpfr_set_d (myce, 0, MPFR_RNDN); // reset the sum to 0 for each plan
    for(j=0; j<myNROW; ++j) {
      for (i = 0; i<myNCOL; ++i) {
        mpfr_set_d (working,-psi*(u(j,i)-exPayoff(k)),MPFR_RNDN); // set working to the exponent as a double
        mpfr_exp(working,working,GMP_RNDN); // exponentiate the working variable
        mpfr_mul_d(working,working,mywt(j,i),GMP_RNDN); // convert working to the product (as an integer)
        mpfr_add(myce,myce,working,GMP_RNDN); // add working to the sum in progress
      }
    }
    mpfr_log(myce,myce,GMP_RNDN); // take log of sum once finished
    myce_out = mpfr_get_d(myce,GMP_RNDN); // convert back to double
    //Rcout << myce_out << "\n";
    
    // Assign switching costs if the plan being examined is not the given choice_noswitch
    // NOTE: require that choice_noswitch is in the choice set before imposing switching costs
    if (isnan(choice_noswitch) || choice_noswitch == k) {
          exVal(k) = exPayoff(k) - (1/psi)*myce_out - psi*allp_prem(k)/1000 - psi*c_mch/1000;
      } else {
          exVal(k) = exPayoff(k) - (1/psi)*myce_out - psi*allp_prem(k)/1000 - psi*c_mch/1000 - sc;
      }
  } // Finishes loop through plans 
  
  // Calculate choice probabilities
  mpfr_set_d(cp1,exVal(0)/sigma_E,GMP_RNDN);
  mpfr_exp(cp1,cp1,GMP_RNDN);
  mpfr_set_d(cp2,exVal(1)/sigma_E,GMP_RNDN);
  mpfr_exp(cp2,cp2,GMP_RNDN);
  mpfr_set_d(cp3,exVal(2)/sigma_E,GMP_RNDN);
  mpfr_exp(cp3,cp3,GMP_RNDN);
  mpfr_set_d(cp4,exVal(3)/sigma_E,GMP_RNDN);
  mpfr_exp(cp4,cp4,GMP_RNDN);
  mpfr_set_d(cp5,exVal(4)/sigma_E,GMP_RNDN);
  mpfr_exp(cp5,cp5,GMP_RNDN);
  mpfr_set_d(cp6,exVal(5)/sigma_E,GMP_RNDN);
  mpfr_exp(cp6,cp6,GMP_RNDN);
  
  // denominator is sum (all non-NANs)
  if (!isnan(mpfr_get_d(cp1,GMP_RNDN))) {
    mpfr_add(sum,sum,cp1,GMP_RNDN);
  }
  if (!isnan(mpfr_get_d(cp2,GMP_RNDN))) {
    mpfr_add(sum,sum,cp2,GMP_RNDN);
  }
  if (!isnan(mpfr_get_d(cp3,GMP_RNDN))) {
    mpfr_add(sum,sum,cp3,GMP_RNDN);
  }
  if (!isnan(mpfr_get_d(cp4,GMP_RNDN))) {
    mpfr_add(sum,sum,cp4,GMP_RNDN);
  }
  if (!isnan(mpfr_get_d(cp5,GMP_RNDN))) {
    mpfr_add(sum,sum,cp5,GMP_RNDN);
  }
  if (!isnan(mpfr_get_d(cp6,GMP_RNDN))) {
    mpfr_add(sum,sum,cp6,GMP_RNDN);
  }
  
  mpfr_div(cp1,cp1,sum,GMP_RNDN);
  mpfr_div(cp2,cp2,sum,GMP_RNDN);
  mpfr_div(cp3,cp3,sum,GMP_RNDN);
  mpfr_div(cp4,cp4,sum,GMP_RNDN);
  mpfr_div(cp5,cp5,sum,GMP_RNDN);
  mpfr_div(cp6,cp6,sum,GMP_RNDN);
  
  choiceProbs(0) = mpfr_get_d(cp1,GMP_RNDN);
  choiceProbs(1) = mpfr_get_d(cp2,GMP_RNDN);
  choiceProbs(2) = mpfr_get_d(cp3,GMP_RNDN);
  choiceProbs(3) = mpfr_get_d(cp4,GMP_RNDN);
  choiceProbs(4) = mpfr_get_d(cp5,GMP_RNDN);
  choiceProbs(5) = mpfr_get_d(cp6,GMP_RNDN);
    
  return(choiceProbs);
}
