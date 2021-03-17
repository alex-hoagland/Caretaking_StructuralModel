# Caretaking_StructuralModel
This code estimates the structural model as seen here (link to be added). The goal of this code is to optimize a likelihood function of about 44 parameters. The model is one of health care spending, in which families first choose an insurance plan from a discrete set of choices, and then choose a continuous amount of medical spending throughout the year. The model includes various degrees of unobserved heterogeneity, which are integrated numerically using Gaussian quadrature. 

The file organization is as follows: 
1. mainEstimation.R loads all data and functions, sets the initial parameters, and runs a nonlinear optimizer over the calculateLikelihood() function. 
2. calculateLikelihood.R contains the main function, which calculates the likelihood function for each household in the sample. This includes the initialization of the Gaussian quadrature types across households (denoted by the "si" variable in loops). The function then executes a parallel loop (%dopar%) across all households. Returns a vector of log(likelihood) for each household. 
3. spendingDensities.R is a subfunction in calculateLikelihood() that calculates the conditional likelihood of a household's medical spending choices given a plan choice. This is done at the individual level and then aggregated to the household level by assuming conditional independence. 
4. findChoices_cpp.cpp is an Rcpp file that is a subfunction in calculateLikelihood(). This file calculates the choice probabilities for each plan in a household's choice set. 
