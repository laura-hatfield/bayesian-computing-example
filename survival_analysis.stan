//bayesian survival analysis

data{
  int<lower=0> p; //number of covariates
  int<lower=0> k; //number of treatments (not including baseline (BMS))
  int<lower=0> Nobs; //number of non-censored observations
  int<lower=0> Ncen; //number of censored observations
  vector<lower=0>[Nobs] Yobs; //survival time for non-censored observations
  vector<lower=0>[Ncen] Ycen; //survival time for censored observations
  matrix[Nobs, p] Xobs; //covariates for non-censored obs
  matrix[Ncen, p] Xcen; //covariates for censored obs
  matrix[Nobs, k] Tobs; //treatments for non-censored
  matrix[Ncen, k] Tcen; //treatments for censored obs 
}


parameters{
  real<lower=0> alpha; //shape parameter
  real beta_0; //intercept (baseline treatment = BMS)
  vector[k] beta_t; //treatment coefficients (contrasts)
  vector[p] beta_x; //predictor coefficients 
}


model{
  alpha ~ normal(0, 1);
  beta_0 ~ normal(0,10);
  beta_t ~ normal(0,10);
  beta_x ~ normal(0,10);
  
  Yobs ~ weibull(alpha, exp(-(beta_0 + Tobs * beta_t + Xobs * beta_x) / alpha));
  target += weibull_lccdf(Ycen | alpha, exp(-(beta_0 + Tcen * beta_t + Xcen * beta_x) / alpha));
}


