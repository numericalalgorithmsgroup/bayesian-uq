
functions {
  vector spec_fun(vector p, vector theta, int dataset);
}

data {
  int<lower=1> M;        // number of frequencies / length of p in data
  int<lower=1> M_pred;   // number of frequencies / length of p in generated data
  vector[M] p;           // positions / frequencies
  vector[M_pred] p_pred; // positions / frequencies
  vector[M] y1;           // periodogram / observations
  vector[M] y2;           // periodogram / observations
  // means and standard deviations of prior parameters
  int<lower=1> K;
}

parameters {
  vector[K] theta;
}

transformed parameters {
  real omega0_c1 = exp(theta[1]);
  real omega0_c2 = exp(theta[2]);
  real sd_in_c1 = exp(theta[3]);
  real sd_in_c2 = exp(theta[4]);
  real zeta = exp(theta[5]);
}

model {
  // Likelihood
  vector[M] spec1;
  vector[M] spec2;
  spec1 = spec_fun(p, theta, 1);
  for (m in 1:M){
    y1[m] ~ exponential( 1 / spec1[m]);
  }

  spec2 = spec_fun(p, theta, 2);
  for (m in 1:M){
    y2[m] ~ exponential( 1 / spec2[m]);
  }
}
