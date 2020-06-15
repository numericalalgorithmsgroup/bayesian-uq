
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
  vector[K] m_theta;
  vector[K] s_theta;
}

parameters {
  vector[K] theta;
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

  // Prior
  for (k in 1:K){
    theta[k] ~ normal(m_theta[k], s_theta[k]);
  }
}
