#include <vector>
#include <fstream>
#include <limits>
#include <iostream>
#include "float.h"
#include "Eigen/Dense"
#include "scalar-typedef.hpp"

#include "stable-sde.hpp"
#include "data-vector.hpp"
#include "vector-io.hpp"
#include "computation.hpp"

using namespace std;
using namespace Eigen;

// Explicit template instantiation(s)
template class Computation<double>;
template double log_whittle_like(Likelihood_mode LikelihoodMode,
                                 Stable_sde<double> & m,
                                 const Data_vector<double> & y,
                                 unsigned int & ifail);

#ifdef USE_DCO_TYPES
template class Computation<gt1s_scalar>;
template gt1s_scalar log_whittle_like(Likelihood_mode LikelihoodMode,
                                 Stable_sde<gt1s_scalar> & m,
                                 const Data_vector<gt1s_scalar> & y,
                                 unsigned int & ifail);

template class Computation<ga1s_scalar>;
template ga1s_scalar log_whittle_like(Likelihood_mode LikelihoodMode,
                                 Stable_sde<ga1s_scalar> & m,
                                 const Data_vector<ga1s_scalar> & y,
                                 unsigned int & ifail);

template class Computation<gt2s_ga1s_scalar>;
template gt2s_ga1s_scalar log_whittle_like(Likelihood_mode LikelihoodMode,
                                 Stable_sde<gt2s_ga1s_scalar> & m,
                                 const Data_vector<gt2s_ga1s_scalar> & y,
                                 unsigned int & ifail);
#endif

template<typename SCALAR_T>
SCALAR_T log_whittle_like(Likelihood_mode LikelihoodMode,
                        Stable_sde<SCALAR_T> & m,
                        const Data_vector<SCALAR_T> & y,
                        unsigned int & ifail){
  int print_flag = 0;
  double inf = std::numeric_limits<double>::infinity();
  // extract variables from y
  Matrix<SCALAR_T, Dynamic, 1> xf = y.get_freq();
  Matrix<SCALAR_T, Dynamic, 1> I = y.get_periodogram();
  int n_freq = y.get_nfreq();
  SCALAR_T tstep = y.get_timestep();

  // evaluate spectral density
  Matrix<SCALAR_T, Dynamic, 1> spec_dens = m.spec_fun(xf, tstep, ifail);

  if (print_flag) cout << "frequencies:" << endl;
  if (print_flag) cout << xf << endl;

  if (print_flag) m.print_members();

  if (print_flag) cout << "periodogram:" << endl;
  if (print_flag) cout << I << endl;

  if (print_flag) cout << "spectral density:" << endl;
  if (print_flag) cout << spec_dens << endl;

  if ( ifail){
    // fixed point is unstable
    return -inf;
  }

  SCALAR_T ll = 0.0;
  switch(LikelihoodMode){
  case EvalWhittle:
    // evaluate Whittle likelihood from periodogram and spectral density
    for (int k = 0; k < n_freq; ++k){
      ll -= log( spec_dens[k] ) + I[k] / spec_dens[k];
    }
    break;
  }
  // when conc=0, this function just checks for stability and returns ll=0 when fixed point is stable
  return ll;
}

#ifdef USE_FFD
template<>
Matrix<double, Dynamic, 1> Computation<double>::grad(Matrix<double, Dynamic, 1>& theta0,
                                                     const Eigen::Matrix<double, Eigen::Dynamic, 1>& xstar_c0_in,
                                                     unsigned int & ifail) {
  double ll0 = eval(theta0, xstar_c0_in, ifail);
  const int n_theta = theta0.size();
  Matrix<double, Dynamic, 1> grad(n_theta);
  for (int i=0; i < n_theta; i++){
    Matrix<double, Dynamic, 1> theta = theta0;
    // if theta(i) is zero, set h to be sqrt of machine precision
    // otherwise set h to be (sqrt of machine precision) * ( absolute value of theta(i) )
    double h = theta0(i) == 0.0 ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON) * abs( theta0(i)  );
    theta(i) += h;
    double ll = eval(theta, xstar_c0_in, ifail);
    grad(i) = (ll - ll0) / h;
  }
  return grad;
}
#endif

#ifdef USE_CFD
template<>
Matrix<double, Dynamic, 1> Computation<double>::grad(Matrix<double, Dynamic, 1>& theta0,
                                                     const Eigen::Matrix<double, Eigen::Dynamic, 1>& xstar_c0_in,
                                                     unsigned int & ifail) {
  const int n_theta = theta0.size();
  Matrix<double, Dynamic, 1> grad(n_theta);
  for (int i=0; i < n_theta; i++){
    Matrix<double, Dynamic, 1> theta_p = theta0, theta_n = theta0;
    // if theta(i) is zero, set h to be sqrt of machine precision
    // otherwise set h to be (sqrt of machine precision) * ( absolute value of theta(i) )
    double h = theta0(i) == 0.0 ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON) * abs( theta0(i)  );
    theta_p(i) += h;
    theta_n(i) -= h;
    double ll_p = eval(theta_p, xstar_c0_in, ifail);
    double ll_n = eval(theta_n, xstar_c0_in, ifail);
    grad(i) = (ll_p - ll_n) / (2*h);
  }
  return grad;
}
#endif

template<>
Matrix<double, Dynamic, Dynamic> Computation<double>::hess(Matrix<double, Dynamic, 1>& theta0,
                                                           const Eigen::Matrix<double, Eigen::Dynamic, 1>& xstar_c0_in,
                                                           unsigned int & ifail) {
  const int n_theta = theta0.size();
  Matrix<double, Dynamic, Dynamic> hess(n_theta, n_theta);
  for (int i=0; i < n_theta; i++){
    for (int j=0; j < n_theta; j++){
      hess(i,j) = 0.0;
    }
  }
  for (int i=0; i < n_theta; i++){
    for (int j=0; j <= i; j++){
      // if theta(i) is zero, set h to be sqrt of machine precision
      // otherwise set h to be (sqrt of machine precision) * ( absolute value of theta(i) )
      double h_i = theta0(i) == 0.0 ? pow(DBL_EPSILON, 1.0/3.0) : pow(DBL_EPSILON, 1.0/3.0) * abs( theta0(i)  );
      double h_j = theta0(j) == 0.0 ? pow(DBL_EPSILON, 1.0/3.0) : pow(DBL_EPSILON, 1.0/3.0) * abs( theta0(j)  );
      double h = min(h_i, h_j);
      Matrix<double, Dynamic, 1> theta_pp = theta0, theta_np = theta0, theta_pn = theta0, theta_nn = theta0;
      theta_pp(i) += h; theta_pp(j) += h;
      theta_np(i) -= h; theta_np(j) += h;
      theta_pn(i) += h; theta_pn(j) -= h;
      theta_nn(i) -= h; theta_nn(j) -= h;
      double ll_pp = eval(theta_pp, xstar_c0_in, ifail);
      double ll_np = eval(theta_np, xstar_c0_in, ifail);
      double ll_pn = eval(theta_pn, xstar_c0_in, ifail);
      double ll_nn = eval(theta_nn, xstar_c0_in, ifail);
      if (ifail){
        hess(i,j) = NAN;
        return hess;
      }
      hess(i,j) =  - (ll_pp - ll_np - ll_pn + ll_nn) / ( 4 * pow(h,2.0) );
    }
  }
  for (int i=0; i < n_theta; i++){
    for (int j=0; j <= i; j++){
      hess(j,i) = hess(i,j);
    }
  }
  return hess;
}

template<>
void Computation<double>::derivatives(Eigen::Matrix<double, Eigen::Dynamic, 1>& theta,
                                      const Eigen::Matrix<double, Eigen::Dynamic, 1>& x_star,
                                      Eigen::Matrix<double, Eigen::Dynamic, 1> & grad,
                                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & hess,
                                      unsigned int & ifail)
{
  // unsigned int ifail=0;
  const int theta_size = theta.size();
  // Matrix<double, Dynamic, Dynamic> hess(theta_size,theta_size);
  // Matrix<double, Dynamic, 1> grad(theta_size,1);
  grad = Computation<double>::grad(theta, x_star, ifail);
  hess = Computation<double>::hess(theta, x_star, ifail);
  return;
}

#ifdef USE_DCO_TYPES
template<>
Matrix<double,Dynamic,1> Computation<ga1s_scalar>::grad(Matrix<ga1s_scalar, Dynamic, 1>& theta,
                                                   const Eigen::Matrix<ga1s_scalar, Eigen::Dynamic, 1>& xstar_c0_in,
                                                   unsigned int & ifail) {
  const int n_theta = theta.size();
  Matrix<double,Dynamic,1>  grad(n_theta);

  // Create tape
  DCO_GA1S_MODE::global_tape = DCO_GA1S_MODE::tape_t::create();

  for (int i=0; i < n_theta; i++){
    DCO_GA1S_MODE::global_tape->register_variable( theta(i) );
  }

  ga1s_scalar ll = eval(theta, xstar_c0_in, ifail);
  dco::derivative( ll ) = 1.0;

  DCO_GA1S_MODE::global_tape->interpret_adjoint();

  for (int i=0; i < n_theta; i++){
    grad(i) = dco::value( dco::derivative( theta(i) ) );
  }
  DCO_GA1S_MODE::global_tape->reset();
  return grad;
}
#endif

#undef GRAD_GT2S_GA1S
#ifdef GRAD_GT2S_GA1S
template<>
Matrix<double,Dynamic,1> Computation<Scalar>::grad(Matrix<Scalar, Dynamic, 1>& theta,
                                                   const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x_star,
                                                   unsigned int & ifail) {
  const int theta_size = theta.size();
  Matrix<double, Dynamic, 1> grad(theta_size);
  // Create tape
  DCO_MODE::global_tape = DCO_MODE::tape_t::create();

  for(int i=0; i<theta_size; i++) {
    // Register inputs
    for(int j=0; j<theta_size; j++) {
      DCO_MODE::global_tape->register_variable(theta(j));
    }
    dco::derivative(dco::value(theta(i))) = 1.0;

    // Actual computation to be differentiated
    Scalar ll = eval(theta, x_star, ifail);

    DCO_MODE::global_tape->register_output_variable(ll);
    dco::derivative(ll) = 1.0;

    DCO_MODE::global_tape->interpret_adjoint();

    // Fill Jacobian
    grad(i) = dco::value(dco::derivative(theta(i)));

    dco::derivative(dco::value(theta(i))) = 0.0;
    DCO_MODE::global_tape->reset();
  }
  return grad;
}
#endif

#undef HESS_GT2S_GA1S
#ifdef HESS_GT2S_GA1S
template<>
Matrix<double,Dynamic,Dynamic> Computation<Scalar>::hess(Matrix<Scalar, Dynamic, 1>& theta,
                                                         const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x_star,
                                                         unsigned int & ifail) {
  const int theta_size = theta.size();
  Matrix<double, Dynamic, Dynamic> hess(theta_size,theta_size);
  // Create tape
  DCO_MODE::global_tape = DCO_MODE::tape_t::create();

  for(int i=0; i<theta_size; i++) {
    // Register inputs
    for(int j=0; j<theta_size; j++) {
      DCO_MODE::global_tape->register_variable(theta(j));
    }
    dco::derivative(dco::value(theta(i))) = 1.0;

    // Actual computation to be differentiated
    Scalar ll = eval(theta, x_star, ifail);

    DCO_MODE::global_tape->register_output_variable(ll);
    dco::derivative(ll) = 1.0;

    DCO_MODE::global_tape->interpret_adjoint();

    // Fill Hessian
    for(int j=0; j<theta_size; j++) hess(j, i) = dco::derivative(dco::derivative(theta(j)));

    dco::derivative(dco::value(theta(i))) = 0.0;
    DCO_MODE::global_tape->reset();
  }

  return hess;
}
#endif

template<typename SCALAR_T>
std::string Computation<SCALAR_T>::report(const unsigned int & ifail) const
{
  std::string message;
  switch(ifail) {
  case 0 :
    message = "SDE has stable fixed points for c=0, c_light, c_deep";
    break;
  case 1 :
    message = "Steady state solver for c=0 did not converge";
    break;
  case 2 :
    message = "Steady state is unstable for c=0";
    break;
  case 3 :
    message = "Steady state solver for c_light did not converge";
    break;
  case 4 :
    message = "Steady state is unstable for c_light";
    break;
  case 5 :
    message = "Steady state solver for c_deep did not converge";
    break;
  case 6 :
    message = "Steady state is unstable for c_deep";
    break;
  case 7 :
    message = "Proposed parameters outside support of prior distribution";
    break;
  default :
    message = "ifail != 0";
    break;
  }
  return message;
}
