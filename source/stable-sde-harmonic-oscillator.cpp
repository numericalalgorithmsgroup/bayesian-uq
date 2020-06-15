#include <vector>
#include <iostream>
#include "Eigen/Dense"

#include "scalar-typedef.hpp"
#include "stable-sde.hpp"
#include "stable-sde-harmonic-oscillator.hpp"

// Explicit template instantiation
template class Stable_sde_harmonic_oscillator<double>;
#ifdef USE_DCO_TYPES
template class Stable_sde_harmonic_oscillator<gt1s_scalar>;
template class Stable_sde_harmonic_oscillator<ga1s_scalar>;
template class Stable_sde_harmonic_oscillator<gt2s_ga1s_scalar>;
#endif

enum p_harmonic_oscillator { OMEGA0, ZETA, SD_IN, SD_OBS};
// enum fp_harmonic_oscillator { SD_OBS};

// note underscore at end of each element of enumeration
enum theta_post_harmonic_oscillator { LOG_OMEGA0_C1_, LOG_OMEGA0_C2_, LOG_SD_IN_C1_, LOG_SD_IN_C2_, LOG_ZETA_};
enum theta_prior_harmonic_oscillator { OMEGA0_C1_, OMEGA0_C2_, SD_IN_C1_, SD_IN_C2_, ZETA_};


using namespace Eigen;
using namespace std;

template<typename SCALAR_T>
Matrix<SCALAR_T, Dynamic, 1> Stable_sde_harmonic_oscillator<SCALAR_T>::dx_dt(const Matrix<SCALAR_T, Dynamic, 1> & x,
                                                                             const Matrix<SCALAR_T, Dynamic, 1> & p) {
  Matrix<SCALAR_T, Dynamic, 1> F( x.size() );
  F[0] =  x[1];
  F[1] = -  pow( p[OMEGA0], 2.0) * x[0] - 2 * p[ZETA] * p[OMEGA0] * x[1] ;
  return(F);
}

template<typename SCALAR_T>
Matrix<SCALAR_T, Dynamic, Dynamic> Stable_sde_harmonic_oscillator<SCALAR_T>::sde_jacobian(const Matrix<SCALAR_T, Dynamic, 1> & x,
                                                                                          const Matrix<SCALAR_T, Dynamic, 1> & p) {
  J = Matrix<SCALAR_T, Dynamic, Dynamic>::Zero(d, d);
  // define Jacobian for harmonic oscillator
  SCALAR_T omega0 = p[0];
  SCALAR_T zeta = p[1];
  J(0,0) = 0; J(0,1) = 1.0 ;
  J(1,0) = - pow( p[OMEGA0], 2.0); J(1,1) = -2 * p[ZETA] * p[OMEGA0];

  return J;
}

template<typename SCALAR_T>
SCALAR_T Stable_sde_harmonic_oscillator<SCALAR_T>::eval_sd_in() {
  return p[SD_IN];
}

template<typename SCALAR_T>
SCALAR_T Stable_sde_harmonic_oscillator<SCALAR_T>::eval_sd_obs() {
  return p[SD_OBS];
}

template<typename SCALAR_T>
void Stable_sde_harmonic_oscillator<SCALAR_T>::set_model_vars_c1(const Matrix<SCALAR_T, Dynamic, 1> & theta) {
  // elements that go into x
  x_star[0] = 0.0;
  x_star[1] = 0.0;

  // elements that go into p
  // p.resize( 4 );
  p[OMEGA0] = exp( theta[LOG_OMEGA0_C1_] );
  p[ZETA]   = exp( theta[LOG_ZETA_] );
  p[SD_IN]  = exp( theta[LOG_SD_IN_C1_] );
}

template<typename SCALAR_T>
void Stable_sde_harmonic_oscillator<SCALAR_T>::set_model_vars_c2(const Matrix<SCALAR_T, Dynamic, 1> & theta) {
  // elements that go into x
  x_star[0] = 0.0;
  x_star[1] = 0.0;

  // elements that go into p
  // p.resize( 4 );
  p[OMEGA0] = exp( theta[LOG_OMEGA0_C2_] );
  p[ZETA]   = exp( theta[LOG_ZETA_] );
  p[SD_IN] =  exp( theta[LOG_SD_IN_C2_] );
}

template<typename SCALAR_T>
Matrix<SCALAR_T,Dynamic,1> Stable_sde_harmonic_oscillator<SCALAR_T>::g(const Matrix<SCALAR_T,Dynamic,1> & theta_post) const
{
  Matrix<SCALAR_T,Dynamic,1> theta_prior( theta_post.size() );
  theta_prior[OMEGA0_C1_] = exp(theta_post[LOG_OMEGA0_C1_]);
  theta_prior[OMEGA0_C2_] = exp(theta_post[LOG_OMEGA0_C2_]);
  theta_prior[SD_IN_C1_]  = exp(theta_post[LOG_SD_IN_C1_]);
  theta_prior[SD_IN_C2_]  = exp(theta_post[LOG_SD_IN_C2_]);
  theta_prior[ZETA_]      = exp(theta_post[LOG_ZETA_]);
  return theta_prior;
}

template<typename SCALAR_T>
Matrix<SCALAR_T,Dynamic,Dynamic> Stable_sde_harmonic_oscillator<SCALAR_T>::g_jacobian(const Matrix<SCALAR_T,Dynamic,1> & theta) const
{
  Matrix<SCALAR_T,Dynamic,Dynamic> J = MatrixXd::Zero( theta.size(), theta.size() );
  J(0,0) = exp(theta[LOG_OMEGA0_C1_]);
  J(1,1) = exp(theta[LOG_OMEGA0_C2_]);
  J(2,2) = exp(theta[LOG_SD_IN_C1_]);
  J(3,3) = exp(theta[LOG_SD_IN_C2_]);
  J(4,4) = exp(theta[LOG_ZETA_]);
  return J;
}
