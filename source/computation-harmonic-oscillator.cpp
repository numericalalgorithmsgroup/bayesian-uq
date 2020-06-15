#include <vector>
#include "Eigen/Dense"

#include "scalar-typedef.hpp"

#include "stable-sde.hpp"
#include "stable-sde-harmonic-oscillator.hpp"
#include "data-vector.hpp"

#include "computation.hpp"
#include "computation-harmonic-oscillator.hpp"

using namespace std;
using namespace Eigen;

// Explicit template instantiation(s)
template class Computation_harmonic_oscillator<double>;
#ifdef USE_DCO_TYPES
template class Computation_harmonic_oscillator<ga1s_scalar>;
template class Computation_harmonic_oscillator<gt1s_scalar>;
template class Computation_harmonic_oscillator<gt2s_ga1s_scalar>;
#endif

template<typename SCALAR_T>
SCALAR_T Computation_harmonic_oscillator<SCALAR_T>::eval(Matrix<SCALAR_T, Dynamic, 1>& theta,
                         const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                         unsigned int & ifail){
  unsigned int  ifail_m = 0, print_flag=0;
  SCALAR_T ll = 0;
  double inf = std::numeric_limits<double>::infinity();
  // initialize Model class
  Stable_sde_harmonic_oscillator<SCALAR_T> m;
  // these members only need to be initialized once
  m.set_d( 2 );           // dimension of state variables
  m.set_input_ind( 1 );    // index for noise input
  m.set_output_ind( 0 );   // index for observation state variable
  // set parameters - needed for fixed parameters used in set_model_vars
  m.set_param_unknown(p);
  // evaluate likelihood at c1
  m.set_model_vars_c1(theta);
  ll += log_whittle_like(EvalWhittle, m, y_light, ifail_m);
  if (ifail_m){
    ifail = 4;
    return -inf;
  }

  // evaluate likelihood at c2
  m.set_model_vars_c2(theta);
  ll += log_whittle_like(EvalWhittle, m, y_deep, ifail_m);
  if (ifail_m){
    ifail = 6;
    return -inf;
  }
  return ll;
}

template<typename SCALAR_T>
Matrix<SCALAR_T, Dynamic, 1> Computation_harmonic_oscillator<SCALAR_T>::find_ss_c0(const Matrix<SCALAR_T, Dynamic, 1>& theta){
  Matrix<SCALAR_T, Dynamic, 1> x_star(2);
  x_star << 0.0, 0.0;
  return x_star;
}
