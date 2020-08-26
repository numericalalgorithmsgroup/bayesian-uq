#include <vector>
#include "Eigen/Dense"
#include "scalar-typedef.hpp"

#include "stable-sde.hpp"
#include "data-vector.hpp"
#include "vector-io.hpp"
#include "computation.hpp"

using namespace std;
using namespace Eigen;

#ifdef USE_DCO_TYPES

template<>
void Computation<gt2s_ga1s_scalar>::derivatives(Matrix<gt2s_ga1s_scalar, Dynamic, 1>& theta,
                                      const Eigen::Matrix<gt2s_ga1s_scalar, Eigen::Dynamic, 1>& x_star,
                                      Matrix<double,Dynamic,1> & grad,
                                      Matrix<double,Dynamic,Dynamic> & hess,
                                      unsigned int & ifail){
  const int theta_size = theta.size();
  grad.resize(theta_size);
  hess.resize(theta_size,theta_size);
  // Create tape
  if( DCO_GA1S_MODE::global_tape==nullptr) {
    DCO_GT2S_GA1S_MODE::global_tape = DCO_GT2S_GA1S_MODE::tape_t::create();
  }
  auto pos = DCO_GT2S_GA1S_MODE::global_tape->get_position();

  for(int i=0; i<theta_size; i++) {
    // Register inputs
    for(int j=0; j<theta_size; j++) {
      DCO_GT2S_GA1S_MODE::global_tape->register_variable(theta(j));
    }
    dco::derivative(dco::value(theta(i))) = 1.0;

    // Actual computation to be differentiated
    gt2s_ga1s_scalar ll = eval(theta, x_star, ifail);

    DCO_GT2S_GA1S_MODE::global_tape->register_output_variable(ll);
    dco::derivative(ll) = 1.0;

    DCO_GT2S_GA1S_MODE::global_tape->interpret_adjoint();

    // Fill Jacobian
    grad(i) = dco::value(dco::derivative(theta(i)));
    // Fill Hessian
    for(int j=0; j<theta_size; j++) hess(j, i) = dco::derivative(dco::derivative(theta(j)));

    dco::derivative(dco::value(theta(i))) = 0.0;
    DCO_GT2S_GA1S_MODE::global_tape->reset();
  }
  DCO_GT2S_GA1S_MODE::global_tape->reset_to( pos );
  DCO_GT2S_GA1S_MODE::tape_t::remove(DCO_GT2S_GA1S_MODE::global_tape);
  return;
}
#endif
