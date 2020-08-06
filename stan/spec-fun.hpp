#ifndef STAN_SPEC_FUN_H_
#define STAN_SPEC_FUN_H_

#include <stan/math.hpp>
// from 3rd party libraries
#define _USE_MATH_DEFINES
#include <math.h>
// from dco Eigen branch
#include "Eigen/Dense"
#include "unsupported/Eigen/dco"
// from bayesian-uq code base
#include "scalar-typedef.hpp"
#include "stable-sde.hpp"
#include "stable-sde-harmonic-oscillator.hpp"

using namespace std;
using namespace Eigen;

template<typename SCALAR_T>
Matrix<SCALAR_T, Dynamic, 1> computation(const Eigen::Matrix<double, Eigen::Dynamic, 1>& xf,
                                         const Matrix<SCALAR_T, Dynamic, 1>& theta,
                                         const int& dataset) {
  int overflow_flag = 0;
  for (int i = 0; i < theta.size(); ++i){
    // if (stan::math::is_inf( exp(theta(i)) )){
    if ( exp(theta(i)) > 1e16 ){
      cout << "Warning Message generated in spec-fun.hpp: exp(theta) is diverging." << endl;
      cout << "  This may happen frequently during Warmup.  It should not happen frequently during Sampling." << endl;
      overflow_flag = 1;
    }
  }
  if (overflow_flag)
    return Matrix<SCALAR_T, Dynamic, 1>::Constant(xf.size(), stan::math::positive_infinity () );

  // initialize Model class
  Stable_sde_harmonic_oscillator<SCALAR_T> m;
  m.set_d( 2 );
  m.set_input_ind( 1 );
  m.set_output_ind( 0 );

  // set parameters - needed for fixed parameters used in set_model_vars
  Matrix<SCALAR_T, Dynamic, 1> p(4);
  p[3] = 0.05;
  m.set_param_unknown(p);  // note that first three elements should get overwritten by set_model_vars

  if (dataset == 1)
    m.set_model_vars_c1(theta);
  else if (dataset == 2)
    m.set_model_vars_c2(theta);
  else
    throw std::logic_error("not implemented");  // this should never be called

  SCALAR_T tstep = 0.01;

  // these are (intermediate) outputs that we need to get derivatives of w.r.t. theta
  unsigned int ifail = 0;
  Matrix<SCALAR_T, Dynamic, 1> spec_dens = m.spec_fun(xf, tstep, ifail);

  return spec_dens;
}

namespace spectral_inference_model_namespace {

#ifdef USE_BASIC_TYPES
  template<typename T1, typename T2>
  Matrix<typename boost::math::tools::promote_args<T1, T2>::type, Eigen::Dynamic, 1>
  spec_fun(const Eigen::Matrix<T1, Eigen::Dynamic, 1> & xf,
           const Eigen::Matrix<T2,Eigen::Dynamic,1>& x,
           const int& dataset, std::ostream* pstream__) {
    throw std::logic_error("not implemented");  // this should never be called
  }

  // this is a template specialization for <double, double>
  template<>
  Matrix<double,Eigen::Dynamic, 1>
  spec_fun<double, double>(const Eigen::Matrix<double, Eigen::Dynamic, 1>& xf,
                           const Eigen::Matrix<double,Eigen::Dynamic,1>& x,
                           const int& dataset, std::ostream* pstream__) {
    return computation(xf, x, dataset);
  }

  // this is a template specialization for <double, stan::math::var>
  template<>
  Matrix<stan::math::var,Eigen::Dynamic, 1>
  spec_fun<double, stan::math::var>(const Eigen::Matrix<double, Eigen::Dynamic, 1>& xf,
                                    const Eigen::Matrix<stan::math::var,Eigen::Dynamic,1>& x_stan,
                                    const int& dataset, std::ostream* pstream__) {
    Eigen::Matrix<double,Eigen::Dynamic,1> a = value_of(x_stan);
    Matrix<double, Dynamic, 1> x0_double = a;
    Matrix<double, Dynamic, 1> spec_dens0 = computation(xf, x0_double, dataset); // Actual computation to be differentiated
    const int num_inputs = x0_double.rows();
    int num_outputs = spec_dens0.rows();
    Eigen::Matrix<double,Eigen::Dynamic,1> fa = spec_dens0;
    vector<vector<double>> grad_fa(num_outputs);
    for(int j=0; j<num_outputs; j++) {
      grad_fa[j].resize(num_inputs);
    }
#ifdef USE_FFD
    for(int i=0; i<num_inputs; i++) {
      Matrix<double, Dynamic, 1> x_double = x0_double;
      double h = x_double(i,0) == 0.0 ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON) * abs( x_double(i,0)  );
      x_double(i,0) += h;
      Eigen::Matrix<double,Eigen::Dynamic,1> spec_dens = computation(xf, x_double, dataset); // Actual computation to be differentiated
      for(int j=0; j<num_outputs; j++) grad_fa[j][i] = (spec_dens(j, 0) - spec_dens0(j, 0)) / h;
    }
#endif
#ifdef USE_CFD
    for(int i=0; i<num_inputs; i++) {
      Matrix<double, Dynamic, 1> x_double_p = x0_double, x_double_n = x0_double;
      double h = x_double_p(i,0) == 0.0 ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON) * abs( x_double_p(i,0)  );
      x_double_p(i,0) += h;
      x_double_n(i,0) -= h;
      Eigen::Matrix<double,Eigen::Dynamic,1> spec_dens_p = computation(xf, x_double_p, dataset); // Actual computation to be differentiated
      Eigen::Matrix<double,Eigen::Dynamic,1> spec_dens_n = computation(xf, x_double_n, dataset);
      for(int j=0; j<num_outputs; j++) grad_fa[j][i] = (spec_dens_p(j, 0) - spec_dens_n(j, 0)) / (2*h);
    }
#endif
    vector<stan::math::var> x_std(x_stan.data(), x_stan.data() + x_stan.rows() * x_stan.cols());
    Eigen::Matrix<stan::math::var,Eigen::Dynamic,1> f_var_jacobian( grad_fa.size() );
    for (int i=0; i < grad_fa.size(); ++i){
      f_var_jacobian(i) = precomputed_gradients(fa(i), x_std, grad_fa[i]);
    }
    return f_var_jacobian;
  }
#endif

#ifdef USE_DCO_TYPES
  template<typename T1, typename T2>
  Matrix<typename boost::math::tools::promote_args<T1, T2>::type, Eigen::Dynamic, 1>
  spec_fun(const Eigen::Matrix<T1, Eigen::Dynamic, 1> & xf,
           const Eigen::Matrix<T2,Eigen::Dynamic,1>& x,
           const int& dataset, std::ostream* pstream__) {
    throw std::logic_error("not implemented");  // this should never be called
  }

  // this is a template specialization for <double, double>
  template<>
  Matrix<double,Eigen::Dynamic, 1>
  spec_fun<double, double>(const Eigen::Matrix<double, Eigen::Dynamic, 1>& xf,
                           const Eigen::Matrix<double,Eigen::Dynamic,1>& x,
                           const int& dataset, std::ostream* pstream__) {
    return computation(xf, x, dataset);
  }

  // this is a template specialization for <double, stan::math::var>
  template<>
  Matrix<stan::math::var,Eigen::Dynamic, 1>
  spec_fun<double, stan::math::var>(const Eigen::Matrix<double, Eigen::Dynamic, 1>& xf,
                                    const Eigen::Matrix<stan::math::var,Eigen::Dynamic,1>& x_stan,
                                    const int& dataset, std::ostream* pstream__) {
    Eigen::Matrix<double,Eigen::Dynamic,1> a = value_of(x_stan);
    Eigen::Matrix<double,Eigen::Dynamic,1> fa;
    Matrix<gt1s_scalar, Dynamic, 1> x_dco = a;
    Matrix<gt1s_scalar, Dynamic, 1> spec_dens;

    const int num_inputs = x_dco.rows();
    int num_outputs;
    vector<vector<dco::mode<gt1s_scalar>::value_t>> grad_fa;
    for(int i=0; i<num_inputs; i++) {
      dco::derivative(x_dco(i, 0)) = 1.0;
      spec_dens = computation(xf, x_dco, dataset); // Actual computation to be differentiated
      dco::derivative(x_dco(i, 0)) = 0.0;
      if(i == 0) {
        num_outputs = spec_dens.rows();
        fa.resize(num_outputs);
        grad_fa.resize(num_outputs);
        for(int j=0; j<num_outputs; j++) {
          fa(j) = dco::value(spec_dens(j));
          grad_fa[j].resize(num_inputs);
        }
      }
      for(int j=0; j<num_outputs; j++) grad_fa[j][i] = dco::derivative(spec_dens(j, 0));
    }

    vector<stan::math::var> x_std(x_stan.data(), x_stan.data() + x_stan.rows() * x_stan.cols());
    Eigen::Matrix<stan::math::var,Eigen::Dynamic,1> f_var_jacobian( grad_fa.size() );
    for (int i=0; i < grad_fa.size(); ++i){
      f_var_jacobian(i) = precomputed_gradients(fa(i), x_std, grad_fa[i]);
    }
    return f_var_jacobian;
  }
#endif
}

#endif
