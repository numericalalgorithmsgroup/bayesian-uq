#include <vector>
#include <random>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw
#include <math.h>
#include "Eigen/Dense"

#include "scalar-typedef.hpp"
#include "data-vector.hpp"
#include "computation.hpp"
#include "shared-probability-distributions.hpp"

#include "mh-sm-mala.hpp"

using namespace std;
using namespace Eigen;

template<typename SCALAR_T_> std::ostream & operator<<(std::ostream & os, const MH_move_smMALA<SCALAR_T_> & mh_move );

template class MH_move_smMALA<double>;
template std::ostream & operator<<(std::ostream & os, const MH_move_smMALA<double> & mh_move );
#ifdef USE_DCO_TYPES
template class MH_move_smMALA<gt2s_ga1s_scalar>;
template std::ostream & operator<<(std::ostream & os, const MH_move_smMALA<gt2s_ga1s_scalar> & mh_move );
#endif

static inline double square(double x){
  return x*x;
}

MatrixXd cov(const MatrixXd & hess0, const double & h, const double & alpha){
  EigenSolver<MatrixXd> es(hess0);
  VectorXd lambda0 = es.eigenvalues().real();
  MatrixXd V = es.eigenvectors().real();
  MatrixXd hess, reg_mat, C;
  ArrayXd c, lambda;

  double min_hess_eigval = lambda0.minCoeff();
  /* cout << "min_hess_eigval:" << endl; */
  /* cout << min_hess_eigval << endl; */

  c = alpha * lambda0;
  lambda = lambda0.array() / ( alpha * lambda0.array() ).tanh();
  C = pow(h, 2.0) * V * lambda.matrix().asDiagonal().inverse() * V.transpose();

  return C;
}


template<typename SCALAR_T>
MH_move_smMALA<SCALAR_T>::MH_move_smMALA(Computation<SCALAR_T> & computation_in,
                                         Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta0,
                                         const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                                         const unsigned int seed,
                                         const double alpha_in,
                                         const double h_in) :
  computation(computation_in),
  gen(seed), h(h_in), alpha(alpha_in)
{
  int print_flag = 0;
  const int theta_size = theta0.size();
  Matrix<double, Dynamic, 1> theta0_double( theta_size );
#ifdef USE_DCO_TYPES
  for (int i = 0; i < theta_size; ++i) theta0_double(i) = dco::value( dco::value( theta0(i) ) );
#else
  for (int i = 0; i < theta_size; ++i) theta0_double(i) = theta0(i);
#endif

  if (print_flag) cout << "theta0:" << endl;
  if (print_flag)   cout << theta0 << endl;

  xstar_c0_theta0 = computation.find_ss_c0(theta0);
  if (ifail){
    call_string = "find_ss_c0 : ";
  }else{
    // evaluate likelihood
    lt0 = computation.eval(theta0, xstar_c0_theta0, ifail);
    if (ifail) call_string = "eval : ";
  }
  Matrix<double, Dynamic, 1> lt0_grad = VectorXd::Zero(theta_size);
  Matrix<double, Dynamic, Dynamic> lt0_hess = MatrixXd::Zero(theta_size,theta_size);
  if (!ifail){
    // evaluate derivatives of likelihood
     computation.derivatives(theta0, xstar_c0_theta0, lt0_grad, lt0_hess, ifail);
     if (ifail) call_string = "derivatives : ";
    // regularize hessian and invert to get covariance matrix
    C0 = cov(lt0_hess, h, alpha);
  }else{
    mu0 = theta0_double;
    C0 = MatrixXd::Identity(theta0.size(), theta0.size());
  }

  mu0 = theta0_double + 0.5 * C0 * lt0_grad;

  if (print_flag) cout << "lt0_grad = " << endl;
  if (print_flag) cout << lt0_grad << endl;
  if (print_flag) cout << "lt0_hess = " << endl;
  if (print_flag) cout << lt0_hess << endl;
  if (print_flag) cout << "C0 = " << endl;
  if (print_flag) cout << C0 << endl;
}

template<typename SCALAR_T>
MH_move_smMALA<SCALAR_T> & MH_move_smMALA<SCALAR_T>::operator()(Matrix<SCALAR_T, Dynamic, 1> & theta0,
                                                                Matrix<SCALAR_T, Dynamic, 1> & theta)
{
  int print_flag = 0;
  ifail = 0; call_string="";
  // sample parameters
  theta = rand_multi_normal(mu0, C0, gen);

  // evaluate steady state at c=0
  const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> xstar_c0_theta = computation.find_ss_c0(theta);
  if (ifail){
    call_string = "find_ss_c0 : ";
    lt_diff = 0.0; q_diff = 0.0;
    ++counter;
    return *this;
  }

  // evaluate likelihood
  lt = computation.eval(theta, xstar_c0_theta, ifail);
  if (ifail){
    call_string = "eval : ";
    lt_diff = 0.0; q_diff = 0.0;
    ++counter;
    return *this;
  }
  if (isnan(lt)){
    call_string = "eval : ";
    ++counter;
    return *this;
  }

  // evaluate derivatives using template specializations
  Matrix<double, Dynamic, 1> lt_grad;
  Matrix<double, Dynamic, Dynamic> lt_hess;
  computation.derivatives(theta, xstar_c0_theta, lt_grad, lt_hess, ifail);
  if (ifail){
    call_string = "grad : ";
    lt_diff = 0.0; q_diff = 0.0;
    ++counter;
    return *this;
  }

  // compute mean and covariance of proposal
  const int theta_size = theta.size();
  Matrix<double, Dynamic, 1> theta_double( theta_size ), theta0_double( theta_size );
#ifdef USE_DCO_TYPES
  for (int i = 0; i < theta_size; ++i) theta_double(i) = dco::value( dco::value( theta(i) ) );
  for (int i = 0; i < theta_size; ++i) theta0_double(i) = dco::value( dco::value( theta0(i) ) );
  lt_diff = dco::value( dco::value(lt) ) - dco::value( dco::value(lt0) );
#else
  for (int i = 0; i < theta_size; ++i) theta_double(i) = theta(i);
  for (int i = 0; i < theta_size; ++i) theta0_double(i) = theta0(i);
  lt_diff = lt - lt0;
#endif

  Matrix<double, Dynamic, Dynamic> C = cov(lt_hess, h, alpha);
  Matrix<double, Dynamic, 1> mu = theta_double + 0.5 * C * lt_grad;


  // printing
  if (print_flag)   cout << "theta0:" << endl;
  if (print_flag)   cout << theta0_double << endl;
  if (print_flag)   cout << "mu:" << endl;
  if (print_flag)   cout << mu << endl;
  if (print_flag)   cout << "C:" << endl;
  if (print_flag)   cout << C << endl;

  if (print_flag)   cout << "theta:" << endl;
  if (print_flag)   cout << theta_double << endl;
  if (print_flag)   cout << "mu0:" << endl;
  if (print_flag)   cout << mu0 << endl;
  if (print_flag)   cout << "C0:" << endl;
  if (print_flag)   cout << C0 << endl;

  // evaluate log proposal densities
  double q_theta0 = lpdf_multi_normal(theta0_double, mu, C);
  double q_theta  = lpdf_multi_normal(theta_double, mu0, C0);

  q_diff = q_theta0 - q_theta;
  aprob = fmin(1.0, exp( lt_diff + q_diff ) );

  // evaluate acceptance probability and accept / reject
  std::uniform_real_distribution<double> std_uniform_rng(0.0, 1.0);
  if (print_flag) cout << "u = " << std_uniform_rng(gen) << endl;
  if (aprob >  std_uniform_rng(gen)){
    theta0 = theta;
    lt0 = lt; mu0 = mu; C0 = C;
    n_acc += 1;
  }

  ++counter;
  return *this;
}

template<typename SCALAR_T_>
std::ostream & operator<<(std::ostream & os, const MH_move_smMALA<SCALAR_T_> & mh_move )
{
  if (mh_move.counter == 0){
    os << left << setw(16) << "it";
    os << left << setw(16) << "lt0";
    os << left << setw(16) << "lt_diff";
    os << left << setw(16) << "q_diff";
    os << endl;

    os << left << setw(16) << mh_move.counter;
    os << left << setw(16) << mh_move.lt0;
    os << left << setw(16) << "-";
    os << left << setw(16) << "-";
    os << left << setw(50) << mh_move.call_string + mh_move.computation.report(mh_move.ifail);
  }else{
    os << left << setw(16) << mh_move.counter;
    os << left << setw(16) << mh_move.lt0;
    os << left << setw(16) << mh_move.lt_diff;
    os << left << setw(16) << mh_move.q_diff;
    os << left << setw(50) << mh_move.call_string + mh_move.computation.report(mh_move.ifail);
  }
  return os;
}

template<typename SCALAR_T>
Matrix<double, Dynamic, Dynamic> MH_move_smMALA<SCALAR_T>::hess_to_prec(const Matrix<double, Dynamic, Dynamic> & hess0){
  int print_flag = 0;
  EigenSolver<MatrixXd> es(hess0);
  VectorXd lambda0 = es.eigenvalues().real();
  MatrixXd V = es.eigenvectors().real();
  MatrixXd hess, reg_mat, C;
  ArrayXd c, lambda;
  min_hess_eigval = lambda0.minCoeff();
  if (print_flag) cout << "alpha = " << endl;
  if (print_flag) cout << alpha << endl;
  lambda = lambda0.array() / ( alpha * lambda0.array() ).tanh();
  if (print_flag) cout << "lambda = " << endl;
  if (print_flag) cout << lambda << endl;
  return V * lambda.matrix().asDiagonal() * V.transpose();
}
