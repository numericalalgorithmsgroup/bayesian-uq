#ifndef STABLE_SDE_H_
#define STABLE_SDE_H_

#include "Eigen/Dense"

template<typename SCALAR_T>
class Stable_sde
{
public:
  virtual Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> dx_dt(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x,
                                                           const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & p) = 0;
  virtual Eigen::Matrix<SCALAR_T, Eigen::Dynamic, Eigen::Dynamic> sde_jacobian(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x,
                                                                               const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & p) = 0;
  virtual SCALAR_T eval_sd_in() = 0;
  virtual SCALAR_T eval_sd_obs() = 0;
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> spec_fun(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic, 1> & xf,
                                                      const SCALAR_T & tstep,
                                                      unsigned int & ifail);
  void eigen_dec();
  void eigen_sort(Eigen::Matrix<std::complex<SCALAR_T>, Eigen::Dynamic, Eigen::Dynamic> & V,
                  Eigen::Matrix<std::complex<SCALAR_T>, Eigen::Dynamic, 1> & eigen_vals);
  void print_members();
  void set_d(const int & d_in){
    d = d_in;
    x_star.resize(d);
  }
  void set_x_star(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x_star_in){
    x_star = x_star_in;
  }
  const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> get_x_star() const{
    return x_star;
  }
  const Eigen::Matrix<std::complex<SCALAR_T>, Eigen::Dynamic, 1> get_lambda_vec() const{
    return lambda_vec;
  }
  const Eigen::Matrix<std::complex<SCALAR_T>, Eigen::Dynamic, Eigen::Dynamic> get_jacobian() const{
    return J;
  }
  void set_input_ind(const int & input_ind_in){
    l_ind = input_ind_in;
  }
  void set_output_ind(const int & output_ind_in){
    r_ind = output_ind_in;
  }
  void set_param_unknown(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & p_in){
    p = p_in;
  }
  const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> get_param_unknown() const{
    return p;
  }
  void set_param_fixed(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & fp_in){
    fp = fp_in;
  }
  // solve for steady-state using model parameters
  // calls to these function should be preceded by a calls to set_param_unknown and set_param_fixed
  void find_ss_c0(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x0, unsigned int & ifail);
  void find_ss(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x0, unsigned int & ifail);
protected:
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> p;  // unknown parameters
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> fp;  // fixed parameters
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> x_star;  // equilibrium value
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, Eigen::Dynamic> J; // matrix for jacobian of differential equation
  // matrices for right and left eigenvectors of J
  Eigen::Matrix<std::complex<SCALAR_T>, Eigen::Dynamic, Eigen::Dynamic> R, L;
  // matrix for eigenvalues of J
  Eigen::Matrix<std::complex<SCALAR_T>, Eigen::Dynamic, 1> lambda_vec;
  // dimension of state space / number of variables
  int d = 0;
  // indices for row / column of right / left eigenvector matrix that is used in spectral density calculation
  int r_ind=0, l_ind=0;
  // dimensions for parameter vector
  int n_param_zeroconc=0;
};

#endif
