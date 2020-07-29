#ifndef MH_SM_MALA_H_
#define MH_SM_MALA_H_

template<typename SCALAR_T>
class MH_move_smMALA
{
  template<typename SCALAR_T_>
  friend std::ostream & operator<<(std::ostream & os, const MH_move_smMALA<SCALAR_T_> & mh_move );
public:
  MH_move_smMALA(Computation<SCALAR_T> & computation_in,
                 Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta0,
                 const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                 const unsigned int seed,
                 const double alpha_in,
                 const double h_in);
  void set_h(const double & h_in){
    h = h_in;
  }
  void set_alpha(const double & alpha_in){
    alpha = alpha_in;
  }
  MH_move_smMALA & operator()(Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & theta0,
                              Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & theta);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_to_prec(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & hess0);
  std::string report(const unsigned int & ifail) const;
  int get_n_acc() const { return n_acc; };
  SCALAR_T get_lt0() const { return lt0; };
protected:
  unsigned int ifail = 0;
  Computation<SCALAR_T> & computation;
  // variables for printing
  unsigned int counter=0;
  std::string call_string="";
  SCALAR_T lt0, lt;
  // equilibrium state at c=0
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> xstar_c0_theta0;

  // derivatives of posterior
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C0;
  // RNG variables for multivariate normal sampling
  std::mt19937 gen;
  double h;
  double alpha, min_hess_eigval;
  double lt_diff, q_diff, aprob;
  int n_acc = 1;
};

#endif
