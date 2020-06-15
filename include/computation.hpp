#ifndef SHARED_FUNCTIONS_LIKELIHOOD_H_
#define SHARED_FUNCTIONS_LIKELIHOOD_H_

#include "stable-sde.hpp"

enum Likelihood_mode {CheckStability = 0,
                      EvalWhittle = 1};

template<typename SCALAR_T>
SCALAR_T log_whittle_like(Likelihood_mode LikelihoodMode,
                        Stable_sde<SCALAR_T> & m,
                        const Data_vector<SCALAR_T> & y,
                        unsigned int & ifail);

template<typename SCALAR_T>
class Computation
{
public:
  Computation(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& p_in,
              const Data_vector<SCALAR_T> & y_in) :
    p(p_in),
    y_c0(y_in)
  { };
  Computation(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& p_in,
              const double & conc_light_in,
              const Data_vector<SCALAR_T> & y_light_in,
              const double & conc_deep_in,
              const Data_vector<SCALAR_T> & y_deep_in) :
    p(p_in),
    conc_light(conc_light_in), y_light(y_light_in),
    conc_deep(conc_deep_in), y_deep(y_deep_in)
  { };

  virtual Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> find_ss_c0(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta) = 0;

  virtual SCALAR_T eval(Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta,
                        const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                        unsigned int & ifail)=0;


  Eigen::Matrix<double, Eigen::Dynamic, 1> grad(Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta0,
                                                const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                                                unsigned int & ifail);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess(Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta0,
                                                             const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                                                             unsigned int & ifail);

  void derivatives(Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta0,
                   const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                   Eigen::Matrix<double,Eigen::Dynamic,1> & grad,
                   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> & hess,
                   unsigned int & ifail);
  std::string report(const unsigned int & ifail) const;
protected:
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> p;
  Data_vector<SCALAR_T> y_c0;
  double conc_light;
  Data_vector<SCALAR_T> y_light;
  double conc_deep;
  Data_vector<SCALAR_T> y_deep;
};

#endif
