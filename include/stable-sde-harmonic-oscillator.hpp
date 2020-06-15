// Derived class for harmonic oscillator model - (a) single dataset, no transformation
#ifndef STABLE_SDE_HARMONIC_OSCILLATOR_A_H_
#define STABLE_SDE_HARMONIC_OSCILLATOR_A_H_

template<typename SCALAR_T>
class Stable_sde_harmonic_oscillator : public Stable_sde<SCALAR_T>
{
public:
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, Eigen::Dynamic> sde_jacobian(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x,
                                                                       const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & p) override;
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> dx_dt(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x,
                                                   const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & p) override;
  SCALAR_T eval_sd_in() override;
  SCALAR_T eval_sd_obs() override;
  void set_model_vars_c1(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & theta);
  void set_model_vars_c2(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & theta);
  Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> g(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> & theta_post) const;
  Eigen::Matrix<SCALAR_T,Eigen::Dynamic,Eigen::Dynamic> g_jacobian(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> & theta) const;

  using Stable_sde<SCALAR_T>::J;
  using Stable_sde<SCALAR_T>::d;
  using Stable_sde<SCALAR_T>::p;
  using Stable_sde<SCALAR_T>::fp;
  using Stable_sde<SCALAR_T>::x_star;
};

#endif
