// Derived class
#ifndef COMPUTATION_A_H_
#define COMPUTATION_A_H_

template<typename SCALAR_T>
class Computation_harmonic_oscillator : public Computation<SCALAR_T>
{
public:
  Computation_harmonic_oscillator(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& fixed_p_in,
                const double & conc_light_in,
                const Data_vector<SCALAR_T> & y_light_in,
              const double & conc_deep_in,
                const Data_vector<SCALAR_T> & y_deep_in) :
    Computation<SCALAR_T>(fixed_p_in,
                          conc_light_in,
                          y_light_in,
                          conc_deep_in,
                          y_deep_in) { };
  SCALAR_T eval(Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta,
                const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& xstar_c0_in,
                unsigned int & ifail) override;
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> find_ss_c0(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1>& theta) override;
protected:
  using Computation<SCALAR_T>::p;
  using Computation<SCALAR_T>::conc_light;
  using Computation<SCALAR_T>::y_light;
  using Computation<SCALAR_T>::conc_deep;
  using Computation<SCALAR_T>::y_deep;
};


#endif
