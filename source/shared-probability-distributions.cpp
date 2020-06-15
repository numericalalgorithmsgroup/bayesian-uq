#include <random>
#include "Eigen/Dense"
#define _USE_MATH_DEFINES
#include <math.h>

#include "scalar-typedef.hpp"

#include "shared-probability-distributions.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace Eigen;

// Explicit template instantiation(s)
template double lpdf_imulti_normal(const Eigen::Matrix<double,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<double,Eigen::Dynamic,1> &meanVec,
                            const Eigen::Matrix<double,Eigen::Dynamic,1> &sdVec);

template double lpdf_imulti_lognormal(const Eigen::Matrix<double,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<double,Eigen::Dynamic,1> &meanVec,
                                      const Eigen::Matrix<double,Eigen::Dynamic,1> &sdVec,
                                      unsigned int & ifail);

#ifdef USE_DCO_TYPES
template gt1s_scalar lpdf_imulti_normal(const Eigen::Matrix<gt1s_scalar,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<gt1s_scalar,Eigen::Dynamic,1> &meanVec,
                            const Eigen::Matrix<gt1s_scalar,Eigen::Dynamic,1> &sdVec);
template ga1s_scalar lpdf_imulti_normal(const Eigen::Matrix<ga1s_scalar,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<ga1s_scalar,Eigen::Dynamic,1> &meanVec,
                            const Eigen::Matrix<ga1s_scalar,Eigen::Dynamic,1> &sdVec);
template gt2s_ga1s_scalar lpdf_imulti_normal(const Eigen::Matrix<gt2s_ga1s_scalar,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<gt2s_ga1s_scalar,Eigen::Dynamic,1> &meanVec,
                            const Eigen::Matrix<gt2s_ga1s_scalar,Eigen::Dynamic,1> &sdVec);

template gt1s_scalar lpdf_imulti_lognormal(const Eigen::Matrix<gt1s_scalar,Eigen::Dynamic,1> &x,
                                           const Eigen::Matrix<gt1s_scalar,Eigen::Dynamic,1> &meanVec,
                                           const Eigen::Matrix<gt1s_scalar,Eigen::Dynamic,1> &sdVec,
                                           unsigned int & ifail);
template ga1s_scalar lpdf_imulti_lognormal(const Eigen::Matrix<ga1s_scalar,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<ga1s_scalar,Eigen::Dynamic,1> &meanVec,
                                           const Eigen::Matrix<ga1s_scalar,Eigen::Dynamic,1> &sdVec,
                                           unsigned int & ifail);
template gt2s_ga1s_scalar lpdf_imulti_lognormal(const Eigen::Matrix<gt2s_ga1s_scalar,Eigen::Dynamic,1> &x,
                                                const Eigen::Matrix<gt2s_ga1s_scalar,Eigen::Dynamic,1> &meanVec,
                                                const Eigen::Matrix<gt2s_ga1s_scalar,Eigen::Dynamic,1> &sdVec,
                                                unsigned int & ifail);
#endif

template<typename SCALAR_T>
SCALAR_T lpdf_univariate_normal(const SCALAR_T & x, const SCALAR_T & mu, const SCALAR_T & sd)
{
  const SCALAR_T logSqrt2Pi = 0.5*std::log(2*M_PI);
  const SCALAR_T z = x - mu;
  return - logSqrt2Pi - log(sd) - 0.5 * (z * z) / (sd * sd);
}

template<typename SCALAR_T>
SCALAR_T lpdf_univariate_lognormal(const SCALAR_T & x, const SCALAR_T & mu, const SCALAR_T & sd)
{
  const SCALAR_T logSqrt2Pi = 0.5*std::log(2*M_PI);
  const SCALAR_T z = log(x) - mu;
  return -log(x) - logSqrt2Pi - log(sd) - 0.5 * (z * z) / (sd * sd);
}




template<typename SCALAR_T>
SCALAR_T lpdf_imulti_lognormal(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &x,
                               const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &meanVec,
                               const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &sdVec,
                               unsigned int & ifail)
{
  SCALAR_T lp = 0;
  for (int i=0; i < x.size(); ++i){
    if (x[i] < 0){
      ifail = 7;
      return - std::numeric_limits<SCALAR_T>::infinity();;
    }
    SCALAR_T mu = log( meanVec[i] / sqrt(1 + pow( sdVec[i] / meanVec[i] , 2.0) ) );
    SCALAR_T log_sd = sqrt( log( 1 + pow( sdVec[i] / meanVec[i] , 2.0)  ) );
    lp += lpdf_univariate_lognormal(x[i], mu, log_sd);
  }
  return lp;
}


template<typename SCALAR_T>
SCALAR_T lpdf_imulti_normal(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &meanVec,
                            const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &sdVec)
{
  SCALAR_T lp = 0;
  for (int i=0; i < x.size(); ++i){
    lp += lpdf_univariate_normal(x[i], meanVec[i], sdVec[i]);
  }
  return lp;
}

// downloaded from https://stackoverflow.com/questions/41538095/evaluate-multivariate-normal-gaussian-density-in-c
double lpdf_multi_normal(const Eigen::VectorXd &x, const Eigen::VectorXd &meanVec,
                         const Eigen::MatrixXd &covMat)
{
    // avoid magic numbers in your code. Compilers will be able to compute this at compile time:
    const double logSqrt2Pi = 0.5*std::log(2*M_PI);
    typedef Eigen::LLT<Eigen::MatrixXd> Chol;
    Chol chol(covMat);
    // Handle non positive definite covariance somehow:
    if(chol.info()!=Eigen::Success) throw "decomposition failed!";
    const Chol::Traits::MatrixL& L = chol.matrixL();
    Eigen::VectorXd centredVec = (x - meanVec);
    // solve  L y = (x - mu) for y
    double quadform = (L.solve(centredVec)).squaredNorm();
    // log(|Sigma|^{-1/2} = - 0.5 log( |Sigma| )
    return -x.rows()*logSqrt2Pi - 0.5*quadform - 0.5 * covMat.colPivHouseholderQr().logAbsDeterminant() ;
}

// Note - this is not an optimal implementation if generating many samples from dist with the same covariance matrix
Matrix<double,Dynamic,1> rand_multi_normal(const Matrix<double,Dynamic,1> mean,
                                           Matrix<double,Dynamic,Dynamic> cov,
                                           std::mt19937 & gen)
{
  std::normal_distribution<> std_normal_rng;
  Eigen::LLT<Eigen::MatrixXd> llt_of_cov = cov.llt();
  Eigen::VectorXd z(cov.cols());
  for (int i = 0; i < cov.cols(); ++i){
    z(i) = std_normal_rng(gen);
  }
  MatrixXd L = llt_of_cov.matrixL();
  return mean + L  * z;
}

// Adapted from https://stackoverflow.com/questions/41538095/evaluate-multivariate-normal-gaussian-density-in-c
#undef MYPREC_VERSION
#ifdef MYPREC_VERSION
double lpdf_multi_normal(const Eigen::VectorXd &x, const Eigen::VectorXd &meanVec, const Eigen::MatrixXd &precMat)
{
    // avoid magic numbers in your code. Compilers will be able to compute this at compile time:
    // const double logSqrt2Pi = 0.5*std::log(2*M_PI);
    typedef Eigen::LLT<Eigen::MatrixXd> llt_of_prec;
    llt_of_prec chol(precMat);
    // Handle non positive definite covariance somehow:
    if(chol.info()!=Eigen::Success) throw "decomposition failed!";
    const llt_of_prec::Traits::MatrixL& L = chol.matrixL();
    double quadform = (L * (x - meanVec)).squaredNorm();
    // return std::exp(-x.rows()*logSqrt2Pi - 0.5*quadform) / L.determinant();
    return - 0.5*quadform +  0.5 * precMat.colPivHouseholderQr().logAbsDeterminant();
}
#endif


#undef STACKOVERFLOW_VERSION
#ifdef STACKOVERFLOW_VERSION
// From https://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
// October 2019
struct normal_random_variable
{
    normal_random_variable(Eigen::MatrixXd const& covar)
        : normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
    {}

    normal_random_variable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
        : mean(mean)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    Eigen::VectorXd operator()() const
    {
        static std::mt19937 gen{ std::random_device{}() };
        static std::normal_distribution<> dist;

        return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
    }
};
#endif

#undef MY_COV_VERSION
#ifdef MY_COV_VERSION
// Note - this is less efficient than stackoverflow version if generating many samples from dist with the same covariance matrix
Matrix<double,Dynamic,1> rand_multi_normal(const Matrix<double,Dynamic,1> mean,
                                           Matrix<double,Dynamic,Dynamic> covar,
                                           std::mt19937 & gen,
                                           std::normal_distribution<> & dist)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
  Matrix<double,Dynamic,Dynamic> transform m= eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](double x) { return dist(gen); });
}
#endif

#undef STAN_PREC_VERSION
#ifdef STAN_PREC_VERSION
// From https://mc-stan.org/math/d4/d28/multi__normal__prec__rng_8hpp_source.html
// October 2019

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
  namespace math {

    template <typename T_loc, class RNG>
    inline typename StdVectorBuilder<true, Eigen::VectorXd, T_loc>::type
    multi_normal_prec_rng(const T_loc &mu, const Eigen::MatrixXd &S, RNG &rng) {
      using boost::normal_distribution;
      using boost::variate_generator;

      static const char *function = "multi_normal_prec_rng";

      check_positive(function, "Precision matrix rows", S.rows());
      check_finite(function, "Precision matrix", S);
      check_symmetric(function, "Precision matrix", S);

      Eigen::LLT<Eigen::MatrixXd> llt_of_S = S.llt();
      check_pos_definite(function, "precision matrix argument", llt_of_S);

      vector_seq_view<T_loc> mu_vec(mu);
      check_positive(function, "number of location parameter vectors",
                     mu_vec.size());
      size_t size_mu = mu_vec[0].size();

      size_t N = mu_vec.size();

      for (size_t i = 1; i < N; i++) {
        int size_mu_new = mu_vec[i].size();
        check_size_match(function,
                         "Size of one of the vectors of "
                         "the location variable",
                         size_mu_new,
                         "Size of another vector of the "
                         "location variable",
                         size_mu);
      }

      for (size_t i = 0; i < N; i++) {
        check_finite(function, "Location parameter", mu_vec[i]);
      }

      check_size_match(function, "Rows of location parameter", size_mu, "Rows of S",
                       S.rows());

      StdVectorBuilder<true, Eigen::VectorXd, T_loc> output(N);

      variate_generator<RNG &, normal_distribution<>> std_normal_rng(
                                                                     rng, normal_distribution<>(0, 1));

      for (size_t n = 0; n < N; ++n) {
        Eigen::VectorXd z(S.cols());
        for (int i = 0; i < S.cols(); i++)
          z(i) = std_normal_rng();

        output[n] = Eigen::VectorXd(mu_vec[n]) + llt_of_S.matrixU().solve(z);
      }

      return output.data();
    }

  }  // namespace math
}  // namespace stan
#endif
