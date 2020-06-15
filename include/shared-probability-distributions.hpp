#ifndef SHARED_FUNCTIONS_DISTRIBUTIONS_H_
#define SHARED_FUNCTIONS_DISTRIBUTIONS_H_

template<typename SCALAR_T>
SCALAR_T lpdf_imulti_lognormal(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &x,
                               const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &meanVec,
                               const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &sdVec,
                               unsigned int & ifail);

template<typename SCALAR_T>
SCALAR_T lpdf_imulti_normal(const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &x,
                            const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &meanVec,
                            const Eigen::Matrix<SCALAR_T,Eigen::Dynamic,1> &sdVec);

double lpdf_multi_normal(const Eigen::VectorXd &x, const Eigen::VectorXd &meanVec,
                         const Eigen::MatrixXd &covMat);

Eigen::Matrix<double,Eigen::Dynamic,1> rand_multi_normal(const Eigen::Matrix<double,Eigen::Dynamic,1> mean,
                                                         Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> prec,
                                                         std::mt19937 & gen);

#endif
