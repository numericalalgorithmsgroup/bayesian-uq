#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <iostream>
#include <numeric>
#include "Eigen/Dense"

#include "scalar-typedef.hpp"
#include "stable-sde.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Eigen;

// Explicit template instantiation
template class Stable_sde<double>;
#ifdef USE_DCO_TYPES
template class Stable_sde<gt1s_scalar>;
template class Stable_sde<ga1s_scalar>;
template class Stable_sde<gt2s_ga1s_scalar>;
#endif

template<typename SCALAR_T>
void Stable_sde<SCALAR_T>::eigen_dec()
{
  // evaluate Jacobian of SDE drift term,
  sde_jacobian(x_star, p);

	// calculate eigendecomposition for right eigenvectors
  EigenSolver<Matrix<SCALAR_T, Dynamic, Dynamic> > es_r(J);
	R = es_r.eigenvectors();
	lambda_vec = es_r.eigenvalues();
  eigen_sort(R, lambda_vec);

  // normalize eigenvectors so that r ^ T * r = 1
  for (int j=0; j < d; j++){
    complex<SCALAR_T> rprod = R.col(j).transpose() * R.col(j);
    R.col(j) = R.col(j) / sqrt(rprod);
  }

	// calculate eigendecomposition for left eigenvectors
  EigenSolver<Matrix<SCALAR_T, Dynamic, Dynamic> > es_l(J.transpose());
	L = es_l.eigenvectors();
	lambda_vec = es_l.eigenvalues();
  eigen_sort(L, lambda_vec);

  // normalize eigenvectors so that L * R = I
  for (int j=0; j < d; j++){
    complex<SCALAR_T> rprod = L.col(j).transpose() * R.col(j);
    L.col(j) = L.col(j) / rprod;
  }
  L.transposeInPlace();
}

// make sure that right and left eigenvectors have the same ordering
template<typename SCALAR_T>
void Stable_sde<SCALAR_T>::eigen_sort(Matrix<std::complex<SCALAR_T>, Dynamic, Dynamic> & V,
                                      Matrix<std::complex<SCALAR_T>, Dynamic, 1> & eigen_vals)
{
  VectorXi idx(eigen_vals.size());
  iota(idx.data(), idx.data()+idx.size(), 0);
  sort(idx.data(), idx.data()+idx.size(),
       [&eigen_vals](int i1, int i2)
       {
         SCALAR_T r1 = abs(eigen_vals(i1)), r2 = abs(eigen_vals(i2));
         SCALAR_T im1 = imag( eigen_vals(i1) ), im2 = imag(eigen_vals(i2));
         SCALAR_T re1 = real( eigen_vals(i1) ), re2 = real(eigen_vals(i2));
         if ( abs( r1 - r2) > 1E-6 ){
           // if sizes of complex numbers are different, sort from larger to smaller size
           return r1 > r2;
         }else if(  abs( im1 - im2) > 1E-6  ){
           // if sizes of complex numbers are equal, sort from larger to smaller imaginary part
           return imag( eigen_vals(i1) ) > imag( eigen_vals(i2) );
         }else{
           // if sizes of complex numbers are equal and the imaginary parts are equal,
           // sort from larger to smaller real part
           return real( eigen_vals(i1) ) > real( eigen_vals(i2) );
         }
       }
       );
  PermutationMatrix<Dynamic,Dynamic> P( idx );
  V *= P;
  eigen_vals = eigen_vals.transpose() * P;

}


template<typename SCALAR_T>
Matrix<SCALAR_T, Dynamic, 1> Stable_sde<SCALAR_T>::spec_fun(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & xf,
                                                            const SCALAR_T & tstep,
                                                            unsigned int & ifail)
{
  int print_flag=0;
  // initialize variables for spectral density calculation
  int n_freq = xf.size();
  Matrix<SCALAR_T, Dynamic, 1> spec_dens(n_freq);
  // compute eigenvalues and eigenvectors
  eigen_dec();
  if ( lambda_vec.real().maxCoeff() > 0){
    if (print_flag) cout << lambda_vec.real();
    // fixed point is unstable
    ifail = 1;
    return spec_dens;
  }
  // define variables for components of eigenvectors that contribute to spectral density
  Matrix<complex<SCALAR_T>, Dynamic, 1> r = R.row(r_ind), l = L.col(l_ind);
  for (int i=0; i < n_freq; i++){
    complex<SCALAR_T> c = static_cast<SCALAR_T>(2*M_PI * xf[i]) * 1i;
    Matrix<complex<SCALAR_T>, Dynamic, 1> diag_vec = (Matrix<complex<SCALAR_T>, Dynamic, 1>::Constant(d,c) - lambda_vec).cwiseInverse();
    complex<SCALAR_T> trans_func = (r.transpose() * diag_vec.asDiagonal() * l)(0);

    // evaluate spectral density
    spec_dens(i) = pow(abs(trans_func),2.0) * pow( eval_sd_in(), 2.0 );
    spec_dens(i) += pow( eval_sd_obs(), 2.0) * tstep  ;
    // spec_dens(i) *=2;
  }
  return spec_dens;
}

template<typename SCALAR_T>
void Stable_sde<SCALAR_T>::find_ss(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x0,
                                   unsigned int & ifail){
  int print_flag = 0;
  // evaluate RHS of differential equation model dx / dt = F(x)
  Matrix<SCALAR_T, Dynamic, 1> x = x0;
  Matrix<SCALAR_T, Dynamic, 1> e = VectorXd::Ones(x.size() );
  int i = 0;
  while (e.norm() > 1e-6 && i < 100){
    Matrix<SCALAR_T, Dynamic, 1> F_x = dx_dt(x, p) ;
    Matrix<SCALAR_T, Dynamic, Dynamic> J = sde_jacobian(x, p);
    if (print_flag) cout << "F_x = " << endl;
    if (print_flag) cout << F_x << endl;
    e = J.colPivHouseholderQr().solve(-F_x);
    x += e;
    ++i;
  }
  if ((e.norm() > 1e-6 ) || ( isnan( e.norm() ) ) ){
    ifail = 2;
  } else {
    // only update x_star if solver has converged
    x_star = x;
  }
}

/* template<typename SCALAR_T> */
/* void Stable_sde<SCALAR_T>::find_ss_c0(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & x0, */
/*                                    unsigned int & ifail){ */
/*   // evaluate RHS of differential equation model dx / dt = F(x) */
/*   Matrix<SCALAR_T, Dynamic, 1> x = x0; */
/*   Matrix<SCALAR_T, Dynamic, 1> e = VectorXd::Ones(x.size() ); */
/*   int i = 0; */
/*   while (e.norm() > 1e-6 && i < 100){ */
/*     Matrix<SCALAR_T, Dynamic, 1> F_x = dx_dt_c0(x) ; */
/*     Matrix<SCALAR_T, Dynamic, Dynamic> J = sde_jacobian_c0(x); */
/*     e = J.colPivHouseholderQr().solve(-F_x); */
/*     x += e; */
/*     ++i; */
/*   } */
/*   if ((e.norm() > 1e-6 ) || ( isnan( e.norm() ) ) ){ */
/*     ifail = 2; */
/*   } else { */
/*     // only update x_star if solver has converged */
/*     x_star = x; */
/*   } */
/* } */


template<typename SCALAR_T>
void Stable_sde<SCALAR_T>::print_members(){
  cout << "Jacobian:" << endl;
  cout << J << endl;
  cout << "Eigenvalues:" << endl;
  cout << lambda_vec << endl << endl;
  cout << "Right eigenvector:" << endl;
  cout << R.row(r_ind).transpose() << endl << endl;
  cout << "Left eigenvector:" << endl;
  cout << L.col(l_ind) << endl << endl;
  cout << "L * R " << endl;
  cout << (L * R).real() << endl << endl;
}
