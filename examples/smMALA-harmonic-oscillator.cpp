#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <float.h>
#include <ctime>

#include "Eigen/Dense"
#include "scalar-typedef.hpp"
#include "stable-sde.hpp"
#include "stable-sde-harmonic-oscillator.hpp"
#include "data-vector.hpp"
#include "vector-io.hpp"
#include "computation.hpp"
#include "computation-harmonic-oscillator.hpp"
#include "mh-sm-mala.hpp"

using namespace std;
using namespace Eigen;

// define functions that do likelihood evaluation and reading data
template<typename SCALAR_T>
void read_in_data(Matrix<SCALAR_T, Dynamic, 1> & theta0,
                  Matrix<SCALAR_T, Dynamic, 1> & p,
                  Matrix<SCALAR_T, Dynamic, 1>& x0)
{
  string data_dir = "examples/data/harmonic-oscillator/";
  ifstream ifs1(data_dir + "initial-theta.dat" );
  read_vec(ifs1, theta0);
  ifstream ifs2(data_dir + "initial-param.dat" );
  read_vec(ifs2, p);
  ifstream ifs3(data_dir + "initial-xstar.dat" );
  read_vec(ifs3, x0);
}

void whittle_mcmc_harmonic_oscillator()
{
  int print_flag = 1;
  unsigned int ifail = 0;
  string data_dir = "examples/data/harmonic-oscillator/";

  unsigned int seed = 1733;
  double alpha = 1e6;
  double h = 1.0;

#ifdef USE_BASIC_TYPES
  // read in data
  Data_vector<double> y_light(data_dir + "data1_input.dat");
  Data_vector<double> y_deep(data_dir + "data2_input.dat");
  Matrix<double, Dynamic, 1> theta0, theta, p, x_star;
  double conc_light, conc_deep;
  read_in_data(theta0, p, x_star);
  // initialize computation object with data that is fixed / constant
  Computation_harmonic_oscillator<double> computation(p, 0.0, y_light, 0.0, y_deep);
  MH_move_smMALA<double> sampler_smMALA(computation, theta0, x_star, seed, alpha, h);
  ofstream ofs("examples/samples/harmonic-oscillator/output_basic.csv" );
#elif USE_DCO_TYPES
  // read in data
  Data_vector<gt2s_ga1s_scalar> y_light(data_dir + "data1_input.dat");
  Data_vector<gt2s_ga1s_scalar> y_deep(data_dir + "data2_input.dat");
  Matrix<gt2s_ga1s_scalar, Dynamic, 1> theta0, theta, p, x_star;
  gt2s_ga1s_scalar conc_light, conc_deep;
  read_in_data(theta0, p, x_star);
  // initialize computation object with data that is fixed / constant
  Computation_harmonic_oscillator<gt2s_ga1s_scalar> computation(p, 0.0, y_light, 0.0, y_deep);
  MH_move_smMALA<gt2s_ga1s_scalar> sampler_smMALA(computation, theta0, x_star, seed, alpha, h);
  ofstream ofs("examples/samples/harmonic-oscillator/output_dco.csv" );
#endif

  // Record time
  clock_t begin = clock();

  // MCMC sampling
  cout << sampler_smMALA << endl;
  ofs << "#     algorithm = lmc" << endl;
  ofs << "#       lmc" << endl;
  ofs << "#         engine = smmala" << endl;

  ofs << "lp__,omega0_c1,omega0_c2,sd_in_c1,sd_in_c2,zeta" << endl;

  int n_iter = 1000;
  int n_stride = 100;
  for (int i=0; i < n_iter; ++i){
    ofs << sampler_smMALA.get_lt0() << ",";
    for (int j=0; j < (theta0.size()-1); ++j){
      ofs << exp( theta0(j) ) << "," ;
    }
    ofs << exp( theta0( theta0.size()-1 ) ) << endl;
    if ( (i+1) % n_stride == 0){
      cout << sampler_smMALA(theta0, theta) << endl;
    }else{
      sampler_smMALA(theta0, theta);
    }
  }
  cout << endl;
  cout << "Acceptance rate : " << (double)sampler_smMALA.get_n_acc() / (double)n_iter << endl;
  cout << endl;

  // record time and evaluate time elapsed
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  ofs << "#" << endl;
  ofs << "#  Elapsed Time: " << 0.0 << " seconds (Warm-up)" << endl;
  ofs << "#                " << elapsed_secs << " seconds (Sampling)" << endl;
  ofs << "#" << endl;
}

int main()
{
  try{
    whittle_mcmc_harmonic_oscillator();
  }catch(runtime_error err) {
    cout << err.what() << endl;
  }catch (const std::exception& err) {
    std::cout << " a standard exception was caught, with message '"
              << err.what() << "'\n";
  }

  return 0;
}
