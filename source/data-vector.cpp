// core headers for STL
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>  // shared_ptr
#include <fstream>
// Eigen library for linear algebra
#include "Eigen/Dense"

#include <algorithm>    // std::min

// fftw library for DFT
#include "fftw3.h"

#include "scalar-typedef.hpp"

#include "vector-io.hpp"
#include "data-vector.hpp"

using namespace std;
using namespace Eigen;

// Explicit template instantiation
template class Data_vector<double>;
#ifdef USE_DCO_TYPES
template class Data_vector<gt1s_scalar>;
template class Data_vector<ga1s_scalar>;
template class Data_vector<gt2s_ga1s_scalar>;
#endif

template<typename SCALAR_T>
Matrix<SCALAR_T, Dynamic, 1> linseq(const double & start, const double & end, const double & step)
{
  double x_i = start;
  int nsteps = (end - start) / step;
  Matrix<SCALAR_T, Dynamic, 1> x(nsteps);
  for (int i = 0; i < nsteps; ++i){
    x[i] = x_i;
    x_i += step;
  }
  return x;
}

template<typename SCALAR_T>
Data_vector<SCALAR_T>::Data_vector(const string & fstr){

  ifstream ifs(fstr);
  if ( !ifs ){
    throw runtime_error("Error: Data vector input file not found");
  }

  // read in data from file
  ifs >> *this ;

  // evaluate periodogram_scale from time_units
  if (time_units == "seconds"){
    periodogram_scale = 1.0;
  }else if (time_units == "milliseconds"){
    periodogram_scale = 1.0E3;
  }else{
    throw runtime_error("Error: Unrecognized time units");
  }

  // apply DFT to time-series
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> I0 = periodogram();

  // evaluate fstep
  int ts_length = ts.size();
  fstep = 1 / ( (ts_length-1) * time_step);

  // create vector of Fourier frequencies and
  xf = linseq<SCALAR_T>(freq_lim[0], freq_lim[1], fstep);

  // subset periodogram using freq_lim
  n_freq = min(xf.size(), I0.size() );
  I.resize(n_freq);
  int i = 0, k=0;
  while ( (i < xf.size() ) && (i < I0.size()) ){
    if (xf[i] >= freq_lim[0]){
      I[k] = I0[i]; // add periodogram value to Data_vector
      ++k; // increment index for I
    }
    ++i; // increment index for I0
  }
}

template<typename SCALAR_T>
ifstream & operator>>(ifstream & ifs, Data_vector<SCALAR_T> & data)
{
  string str, line, name, unit; // declare variables that stream reads into
  double val, val1, val2;

  // process line for time_step
  getline(ifs,line);
  istringstream ss_tstep(line); // bind record to the line we have just read
  ss_tstep >> str >> name >> val;
  if (name != "time_step"){
    throw runtime_error("Error: Mismatch between name in file and name in class - time_step");
  }else{
    data.time_step = val;
  }

  // process line for time_units
  getline(ifs,line);
  istringstream ss_tunits(line); // bind record to the line we have just read
  ss_tunits >> str >> name >> unit;
  if (name != "time_units"){
    throw runtime_error("Error: Expected time_units string, observed" + name);
  }else{
    data.time_units = unit;
  }

  // process line for freq_lim
  getline(ifs,line);
  istringstream ss_flim(line); // bind record to the line we have just read
  ss_flim >> str >> name >> data.freq_lim[0] >> data.freq_lim[1];
  if (name != "freq_lim"){
    throw runtime_error("Error: Expected freq_lim string, observed" + name);
  }

	data.ts.clear();
	while (ifs >> val){
		data.ts.push_back(val);
  }

	return ifs;

}

template<typename SCALAR_T>
Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> Data_vector<SCALAR_T>::periodogram()
{
  fftw_complex *dft;
  fftw_plan p;

  int N_FFT = ts.size();
  // get pointer to time-series data
  double* p_ts = reinterpret_cast<double*>(&ts[0]);
  // initialize vector that will contain periodogram
  Matrix<SCALAR_T, Dynamic, 1> I = Matrix<SCALAR_T, Dynamic, 1>::Zero(N_FFT / 2);
  // allocate memory for DFT
  dft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_FFT);

  // create a plan
  p = fftw_plan_dft_r2c_1d(N_FFT, p_ts, dft, FFTW_ESTIMATE);

  // execute the FFT plan
  fftw_execute(p);

  // regular f are complex and neg. doubled
  for (int i = 1; i < N_FFT / 2; i++){
    I[i-1] += (dft[i][0] * dft[i][0] + dft[i][1] * dft[i][1]);
  }

  // norm PSD - multiply by observation time-step, divide by # of obs.
  for (int i = 0; i < N_FFT / 2; i++){
		I[i] *= periodogram_scale * time_step / N_FFT;
  }

  // execute
  fftw_execute(p); /* repeat as needed */

  fftw_destroy_plan(p);
  fftw_free(dft);

  return I;
}

template<typename SCALAR_T>
vector<double> Data_vector<SCALAR_T>::pwelch(const vector<double> & win, double overlap)
{
	// determine data windowing
	if (win.size() > ts.size())
		throw runtime_error("Specified window is larger than data.");

  auto N_step = static_cast<unsigned int>(round((1.0 - overlap) * win.size())); // stepping to next window

  if (N_step < 1)
		throw runtime_error("Window size times overlap is too small.");

	auto N_FFT = static_cast<unsigned int>(pow(2.0, ceil(log2(win.size())))); // actual length of zero-padded window
	// FFTW variables for real data size N_FFT -> complex DFT size N_FFT/2 + 1
	// inspace transform, hence reserve real space twice the complex: 2*(NFFT/2+1)
	auto* data_win = static_cast<double *>(fftw_malloc(sizeof(double) * (N_FFT + 2))); // real data

	psd_welch.assign(N_FFT / 2 + 1, 0.0);
	unsigned int L_win = 0; // number of windows
	for (unsigned long pos = 0; pos + win.size() - 1 < ts.size(); pos += N_step) {
		for (unsigned int i = 0; i < win.size(); i++){ // window part of the data
      data_win[i] = ts[pos + i] * win[i];
    }
  }
  return psd_welch;

}

vector<double> gcos_win(unsigned int N, vector<double> a)
{
	const auto pi = 3.14159265358979323846;

	if (a.size() < 1 || N < 2)
		throw runtime_error("Generalized cosine window: bad function call.");
	vector<double> win;
	win.assign(N, a[0]);
	for (unsigned int k = 1; k < a.size(); k++)
		for (unsigned int n = 0; n < N; n++)
			win[n] += a[k] * cos(2.0 * pi * k * n / (N - 1));
	return win;
}
