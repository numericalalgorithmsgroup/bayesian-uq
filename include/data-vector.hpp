#ifndef DATA_VECTOR_H_
#define DATA_VECTOR_H_

template<typename SCALAR_T>
Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> linseq(const double & start, const double & end, const double & step);

template<typename SCALAR_T>
class Data_vector
{
  template<typename Scalar>
  friend std::ifstream & operator>>(std::ifstream & ifs, Data_vector<Scalar> & p);
  friend std::ofstream & operator<<(std::ofstream & ofs, const Data_vector<double> & p);
public:
  Data_vector() = default;
  Data_vector(const std::string & fstr);
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> periodogram();
  std::vector<double> pwelch(const std::vector<double> & win, double overlap);
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> get_freq() const {
    return xf;
  }
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> get_periodogram() const {
    return I;
  }
  int get_nfreq() const {
    return n_freq;
  }
  double get_timestep() const {
    return time_step;
  }
private:
  int n_freq = 0;
  std::vector<double> ts; // time-series
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> I;  // periodogram
  std::vector<double> psd_welch;  // Welch estimate of spectral density
  double time_step;
  std::string time_units;
  double fstep;
  Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> xf;
  double periodogram_scale;
  std::vector<double> freq_lim = std::vector<double>(2);
};

std::vector<double> gcos_win(unsigned int N, std::vector<double> a);

#endif
