#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip> // std::setw
#include "Eigen/Dense"

#include "scalar-typedef.hpp"
#include "vector-io.hpp"

// Explicit template instantiation
template void read_vec<double>(std::ifstream & ifs, Eigen::Matrix<double, Eigen::Dynamic, 1> & p);
#ifdef USE_DCO_TYPES
template void read_vec<gt1s_scalar>(std::ifstream & ifs, Eigen::Matrix<gt1s_scalar, Eigen::Dynamic, 1> & p);
template void read_vec<ga1s_scalar>(std::ifstream & ifs, Eigen::Matrix<ga1s_scalar, Eigen::Dynamic, 1> & p);
template void read_vec<gt2s_ga1s_scalar>(std::ifstream & ifs, Eigen::Matrix<gt2s_ga1s_scalar, Eigen::Dynamic, 1> & p);
#endif

using namespace std;
using namespace Eigen;

// function for inputting parameters from a file
template<typename SCALAR_T>
void read_vec(ifstream & ifs, Matrix<SCALAR_T, Dynamic, 1> & p)
{
  string line, name; // declare variables that stream reads into
  SCALAR_T val;

  // count number of lines in the file
  int n_lines = 0;
  while (getline(ifs,line)){
    ++n_lines; // increment line counter
  }
  // resize vector to match number of lines in file
  p.resize(n_lines,1);

  ifs.clear(); // reset the state of the file stream
  ifs.seekg(0, std::ios::beg); // go back to the beginning of the file
  int i = 0; // start a counter for indexing value vector
  cout << "Reading into Eigen vector..." << endl;
  while (getline(ifs,line)){
    if (!line.empty() ){
      istringstream record(line); // bind record to the line we have just read
      record >> name >> val; // read in parameter name and value
      p[i] = val; // append value to Vector
      cout << left << setw(20) << name << val << endl;
    }
    ++i; // increment counter
  }
  cout << endl;
  if (i == n_lines){
    // cout << "Number of parameters read in = " << p.vals->size() << endl;
  }else{
    throw runtime_error("Error: Number of parameters read in not equal to number of lines in file");
  }
}

ofstream & write_vec(ofstream & ofs, const vector<double> & data)
{
	for (auto i = 0; i < data.size(); i++)
		ofs << data[i] << endl;
	return ofs;
}
