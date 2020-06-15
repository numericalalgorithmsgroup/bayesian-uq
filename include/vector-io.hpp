#ifndef VECTOR_IO_H_
#define VECTOR_IO_H_

template<typename SCALAR_T>
void read_vec(std::ifstream & ifs, Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> & p);

std::ofstream & write_vec(std::ofstream & ofs, const std::vector<double>& data);

#endif
