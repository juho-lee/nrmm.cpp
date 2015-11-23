#ifndef VEC_H_
#define VEC_H_

#include <fstream>
#include <iterator>
#include <vector>
#include <Eigen/Dense>

namespace npbayes
{
	// real valued vectors and Matrices
	// vectors are column vectors
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat;

	int ReadMat(const char *filename, Mat &X);
	int WriteMat(const char *filename, const Mat &X);
}

#endif

