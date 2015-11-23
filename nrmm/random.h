#ifndef RANDOM_H_
#define RANDOM_H_

#include <ctime>
#include <cstdlib>
#include <cMath>
#include <vector>
#include <map>

namespace npbayes
{
	// random uniform
	double Randu(void);

	// random unit normal
	double Randn(void);
	
	// random exponential
	double RandExp(double lambda);

	// random multinomial
	int RandMult(const std::vector<double> &p);

	// random gamma
	double RandGamma(double a, double b);

	// random beta
	double RandBeta(double a, double b);
    
	//// random wishart
	//void rand_wishart(int df, const Mat &chol_Scale, Mat &X_chol, Mat &X);
}


#endif