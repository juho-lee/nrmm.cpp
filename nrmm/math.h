#ifndef MATH_H_
#define MATH_H_

#include <vector>
#include "def.h"

namespace npbayes {
	// log(Gamma(x))
	// Written by Tom Minka, in lightspeed MATLAB toolbox 
	// http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
	double LogGamma(double x);

	// Regularized upper incomplete gamma function
	// Q(a, x) = Gamma(a, x) / Gamma(x)
	// Written by Kai Zhang, copied from
	// https://sites.google.com/site/kaizhangstatmech/code/special-functions
	double Q(double a, double x);

	// log(exp(a) + exp(b))
	double LogSumExp(double a, double b);	
}

#endif