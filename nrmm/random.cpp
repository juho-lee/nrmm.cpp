#include "random.h"

namespace npbayes
{
	double Randu(void)
	{
		return ((double)rand() / RAND_MAX);
	}

	int RandMult(const std::vector<double> &p)
	{
        std::vector<double> cp(p.size(), 0);
		cp[0] = p[0];
		int i = 1;
		for (i = 1; i < p.size(); ++i)
			cp[i] = cp[i - 1] + p[i];
		double r = cp[p.size() - 1] * Randu();
		i = -1;
		while (r > cp[++i]);
		return i;
	}

	// http://stackoverflow.com/a/10645091
	double Randn(void)
	{
		double r = 0, u, v;
		while (r == 0 || r > 1) {
			u = 2 * Randu() - 1;
			v = 2 * Randu() - 1;
			r = u*u + v*v;
		}
		return (u * sqrt(-2 * log(r) / r));
	}

	double RandExp(double lambda)
	{
		double u;
		do { u = Randu(); } while (u == 1);
		return -log(u) / lambda;
	}

	/* Algorithm:
	* G. Marsaglia and W.W. Tsang, A simple method for generating gamma
	* variables, ACM Transactions on MatheMatical Software, Vol. 26, No. 3,
	* Pages 363-372, September, 2000.
	* http://portal.acm.org/citation.cfm?id=358414
	*/
	double RandGamma(double a, double b)
	{
		double boost, d, c, v;
		if (a < 1) {
			/* boost using Marsaglia's (1961) method: gam(a) = gam(a+1)*U^(1/a) */
			boost = exp(log(Randu()) / a);
			a++;
		}
		else boost = 1;

		d = a - 1.0 / 3; c = 1.0 / sqrt(9 * d);
		while (1) {
			double x, u;
			do {
				x = Randn();
				v = 1 + c * x;
			} while (v <= 0);
			v = v * v * v;
			x = x * x;
			u = Randu();
			if ((u < 1 - .0331*x*x) || (log(u) < 0.5*x + d*(1 - v + log(v))))
				break;
		}
		return (boost*d*v / b);
	}

	double RandBeta(double a, double b)
	{
		double x = RandGamma(a, 1), y = RandGamma(b, 1);
		return (x / (x + y));
	}
    
	//// sampling from Wishart distribution using Bartlett decomposition
	//// http://en.wikipedia.org/wiki/Wishart_distribution
	//void rand_wishart(int df, const Mat &chol_Scale, Mat &X_chol, Mat &X)
	//{
	//	int d = chol_Scale.rows();
	//	Mat A = Mat::Zero(d, d);
	//	for (int i = 0; i < d; ++i) {
	//		A(i, i) = sqrt(rand_gamma(0.5*(df - i), 0.5));
	//		for (int j = 0; j < i; ++j)
	//			A(i, j) = randn();
	//	}
	//	X_chol = chol_Scale * A;
	//	X = X_chol * X_chol.transpose();
	//}
}