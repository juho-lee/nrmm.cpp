#include "Math.h"

namespace npbayes
{
	// log(Gamma(x))
	// Written by Tom Minka, in lightspeed MATLAB toolbox 
	// http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
	double LogGamma(double x)
	{
		const static double ln_sqrt_2PI = 0.91893853320467274178;
		static double gamma_series[] = {
			76.18009172947146,
			-86.50532032941677,
			24.01409824083091,
			-1.231739572450155,
			0.1208650973866179e-2,
			-0.5395239384953e-5
		};

		int i;
		double denom, x1, series;
		if (x <= 0) return INF;
		if (x == 1 || x == 2) return 0;

		/* Lanczos method */
		denom = x + 1;
		x1 = x + 5.5;
		series = 1.000000000190015;
		for (i = 0; i < 6; i++) {
			series += gamma_series[i] / denom;
			denom += 1.0;
		}

		return (ln_sqrt_2PI + (x + 0.5) * log(x1) - x1 + log(series / x));
	}

	// Regularized upper incomplete gamma function
	// Q(a, z) = Gamma(a, z) / Gamma(z)
	// Written by Kai Zhang, copied from
	// https://sites.google.com/site/kaizhangstatmech/code/special-functions

	////////////////////////////////////////////////////////////////////////////////////
	// Returns the incomplete gamma function P(a; x) 
	// evaluated by its series representation as gamser.
	// Also returns ln Gamma(a) as gln.
	void gser(double *gamser, double a, double x, double *gln)
	{		
		const int ITMAX = 100;
		const double EPS = 3.0e-7;
		void nrerror(char error_text[]);
		int n;
		double sum, del, ap;
		*gln = LogGamma(a);
		if (x <= 0.0) {
			// if (x < 0.0) printf("x less than 0 in routine gser!!!\n");
			*gamser = 0.0;
			return;
		}
		else {
			ap = a;
			del = sum = 1.0 / a;
			for (n = 1; n <= ITMAX; n++) {
				++ap;
				del *= x / ap;
				sum += del;
				if (fabs(del) < fabs(sum)*EPS) {
					*gamser = sum*exp(-x + a*log(x) - (*gln));
					return;
				}
			}
		}
		return;
	}

	// Returns the incomplete gamma function Q(a; x) 
	// evaluated by its continued fraction representation as gammcf. 
	// Also returns ln Gamma(a) as gln.
	void gcf(double *gammcf, double a, double x, double *gln)
	{
		const int ITMAX = 100;
		const double EPS = 3.0e-7;
		const double FPMIN = 1.0e-30;		

		int i;
		double an, b, c, d, del, h;

		*gln = LogGamma(a);
		//SuffStats up for evaluating continued fraction by modified Lentz's method (x5.2) with b0 = 0.
		b = x + 1.0 - a; 
		c = 1.0 / FPMIN;
		d = 1.0 / b;
		h = d;
		//Iterate to convergence.
		for (i = 1; i <= ITMAX; i++) {   
			an = -i*(i - a);
			b += 2.0;
			d = an*d + b;
			if (fabs(d) < FPMIN) d = FPMIN;
			c = b + an / c;
			if (fabs(c) < FPMIN) c = FPMIN;
			d = 1.0 / d;
			del = d*c;
			h *= del;
			if (fabs(del - 1.0) < EPS) break;
		}
		// if (i > ITMAX) printf("a too large, ITMAX too small in gcf!!!\n");
		//Put factors in front.
		*gammcf = exp(-x + a*log(x) - (*gln))*h; 
	}

	double Q(double a, double x)
	{
		if (a < 0) {
			return (Q(1 + a, x) - (pow(x, a) * exp(-x - LogGamma(1 + a))));
		}
		double gamser, gammcf, gln;
		if (x < a + 1.0) {
			gser(&gamser, a, x, &gln);
			return 1.0 - gamser;
		}
		else {
			gcf(&gammcf, a, x, &gln);
			return gammcf;
		}		
	}
	////////////////////////////////////////////////////////////////////////////////////

	double LogSumExp(double a, double b)
	{
        if (a == -INF && b == -INF) return -INF;
		if (a > b) return a + log(1 + exp(b - a));
		else return b + log(1 + exp(a - b));
	}
}