#include "slice_sample.h"

namespace npbayes
{
	double SliceSample(fun log_f, double x0, 
		double width, double ub, double lb, int miter)
	{
		double y, xl, xr, x;
		y = log_f(x0) + log(Randu());
		xl = x0 - width * Randu();
		xr = xl + width;
		for (int j = 0; j < miter; ++j) {
			if (xl < lb) { xl = lb; break; }
			else if (log_f(xl) > y) xl -= width;
			else break;
		}
		for (int j = 0; j < miter; ++j) {
			if (xr > ub) { xr = ub; break; }
			else if (log_f(xr) > y) xr += width;
			else break;
		}
		x = Randu() * (xr - xl) + xl;
		for (int j = 0; j < miter; ++j) {
			if (log_f(x) > y) break;
			else {
				if (x > x0) xr = x;
				else xl = x;
				x = Randu() * (xr - xl) + xl;
				if (x > ub) x = ub;
				else if (x < lb) x = lb;
			}
		}
		return x;
	}
}