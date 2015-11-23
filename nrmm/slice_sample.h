#ifndef SLICE_SAMPLE_H_
#define SLICE_SAMPLE_H_

// implementation of slice sampler (Neal 2003)

#include <functional>
#include "def.h"
#include "random.h"

namespace npbayes
{
	typedef std::function<double(double)> fun;
	using namespace std::placeholders;
	double SliceSample(fun log_f, double x0, double width = 1,
		double ub = INF, double lb = -INF, int miter = 200);
}


#endif