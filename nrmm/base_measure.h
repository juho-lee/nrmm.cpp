#ifndef BASE_MEASURE_H_
#define BASE_MEASURE_H_

#include "def.h"
#include "suff_stats.h"

namespace npbayes
{
	class BaseMeasure
	{
	public:
		BaseMeasure(void) { }
		virtual ~BaseMeasure(void) { }
		virtual double LogMarginal(SuffStats *ss) const = 0;
		virtual double LogMarginal(SuffStats *ss0, SuffStats *ss1) const = 0;
		double LogPred(SuffStats *ss0, SuffStats *ss1) const 
		{
			if (ss0 == 0) return -INF;
			else if (ss1 == 0) return LogMarginal(ss0);
			else return LogMarginal(ss0, ss1) - LogMarginal(ss1);
		}
	};
}

#endif