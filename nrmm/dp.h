#ifndef DP_H_
#define DP_H_

#include "Math.h"
#include "nrm.h"

namespace npbayes
{
	class DP : public NRM
	{
	public:
		DP(BaseMeasure *H) : NRM(H), alpha_(1), log_alpha_(0) { }
		~DP(void) { }
		double LogKappa(int n) const { return log_alpha_ + LogGamma(n); }
		double LogKappaJoin(int n) const { return log(n); }
		double LogKappaNew(void) const { return log_alpha_; }

		double LogJointConst(const ClusterSet &clset) const
		{
			int n = 0;
			foreach(it, clset) n += it->n;
			return LogGamma(alpha_) - LogGamma(alpha_ + n);
		}
	private:
		double alpha_;
		double log_alpha_;
	};
}

#endif