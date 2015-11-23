#ifndef NRM_H_
#define NRM_H_

#include "cluster.h"
#include "random.h"
#include "base_measure.h"

namespace npbayes
{
	class NRM
	{
	public:
		NRM(BaseMeasure *H) : u_(RandGamma(1.0, 1.0)), H_(H) { }
		~NRM(void) { }

		virtual double LogKappa(int n) const = 0;
		virtual double LogKappaJoin(int n) const = 0;
		virtual double LogKappaNew(void) const = 0;		
		double LogMarginal(SuffStats *ss) const
		{
			return H_->LogMarginal(ss); 
		}
		double LogMarginal(SuffStats *ss0, SuffStats *ss1) const
		{
			return H_->LogMarginal(ss0, ss1);
		}
		double LogPred(SuffStats *ss0, SuffStats *ss1) const
		{
			return H_->LogPred(ss0, ss1); 
		}

		virtual double LogJointConst(const ClusterSet &clset) const { return 0; }
		double LogJoint(const ClusterSet &clset) const {
			double y = LogJointConst(clset);
			foreach(it, clset) y += LogKappa(it->n) + H_->LogMarginal(it->ss);
			return y;
		}

		virtual void SampleU(const ClusterSet &clset) { }
		virtual void SampleHyperparams(const ClusterSet &clset) { }

		////virtual void SampleU(double L, const StickSet &stkset) { }
		////virtual void SampleHyperparams(double L, const StickSet &stkset) { }
		////virtual void SampleJumps(StickSet &stkset) { }
		////virtual void SampleJumps(double L, StickSet &stkset, int maxlen = 1000) { }		
	protected:
		double u_;
		BaseMeasure *H_;
	};
}

#endif