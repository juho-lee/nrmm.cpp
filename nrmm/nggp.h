#ifndef NGGP_H_
#define NGGP_H_

#include "math.h"
#include "nrm.h"
#include "slice_sample.h"

namespace npbayes
{
	class NGGP : public NRM
	{
	public:
		NGGP(BaseMeasure *H) : NRM(H) 
		{ 
			alpha_ = RandGamma(a_alpha_, b_alpha_);
			sigma_ = RandBeta(a_sigma_, b_sigma_);
		}
		~NGGP(void) { }

		double LogKappa(int n) const;
		double LogKappaJoin(int n) const;
		double LogKappaNew(void) const;

		double LogJointConst(const ClusterSet &clset) const;

		void SampleU(const ClusterSet &clset);
		void SampleHyperparams(const ClusterSet &clset);

		//void SampleU(double L, const StickSet &stkset);
		//void SampleHyperParams(double L, const StickSet &stkset);
		//void SampleJumps(StickSet &stkset);
		//void SampleJumps(double L, StickSet &stkset, int maxlen = 1000);
		
		static const double a_alpha_, b_alpha_;
		static const double a_sigma_, b_sigma_;
	private:		
		double alpha_;		
		double sigma_;
	};
}

#endif